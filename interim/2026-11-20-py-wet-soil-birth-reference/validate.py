"""Explicit independent checks plus unchanged existing CI Python inventories."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import unittest
import checks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GROUPS = ('python-core-0', 'python-core-1', 'python-mean-evidence', 'python-spatial-evidence')


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def sources():
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    paths = [ROOT/p for p in tracked if p and (Path(p).suffix in ('.py', '.rs', '.toml', '.yml', '.yaml')
                                              or p.endswith('Cargo.lock'))]
    paths += list(HERE.glob('*.py'))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(paths))}


def flatten(suite):
    for test in suite:
        if isinstance(test, unittest.TestSuite):
            yield from flatten(test)
        else:
            yield test


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    frozen = sources()
    write(out/'sources-before.json', frozen)
    suite = unittest.defaultTestLoader.loadTestsFromModule(checks)
    ids = sorted(test.id() for test in flatten(suite))
    with (out/'new-checks.log').open('w', encoding='utf-8') as log:
        result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    new = dict(tests_run=result.testsRun, test_ids=ids, failures=len(result.failures), errors=len(result.errors),
               skipped=len(result.skipped), expected_failures=len(result.expectedFailures),
               unexpected_successes=len(result.unexpectedSuccesses))
    valid = result.wasSuccessful() and result.testsRun == len(ids) and not (
        result.skipped or result.expectedFailures or result.unexpectedSuccesses)
    reports = {}
    for group in GROUPS:
        folder = out/group
        with (out/(group+'.log')).open('w', encoding='utf-8') as log:
            code = subprocess.call([sys.executable, str(ROOT/'tools/run_ci_python_tests.py'),
                                    '--group', group, '--output', str(folder)], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        report = json.loads((folder/'python-test-inventory.json').read_text(encoding='utf-8'))
        reports[group] = report
        valid &= code == 0 and report['exit_code'] == 0
    inventories = [report['discovered_test_ids'] for report in reports.values()]
    valid &= all(x == inventories[0] for x in inventories)
    final_inventory = out/'final-inventory'
    code = subprocess.call([sys.executable, str(ROOT/'tools/run_ci_python_tests.py'), '--group', GROUPS[0],
                            '--output', str(final_inventory), '--inventory-only'], cwd=ROOT)
    after_ids = json.loads((final_inventory/'python-test-inventory.json').read_text())['discovered_test_ids']
    valid &= code == 0 and after_ids == inventories[0]
    after = sources()
    write(out/'sources-after.json', after)
    valid &= frozen == after
    total = sum(report.get('tests_run', 0) for report in reports.values())
    executed_ids = [identity for r in reports.values() for identity in r['all_group_test_ids'][r['group']]]
    valid &= sorted(executed_ids) == inventories[0] and len(executed_ids) == len(set(executed_ids))
    summary = dict(validation_success=bool(valid), a03_complete=False, accepted=None,
                   new_checks=new, existing_groups={g: {k: v for k, v in r.items() if k not in
                       ('discovered_test_ids', 'all_group_test_ids')} for g, r in reports.items()},
                   existing_tests=total, total_tests=total+result.testsRun,
                   existing_collection_unchanged=after_ids == inventories[0],
                   existing_ids_partitioned_exactly_once=sorted(executed_ids) == inventories[0],
                   frozen_sources=len(frozen), source_hashes_unchanged=frozen == after,
                   python=sys.version, elapsed_seconds=time.monotonic()-start,
                   windows_ci='未確認')
    write(out/'validation.json', summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return int(not valid)


if __name__ == '__main__':
    raise SystemExit(main())
