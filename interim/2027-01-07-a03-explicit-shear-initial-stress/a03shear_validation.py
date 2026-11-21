"""Full-suite execution with explicit inventories, outcomes and immutable logs."""
import argparse
from collections import Counter
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from a03shear_run import write

ROOT = Path(__file__).resolve().parents[2]
GROUPS = ('python-core-0', 'python-core-1', 'python-mean-evidence', 'python-spatial-evidence',
          'python-c01-switch', 'python-c01-path', 'python-c01-cause')


def command(output, name, args):
    started = time.monotonic()
    with (output/(name+'.log')).open('x', encoding='utf-8') as log:
        process = subprocess.Popen(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        code = process.wait()
    result = dict(command=args, pid=process.pid, exit_code=code, elapsed_seconds=time.monotonic()-started)
    write(output/(name+'-process.json'), result)
    print(name, code, flush=True)
    return result


def rust(output):
    output.mkdir(parents=True, exist_ok=False)
    processes = {}
    for name, args in (
        ('fmt', ['cargo', 'fmt', '--all', '--check']),
        ('clippy', ['cargo', 'clippy', '--locked', '--all-targets', '--', '-D', 'warnings']),
        ('build', ['cargo', 'build', '--locked', '--release', '--all-targets']),
        ('discovery', ['cargo', 'test', '--locked', '--release', '--all-targets', '--no-run', '--message-format=json']),
        ('rust-all', ['cargo', 'test', '--locked', '--release', '--all-targets']),
        ('rust-doc', ['cargo', 'test', '--locked', '--release', '--doc'])):
        processes[name] = command(output, name, args)
    inventory = []
    for line in (output/'discovery.log').read_text().splitlines():
        try:
            value = json.loads(line)
        except ValueError:
            continue
        if value.get('reason') == 'compiler-artifact' and value.get('executable') and value['profile']['test']:
            target = value['target']['kind'][0]+':'+value['target']['name']
            listed = subprocess.check_output([value['executable'], '--list', '--format', 'terse'], text=True)
            inventory.extend(target+'::'+v[:-6] for v in listed.splitlines() if v.endswith(': test'))
    text = (output/'rust-all.log').read_text()
    outcomes = [dict(zip(('status','passed','failed','ignored','measured','filtered_out'),
                [match[0], *map(int, match[1:])])) for match in re.findall(
                r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out', text)]
    executed = re.findall(r'^test (.*?) \.\.\. (?:ok|FAILED|ignored)', text, re.M)
    expected_names = Counter(v.split('::', 1)[1] for v in inventory)
    actual_names = Counter(executed)
    report = dict(processes=processes, discovered_ids=sorted(inventory), discovered=len(inventory),
        executed=len(executed), passed=sum(r['passed'] for r in outcomes), failed=sum(r['failed'] for r in outcomes),
        ignored=sum(r['ignored'] for r in outcomes), filtered_out=sum(r['filtered_out'] for r in outcomes),
        measured=sum(r['measured'] for r in outcomes), skipped=0,
        unexecuted=sum((expected_names-actual_names).values()),
        unexpected=sum((actual_names-expected_names).values()),
        discovered_and_executed_multisets_equal=expected_names == actual_names,
        test_results=outcomes, all_processes_exit_zero=all(p['exit_code']==0 for p in processes.values()))
    write(output/'summary.json', report)


def source_binding():
    files = {p.decode('utf-8') for p in subprocess.check_output(['git','ls-files','-z'], cwd=ROOT).split(b'\0') if p}
    files.update(p.relative_to(ROOT).as_posix() for p in (ROOT/'interim/2027-01-07-a03-explicit-shear-initial-stress').iterdir() if p.suffix in ('.py', '.rs'))
    files.add('tests/test_a03_explicit_shear_initial_stress.py')
    paths = sorted(p for p in files if Path(p).suffix in ('.rs','.py','.toml','.ps1') or Path(p).name == 'Cargo.lock')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def python(output):
    output.mkdir(parents=True, exist_ok=False)
    before = source_binding()
    write(output/'source-binding-before.json', before)
    def run(group):
        return command(output, group, [sys.executable, 'tools/run_ci_python_evidence_shards.py',
                         '--group', group, '--output', str(output/group)])
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        processes = list(pool.map(run, GROUPS))
    after = source_binding()
    reports = [json.loads((output/g/'python-test-inventory.json').read_text()) for g in GROUPS]
    ids = [test for r in reports for test in r.get('test_elapsed_seconds', {})]
    discovered = sorted(set(test for r in reports for test in r.get('discovered_test_ids', [])))
    assigned = [test for r in reports for test in r.get('all_group_test_ids', {}).get(r['group'], [])]
    report = dict(processes=processes, discovered=len(discovered), executed=len(ids),
        executed_ids=sorted(ids), discovered_ids=discovered, assigned_ids=sorted(assigned),
        passed=sum(r.get('tests_run',0)-r.get('failures',0)-r.get('errors',0)-r.get('skipped',0)-r.get('expected_failures',0)-r.get('unexpected_successes',0) for r in reports),
        **{key: sum(r.get(key,0) for r in reports) for key in ('failures','errors','skipped','expected_failures','unexpected_successes')},
        unexecuted=len(set(discovered)-set(ids)), duplicate_ids=len(ids)-len(set(ids)),
        all_discovered_assigned_executed_exactly_once=discovered == sorted(ids) == sorted(assigned),
        source_count=len(before), sources_unchanged=before == after,
        changed_sources=[p for p in before if before.get(p)!=after.get(p)],
        all_processes_exit_zero=all(p['exit_code']==0 for p in processes))
    write(output/'source-binding-after.json', after)
    write(output/'summary.json', report)


def collect(output):
    r = json.loads((output/'rust/summary.json').read_text())
    p = json.loads((output/'python/summary.json').read_text())
    ok = (r['all_processes_exit_zero'] and r['discovered_and_executed_multisets_equal']
        and not any(r[k] for k in ('failed','ignored','filtered_out','unexecuted','unexpected'))
        and p['all_processes_exit_zero'] and p['all_discovered_assigned_executed_exactly_once'] and p['sources_unchanged']
        and not any(p[k] for k in ('failures','errors','skipped','expected_failures','unexpected_successes','unexecuted','duplicate_ids')))
    write(output/'validation.json', dict(passed=ok, software_tests_passed=ok, a03_complete=False,
        diagnostic_only=True, source_commit_base='cdaf850c0347c8d841decce7ddff608d03ae7933',
        rust={k:r[k] for k in ('discovered','executed','passed','failed','ignored','filtered_out','measured','skipped','unexecuted','unexpected')},
        python={k:p[k] for k in ('discovered','executed','passed','failures','errors','skipped','expected_failures','unexpected_successes','unexecuted','duplicate_ids')},
        rust_processes=r['processes'], python_processes=p['processes'],
        source_count=p['source_count'], source_unchanged_during_full_python=p['sources_unchanged'],
        inventories_match=True if ok else False, windows_ci_after_push='未確定',
        acceptance_tolerance_changed=False, these_counts_do_not_override_acceptance=True))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('rust','python','collect'))
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    globals()[a.mode](a.output.resolve())
