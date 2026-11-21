"""Run only the related software tests; numerical acceptance stays false."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import study

TESTS = ('test_wet_joint_finer_full_history.py', 'test_wet_joint_simultaneous.py',
         'test_cut_volume_history.py', 'test_hydrostatic_cut_volume.py', 'test_nc_midpoint_resume.py')


def preservation():
    baseline = json.loads((study.HERE/'protected-sha256.json').read_text())
    changed = []; missing = []
    for name, expected in baseline['files'].items():
        path = study.ROOT/name
        if not path.is_file():
            missing.append(name)
        elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            changed.append(name)
    return dict(starting_commit=baseline['starting_commit'], parent=baseline['parent'],
                protected_files=len(baseline['files']), unchanged=len(baseline['files'])-len(changed)-len(missing),
                changed=changed, missing=missing, all_existing_files_unchanged=not changed and not missing)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    before = preservation(); rows = []
    source_paths = [study.HERE/n for n in ('study.py', 'reference_history.py', 'compact.py', 'validate.py')]
    source_paths += [study.ROOT/'tests'/n for n in TESTS]
    source_hashes = {str(p.relative_to(study.ROOT)): study.sha(p.read_bytes()) for p in source_paths}
    study.write(output/'tested-source-sha256.json', source_hashes)
    for name in TESTS:
        command = [sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-p', name, '-v']
        result = subprocess.run(command, cwd=study.ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        log = name.replace('.py', '.log')
        (output/log).write_text(result.stdout, encoding='utf-8')
        match = re.search(r'Ran (\d+) tests? in ([0-9.]+)s', result.stdout)
        assert match, ('missing unittest count', name)
        rows.append(dict(log=log, tests=int(match[1]), unittest_reported_seconds=float(match[2]),
                         exit_code=result.returncode, software_check_passed=result.returncode == 0,
                         skips=len(re.findall(r'\.\.\. skipped ', result.stdout))))
    after = preservation()
    same_sources = all(study.sha((study.ROOT/p).read_bytes()) == h for p, h in source_hashes.items())
    report = dict(schema='soil-fem-wet-joint-finer-validation/1', checks=rows,
                  python_tests=sum(r['tests'] for r in rows),
                  software_checks_passed=all(r['software_check_passed'] for r in rows) and same_sources,
                  numerical_acceptance_passed=False, tested_sources_unchanged=same_sources,
                  all_existing_files_unchanged=before['all_existing_files_unchanged'] and after['all_existing_files_unchanged'],
                  rust_build_rerun=False, rust_tests_rerun=False, python_full_suite_rerun=False,
                  windows_ci='unverified', windows_numerical_reproduction='unverified')
    study.write(output/'validation.json', report); study.write(output/'preservation.json', after)
    assert report['software_checks_passed'] and report['all_existing_files_unchanged'], report
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output.resolve()), ensure_ascii=False))
