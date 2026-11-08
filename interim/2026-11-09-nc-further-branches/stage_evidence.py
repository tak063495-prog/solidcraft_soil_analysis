"""Stage auditable evidence, retaining full fields and actual solved systems.

Material process IO is hashed in full but excluded from the distribution, as in
prior diagnostics. Raw IO stays in the execution workspace. No solver is run.
"""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        for p in sorted(source.rglob('*')):
            if p.is_file():
                copy(p, target/p.relative_to(source))
    elif source.is_file():
        # Copies avoid coupling staged evidence to subsequent workspace changes.
        shutil.copyfile(source, target)


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def stage(work, out):
    out.mkdir(parents=True, exist_ok=False)
    old = work/'beyond-verified/evidence'
    copy(old/'origin', out/'origin')
    copy(work/'prior-seed/independent-arc-history400', out/'original-independent-arc')
    for kind in ['native', 'independent']:
        for label in ['whole', 'fine-tail']:
            copy(old/f'extended-{kind}-{label}/restart.json', out/f'parents/{kind}-{label}.json')
    processes = []
    for p in sorted(work.iterdir()):
        if not p.is_dir() or not p.name.startswith(('extended-', 'roots-')):
            continue
        if p.name in ['roots-native-anchor1024', 'roots-independent-anchor1024']:
            dest = out/'failed-executions'/p.name
        else:
            dest = out/p.name
            report = json.loads((p/'report.json').read_text())
            # Each successful CLI prints exactly its final report; traceback/error logs fail this check.
            labels = [p.name.removeprefix('extended-'), p.name]
            if p.name.endswith('-reverse'):
                labels = [p.name.removeprefix('extended-')+'-fixed']
            log = next((work/(name+'.log') for name in labels if (work/(name+'.log')).exists()), None)
            if log is None or json.loads(log.read_text()) != report:
                raise ValueError('successful process log/report mismatch: '+p.name)
            processes.append(dict(folder=p.name, completed_cli_report_exact=True, exit_code=0,
                                  elapsed_seconds=report['elapsed_seconds'], passed=report['passed'], stop_reason=report['stop_reason']))
        dest.mkdir(parents=True)
        for item in sorted(p.iterdir()):
            if item.is_file():
                copy(item, dest/item.name)
        if p.name.startswith('roots-'):
            kind = 'native' if 'native' in p.name else 'independent'
            seed = old/f'extended-{kind}-whole/restart.json' if 'anchor1024' in p.name else work/f'extended-{kind}-long/restart.json'
            controls = json.loads((p/'controls.json').read_text())
            if sha(seed) != controls['binding']['seed_sha256']:
                raise ValueError('root seed changed')
            copy(seed, dest/'seed.json')
        if (p/'material').exists():
            files = {x.relative_to(p/'material').as_posix():dict(sha256=sha(x), bytes=x.stat().st_size) for x in sorted((p/'material').rglob('*')) if x.is_file()}
            write(dest/'native-material-trial-file-sha256.json', dict(files=files, files_bundled=False,
                files_retained_in_workspace=True, note='All raw request/result/log bytes hashed; not distributed'))
    for p in sorted(work.glob('*.log')):
        copy(p, out/'logs'/p.name)
    for name in ['mid-checkpoint-captures.json', 'native-actual-mid-checkpoint.json', 'independent-actual-mid-checkpoint.json']:
        copy(work/name, out/name)
    copy(work/'failed-sources', out/'failed-executions/source')
    copy(work/'prepared/isolated/src', out/'observer-source/src')
    copy(work/'prepared/isolated/examples', out/'observer-source/examples')
    for name in ['Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml']:
        copy(work/'prepared/isolated'/name, out/'observer-source'/name)
    copy(work/'prepared/source-sha256.json', out/'observer-source-sha256.json')
    binary = work/'prepared/target/release/examples/material_diagnostic'
    write(out/'binary-sha256.json', dict(material_diagnostic=sha(binary), binary_bundled=False,
        exact_prior_binary=sha(binary) == json.loads((old/'binary-sha256.json').read_text())['material_diagnostic']))
    copy(work/'prior-seed/verification.json', out/'prior-seed-verification.json')
    copy(work/'beyond-verified/remeasurement.json', out/'beyond-remeasurement.json')
    write(out/'numerical-processes.json', processes)
    write(out/'execution-issues.json', dict(excluded_from_successful_processes=[
        dict(log='logs/native-reverse.log', exit_code=2, numerical_samples=0, reason='missing --output parser definition'),
        dict(log='logs/independent-reverse.log', exit_code=2, numerical_samples=0, reason='missing --output parser definition'),
        dict(folder='failed-executions/roots-native-anchor1024', exit_code=1, reason='numpy int64 in final JSON report; raw trials and old driver retained'),
        dict(folder='failed-executions/roots-independent-anchor1024', exit_code=1, reason='numpy int64 in final JSON report; raw trials and old driver retained')],
        numerical_criteria_changed=False))
    source = {p.name:sha(p) for p in HERE.glob('*.py')}
    for name in source:
        copy(HERE/name, out/'diagnostic-source'/name)
    write(out/'diagnostic-source-sha256.json', source)
    return dict(processes=len(processes), passed=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(stage(a.work.resolve(), a.output.resolve()), indent=2))
