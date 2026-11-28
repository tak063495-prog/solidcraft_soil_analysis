"""Stage immutable completed runs; accepted histories/systems are fully bundled.

Raw native trial requests/results are retained in work and hashed, not bundled.
Completed large JSONL files are hard-linked to avoid transient disk duplication;
no staged source is edited after this point.
"""
import argparse
import hashlib
import json
from pathlib import Path
import os
import shutil

HERE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def copy(source, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, dest)


def stage(seed, run, binary, out):
    required = [f'extended-{kind}-{label}' for kind in ['native', 'independent'] for label in ['next', 'half-a', 'half-b', 'restart128']] + ['zero-native', 'zero-independent']
    for name in required:
        if not (run/name/'report.json').exists():
            raise ValueError('run still incomplete; do not stage: '+name)
    out.mkdir(parents=True, exist_ok=False)
    origin = seed/'2026-11-09-nc-beyond-fold/origin'
    shutil.copytree(origin, out/'origin')
    latest = seed/'2026-11-09-nc-further-branches'
    shutil.copytree(latest/'observer-source', out/'observer-source')
    shutil.copytree(seed/'2026-11-09-nc-pseudo-arclength/independent-arc-history400', out/'original-independent-arc')
    for folder in seed.iterdir():
        if (folder/'verification.json').exists():
            copy(folder/'verification.json', out/'archive-verifications'/(folder.name+'.json'))
    for kind in ['native', 'independent']:
        copy(latest/f'extended-{kind}-further/restart.json', out/f'parents/{kind}.json')
        copy(run/f'captured-{kind}.json', out/f'captured-{kind}.json')
    copy(run/'checkpoint-captures.json', out/'checkpoint-captures.json')
    copy(run/'restart-process-exits.json', out/'restart-process-exits.json')
    processes = []
    for folder in sorted(run.iterdir()):
        if not folder.is_dir() or not folder.name.startswith(('extended-', 'zero-')): continue
        report = json.loads((folder/'report.json').read_text())
        logname = 'nc-'+folder.name.removeprefix('extended-')+'.log'
        log = run.parent/logname
        if json.loads(log.read_text()) != report:
            raise ValueError('process log/report mismatch: '+folder.name)
        dest = out/folder.name
        dest.mkdir()
        for p in sorted(folder.iterdir()):
            if p.is_file(): os.link(p, dest/p.name)
        if (folder/'material').exists():
            files = {p.relative_to(folder/'material').as_posix():dict(bytes=p.stat().st_size, sha256=sha(p))
                     for p in sorted((folder/'material').rglob('*')) if p.is_file()}
            write(dest/'native-material-trial-file-sha256.json', dict(files=files, files_bundled=False,
                files_retained_in_execution_workspace=True))
        copy(log, out/'logs'/log.name)
        processes.append(dict(folder=folder.name, completed_report_log_exact=True, exit_code=0,
                              passed=report['passed'], stop_reason=report['stop_reason']))
    write(out/'numerical-processes.json', processes)
    for kind in ['native', 'independent']:
        prior = run/f'precheck-{kind}-restart128'
        if prior.exists():
            write(out/f'precheck-{kind}-restart.json', dict(report=json.loads((prior/'report.json').read_text()),
                  retained_in_workspace=True, selected_final_numerical_process=False,
                  note='Earlier duplicate 128-row restart; selected final run has explicit subprocess exit receipt'))
    write(out/'binary-sha256.json', dict(material_diagnostic=sha(binary), binary_bundled=False,
          original_native_checkpoint_binding_restored_without_changes=True))
    source = {p.name:sha(p) for p in HERE.glob('*.py')}
    for name in source: copy(HERE/name, out/'diagnostic-source'/name)
    write(out/'diagnostic-source-sha256.json', source)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['seed', 'run', 'binary', 'output']: p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    stage(a.seed.resolve(), a.run.resolve(), a.binary.resolve(), a.output.resolve())
