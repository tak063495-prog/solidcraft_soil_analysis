"""Prepare verified existing seeds and list/run only the existing small case."""
import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = HERE.parent/'2026-11-09-nc-pseudo-arclength'


def commands(out):
    py = sys.executable
    verified = out/'prior-verified'
    prepared = out/'prepared'
    binary = prepared/'target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    items = [
        ('verify-prior', [py, str(PRIOR/'remeasure.py'), '--output', str(verified)]),
        ('prepare-source', [py, str(PRIOR/'prepare_evidence.py'), '--output', str(prepared), '--source-only']),
        ('build-observer', ['cargo', 'build', '--locked', '--offline', '--release', '--manifest-path', str(prepared/'isolated/Cargo.toml'),
                            '--target-dir', str(prepared/'target'), '--example', 'material_diagnostic'])]
    for kind in ['native', 'independent']:
        args = ['--map', kind]
        if kind == 'native':
            args += ['--evidence', str(verified/'evidence/prior'), '--binary', str(binary)]
        else:
            args += ['--reference', str(verified/'evidence/independent-input-rk4-fine'),
                     '--arc', str(verified/'evidence/independent-arc-history400')]
        for label, ds, count, resume in [('coarse', 1e-4, 512, None), ('fine', 5e-5, 1024, None),
                                         ('whole', 1e-4, 1024, None), ('resumed', 1e-4, 512, 'coarse'),
                                         ('fine-tail', 5e-5, 1024, 'fine')]:
            command = [py, str(HERE/'continue_path.py'), *args, '--output', str(out/f'extended-{kind}-{label}'),
                       '--step-size', str(ds), '--count', str(count), '--budget-seconds', '600']
            if resume:
                command += ['--resume', str(out/f'extended-{kind}-{resume}/restart.json')]
            items.append((f'{kind}-{label}', command))
    return items


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    records = []
    for label, command in commands(out):
        started = time.monotonic()
        with (out/(label+'.log')).open('x') as log:
            result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        records.append(dict(label=label, command=command, exit_code=result.returncode, elapsed_seconds=time.monotonic()-started))
        (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
        result.check_returncode()
        report = out/f'extended-{label}'/'report.json'
        if report.exists() and json.loads(report.read_text())['equilibrium_samples'] == 0:
            raise ValueError('diagnostic stopped before any sample; retain report, never infer failure or acceptance')
    # Reuse only immutable original data needed for force reassembly.
    origin = out/'origin'; origin.mkdir()
    shutil.copytree(out/'prior-verified/evidence/prior', origin/'prior')
    shutil.copytree(out/'prior-verified/evidence/independent-input-rk4-fine', origin/'reference')
    subprocess.run([sys.executable, str(HERE/'audit.py'), '--evidence', str(out), '--output', str(out/'analysis')], cwd=ROOT, check=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, required=True); p.add_argument('--list-only', action='store_true')
    a = p.parse_args(); out = a.output.resolve()
    if a.list_only:
        print(json.dumps(commands(out), indent=2))
    else:
        run(out)
