"""Freeze compatible reference sources, inputs and accepted equilibrium frames."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import freeze_sources


def run(base, output, steps):
    output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    freeze_sources.freeze(output)
    sources, inputs = output/'reference-sources', output/'inputs'
    sources.mkdir()
    inputs.mkdir()
    for name in ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json'):
        shutil.copyfile(base/'inputs'/name, inputs/name)
    for name in ('mcc_compatible_sparse_run.py', 'mcc_compatible_sparse_reference.py',
                 'mcc_rate_fem.py', 'mcc_rate_midpoint.py', 'fem_reference.py', 'verify_high_order.py', 'high_order_io.py'):
        shutil.copyfile(root/'tools'/name, sources/name)
    for name in ('requirements.lock.txt', 'requirements-reference.lock.txt'):
        shutil.copyfile(root/name, sources/name)
    command = [sys.executable, '-X', 'utf8', str(sources/'mcc_compatible_sparse_run.py'),
               *[str(inputs/name) for name in ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json')],
               str(output/'cycle.json'), '--steps', str(steps),
               '--dependency-lock', str(sources/'requirements-reference.lock.txt'),
               '--frame-directory', str(output/'frames')]
    start = time.monotonic()
    status = dict(schema='soil-fem-compatible-sparse-study-progress/1', status='running', reference_steps=steps,
                  goal_complete=False, p6_complete=False, command=command)
    (output/'progress.json').write_text(json.dumps(status, indent=2), encoding='utf-8')
    code = subprocess.call(command)
    status.update(status='completed' if code == 0 else 'failed', exit_code=code, elapsed_seconds=time.monotonic()-start)
    temporary = output/'progress.tmp'
    temporary.write_text(json.dumps(status, indent=2), encoding='utf-8')
    temporary.replace(output/'progress.json')
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--steps', type=int, default=100)
    args = parser.parse_args()
    sys.exit(run(args.base.resolve(), args.output.resolve(), args.steps))
