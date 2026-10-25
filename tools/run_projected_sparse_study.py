"""Freeze an independent sparse initial/construction study before execution."""
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
    inputs, sources = output/'inputs', output/'reference-sources'
    inputs.mkdir()
    sources.mkdir()
    for name in ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json'):
        shutil.copyfile(base/'inputs'/name, inputs/name)
    for name in ('mcc_projected_sparse_run.py', 'mcc_projected_sparse_reference.py',
                 'mcc_projected_volume_diagnostic.py', 'mcc_projected_cycle_diagnostic.py',
                 'mcc_rate_fem.py', 'mcc_rate_midpoint.py', 'fem_reference.py',
                 'verify_high_order.py', 'high_order_io.py'):
        shutil.copyfile(root/'tools'/name, sources/name)
    for name in ('requirements.lock.txt', 'requirements-reference.lock.txt'):
        shutil.copyfile(root/name, sources/name)
    started = time.monotonic()
    progress = dict(schema='soil-fem-sparse-reference-study-progress/1', status='prepared',
                    steps=steps, goal_complete=False, p6_complete=False, phases=[])

    def save():
        progress['elapsed_seconds'] = time.monotonic()-started
        temporary = output/'progress.tmp'
        temporary.write_text(json.dumps(progress, indent=2), encoding='utf-8')
        temporary.replace(output/'progress.json')

    for phase in ('initial', 'cycle'):
        progress['status'] = phase+'_running'
        save()
        command = [sys.executable, '-X', 'utf8', str(sources/'mcc_projected_sparse_run.py'), phase,
                   *[str(inputs/name) for name in ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json')],
                   str(output/(phase+'.json')), '--steps', str(steps),
                   '--dependency-lock', str(sources/'requirements-reference.lock.txt')]
        if phase == 'cycle':
            command += ['--initial', str(output/'initial.json')]
        code = subprocess.call(command)
        progress['phases'].append(dict(phase=phase, exit_code=code))
        if code:
            progress['status'] = 'failed'
            save()
            return code
        progress['status'] = phase+'_saved'
        save()
    progress['status'] = 'completed'
    save()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--steps', type=int, default=100)
    args = parser.parse_args()
    sys.exit(run(args.base.resolve(), args.output.resolve(), args.steps))
