"""Freeze an independent sparse initial/construction study before execution."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import freeze_sources


def run(base, output, steps, factor_reuse=False, matrix_free_factor_reuse=False,coercive_factor_reuse=False):
    if sum(map(bool,(factor_reuse,matrix_free_factor_reuse,coercive_factor_reuse)))>1:raise ValueError('choose one experimental reference backend')
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
    if factor_reuse or matrix_free_factor_reuse or coercive_factor_reuse:
        shutil.copyfile(root/'tools/mcc_projected_factor_reference.py',sources/'mcc_projected_factor_reference.py')
    if matrix_free_factor_reuse or coercive_factor_reuse:
        shutil.copyfile(root/'tools/mcc_projected_matrixfree_reference.py',sources/'mcc_projected_matrixfree_reference.py')
    if coercive_factor_reuse:
        shutil.copyfile(root/'tools/mcc_projected_coercive_reference.py',sources/'mcc_projected_coercive_reference.py')
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
        if factor_reuse:command.append('--factor-reuse')
        if matrix_free_factor_reuse:command.append('--matrix-free-factor-reuse')
        if coercive_factor_reuse:command.append('--coercive-factor-reuse')
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
    backend = parser.add_mutually_exclusive_group()
    backend.add_argument('--factor-reuse', action='store_true',
        help='experimental corrected-factor backend; keeps original equations and residual gates')
    backend.add_argument('--matrix-free-factor-reuse', action='store_true',
        help='experimental original-operator backend; original assembled/fresh fallback retained')
    backend.add_argument('--coercive-factor-reuse', action='store_true',
        help='experimental coercive active-set shared factor; fixed-space original equations and gates retained')
    args = parser.parse_args()
    sys.exit(run(args.base.resolve(),args.output.resolve(),args.steps,args.factor_reuse,args.matrix_free_factor_reuse,args.coercive_factor_reuse))
