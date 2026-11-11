"""Reproduce the small specimen; larger meshes/longer paths need --user-pc."""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import json
from pathlib import Path
import subprocess
import shutil
import time
import traceback
import numpy as np
from .geometry import refine
from .solver import Model, run
from .compare import read, write, digest, save, pair, load_series, spatial_series, diagnose, common, metric
from .native import configuration, convert

ROOT = Path(__file__).resolve().parents[3]
BASE = '5f564621c7396c2a605edbfe5451cd800cf318d3'
INPUTS = ('joint-water-column-tet10-mesh.json', 'ground-joint-water-mc.json', 'stages-joint-water-cycle.json')


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def differences_by_stage(a, b):
    keys, (ia, ib) = common([a, b]); rows = []
    for name in dict.fromkeys(k[0] for k in keys):
        mask = [i for i, key in enumerate(keys) if key[0] == name]
        rows.append(dict(stage=name, plastic_slip=metric(a['fields']['plastic_slip'][ia[mask]], b['fields']['plastic_slip'][ib[mask]])))
    return rows


def execute(args):
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    args._output_created = True
    write(output/'request.json', dict(starting_commit=BASE, steps=args.steps, levels=args.levels,
                                    user_pc=args.user_pc, numerical_acceptance=False))
    if (not args.user_pc and (max(args.steps) > 128 or max(args.levels) > 1)):
        raise ValueError('over cloud budget: use the documented Windows --user-pc commands')
    if args.levels != list(range(len(args.levels))) or len(set(args.steps)) != len(args.steps) or args.steps != sorted(args.steps):
        raise ValueError('levels must be 0..n; steps must be strictly increasing')
    if min(args.steps) < 2:
        raise ValueError('at least two global increments required')
    if args.native_solver and args.native_steps and not set(args.native_steps) <= set(args.steps):
        raise ValueError('--native-steps must be a subset of --steps')
    mesh, ground, case = [read(ROOT/'examples'/name) for name in INPUTS]
    for name, value in zip(('mesh-h0', 'ground', 'stages'), (mesh, ground, case)):
        write(output/(name+'.json'), value)
    solver_snapshot = None
    if args.native_solver:
        solver_snapshot = output/('native-probe.exe' if args.native_solver.suffix.lower() == '.exe' else 'native-probe')
        shutil.copy2(args.native_solver, solver_snapshot)
    input_hashes = {name: digest(ROOT/'examples'/name) for name in INPUTS}
    sources = {str(path.relative_to(ROOT)): digest(path) for path in sorted(Path(__file__).parent.glob('*.py'))}
    write(output/'provenance.json', dict(starting_commit=BASE, input_sha256=input_hashes, reference_source_sha256=sources,
                                       native_solver_sha256=digest(solver_snapshot) if solver_snapshot else None))
    triangles = [np.array(mesh['nodes'])[face].tolist() for face in mesh['face_sets'][case['joint_system']['joints'][0]['side_a']]]
    write(output/'parents.json', triangles)
    results = {}; native = {}; diagnostics = {}; timings = {}; models = {}
    for level in args.levels:
        if level:
            mesh = refine(mesh); write(output/f'mesh-h{level}.json', mesh)
        model = Model(mesh, ground, case, pc=args.user_pc); models[level] = model
        for n in args.steps:
            name = f'h{level}-g{n}-l8'; started = time.monotonic()
            # Preserve every committed point even if a later iteration fails.
            with (output/(name+'-partial.jsonl')).open('x', encoding='utf-8') as stream:
                def sink(row, fields):
                    stream.write(json.dumps(dict(row=row, fields=fields), default=json_value, allow_nan=False)+'\n'); stream.flush()
                value = run(model, n, 8, sink=sink)
            save(output/name, value); results[level, n] = value
            diagnostics[name] = diagnose(value, model); timings[name] = time.monotonic()-started
            print('reference complete', name, flush=True)
            if args.native_solver and n in args.native_steps:
                native_name = f'h{level}-g{n}-native'; cfg = configuration(case, n)
                write(output/(native_name+'-stages.json'), cfg)
                folder = output/(native_name+'-raw')
                command = [str(solver_snapshot), str(output/f'mesh-h{level}.json'), str(output/'ground.json'), str(output/(native_name+'-stages.json')), str(folder)]
                if args.user_pc:
                    command.append('--user-pc')
                write(output/(native_name+'-command.json'), dict(argv=command, timeout_seconds=args.native_timeout_seconds))
                started = time.monotonic()
                with (output/(native_name+'.log')).open('x') as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=args.native_timeout_seconds, check=True)
                value = convert(folder, model, n); save(output/native_name, value); native[level, n] = value
                diagnostics[native_name] = diagnose(value, model); timings[native_name] = time.monotonic()-started
                print('native complete', native_name, flush=True)
    coarse_model = models[0]; diagnostic_n = min(64, max(args.steps))
    if diagnostic_n not in args.steps:
        diagnostic_n = args.steps[-1]
    local_values = [results[0, diagnostic_n]]
    for local in (16, 32):
        name = f'h0-g{diagnostic_n}-l{local}'; started = time.monotonic()
        value = run(coarse_model, diagnostic_n, local); save(output/name, value)
        local_values.append(value); timings[name] = time.monotonic()-started
    tighter_n = args.steps[min(1, len(args.steps)-1)]
    tight = run(coarse_model, tighter_n, 8, tolerance=1e-11); save(output/f'h0-g{tighter_n}-tight', tight)
    from .compare import observed
    endpoint_comparisons = {}
    for (h, n), value in native.items():
        selected = (np.arange(7)+1)*n-1; ref = results[h, n]
        rows = [{key: metric(value['fields'][key][i], ref['fields'][key][i]) for key in ('u', 'traction')} for i in selected]
        endpoint_comparisons[f'h{h}-g{n}'] = dict(frames=rows, existing_displacement_gate=.005, existing_traction_gate=.01,
            existing_endpoint_gates_passed=all(row['u']['relative_l2'] < .005 and row['traction']['relative_l2'] < .01 for row in rows),
            full_history_acceptance=False)
    keys, ids = common(local_values)
    slips = [v['fields']['plastic_slip'][ii] for v, ii in zip(local_values, ids)]
    report = dict(schema='py-wet-joint-reference/1', starting_commit=BASE,
                  diagnostic_completed=True, passed=False, numerical_acceptance=False,
                  load_series={f'h{level}': load_series([results[level, n] for n in args.steps]) for level in args.levels},
                  spatial_series={str(n): spatial_series([results[level, n] for level in args.levels], triangles) for n in args.steps} if len(args.levels) > 1 else {},
                  native_to_reference={f'h{h}-g{n}': pair(value, results[h, n]) for (h, n), value in native.items()},
                  existing_endpoint_comparisons=endpoint_comparisons,
                  native_load_series={f'h{h}': load_series([native[h, n] for n in args.native_steps]) for h in args.levels} if len(args.native_steps) >= 2 and args.native_solver else {},
                  native_spatial_series={str(n): spatial_series([native[h, n] for h in args.levels], triangles) for n in args.native_steps} if len(args.levels) > 1 and args.native_solver else {},
                  contact_partition_pairs=[pair(local_values[0], local_values[1]), pair(local_values[1], local_values[2])],
                  contact_partition_order=observed(slips[0]-slips[1], slips[1]-slips[2]),
                  tightened_force_comparison=pair(results[0, tighter_n], tight),
                  tightened_force_steps=tighter_n, stage_differences={}, timings_seconds=timings)
    for n in args.steps:
        if 2*n in args.steps:
            report['stage_differences'][f'{n}/{2*n}'] = differences_by_stage(results[0, n], results[0, 2*n])
    report['maximum_reference_reassembled_residual'] = max(v['maximum_reassembled_force_residual'] for k, v in diagnostics.items() if 'native' not in k)
    report['maximum_native_reassembled_residual'] = max((v['maximum_reassembled_force_residual'] for k, v in diagnostics.items() if 'native' in k), default=None)
    write(output/'diagnostics.json', diagnostics); write(output/'measurements.json', report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', nargs='+', type=int, default=[16, 32, 64, 128]); p.add_argument('--levels', nargs='+', type=int, default=[0, 1])
    p.add_argument('--native-solver', type=Path); p.add_argument('--native-steps', nargs='+', type=int, default=[16, 32, 64])
    p.add_argument('--native-timeout-seconds', type=float, default=900); p.add_argument('--user-pc', action='store_true')
    args = p.parse_args()
    try:
        execute(args)
    except BaseException as error:
        if getattr(args, '_output_created', False) and args.output.exists() and not (args.output/'failure.json').exists():
            write(args.output/'failure.json', dict(completed=False, passed=False, error=repr(error), traceback=traceback.format_exc(),
                                                  reason='timeout' if isinstance(error, subprocess.TimeoutExpired) else 'cancelled' if isinstance(error, KeyboardInterrupt) else 'error'))
        raise


if __name__ == '__main__':
    main()
