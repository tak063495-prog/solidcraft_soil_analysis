"""Extend only the existing seven-stage small wet JOINT series; never certify it."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import audit_wet_joint_bulk_history as bulk
from diagnose_cut_volume_history import validate_path, summary
from wet_joint_cycle_reference import serializable, differences

BASE = 'c71777d2dfd5cde12774981bdd357cf46c69fa06'
ROUTES = (16, 32, 64, 128)


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(serializable(value), stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def measure(folder):
    values = bulk.history.inputs()
    assert tuple(read(folder/(n+'.json')) for n in ('mesh', 'ground', 'stages')) == values, 'changed input'
    paths = {}; results = {}
    for steps in ROUTES:
        route = folder/f'g{steps}-l8'
        result = read(route/'result.json')
        results[steps] = result
        validate_path(result, *values, steps, 8)
        report = bulk.audit(route/'accepted', result, values)
        paths[route.name] = report
    comparisons = {}
    for coarse, fine in zip(ROUTES, ROUTES[1:]):
        comparisons[f'{coarse}__{fine}'] = bulk.history.compare_accepted(
            folder/f'g{coarse}-l8/accepted', folder/f'g{fine}-l8/accepted', ratio=2)
    # Same 112 coordinates for all pairs; do not imply unmatched fine points
    # are compared, or compare maxima on different sampling grids as rates.
    fixed_grid = {}
    for steps in ROUTES[1:]:
        fixed_grid[f'16__{steps}'] = bulk.history.compare_accepted(
            folder/'g16-l8/accepted', folder/f'g{steps}-l8/accepted', ratio=steps//16)
    endpoints = {}
    for coarse, fine in zip(ROUTES, ROUTES[1:]):
        rows = differences(results[coarse], results[fine])
        endpoints[f'{coarse}__{fine}'] = dict(fields=rows, maxima=summary(rows))
    return dict(schema='soil-fem-wet-joint-small-refinement/1', starting_commit=BASE,
                passed=False, diagnostic_completed=True, input_only_reference=True,
                native_rerun=False, independent_bulk_plastic_solution=False,
                full_history_convergence_accepted=False, completion_restart_verified=False,
                p3b_p4_complete=False, structure_water_joint_complete=False, p6_complete=False,
                global_steps=list(ROUTES), contact_steps=8, volume_order=7,
                paths=paths, adjacent_comparisons=comparisons, fixed_16_grid=fixed_grid,
                seven_endpoint_comparisons=endpoints,
                limitations=[
                    'Existing affine-cut hydrostatic GL7 reference, MC elastic bulk branch and Coulomb JOINT only.',
                    'No new all-history acceptance tolerance, extrapolated solution or inferred convergence order.',
                    'Fine unmatched points are audited internally but are not compared across routes.',
                    'No new native comparison, bulk plastic 3D completion or actual solver checkpoint restart.'])


def run(folder):
    folder.mkdir(parents=True, exist_ok=False)
    values = bulk.history.inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        write(folder/(name+'.json'), value)
    timings = []
    for steps in ROUTES:
        route = folder/f'g{steps}-l8'; route.mkdir()
        sink = bulk.history.HistoryWriter(route/'accepted', chunk_size=16)
        started = time.monotonic()
        result = bulk.history.cut_volume_reference(*values, steps, 8, 7, history_sink=sink)
        sink.close()
        write(route/'result.json', result)
        timings.append(dict(route=route.name, elapsed_seconds=time.monotonic()-started))
    write(folder/'timings.json', timings)
    report = measure(folder)
    write(folder/'measurements.json', report)
    return report


def replay(folder):
    expected = read(folder/'measurements.json')
    observed = measure(folder)
    if expected != observed:
        raise ValueError('remeasured fields or acceptance flags differ')
    return dict(all_measurements_exact=True, numerical_acceptance=False,
                accepted_states=sum(p['accepted_states'] for p in observed['paths'].values()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--replay', type=Path)
    args = p.parse_args()
    if args.replay:
        write(args.output, replay(args.replay.resolve()))
    else:
        run(args.output.resolve())
