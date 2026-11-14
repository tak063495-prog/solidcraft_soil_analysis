"""Separate global/contact/volume errors of the existing wet JOINT cycle.

A completed diagnostic is not full-history acceptance. No new history error
threshold is inferred from measurements. References consume inputs only.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import copy
import json
from pathlib import Path
import subprocess
import time
import fem_reference as fem
from verify_wet_structures import read, sha, without_time
from verify_wet_joint_cycle import write, restarts, path_evidence
from wet_joint_cycle_reference import reference, discrete_volume_reference, native_fields, differences, compare
from verify_joint_water import audit

ROOT = Path(__file__).resolve().parents[1]
BASE = 'a4f515636292be61069f5b2abc4279c1263d06be'
NAMES = ('mesh', 'ground', 'stages')
INPUTS = ('joint-water-column-tet10-mesh', 'ground-joint-water-mc', 'stages-joint-water-cycle')


def configuration(original, steps, maximum, rejected=False):
    case = copy.deepcopy(original)
    case['settings']['initial_increment'] = 1/steps
    case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'], 1/(2*steps))
    selected = case['verification_reject_once_stages']
    for stage in case['stages']:
        stage['initial_increment'] = (2 if rejected and stage['name'] in selected else 1)/steps
    if not rejected:
        case['verification_reject_once_stages'] = []
    case['joint_system']['maximum_increment_m'] = maximum
    # Keep the original tangent width; all study maxima are >= 5e-6 m.
    return case


def schedule(run, steps):
    assert run['complete'] and len(run['frames']) == 7, 'incomplete history'
    phases = [('geostatic', run['geostatic']['attempts'])]+[
        (s['name'], [a for a in run['attempts'] if a['stage'] == s['name']]) for s in run['case']['stages']]
    rows = []
    for name, attempts in phases:
        accepted = [a for a in attempts if a['accepted']]
        assert [a['progress_to'] for a in accepted] == [i/steps for i in range(1, steps+1)], ('global schedule changed', name)
        assert all(a['relative_force_residual'] < run['case']['settings']['relative_force_tolerance'] for a in accepted)
        solves = [s for a in accepted for s in a['linear_solves']]
        assert all(s['true_relative_residual'] <= run['case']['settings']['relative_linear_tolerance'] for s in solves)
        assert all(a['committed_state_unchanged_on_rejection'] for a in attempts if not a['accepted'])
        rows.append(dict(stage=name, accepted=len(accepted), rejected=len(attempts)-len(accepted),
                         maximum_relative_force_residual=max(a['relative_force_residual'] for a in accepted),
                         maximum_true_relative_linear_residual=max((s['true_relative_residual'] for s in solves), default=0.)))
    return rows


def all_fields(a, b):
    keys = a['frames'][0].keys()
    return differences(a, dict(frames=[{k: v for k, v in f.items() if k in keys} for f in b['frames']]))


def measure(output, steps, restart_rows):
    n = steps; paths = {}; refs = {}
    for g, local, volume in [(n, 1, 'discrete'), (n, 8, 'discrete'), (n, 32, 'discrete'),
                              (2*n, 8, 'discrete'), (4*n, 8, 'discrete'),
                              (4*n, 8, 'gl5'), (4*n, 8, 'gl7')]:
        key = f'g{g}-l{local}-{volume}'; refs[key] = read(output/(key+'.json'))
        assert len(refs[key]['frames']) == 7 and len(refs[key]['path']) == 7*g
        assert not refs[key]['reference_uses_native_results']
    native = {}
    for g, local in [(n, 5e-4), (n, 5e-5), (n, 5e-6), (2*n, 5e-6), (4*n, 5e-6)]:
        key = f'g{g}-m{local:g}'; r = read(output/key/'stage-results.json')
        paths[key] = schedule(r, g); native[key] = r
    clean = native[f'g{n}-m5e-06']; rejected = read(output/'rejected/stage-results.json')
    rejection_schedule = schedule(rejected, n)
    assert without_time(clean['frames']) == without_time(rejected['frames']), 'rejected trial changed frames'
    assert without_time(clean['last_committed']) == without_time(rejected['last_committed']), 'rejected trial changed state'
    injected = [a for a in rejected['attempts'] if not a['accepted']]
    assert len(injected) == 3 and all(a['reason'] == 'verification_injected_rejection_after_convergence' for a in injected)
    assert any(a.get('rejected_trial_joint_slip_changed') for a in injected)
    fine = refs[f'g{4*n}-l8-discrete']; continuum = refs[f'g{4*n}-l8-gl7']
    comparisons = {}; old_gates = {}
    for key, r in native.items():
        comparisons[key] = dict(discrete_volume=all_fields(native_fields(r), fine),
                                independent_gl7_volume=all_fields(native_fields(r), continuum))
        for label, expected in [('discrete', fine), ('gl7', continuum)]:
            try:
                compare(r, expected); old_gates[key+'/'+label] = True
            except AssertionError:
                old_gates[key+'/'+label] = False
    return dict(schema='soil-fem-wet-joint-separated-convergence/1', starting_commit=BASE, passed=False,
                diagnostic_completed=True, reference_uses_native_results=False,
                solver_sha256=clean['provenance']['solver_binary_sha256'], global_steps=[n, 2*n, 4*n],
                native_schedules=paths, rejection_schedule=rejection_schedule,
                actual_open_slip=path_evidence(rejected),
                stored_field_audits={key: audit(read(output/'mesh.json'), fem.normalize_ground(read(output/'ground.json')), r)
                                     for key, r in native.items()},
                independent_maximum_force_residual={key: r['maximum_relative_force_residual'] for key, r in refs.items()},
                rejected_history_exact_same_accepted_path=True, restart_checks=restart_rows,
                native_vs_reference=comparisons, existing_displacement_traction_gates=old_gates,
                native_global=[all_fields(native_fields(native[f'g{a}-m5e-06']), native_fields(native[f'g{b}-m5e-06']))
                               for a, b in [(n, 2*n), (2*n, 4*n)]],
                native_local=[all_fields(native_fields(native[f'g{n}-m{a:g}']), native_fields(native[f'g{n}-m{b:g}']))
                              for a, b in [(5e-4, 5e-5), (5e-5, 5e-6)]],
                independent_global=[differences(refs[f'g{a}-l8-discrete'], refs[f'g{b}-l8-discrete'])
                                    for a, b in [(n, 2*n), (2*n, 4*n)]],
                independent_local=[differences(refs[f'g{n}-l{a}-discrete'], refs[f'g{n}-l{b}-discrete'])
                                   for a, b in [(1, 8), (8, 32)]],
                volume_discretization=differences(fine, continuum),
                independent_volume_order=differences(refs[f'g{4*n}-l8-gl5'], continuum),
                full_history_acceptance='unresolved: global path and volume errors remain measured; no existing full-history error criterion certifies them',
                numerical_work='Endpoint work and numerical dissipation use different local partitions. Compare absolute magnitudes and the energy identity; relative error of a vanishing numerical quantity is not a continuum convergence certificate.',
                windows_ci='unverified', p3b_p4_complete=False, p6_complete=False)


def run(solver, output, steps=64, restart=True):
    assert steps >= 4 and steps & (steps-1) == 0
    output.mkdir(parents=True, exist_ok=False)
    mesh, ground, original = [read(ROOT/'examples'/(name+'.json')) for name in INPUTS]
    for name, value in zip(NAMES, (mesh, ground, original)):
        write(output/(name+'.json'), value)
    g = fem.normalize_ground(copy.deepcopy(ground)); c = fem.normalize_stage_case(copy.deepcopy(original))
    # Every independent solve precedes reading any native result.
    for count, local, volume in [(steps, 1, 'discrete'), (steps, 8, 'discrete'), (steps, 32, 'discrete'),
                                 (2*steps, 8, 'discrete'), (4*steps, 8, 'discrete'),
                                 (4*steps, 8, 'gl5'), (4*steps, 8, 'gl7')]:
        expected = (discrete_volume_reference(mesh, g, c, count, local) if volume == 'discrete'
                    else reference(mesh, g, c, count, local, int(volume[2:])))
        write(output/f'g{count}-l{local}-{volume}.json', expected)
        print('reference', count, local, volume, flush=True)
    timings = []
    for count, maximum, reject in [(steps, 5e-4, False), (steps, 5e-5, False), (steps, 5e-6, False),
                                   (2*steps, 5e-6, False), (4*steps, 5e-6, False), (steps, 5e-6, True)]:
        key = 'rejected' if reject else f'g{count}-m{maximum:g}'
        folder = output/key; folder.mkdir()
        for name, value in zip(NAMES, (mesh, ground, configuration(original, count, maximum, reject))):
            write(folder/(name+'.json'), value)
        start = time.monotonic()
        with (folder/'solve.log').open('x') as log:
            result = subprocess.run([str(solver), 'solve-stages', *[str(folder/(name+'.json')) for name in NAMES],
                                     str(folder/'stage-results.json')], stdout=log, stderr=subprocess.STDOUT, timeout=300)
        assert result.returncode == 0, ('solver failed', key)
        schedule(read(folder/'stage-results.json'), count)
        timings.append(dict(route=key, elapsed_seconds=time.monotonic()-start, exit_code=result.returncode))
        print('native', key, flush=True)
    key = f'g{4*steps}-m5e-06'; folder = output/key
    restart_rows = restarts(solver, folder, read(folder/'stage-results.json')) if restart else []
    write(output/'restarts.json', restart_rows); write(output/'timings.json', timings)
    report = measure(output, steps, restart_rows); write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=64)
    args = parser.parse_args()
    result = run(args.solver.resolve(), args.output.resolve(), args.steps)
    print(json.dumps(dict(passed=result['passed'], diagnostic_completed=result['diagnostic_completed'])))
