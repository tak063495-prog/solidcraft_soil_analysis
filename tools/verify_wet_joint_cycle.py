"""Reproducible diagnostics for the existing small wet opening/sliding cycle.

Always report passed=false: full-history convergence/acceptance is unresolved.
Existing u/traction gates are reported separately and cannot certify histories.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import copy
import json
from pathlib import Path
from export_stages import export
import fem_reference as fem
from verify_checkpoint import invoke
from verify_joint_water import verify_package
from verify_wet_structures import read, sha, without_time
from wet_joint_cycle_reference import reference, serializable, differences, native_fields, compare

ROOT = Path(__file__).resolve().parents[1]


def write(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(serializable(value), stream, indent=2, allow_nan=False)


def path_evidence(run):
    assert run['complete'] and len(run['frames']) == 7
    rows = []
    for f in run['frames']:
        points = f['state']['joints']['points']; result = f['joint_results']['points']
        rows.append(dict(frame=f['name'], points=len(points),
                         status_counts={s: sum(p['status'] == s for p in points) for s in ('open', 'slip', 'stick')},
                         minimum_gap_m=min(p['gap_m'] for p in result),
                         maximum_gap_m=max(p['gap_m'] for p in result),
                         maximum_accumulated_plastic_slip_m=max(p['accumulated_plastic_slip_m'] for p in points)))
    opened = next(f for f in run['frames'] if f['name'] == 'open')
    assert all(p['status'] == 'open' and p['traction_kpa'] == [0., 0., 0.]
               for p in opened['state']['joints']['points'])
    assert all(p['opening_m'] > 1e-6 and p['hydraulic']['physical_pore_pressure_kpa'] == 2.5
               for p in opened['joint_results']['points'])
    sheared = next(f for f in run['frames'] if f['name'] == 'shear')
    assert all(p['status'] == 'slip' and p['accumulated_plastic_slip_m'] > 0
               for p in sheared['state']['joints']['points'])
    rejects = [a for a in run['geostatic']['attempts']+run['attempts'] if not a['accepted']]
    assert rejects and all(a['committed_state_unchanged_on_rejection'] for a in rejects)
    assert any(a.get('rejected_trial_joint_slip_changed') for a in rejects)
    return dict(frames=rows, rejected_attempts=len(rejects), all_rejections_preserve_committed_state=True)


def restarts(solver, folder, run):
    """One strictly interior accepted point of every phase; no resealing."""
    selected = []; total = 0
    for name, attempts in [('geostatic', run['geostatic']['attempts'])]+[
            (s['name'], [a for a in run['attempts'] if a['stage'] == s['name']]) for s in run['case']['stages']]:
        accepted = [a for a in attempts if a['accepted']]
        interior = next((i for i, a in enumerate(accepted) if 0 < a['progress_to'] < 1), None)
        assert interior is not None, ('no interior checkpoint', name)
        selected.append((name, total+interior+1, accepted[interior]['progress_to']))
        total += len(accepted)
    rows = []
    for name, count, progress in selected:
        cp = folder/f'checkpoint-{name}.json'; out = folder/f'resumed-{name}.json'
        inputs = [str(folder/(n+'.json')) for n in ('mesh', 'ground', 'stages')]
        invoke([str(solver), 'checkpoint-stages', *inputs, str(count), str(cp)], folder/f'checkpoint-{name}.log')
        before = sha(cp)
        invoke([str(solver), 'resume-stages', *inputs, str(cp), str(out)], folder/f'resume-{name}.log')
        assert before == sha(cp), 'checkpoint changed'
        assert without_time(read(out)) == without_time(run), ('whole Run restart mismatch', name)
        rows.append(dict(stage=name, progress=progress, count=count, checkpoint_sha256=before,
                         whole_run_exact_except_elapsed_seconds=True))
    return rows


def run(solver, output, counts=(32, 64), restart=True):
    assert len(counts) == 2 and counts[1] == 2*counts[0] and counts[0] >= 2
    solver = Path(solver).resolve(); output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    mesh = read(ROOT/'examples/joint-water-column-tet10-mesh.json')
    ground = read(ROOT/'examples/ground-joint-water-mc.json')
    original = read(ROOT/'examples/stages-joint-water-cycle.json')
    # Solve references before loading any native result. Only these three
    # inputs and the explicit independent numerical controls are arguments.
    g = fem.normalize_ground(copy.deepcopy(ground))
    c = fem.normalize_stage_case(copy.deepcopy(original))
    refs = []
    for count in counts:
        r = reference(mesh, g, c, count, count, 5); refs.append(r)
        write(output/f'reference-{count}.json', r)
    quadrature = reference(mesh, g, c, counts[1], counts[1], 7)
    write(output/'reference-order7.json', quadrature)
    native = []; checks = {}; routes = {}; audits = {}; comparisons = {}
    for name, count in [('original', None)]+[(f'h{n}', n) for n in counts]:
        folder = output/name; folder.mkdir(); case = copy.deepcopy(original)
        if count is not None:
            case['settings']['initial_increment'] = 1/count
            for s in case['stages']:
                s['initial_increment'] = 1/count
            case['joint_system']['maximum_increment_m'] /= 10
        for key, value in [('mesh', mesh), ('ground', ground), ('stages', case)]:
            write(folder/(key+'.json'), value)
        invoke([str(solver), 'solve-stages', *[str(folder/(n+'.json')) for n in ('mesh', 'ground', 'stages')],
                str(folder/'stage-results.json')], folder/'solve.log')
        n = read(folder/'stage-results.json'); native.append(n); routes[name] = path_evidence(n)
        export(folder/'mesh.json', folder/'stage-results.json', folder/'results')
        audits[name] = verify_package(folder/'mesh.json', folder/'ground.json', folder/'stage-results.json', folder/'results', solver)
        # A failure of existing gates remains evidence, not a raised success.
        try:
            comparisons[name] = compare(n, refs[-1]); checks[name] = True
        except AssertionError as error:
            checks[name] = False
            write(folder/'comparison-rejection.json', dict(reason=str(error)))
            expected = copy.deepcopy(refs[-1]); keys = native_fields(n)['frames'][0].keys()
            for f in expected['frames']:
                for key in list(f.keys()-keys):
                    del f[key]
            comparisons[name] = differences(native_fields(n), expected)
    restart_rows = restarts(solver, output/'original', native[0]) if restart else []
    report = dict(schema='soil-fem-wet-joint-cycle-diagnostics/1', passed=False,
                  starting_commit='b04c80e3e8f6608b82f5c4f8cc936abde62fdc07',
                  independent_reference_completed=True, reference_uses_native_results=False,
                  solver_sha256=sha(solver), input_sha256={str(p.relative_to(ROOT)): sha(p) for p in (
                      ROOT/'examples/joint-water-column-tet10-mesh.json', ROOT/'examples/ground-joint-water-mc.json',
                      ROOT/'examples/stages-joint-water-cycle.json')},
                  reference_global_and_local_subincrements=list(counts), native_routes=routes,
                  native_vs_reference=comparisons, existing_displacement_and_traction_gates=checks,
                  independent_step_sensitivity=differences(refs[0], refs[1]),
                  independent_volume_order_sensitivity=differences(refs[-1], quadrature),
                  native_step_sensitivity=differences(native_fields(native[1]), native_fields(native[2])),
                  restart_checks=restart_rows, native_stored_field_and_output_audits=audits,
                  limitation='Fixed hydrostatic CONTINUOUS water and elastic MC bulk only. Full contact histories are measurements, not certified convergence. Partial volume wetting quadrature and global/local path sensitivity remain explicit. No beams, plates, Darcy flow or general plastic bulk acceptance.',
                  nc=dict(executed=False, passed=False, reason='No evidence-based fundamentally different short experiment; do not repeat precision/local-tolerance tuning.'),
                  p3b_p4_rows_complete=False, structure_water_joint_complete=False, p6_complete=False,
                  windows_ci='unverified')
    write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=32)
    args = parser.parse_args()
    result = run(args.solver, args.output, (args.steps, 2*args.steps))
    print(json.dumps({k: result[k] for k in ('passed', 'independent_reference_completed', 'existing_displacement_and_traction_gates', 'windows_ci')}, indent=2))
