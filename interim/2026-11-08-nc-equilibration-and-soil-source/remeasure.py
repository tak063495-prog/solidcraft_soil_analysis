"""Reassemble saved evidence; this never resumes or reseals a solver state."""
import argparse
import json
import os
from pathlib import Path
import sys
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
import fem_reference as fem
import probe_low_preload_nc as low
import verify_wet_structures as wet
from verify_soil_source_column import audit


def serial(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [serial(v) for v in value]
    return value


def trace_history(path, run):
    origin = run['frames'][2]['state']
    inactive = ~np.asarray(run['last_committed']['active_elements'])
    offsets = np.r_[0, np.cumsum(np.bincount([p['element_index'] for p in run['integration_points']], minlength=len(inactive)))]
    removed = np.concatenate([np.arange(offsets[i], offsets[i+1]) for i in np.flatnonzero(inactive)])
    counts = {}; commits = set(); last = None
    for line in path.open(encoding='utf-8'):
        row = json.loads(line); kind = row['kind']; counts[kind] = counts.get(kind, 0)+1
        if kind == 'iterate':
            old = row['committed']; commits.add(old['stage_progress'])
            assert old['joints'] == origin['joints']
            for i in removed: assert old['material_states'][i] == origin['material_states'][i]
            last = dict(stage=row['stage'], progress=row['progress'], iteration=row['iteration'],
                        relative_force_residual=row['relative_force_residual'],
                        maximum_partition=max(row['partitions']), partitions=row['partitions'])
        if kind in ('iterate', 'line_search') and 'states' in row:
            for i in removed: assert row['states'][i] == origin['material_states'][i]
    return dict(event_counts=counts, observed_committed_progresses=sorted(commits),
                removed_materials_frozen_in_all_saved_iterates_and_candidates=True,
                removed_contact_history_frozen_in_all_saved_committed_states=True,
                last_saved_iterate=last, full_path_history_saved=False)


def decimal_audit(folder):
    from decimal import Decimal, localcontext
    rows = []
    with localcontext() as context:
        context.prec = 70
        for path in sorted(folder.glob('linear-*.json')):
            value = wet.read(path); columns = value['matrix_columns']
            x = [Decimal.from_float(v) for v in value['x']]
            rhs = [Decimal.from_float(v) for v in value['rhs']]
            residual = [sum((Decimal.from_float(columns[j][i])*x[j] for j in range(len(x))), Decimal(0))-rhs[i] for i in range(len(rhs))]
            exact = float((sum(v*v for v in residual)/max(sum(v*v for v in rhs), Decimal('1e-60'))).sqrt())
            observed = value['compensated_residual']
            assert abs(exact-observed) <= 1e-12*max(exact, 1e-30)
            if value['true_residual'] <= 1e-10: assert exact <= 1e-10
            rows.append(dict(file=path.name, decimal70_true_residual=exact,
                             compensated_residual=observed, ordinary_residual=value['ordinary_residual'],
                             linear_gate_passed=value['true_residual'] <= 1e-10))
    assert rows
    accepted = [v for v in rows if v['linear_gate_passed']]
    return dict(all_dumped_linear_solves=len(rows), passing_linear_solves=len(accepted),
                passing_solves_maximum_decimal70_true_residual=max(v['decimal70_true_residual'] for v in accepted),
                passing_solves_with_ordinary_product_above_tolerance=sum(v['ordinary_residual'] > 1e-10 for v in accepted),
                every_compensated_residual_verified_against_decimal70=True, rows=rows,
                full_excavation_acceptance=False)


def measure(folder):
    nc = []
    low.check_inputs(folder/'nc/input', excavation_increment=1/128, minimum_increment=1/4194304,
                     dense_lu=True, tangent_perturbation=1e-9)
    mesh, ground, case = [low.read(folder/'nc/input'/(n+'.json')) for n in low.NAMES]
    for name in ('scaled', 'implicit', 'dd', 'accurate'):
        run = low.read(folder/'nc'/name/'run.json')
        attempts = run['geostatic']['attempts']+run['attempts']
        accepted = [a for a in attempts if a['accepted']]
        linears = [s['true_relative_residual'] for a in accepted for s in a['linear_solves']]
        forces = [a['relative_force_residual'] for a in accepted]
        assert max(linears) <= 1e-10 and max(forces) <= 1e-8
        report = low.inspect_run(mesh, ground, case, run)
        nc.append(dict(name=name, passed=False, production_adopted=False, complete=run['complete'],
                       maximum_accepted_true_linear_residual=max(linears),
                       maximum_accepted_force_residual=max(forces),
                       accepted_linear_reports=len(linears), failure_reason=run['failure_reason'],
                       last_rejected_attempt=run['attempts'][-1],
                       remaining_holding_fraction=1-run['last_committed']['stage_progress'], **report))
    traces = {name: trace_history(folder/'nc'/name/'trace.jsonl', low.read(folder/'nc'/name/'run.json')) for name in ('dd', 'accurate')}
    precision = decimal_audit(folder/'nc/accurate-linear')
    references = []
    for count in (128, 256):
        ref = wet.read(folder/f'nc/reference{count}/reference.json')
        assert not ref['passed']
        references.append({k: v for k, v in ref.items() if k != 'last_accepted'} | dict(last_accepted_progress=ref['last_accepted']['progress'], last_accepted_phase=ref['last_accepted']['name']))
    from measure_low_preload_nc_followup import fields, difference
    sensitivity = []
    for name in ('initial_equilibrated.json', 'load.json', 'unload.json', 'excavate_upper-progress-0.875.json', 'excavate_upper-progress-0.9375.json'):
        a = wet.read(folder/'nc/reference128'/name); b = wet.read(folder/'nc/reference256'/name)
        sensitivity.append(dict(point=name, fields=difference(fields(a, True), fields(b, True), b['active_elements'])))
    soil = folder/'soil-source'
    mesh = wet.read(soil/'mesh.json'); ground = fem.normalize_ground(wet.read(soil/'ground.json'))
    case = fem.normalize_stage_case(wet.read(soil/'stages.json')); native = wet.read(soil/'stage-results.json')
    expected = wet.reference(mesh, ground, case, water_modes=('soil',))
    assert serial([{k: v for k, v in row.items() if k != 'darcy'} for row in expected]) == wet.read(soil/'independent-3d-fields.json')
    rows = wet.compare(mesh, ground, case, native, expected)
    saved = wet.read(soil/'verification.json'); assert rows == saved['frames']
    fine = wet.reference(mesh, ground, case, order=7, water_modes=('soil',))
    refinement = max(float(np.linalg.norm(a['u']-b['u'])/max(np.linalg.norm(b['u']), 1e-12)) for a, b in zip(expected, fine, strict=True))
    assert refinement == saved['independent_gl5_to_gl7_relative_difference']
    analytic = audit(mesh, ground, case, native); assert analytic == saved['analytic_soil_source']
    clean = wet.read(soil/'clean-result.json')
    for key in ('last_committed', 'frames', 'flow_audits'):
        assert wet.without_time(clean[key]) == wet.without_time(native[key])
    for index, _ in enumerate(saved['checkpoints']):
        restored = wet.read(soil/f'resumed-{index}.json'); assert restored['complete']
        for key in ('last_committed', 'frames', 'flow_audits', 'excavation_audits', 'transition_audits', 'attempts'):
            assert wet.without_time(restored.get(key)) == wet.without_time(native.get(key))
    import meshio
    point_soil = np.array([wet.Water(mesh, ground, case).soil[p['element_index']] for p in native['integration_points']])
    for index, ref in enumerate(expected):
        gp = meshio.read(soil/'exports'/f'integration-points-{index:03}.vtu').point_data; ids = gp['integration_point_id']
        np.testing.assert_array_equal(gp['material_stress_reference'], (~point_soil[ids]).astype(np.uint8))
        np.testing.assert_allclose(gp['pore_pressure_kpa'], ref['pore'][ids], rtol=1e-10, atol=1e-8)
        total = ref['stress'][ids].copy(); total[:, :3] += ref['pore'][ids, None]
        np.testing.assert_allclose(gp['total_stress_compression_positive_kpa'], total, rtol=1e-7, atol=1e-7)
    return dict(schema='soil-fem-nc-equilibration-soil-source-measurements/1',
                starting_commit='b020b492ea091288dbfffac0514c570898f8c408',
                low_preload_nc=dict(passed=False, production_adopted=False, experiments=nc,
                    full_release_verified=False, independent_plastic_3d_full_path_verified=False,
                    saved_terminal_history=traces, accurate_linear_decimal70_audit=precision,
                    independent_references=references, partial_independent_step_sensitivity=sensitivity,
                    complete_path_step_sensitivity_verified=False, resume_to_full_release_verified=False),
                soil_source=dict(limited_passed=True, scope=saved['scope'], frames=rows,
                    maximum_relative_displacement_error=max(r['relative_displacement_error'] for r in rows),
                    maximum_relative_stress_error=max(r['relative_stress_error'] for r in rows),
                    gl5_to_gl7_relative_difference=refinement, analytic=analytic,
                    restart_checks=len(saved['checkpoints']), clean_rejected_path_exact=True,
                    independent_3d_fields_exact=True, exported_pore_total_stress_checked=True),
                windows_new_tip_ci='unconfirmed', p3b_p4_rows_complete=False, p6_complete=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', required=True, type=Path); p.add_argument('--output', required=True, type=Path)
    args = p.parse_args(); result = measure(args.evidence)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
