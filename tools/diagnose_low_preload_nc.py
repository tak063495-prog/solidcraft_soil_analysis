"""Measure saved 5 kPa NC fields without certifying excavation or stability.

Native fields are used here only for diagnostics, never to generate an
independent reference solution. The production solver and seals stay intact.
"""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import fem_reference as fem
import mcc_rate_fem as rate
from low_preload_nc_reference import ClosedStickSystem
from measure_low_preload_nc_followup import fields
import probe_low_preload_nc as low
from verify_high_order import physical, prepare


def matrix(model, local):
    ke = np.einsum('eqai,eqab,eqbj,eq->eij', model.rb,
                   local.reshape(*model.weight.shape, 6, 6), model.rb, model.weight, optimize=True)
    values = np.bincount(model.inverse, weights=ke.ravel()[model.mask],
                         minlength=len(model.sparse_rows))
    k = model.kr.copy()
    k[model.sparse_rows, model.sparse_cols] += values
    return k


def affine_load_check(model, state, endload, holding):
    expected = endload+(1-state['stage_progress'])*holding
    actual = np.asarray(state['applied_nodal_load_kn'])
    error = actual-expected
    relative = float(np.linalg.norm(model.reduce(error))/max(np.linalg.norm(actual), 1.))
    # The original force tolerance, not a threshold inferred from this run.
    if relative > 1e-8:
        raise ValueError('saved excavation load does not match original affine release')
    return dict(maximum_absolute_kn=float(np.max(np.abs(error))), relative_free_error=relative,
                remaining_holding_norm_kn=float(np.linalg.norm((1-state['stage_progress'])*holding)),
                remaining_holding_relative_fraction=1-state['stage_progress'])


def state_diagnostics(mesh, ground, case, run):
    state = run['last_committed']
    model = ClosedStickSystem(mesh, ground, state['active_elements'], case)
    f = fields(state)
    material = np.column_stack((f['stress'], np.log(f['pc']), np.log1p(f['void_ratio']),
                                f['plastic_strain'])).reshape(len(mesh['tetrahedra']), -1, 14)
    material = material[state['active_elements']].reshape(-1, 14)
    holding = np.asarray(run['excavation_audits'][0]['initial_holding_load_kn'])
    endload = physical(mesh, ground, np.asarray(state['active_elements']), 0., True,
                       joint_system=model.joint)[0].ravel()
    yf, elastic, normal, dn, den = rate.terms(material, model.parameters)
    # This selects the existing continuum loading/unloading branches using
    # the release direction. It does not solve a new equilibrium state.
    plastic = model.select(material, -holding, yf >= -2e-9)
    ep = elastic.copy()
    ep[plastic] -= np.einsum('qa,qb,q->qab', dn[plastic], dn[plastic], 1/den[plastic])
    k = matrix(model, ep)
    eigenvalues, modes = np.linalg.eigh((k+k.T)/2)
    singular = np.linalg.svd(k, compute_uv=False)
    rhs = model.reduce(-holding)
    du = np.linalg.solve(k, rhs)
    local_min = np.linalg.eigvalsh(ep)[:, 0]
    return dict(progress=state['stage_progress'], active_material_points=len(material),
                active_contact_points=int(model.live.sum()),
                stored_contact_statuses={s: sum(p['status'] == s for p in state['joints']['points'])
                                         for s in sorted({p['status'] for p in state['joints']['points']})},
                normalized_yield_range=[float(yf.min()), float(yf.max())],
                selected_plastic_loading_points=int(plastic.sum()),
                mean_pressure_range_kpa=[float(material[:, :3].mean(1).min()),
                                         float(material[:, :3].mean(1).max())],
                negative_local_continuum_tangents=int(np.sum(local_min < 0)),
                continuum_symmetric_eigenvalues_kn_m=eigenvalues.tolist(),
                continuum_singular_values_kn_m=singular.tolist(),
                continuum_condition_number=float(singular[0]/singular[-1]),
                continuum_asymmetry_relative=float(np.linalg.norm(k-k.T)/np.linalg.norm(k)),
                smallest_mode_release_projection_kn=float(modes[:, 0]@rhs),
                reduced_displacement_rate_norm_m=float(np.linalg.norm(du)),
                continuum_linear_true_relative_residual=float(np.linalg.norm(k@du-rhs)/np.linalg.norm(rhs)),
                affine_release=affine_load_check(model, state, endload, holding),
                history_and_force=low.inspect_run(mesh, ground, case, run),
                interpretation='保存場の連続体接線診断。離散戻り写像接線・独立参照・安定/破壊判定ではない')


def progress_diagnostics(path):
    import json
    events = [json.loads(line)['progress'] for line in path.read_text().splitlines()]
    starts = [i for i, e in enumerate(events) if e['kind'] == 'attempt_started']
    trials = []
    for a, b in zip(starts, starts[1:]+[len(events)]):
        group = events[a:b]
        newton = [e for e in group if e['kind'] == 'newton']
        end = group[-1]
        if not newton:
            continue
        trials.append(dict(progress=group[0]['progress'], stage=group[0]['stage'], outcome=end['kind'],
                           reason=end.get('reason'),
                           residuals=[e['relative_force_residual'] for e in newton],
                           maximum_partitions=[e['maximum_material_partitions'] for e in newton],
                           linear_true_residuals=[e['linear_solve']['true_relative_residual']
                                                  for e in group if e.get('linear_solve')]))
    return trials


def material_measurements(job, mesh, ground, executable, output):
    previous = low.read(job/'recovery-checkpoint-0.json')['run']['last_committed']
    current = low.read(job/'recovery-checkpoint-1.json')['run']['last_committed']
    if previous['stage_progress'] >= current['stage_progress']:
        raise ValueError('expected consecutive generation 258/259 checkpoints')
    q, offsets = prepare(mesh)
    delta = np.asarray(current['displacement_m'])-np.asarray(previous['displacement_m'])
    requests = []
    parameters = {m['id']: m['mcc'] for m in ground['materials']}
    trials = progress_diagnostics(job/'progress.jsonl')
    accepted = next(t for t in reversed(trials) if t['stage'] == 'excavate_upper'
                    and t['progress'] == current['stage_progress'] and t['outcome'] == 'accepted')
    bound = max(accepted['maximum_partitions'])
    # Per-point partitions were not persisted. Sweep the observed upper bound
    # instead of claiming to have reconstructed the hidden inventory.
    if not 1 <= bound <= 32:
        raise ValueError('material diagnostic sweep exceeds its small-case budget')
    for ei, quadrature in enumerate(q):
        if not previous['active_elements'][ei]:
            continue
        strain = -np.einsum('qai,i->qa', quadrature['b'], delta[fem.cells(mesh)[ei]].ravel())
        for qi, de in enumerate(strain):
            for partitions in range(1, bound+1):
                requests.append(dict(parameters=parameters[mesh['material_ids'][ei]],
                                     state=previous['material_states'][offsets[ei]+qi],
                                     strain=de.tolist(), partitions=partitions))
    low.write(output/'material-requests.json', requests)
    result = subprocess.run([str(executable), str(output/'material-requests.json'),
                             str(output/'material-maps.json')], capture_output=True, timeout=60, check=False)
    (output/'material-probe.log').write_bytes(result.stdout+result.stderr)
    if result.returncode:
        raise ValueError('material-map probe failed')
    points = low.read(output/'material-maps.json')['points']
    selected = []
    for i in range(len(points)//bound):
        target = current['material_states'][i]
        def mismatch(p):
            if 'error' in p:
                return float('inf')
            return max(float(np.max(np.abs(np.asarray(p['updated_state'][key])-np.asarray(value))))
                       for key, value in target.items())
        candidates = points[i*bound:(i+1)*bound]
        best = min(candidates, key=mismatch)
        selected.append(dict(point=i, best_partition=candidates.index(best)+1,
                             maximum_absolute_state_component_error=mismatch(best),
                             primal_state_exact=best.get('primal_state_exact', False),
                             input_state_unchanged=best.get('input_state_unchanged', False),
                             finite_difference_relative_errors=[d['relative_to_implicit']
                                                                 for d in best.get('finite_differences', [])],
                             branch_crossing_columns=[sum(a != b for a, b in d['plus_minus_plastic_substeps'])
                                                      for d in best.get('finite_differences', [])]))
    return dict(points=selected, swept_partition_upper_bound=bound,
                local_map_errors=sum('error' in p for p in points),
                widths=[1e-7, 1e-8, 1e-9],
                limitations=['保存された最大分割10以内の候補照合。点ごとの実分割数の証明ではない',
                             '最終受理増分の材料比較のみ。失敗Newton反復の全材料場は未保存'])


def diagnose(previous_job, profiles, executable, output):
    output.mkdir(parents=True, exist_ok=False)
    controls = dict(excavation_increment=1/128, minimum_increment=1/16384, dense_lu=True)
    low.check_inputs(previous_job/'inputs', job=True, **controls)
    pointer = low.read(previous_job/'recovery.json')
    checkpoint = previous_job/pointer['checkpoint']
    if low.sha(checkpoint) != pointer['checkpoint_sha256']:
        raise ValueError('original checkpoint checksum mismatch')
    mesh, ground, case = low.expected_inputs(**controls)
    original = low.read(checkpoint)['run']
    original_solver = low.read(previous_job/'job.json')['binding']['identity']['solver_binary_sha256']
    samples = [dict(name='original-highest', **state_diagnostics(mesh, ground, case, original))]
    reports = []
    for folder in profiles:
        report = low.read(folder/'probe.json')
        if report['passed'] or report['solver_sha256'] != original_solver:
            raise ValueError('diagnostics require the original numerical solver and uncertified reports')
        low.check_inputs(folder/'input', **report['numerical_controls'])
        run = low.read(folder/'observed-run.json')
        samples.append(dict(name=folder.name, **state_diagnostics(mesh, ground, run['case'], run)))
        reports.append(dict(name=folder.name, probe=report,
                            trials=progress_diagnostics(folder/'job/progress.jsonl')))
    maps = material_measurements(previous_job, mesh, ground, executable, output)
    summary = dict(schema='soil-fem-low-preload-nc-cause-diagnosis/1', passed=False,
                   classification='未収束・破壊判定なし', base_commit='8e772568f5b9f773398b984289463b0905156cae',
                   highest_accepted_progress=max(s['progress'] for s in samples),
                   checkpoint_sha256=low.sha(checkpoint), source_job_id=pointer['job_id'],
                   original_solver_sha256=low.read(previous_job/'job.json')['binding']['identity']['solver_binary_sha256'],
                   material_probe_sha256=low.sha(executable), samples=samples, profiles=reports,
                   material_maps=maps, fully_released=False,
                   independent_plastic_3d_full_path_accepted=False, p3b04_complete=False, p6_complete=False,
                   missing=['進捗1.0と全解放・凍結', '独立塑性3D全経路と両実装の刻み感度',
                            '全撤去完了までの再開一致', '失敗反復の全材料場・分割数・線探索候補の計測'],
                   no_new_acceptance_thresholds=True)
    summary['evidence_sha256'] = {p.name: low.sha(p) for p in output.iterdir() if p.is_file()}
    low.write(output/'measurements.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous-job', type=Path, required=True)
    parser.add_argument('--profiles', nargs='+', type=Path, required=True)
    parser.add_argument('--material-probe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = diagnose(args.previous_job, args.profiles, args.material_probe.resolve(), args.output)
    print(result['highest_accepted_progress'], 'passed=False')
