"""Measure partial NC fields and refinement without inventing acceptance limits."""
import argparse
import json
from pathlib import Path

import numpy as np
import fem_reference as fem
import probe_low_preload_nc as low


def fields(state, reference=False):
    if reference:
        rows = np.asarray(state['material_state']).reshape(-1, 14)
        return dict(displacement=np.asarray(state['displacement_m']), stress=rows[:, :6],
                    pc=np.exp(rows[:, 6]), void_ratio=np.exp(rows[:, 7])-1., plastic_strain=rows[:, 8:],
                    **{k: np.asarray(state['contacts'][k]) for k in ('jump_m', 'traction_kpa', 'plastic_slip_m')})
    rows = state['material_states']
    return dict(displacement=np.asarray(state['displacement_m']),
                stress=np.asarray([p['effective_stress_compression_positive_kpa'] for p in rows]),
                pc=np.asarray([p['preconsolidation_pressure_kpa'] for p in rows]),
                void_ratio=np.asarray([p['void_ratio'] for p in rows]),
                plastic_strain=np.asarray([p['plastic_strain_compression_positive'] for p in rows]),
                **{k: np.asarray([p[k] for p in state['joints']['points']])
                   for k in ('jump_m', 'traction_kpa', 'plastic_slip_m')})


def difference(observed, expected, active):
    out = {}
    for key, actual in observed.items():
        target = expected[key]
        if actual.shape != target.shape:
            raise ValueError('field inventory differs: '+key)
        delta = actual-target
        out[key] = dict(relative_l2_all=float(np.linalg.norm(delta)/max(np.linalg.norm(target), 1e-30)),
                        maximum_absolute_all=float(np.max(np.abs(delta))))
        if key in ('stress', 'pc', 'void_ratio', 'plastic_strain'):
            mask = np.repeat(active, len(target)//len(active))
            out[key]['relative_l2_active'] = float(np.linalg.norm(delta[mask])/max(np.linalg.norm(target[mask]), 1e-30))
    return out


def measure(root, output):
    output.mkdir(parents=True, exist_ok=False)
    profiles = []
    solver_digests = set()
    for p in sorted(root.glob('*/probe.json')):
        report = low.read(p)
        controls = report.get('numerical_controls', {})
        if report['observation'] is None and controls.get('dense_lu'):
            # Retain the initial rejected dense setup (the generator's sparse
            # DOF budget exceeded dense's existing limit). No state was solved.
            expected = low.expected_inputs(**controls)
            actual_case = low.read(p.parent/'input/stages.json')
            expected[2]['settings']['maximum_free_dofs'] = actual_case['settings']['maximum_free_dofs']
            for name, value in zip(low.NAMES, expected):
                if low.read(p.parent/'input'/(name+'.json')) != value:
                    raise ValueError('rejected setup also changed physical inputs')
        else:
            low.check_inputs(p.parent/'input', **controls)
        solver_digests.add(report['solver_sha256'])
        if report['passed'] or len(solver_digests) != 1:
            raise ValueError('unexpected certified result or different numerical solver')
        observation = report['observation']
        item = dict(name=p.parent.name, controls=controls, elapsed_seconds=report['elapsed_seconds'],
                    exit_code=report['exit_code'], stop_reason=report['stop_reason'],
                    solver_complete=report['solver_complete'], observation=observation,
                    failure_reason=report.get('solver_failure_reason'))
        if item['failure_reason'] is None:
            item['failure_reason'] = low.read(p.parent/'job/status.json').get('reason')
        if (p.parent/'job/progress.jsonl').exists():
            events = [json.loads(line)['progress'] for line in (p.parent/'job/progress.jsonl').read_text().splitlines()]
            item['new_rejected_attempts'] = [{k: e[k] for k in ('stage', 'progress', 'reason')}
                                             for e in events if e['kind'] == 'rejected']
            item['largest_observed_material_partition'] = max((e.get('maximum_material_partitions') or 0 for e in events), default=0)
        profiles.append(item)
    mesh, ground, case = low.expected_inputs()
    native = {}
    audits = []
    for p in sorted((root/'sampled-checkpoints').glob('*.json')):
        checkpoint = low.read(p); run = checkpoint['run']; state = run['last_committed']
        profile, progress_text = p.stem.rsplit('-', 1)
        native[profile, float(progress_text)] = state
        # Validate geometry/GP order against independent input quadrature.
        index = 0
        for ei in range(len(mesh['tetrahedra'])):
            q = fem.volume(mesh, ei)
            for qi, weight in enumerate(q['weight']):
                item = run['integration_points'][index]; index += 1
                if item['element_index'] != ei or item['local_index'] != qi:
                    raise ValueError('native GP order differs from independent mesh')
                if (not np.allclose(item['position_m'], q['position'][qi], rtol=0., atol=1e-8) or
                        not np.isclose(item['weight_m3'], weight, rtol=1e-10, atol=1e-12)):
                    raise ValueError('native GP geometry differs from independent mesh')
        if index != len(run['integration_points']): raise ValueError('unexpected native GP inventory')
        audits.append(dict(file=p.name, **low.inspect_run(mesh, ground, run['case'], run)))
    comparisons = []; refinement = []
    references = []
    for method, prefix in (('midpoint', 'reference-'), ('rk4', 'reference-rk4-')):
        a_dir, b_dir = root/(prefix+'128'), root/(prefix+'256')
        if not (a_dir/'reference.json').exists() or not (b_dir/'reference.json').exists():
            continue
        for folder in (a_dir, b_dir):
            report = low.read(folder/'reference.json')
            references.append({k: v for k, v in report.items() if k != 'last_accepted'} |
                              dict(name=folder.name, last_accepted_progress=report['last_accepted']['progress'],
                                   last_accepted_phase=report['last_accepted']['name']))
        for filename in ('initial_equilibrated.json', 'load.json', 'unload.json',
                         'excavate_upper-progress-0.875.json', 'excavate_upper-progress-0.9375.json'):
            if not (a_dir/filename).exists() or not (b_dir/filename).exists(): continue
            a, b = low.read(a_dir/filename), low.read(b_dir/filename)
            refinement.append(dict(method=method, point=filename,
                                   fields=difference(fields(a, True), fields(b, True), b['active_elements'])))
            for profile in ('step-1-32-newton60', 'step-1-64-newton60'):
                if filename.startswith('excavate'):
                    actual = native[profile, b['progress']]
                else:
                    run = low.read(root/profile/'observed-run.json')
                    actual = next(f['state'] for f in run['frames'] if f['name']+'.json' == filename)
                if actual['active_elements'] != b['active_elements']: raise ValueError('active set differs')
                comparisons.append(dict(method=method, native_profile=profile, point=filename,
                                         fields=difference(fields(actual), fields(b, True), b['active_elements'])))
    native_refinement = []
    for progress in (.875, .9375):
        a = native['step-1-32-newton60', progress]; b = native['step-1-64-newton60', progress]
        native_refinement.append(dict(progress=progress, fields=difference(fields(a), fields(b), b['active_elements'])))
    highest = max((p for p in profiles if p['observation'] and p['observation']['stage_index'] == 3),
                  key=lambda p: p['observation']['stage_progress'])
    report = dict(schema='soil-fem-low-preload-nc-followup-measurements/1', passed=False,
                  classification='未達・破壊判定なし', base_commit='85556a460179b96ccd95f959cdca7a3b6ff8e156',
                  highest_accepted_progress=highest['observation']['stage_progress'], highest_profile=highest['name'],
                  fully_released=False, independent_plastic_3d_full_path_accepted=False,
                  new_acceptance_thresholds=False, p3b04_complete=False, p4_complete=False, p6_complete=False,
                  profiles=profiles, sampled_native_history_audits=audits,
                  independent_references=references, independent_refinement=refinement,
                  native_refinement=native_refinement, native_vs_independent=comparisons,
                  partial_phase_restarts=low.read(root/'phase-restarts/verification.json'),
                  missing=['進行度1.0・保持荷重ゼロ・全解放後の凍結', '独立塑性3D全撤去経路の完了と刻み収束',
                           '全撤去完了までの直接経路・途中再開一致'],
                  limitations=['全場差は到達済み共通点の測定値で合格値ではない',
                               '接触端点仕事・数値散逸は本体履歴の既存JOINT監査で検査。連続率参照の仕事と同一視しない'])
    low.write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = measure(args.runs, args.output)
    print(result['highest_profile'], result['highest_accepted_progress'], 'passed=False')
