"""Reassemble every new arc state; keep all acceptance claims false."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from diagnose import read, rows, spectrum, state_array, write


def audit_arc(folder, evidence):
    c = read(folder/'controls.json')
    kind = c['map']
    original = read(evidence/(kind+'-complement-refined128')/'controls.json')
    B, R, weights = (np.asarray(original[k]) for k in ('B', 'transform', 'weights'))
    active = np.asarray(original['active_elements'])
    pi = original['active_gp']
    initial = next(rows(evidence/(kind+'-complement-refined128')/'samples.jsonl.gz'))
    samples = list(rows(folder/'samples.jsonl.gz'))
    linear = list(rows(folder/'linear-systems.jsonl.gz'))
    events = read(folder/'events.json')
    final, holding, group = (np.asarray(original[k]) for k in ('final_load_kn', 'holding_kn', 'vertical_holding_kn'))
    maximum_residual = 0.
    maximum_strain = 0.
    trajectory = []
    cursor = 0
    for sample, event in zip(samples, (x for x in events if x['accepted'])):
        if sample['step'] != event['step']:
            raise ValueError('arc event/sample step differs')
        s = sample['coordinate']
        if not .5 <= s <= 1:
            raise ValueError('arc sample outside component path')
        remaining = (2-2*s)*group
        allowance = 5e-14*np.maximum(1., np.abs(holding))
        if np.any(np.abs(remaining) > np.abs(holding)+allowance) or np.any(remaining*holding < -allowance):
            raise ValueError('arc accepted holding bounds exceeded')
        stress = state_array(sample, kind, active, pi)[:, :6]
        force = -np.einsum('qai,qa,q->i', B, stress, weights)
        load = final+remaining
        residual = float(np.linalg.norm(R.T@load-force)/max(np.linalg.norm(load), 1.))
        maximum_residual = max(maximum_residual, residual)
        if residual >= 1e-8:
            raise ValueError('arc independently reassembled force gate exceeded')
        if sample['contacts'] != initial['contacts']:
            raise ValueError('arc frozen contact history changed')
        if kind == 'native':
            frozen = all(sample['states'][i] == initial['states'][i] for i in range(324) if i not in pi)
        else:
            frozen = np.array_equal(np.asarray(sample['states'])[~active], np.asarray(initial['states'])[~active])
        if not frozen:
            raise ValueError('arc removed material history changed')
        strain = float(np.abs(B@np.asarray(sample['reduced_increment_m'])).sum(axis=1).max())
        maximum_strain = max(maximum_strain, strain)
        if strain > .001:
            raise ValueError('existing arc increment strain bound exceeded')
        count = len(event['events'])
        actual = linear[cursor+count-1]
        tangent = linear[cursor+count]
        n = R.shape[1]
        K = np.asarray(actual['matrix'])[:n, :n]
        if not np.array_equal(K, np.asarray(tangent['matrix'])[:n, :n]):
            raise ValueError('arc corrector/tangent matrix order differs')
        spec = spectrum(K)
        trajectory.append(dict(step=sample['step'], coordinate=s, coordinate_increment=sample['coordinate_increment'],
                               smallest_singular=spec['singular_values'][-1], symmetric_min=spec['symmetric_part_eigenvalues'][0],
                               condition_2=spec['condition_2'], relative_force_residual=residual,
                               remaining_holding_norm_kn=float(np.linalg.norm(remaining)),
                               maximum_correction_linear_residual=max(x['true_relative_residual'] for x in linear[cursor:cursor+count+1])))
        cursor += count+1
    if len(samples) != sum(x['accepted'] for x in events):
        raise ValueError('arc accepted field count differs')
    maximum_linear = max((x['true_relative_residual'] for x in linear), default=0.)
    # Every solve of these completed point-count paths must retain the original gate.
    if all(x['accepted'] for x in events) and (cursor != len(linear) or maximum_linear > 1e-10):
        raise ValueError('arc complete linear stream or gate differs')
    peak = max(range(len(trajectory)), key=lambda i: trajectory[i]['coordinate']) if trajectory else None
    report = read(folder/'report.json')
    return dict(report=report, anchor_check=read(folder/'anchor-check.json'), accepted_full_fields=len(samples),
                maximum_reassembled_force_residual=maximum_residual, maximum_true_linear_residual=maximum_linear,
                maximum_total_increment_strain_l1=maximum_strain, removed_material_and_contacts_frozen=True,
                original_component_bounds_verified=True, actual_augmented_systems=len(linear),
                trajectory=trajectory, peak_index=peak,
                peak=trajectory[peak] if peak is not None else None,
                force_coordinate_is_not_proportional_lambda=True, complete_history_accepted=False)


def material_summary(measured, evidence):
    result = {}
    for name, witnesses in measured['material_replays'].items():
        controls = read(evidence/name/'controls.json')
        B, weights = np.asarray(controls['B']), np.asarray(controls['weights'])
        records = []
        for witness in witnesses:
            points = witness.get('points', [])
            row = dict(label=witness['label'], local_errors=witness['local_errors'],
                       relative_to_archived_actual_tangent=witness.get('relative_to_archived_actual_tangent'),
                       newton_event=witness['event'], finite_differences=[], partitions={})
            for j in range(3):
                f = [p['finite_differences'][j] for p in points]
                row['finite_differences'].append(dict(width=f[0]['width'],
                    changed_branch_points=sum(any(a != p['plastic_substeps'] or b != p['plastic_substeps']
                      for a, b in q['plastic_substeps_plus_minus']) for p, q in zip(points, f)),
                    maximum_central_relative_to_ad=max(x['central_relative_to_ad'] for x in f),
                    maximum_forward_relative_to_ad=max(x['forward_relative_to_ad'] for x in f),
                    maximum_backward_relative_to_ad=max(x['backward_relative_to_ad'] for x in f)))
            for count, values in witness.get('partition_sensitivity', {}).items():
                errors = [x for x in values if 'error' in x]
                base = {x['point']: x for x in points}
                stress = []
                tangent = []
                for value in values:
                    if 'error' in value:
                        continue
                    old = base[value['point']]
                    stress.append(float(np.max(np.abs(np.asarray(value['updated_state']['effective_stress_compression_positive_kpa'])-
                                                       np.asarray(old['updated_state']['effective_stress_compression_positive_kpa'])))))
                    tangent.append(float(np.linalg.norm(np.asarray(value['implicit_tangent'])-np.asarray(old['implicit_tangent']))/
                                         np.linalg.norm(old['implicit_tangent'])))
                row['partitions'][count] = dict(errors=errors, maximum_stress_absolute_difference_kpa=max(stress, default=None),
                                                maximum_local_tangent_relative_difference=max(tangent, default=None))
            row.update(material_points=len(points), base_substep_counts=sorted(set(len(x['substeps']) for x in points)),
                       base_plastic_points=sum(x['plastic_substeps'] > 0 for x in points))
            D = np.asarray([x['implicit_tangent'] for x in points])
            mode = witness['replayed_tangent']
            left, right = (np.asarray(mode[k]) for k in ('smallest_left_vector', 'smallest_right_vector'))
            el, er = B@left, B@right
            contributions = weights*np.einsum('qa,qab,qb->q', el, D, er)
            row['actual_algorithmic_mode'] = dict(
                smallest_singular=mode['singular_values'][-1],
                material_local_minimum_singular=float(np.linalg.svd(D, compute_uv=False)[:, -1].min()),
                material_negative_symmetric_points=int((np.linalg.eigvalsh((D+D.transpose(0, 2, 1))/2)[:, 0] < 0).sum()),
                mode_contribution_sum=float(contributions.sum()),
                absolute_mode_contribution_sum=float(np.abs(contributions).sum()),
                gp_contributions=[dict(gp=p['point'], left_mode_strain=a.tolist(), right_mode_strain=b.tolist(),
                                       weighted_bilinear_contribution=float(value))
                                  for p, a, b, value in zip(points, el, er, contributions)])
            records.append(row)
        result[name] = records
    return result


def physical_fields(sample, kind):
    if kind == 'native':
        material = sample['states']
        return dict(displacement_m=np.asarray(sample['displacement_m']),
            stress_kpa=np.asarray([x['effective_stress_compression_positive_kpa'] for x in material]),
            plastic_strain=np.asarray([x['plastic_strain_compression_positive'] for x in material]),
            pc_kpa=np.asarray([x['preconsolidation_pressure_kpa'] for x in material]),
            void_ratio=np.asarray([x['void_ratio'] for x in material]))
    state = np.asarray(sample['states']).reshape(-1, 14)
    return dict(displacement_m=np.asarray(sample['displacement_m']), stress_kpa=state[:, :6],
                plastic_strain=state[:, 8:], pc_kpa=np.exp(state[:, 6]), void_ratio=np.exp(state[:, 7])-1)


def compare_fields(a, b, kind):
    pooled = {}
    contacts_exact = True
    for x, y in zip(a, b):
        fa, fb = physical_fields(x, kind), physical_fields(y, kind)
        contacts_exact &= x['contacts'] == y['contacts']
        for key in fa:
            d2 = float(np.sum((fa[key]-fb[key])**2))
            r2 = float(np.sum(fb[key]**2))
            row = pooled.setdefault(key, dict(difference_squared=0., reference_squared=0., maximum_absolute=0., maximum_state_relative_l2=0.))
            row['difference_squared'] += d2
            row['reference_squared'] += r2
            row['maximum_absolute'] = max(row['maximum_absolute'], float(np.max(np.abs(fa[key]-fb[key]))))
            row['maximum_state_relative_l2'] = max(row['maximum_state_relative_l2'], float(np.sqrt(d2)/max(np.sqrt(r2), 1e-12)))
    for row in pooled.values():
        row['pooled_relative_l2'] = float(np.sqrt(row['difference_squared'])/max(np.sqrt(row['reference_squared']), 1e-12))
    return dict(fields=pooled, all_contacts_exact=contacts_exact, complete_history_converged=False)


def run(source, evidence, output):
    tangent = read(source/'tangent/measurements.json')
    result = dict(starting_commit='113d83b8a37b3b68ce46eb00b6664f042334eae0', passed=False,
                  full_removal_accepted=False, full_release_freeze_verified=False, history_converged=False,
                  independent_plastic_3d_to_completion=False, restart_to_completion_verified=False,
                  static_limit_confirmed=False, no_other_bounded_path_proven=False, production_adopted=False,
                  residual_domain_lambda_zero=read(source/'lambda-zero/measurements.json'),
                  fresh_prefix=read(source/'fresh-prefix/verification.json'), material_replay=material_summary(tangent, evidence),
                  tangent_paths={}, arc_paths={}, restart={}, step_sensitivity={}, material_step_sensitivity={})
    for name, path in tangent['paths'].items():
        frames = []
        for frame in path['frames']:
            gp = frame['material_points']
            frames.append(dict(coordinate=frame['coordinate'], continuum_smallest_singular=frame['singular_values'][-1],
                               continuum_symmetric_min=frame['symmetric_part_eigenvalues'][0],
                               elastic_control_min=frame['elastic_control']['symmetric_part_eigenvalues'][0],
                               plastic_points=frame['plastic_points'], local_negative_symmetric_points=sum(x['local_symmetric_eigenvalues'][0] < 0 for x in gp),
                               minimum_local_singular=min(x['local_singular_values'][-1] for x in gp),
                               minimum_denominator_ratio=min(x['denominator_over_elastic_term'] for x in gp),
                               mode_contribution_sum=frame['mode_contribution_sum'], absolute_mode_contribution_sum=frame['absolute_mode_contribution_sum']))
        result['tangent_paths'][name] = dict(report=path['report'], frames=frames,
            worst_retained_condition_2=path['worst_retained_matrix']['condition_2'],
            all_saved_contacts_identical=path['all_saved_contacts_identical'])
    folders = sorted(p for p in source.glob('arc-*') if p.is_dir())
    for folder in folders:
        result['arc_paths'][folder.name] = audit_arc(folder, evidence)
        if folder.name.endswith('-resumed'):
            direct = source/folder.name.removesuffix('-resumed')
            cp = read(direct/'middle.json')['payload']
            a = [x for x in rows(direct/'samples.jsonl.gz') if x['step'] > cp['step']]
            b = list(rows(folder/'samples.jsonl.gz'))
            if a != b or read(direct/'events.json')[cp['step']+1:] != read(folder/'events.json'):
                raise ValueError('arc resumed fields/events differ')
            if read(direct/'restart.json') != read(folder/'restart.json'):
                raise ValueError('arc final sealed checkpoint differs after resume')
            direct_events = read(direct/'events.json')
            offset = sum(len(event['events'])+1 for event in direct_events if event['accepted'] and event['step'] <= cp['step'])
            tail = list(rows(direct/'linear-systems.jsonl.gz'))[offset:]
            resumed_systems = list(rows(folder/'linear-systems.jsonl.gz'))
            if tail != resumed_systems:
                raise ValueError('arc resumed actual augmented linear systems differ')
            result['restart'][direct.name] = dict(start_coordinate=cp['coordinate'], compared_full_fields=len(a),
                full_fields_events_and_final_checkpoint_exact=True,
                all_actual_augmented_linear_systems_exact=True, compared_actual_systems=len(tail), completion_verified=False)
    for kind in ('native', 'independent'):
        a = list(rows(source/f'arc-{kind}-0.0001-m4/samples.jsonl.gz'))
        b = list(rows(source/f'arc-{kind}-5e-05-m4/samples.jsonl.gz'))[1::2]
        if len(a) != len(b):
            raise ValueError('arc length alignment differs')
        result['step_sensitivity'][kind] = dict(common_arc_lengths=len(a),
            maximum_coordinate_absolute_difference=max(abs(x['coordinate']-y['coordinate']) for x, y in zip(a, b)),
            maximum_displacement_absolute_difference_m=max(float(np.max(np.abs(np.asarray(x['displacement_m'])-np.asarray(y['displacement_m'])))) for x, y in zip(a, b)),
            complete_history_converged=False, comparison_basis='same accumulated arc length, not the same load coordinate',
            all_field_comparison=compare_fields(a, b, kind))
    a = list(rows(source/'arc-independent-0.0001-m4/samples.jsonl.gz'))
    b = list(rows(source/'arc-independent-0.0001-m8/samples.jsonl.gz'))
    result['material_step_sensitivity'] = dict(common_arc_lengths=len(a),
        maximum_coordinate_absolute_difference=max(abs(x['coordinate']-y['coordinate']) for x, y in zip(a, b)),
        maximum_displacement_absolute_difference_m=max(float(np.max(np.abs(np.asarray(x['displacement_m'])-np.asarray(y['displacement_m'])))) for x, y in zip(a, b)),
        complete_history_converged=False, all_field_comparison=compare_fields(a, b, 'independent'))
    write(output, result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'evidence', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    run(a.source.resolve(), a.evidence.resolve(), a.output.resolve())
