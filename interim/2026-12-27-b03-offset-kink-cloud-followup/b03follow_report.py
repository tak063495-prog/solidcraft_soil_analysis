"""Summarize new observations; false acceptance is intentionally independent of tests."""
import argparse
import json
from pathlib import Path
from b03follow_math import write
from b03follow_patch import audit as patch_audit


def report(rotation, boundary, work, precision):
    read = lambda p: json.loads((p/'summary.json').read_text(encoding='utf-8'))
    r, b, w, p = map(read, (rotation, boundary, work, precision))
    maxima = {}
    for side in ('native', 'independent'):
        maxima[side] = dict(
            translation_covariance_max_abs_m=max(c['rotation_covariance'][side]['translation_max_abs_m'] for c in r['cases']),
            rotation_covariance_max_abs_rad=max(c['rotation_covariance'][side]['rotation_max_abs_rad'] for c in r['cases']),
            stiffness_covariance_relative=max(c['rotation_covariance'][side]['relative_stiffness_difference'] for c in r['cases']),
            rigid_normalized_residual=max(max(c[side]['rigid_normalized_residuals']) for c in r['cases']),
            additional_zero_modes=sorted(set(c[side]['additional_zero_modes'] for c in r['cases'])),
            force_residual_l2_kn=max(c[side]['force_residual_l2_kn'] for c in r['cases']),
            moment_residual_l2_kn_m=max(c[side]['moment_residual_l2_kn_m'] for c in r['cases']),
            work_balance_max_abs_kn_m=max(abs(c[side]['work_balance_kn_m']) for c in r['cases']),
            work_covariance_max_abs_kn_m=max(abs(c['rotation_covariance'][side]['external_work_difference_kn_m']) for c in r['cases']),
            reduced_condition_number=max(c[side]['reduced_condition_number'] for c in r['cases']),
            fields={field: dict(
                raw_max_abs={key: max(c['rotation_covariance'][side]['fields'][field]['raw_max_abs'][key] for c in r['cases']) for key in r['cases'][0]['rotation_covariance'][side]['fields'][field]['raw_max_abs']},
                derived_stress_max_abs_kpa=max(c['rotation_covariance'][side]['fields'][field]['derived_top_bottom_inplane_stress_max_abs_kpa'] for c in r['cases']),
                attachment_moment_max_abs_kn=max(c['rotation_covariance'][side]['fields'][field]['attachment_moment_max_abs_kn'] for c in r['cases'])) for field in ('solved', 'polynomial')})
    old_rotation = 1e-9
    exceedances = [dict(id=c['id'], maximum_absolute_rotation_difference_rad=c['native_independent_rotation_max_abs_rad'])
        for c in r['cases'] if c['native_independent_rotation_max_abs_rad'] > old_rotation]
    return dict(schema='b03follow-observations/1', passed=False, b03_complete=False, diagnostic_only=True,
        general_joint_acceptance=None, production_scope_expanded=False,
        base_commit='672824272e2baa1fe6ea776c8fc5889eacfc9a15',
        folder_date_is_request_label=True, execution_date_japan='2026-11-15',
        rotation=dict(native_solves=r['native_solves'], independent_solves=r['independent_solves'],
            max_tri6=r['maximum_tri6'], max_dofs=r['maximum_dofs'], maxima=maxima,
            maximum_native_independent_rotation_difference_rad=max(c['native_independent_rotation_max_abs_rad'] for c in r['cases']),
            maximum_native_independent_translation_difference_m=max(c['native_independent_translation_max_abs_m'] for c in r['cases']),
            original_rotation_absolute_comparison_magnitude_rad=old_rotation,
            exceedances_against_original_absolute_magnitude=exceedances,
            exceedances_not_erased_by_decimal_diagnostic=True,
            interpretation='Observed float64 solve sensitivity amplified by condition number; fixed polynomial fields and stiffness are covariant to arithmetic scale. Selected exact-decimal solves of the same rounded matrices greatly reduce response differences. This is limited evidence, not a new tolerance or a continuum reference.',
            full_patch_results=r, selected_decimal_results=p),
        boundary=b, work=w, polynomial_patch_consistency=patch_audit(rotation),
        contributions={
            'S1': 'No acceptance of nonzero-arm shared-ID folds/steps; free saved patch accounts do not define a joint.',
            'S2': 'Signed-arm/zero/sign-only/e-step refusals unchanged. Coupled accepted cases keep existing distinct IDs.',
            'S3': 'Quantifies face-order and current-thickness boundary intervals; six shared-node list orders unchanged. Local six vertex orders keep the signed arm in this axis-aligned family.',
            'S4': '45+45 permitted distorted flat patches: six physical rigid modes, no additional algebraic zero modes; no folded-joint rotational contract inferred.',
            'S5': 'Installation-reference work explains the saved coupled F*d-2U term; no new installation/restart implementation or tests.',
            'S6': 'New action-surface accounts on 46 saved free states plus 30 saved coupled states with two operators; no new comprehensive work acceptance criterion.',
            'S7': 'Original finite abs(e)<=t/2 guard unchanged; no extrapolation.',
            'S8': 'Only flat straight-midpoint geometry; curved midsides/directors/attachment changes still refused.',
            'S9': 'No compact/Native preparation/display changes; original headless refusals remain.',
            'S10': '2/8/18 patch sizes are different free manufactured-load problems, not a continuum h-convergence series or an extra soil refinement level.',
            'S11': 'Instantaneous work and installation reference only; cumulative endpoint-trapezoid sensitivity neither repeated nor resolved.',
            'S12': 'No JOINT/water/plastic/anchor/activity combination or ownership contract added.'},
        disposition=dict(requested_three_cloud_diagnostics_completed=True,
            unconditional_additional_cloud_diagnostic_identified=False,
            cloud_work_after_specification_decision='Implement selected joint/rotation/action-surface/reference contract, then intervention and independent small-patch equilibrium. This is blocked by undefined specification, not by cloud resources.',
            pc_required_for_current_three=False, pc_optional_small_reproduction=True,
            original_next_soil_refinement='規模上限のため確認不可・未達',
            rejected_junction_single_physical_cause=None,
            user_specification_decisions=['attachment/neutral point and finite intersection region',
                'shared/released rotations, drilling and torsional transfer',
                'symmetric arm comparison reference/thickness/representative policy',
                'work-conjugate load/reaction action surface and installation reference',
                'ownership/generation, authentication/restart, compact/Native interpretation',
                'physical reference problem and original comprehensive criteria']),
        other_rows_complete={key: False for key in ('A-01','A-02','A-03','A-04','B-01','B-02','C-01','C-02','C-03','C-04','NC-01','P6')},
        windows_ci_after_push='未確定')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('rotation', 'boundary', 'work', 'precision', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    write(a.output, report(a.rotation, a.boundary, a.work, a.precision))
