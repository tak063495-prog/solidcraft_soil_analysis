"""Reanalyse saved small assemblies and unchanged old series. No solver calls."""
import argparse
import json
from pathlib import Path
from b03cause_math import observed_order

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def prior_assessment():
    old = ROOT/'interim/2026-12-02-tri6-offset-kink-b03'
    data = json.loads((old/'convergence.json').read_text(encoding='utf-8'))
    rows = []
    for family, series in sorted(data['series'].items()):
        space = series['space']
        loads = series['load']
        row = dict(family=family, spatial_levels=2, spatial_observed_order=None,
                   original_spatial_acceptance=None, native={}, independent={})
        for key in ('native', 'independent'):
            # Preserve original values verbatim; a new norm does not overwrite them.
            row[key]['original_space_difference'] = space[key]['pairs'][0]['absolute_l2_all_points']
            row[key]['original_load_work_trend'] = loads[key]['trends']['cumulative_total_work_kn_m']
        rows.append(row)
    return dict(saved_values_unchanged=True, new_benchmark_solves=0, new_h_series_solves=0,
        original_b03_spatial_levels=2, spatial_observed_order=observed_order([0., 1.]),
        reference_is_finite_flat_discretization=True, continuum_reference_certified=False,
        richardson_or_extrapolation_used=False, literature_values_adopted_as_acceptance=False,
        b02_old_failure_and_negative_orders_replaced=False, series=rows)


def report(patches, controls):
    p = json.loads((patches/'summary.json').read_text(encoding='utf-8'))
    c = json.loads((controls/'summary.json').read_text(encoding='utf-8'))
    energy = []
    for path in sorted((patches/'assemblies').glob('*.json')):
        raw = json.loads(path.read_text(encoding='utf-8'))
        for i, element in enumerate(raw['element_observations']):
            trial = next(f['recovery'] for f in element['kernel_fields'] if f['name']=='elastic-trial')
            energy.append(dict(case=path.stem, element=i,
                energy_kn_m={k:trial[k] for k in ('membrane_energy_kn_m', 'bending_energy_kn_m',
                    'transverse_shear_energy_kn_m', 'drilling_penalty_energy_kn_m')},
                virtual_work_residual_kn_m=element['virtual_work_residual_kn_m'],
                omitted_arm_couple_work_error_kn_m=element['omitted_arm_couple_work_error_kn_m'],
                arm_energy_change_kn_m=element['membrane_bending_arm_energy_difference_kn_m']))
    maxkeys = ('global_rigid_normalized_residual', 'force_residual_l2_kn',
               'moment_residual_l2_kn_m', 'work_balance_residual_kn_m',
               'maximum_omitted_arm_couple_work_error_kn_m')
    unresolved = [
        dict(id='S1', symptom='nonzero shared-ID folds and steps',
             axes=['i','ii','iii','iv'], primary_contribution='signed-arm contract and unspecified joint/rotation ownership',
             unique_numerical_cause=None, status='未達・現行拒否', next='仕様・実装決定'),
        dict(id='S2', symptom='one-zero/e-step/only normal or e reversal',
             axes=['i','iii','iv'], primary_contribution='different arms, neutral points and work-conjugate couples',
             unique_numerical_cause=None, status='未達・現行拒否', next='仕様・実装決定'),
        dict(id='S3', symptom='normal averaging/sign/order/rounding', axes=['i','ii'],
             primary_contribution='constant facet normals; adjacent sequential tolerance and current thickness; actual near-boundary order dependence',
             unique_numerical_cause=None, status='限定観測・一般一意性未達', next='クラウド小型交差対照可能'),
        dict(id='S4', symptom='rigid arms/rotation continuity/drilling/membrane-bending coupling', axes=['iii','iv'],
             primary_contribution='per-facet A transform is consistent in tested patches; joint physical contract remains absent',
             unique_numerical_cause=None, status='限定確認・一般接合未達', next='仕様決定・小型全体対照'),
        dict(id='S5', symptom='installation references/authentication/restart for general joints', axes=['ii','iv'],
             primary_contribution='old limited references do not define a new joint state or generation',
             unique_numerical_cause=None, status='未定義・未達', next='仕様・実装決定'),
        dict(id='S6', symptom='reaction/work/general field acceptance', axes=['i','iv','v'],
             primary_contribution='missing original comprehensive gates; algebraic patch work is not physical joint validation',
             unique_numerical_cause=None, status='未定義・未達', next='対象物理と元基準の充足証拠が必要'),
        dict(id='S7', symptom='out-of-range e', axes=['i','iii'],
             primary_contribution='explicit finite abs(e)<=t/2 model guard; no solve beyond it',
             unique_numerical_cause=None, status='未達・現行拒否', next='仕様・実装決定'),
        dict(id='S8', symptom='curved midsides or changed attachment/unknown nodal e or normal', axes=['i','ii','iii'],
             primary_contribution='affine geometry, declared attachment and strict input types; curved directors absent',
             unique_numerical_cause=None, status='未達・現行拒否', next='B-01/接続仕様の決定'),
        dict(id='S9', symptom='compact/Native preparation and display', axes=['i','ii','iv'],
             primary_contribution='conversion/display/reference/authentication contract unimplemented',
             unique_numerical_cause=None, status='未達・準備拒否、表示未実行', next='仕様・実装待ち、表示はPC'),
        dict(id='S10', symptom='12 to 96 spatial series and reference solution', axes=['v'],
             primary_contribution='only two levels; literature joint/support/load/model definitions differ',
             unique_numerical_cause=None, status='次数未認定・未達', next='小型h系列と仕様判断。上限超は確認不可'),
        dict(id='S11', symptom='load-step cumulative work sensitivity', axes=['i','v'],
             primary_contribution='old endpoint quadrature with installation/removal switches; near-first-order finite differences',
             unique_numerical_cause=None, status='限定観測・物理仕事未達', next='クラウド小型交差対照可能'),
        dict(id='S12', symptom='offset plus JOINT/water/plasticity/anchors/soil activity', axes=['ii','iii','iv'],
             primary_contribution='ownership, material/water generations and initial-state/work contracts remain undefined',
             unique_numerical_cause=None, status='B-02/A群の未達・拒否を保持', next='仕様・実装決定'),
    ]
    return dict(schema='b03cause-measurements/1', passed=False, b03_complete=False, diagnostic_only=True,
        production_scope_expanded=False, general_joint_acceptance=None,
        original_all_field_acceptance=None, all_other_remainders_unmet=True,
        cloud_diagnostics_remaining=True, specification_decisions_needed=True, pc_executed=False,
        maximum_tri6=2, maximum_structural_dofs=72,
        formal_patch_observations=p, production_contract_observations=c,
        local_trial_energy_observations=energy,
        maxima={k:max(abs(r[k]) for r in p['cases']) for k in maxkeys},
        maximum_input_only_relative_stiffness_difference=max(r['input_only_independent']['relative_stiffness_difference'] for r in p['cases']),
        maximum_node_permutation_stiffness_difference=max(r['physical_id_stiffness_relative_difference'] for r in p['node_order']),
        maximum_local_virtual_work_residual_kn_m=max(abs(r['virtual_work_residual_kn_m']) for r in energy),
        prior_reference_assessment=prior_assessment(), unresolved_symptoms=unresolved,
        limitations=['free structural patch only; no soil-coupled rejected production path',
            'input-only independent polynomial reference shares flat Mindlin and drilling model, not an exact shell continuum',
            'trial energy differences are not causal percentages; no coefficient/support/tolerance tuning',
            'rotated force seed differs; only stiffness covariance is compared',
            'no general corner attachment, installation, reaction/work acceptance or new restart defined'],
        numerical_tolerances_relaxed=False, old_values_rewritten=False,
        nc_code_or_settings_changed=False, production_rust_changed=False,
        windows_ci_after_push='未確定')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--patches', type=Path, required=True)
    p.add_argument('--controls', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); result = report(a.patches, a.controls)
    a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result['maxima']))
