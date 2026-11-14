"""Summarize measured observations without changing any acceptance definition."""
import argparse
import math
import numpy as np
from pathlib import Path
from common import HERE,BASE,read,write,sha

def summarize(formal,analysis,localization,controls,output):
    fresh=read(formal/'fresh-solves.json');details=read(analysis/'mesh-2-residual.json')['details']
    peak={kind:max(details[kind]['observations'],key=lambda r:r['residuals']['original']['relative_l2']) for kind in ('native','independent')}
    linear=next(r for r in details['independent']['all_direct_linear_observations'] if r['accepted_index']==peak['independent']['accepted_index'])
    n=peak['native'];stages=read(formal/'mesh-2/stages.json')['stages'];stage_name=stages[n['stage_index']-1]['name']
    attempt=next(a for a in details['native_solver_attempts']['attempts'] if a['accepted'] and a['stage']==stage_name and a['progress_to']==n['progress'])
    f=read(analysis/'field-refinement.json');observations={};fit=[];raw=[];direct=[];conditions=[];owners=[];work_errors=[]
    for case in fresh['cases']:
        observations[case['name']]=read(analysis/(case['name']+'-fields.json'))
        for kind,rows in observations[case['name']].items():
            for r in rows:
                work_errors.append(abs(r['work_recomposition_error_kn_m']))
                if 'recovery_fit_maximum_kn_m' not in r:continue
                fit.append(r['recovery_fit_maximum_kn_m']);raw.append(r['native_raw_vs_kinematic_maximum_kn_m']);conditions.append(r['fit_condition_maximum']);owners+=r['owner_counts']
                direct.append(float(np.max(np.abs(np.asarray(r['shear_recovered_global_kn_m'])-r['shear_direct_global_kn_m']))))
    loc=read(localization)['observations'];source=[]
    for r in loc:
        if r['kind']=='independent' and r['steps']==4:
            source.append(dict(tet10=12*8**r['level'],case=r['case'],nodal_source_tip_work_kn_m=r['nodal_source_tip_work_kn_m'],largest_raw_shear=r['largest_raw_shear_samples'][0],plate_energy_kn_m=r['plate_energy_kn_m']))
    source.sort(key=lambda r:r['tet10']);growth=[]
    for a,b in zip(source,source[1:]):
        da=a['nodal_source_tip_work_kn_m'];db=b['nodal_source_tip_work_kn_m'];delta={k:db[k]-da[k] for k in da}
        growth.append(dict(tet10=[a['tet10'],b['tet10']],signed_work_changes_kn_m=delta,attachment_rotation_fraction_of_nodal_work_change=delta['attachment_rotation']/delta['sum']))
    pairs=[p for p in f['pairs'] if p['series']=='space' and p['steps']==4 and p['kind']=='native']
    orders={key:math.log2(pairs[0]['absolute_l2_all_points'][key]/pairs[1]['absolute_l2_all_points'][key]) for key in ('shear_recovered_global_kn_m','work_sum_kn_m')}
    result=dict(passed=False,b02_complete=False,base_commit='1885f1f374ca417f6ceede181db2621015ccc0b8',
        historical_fixed=dict(path=str(BASE.name),record_sha256=sha(BASE/'measurements.json'),
            native_reassembly_relative=3.871041030426841e-10,independent_reassembly_relative=1.0508229745493171e-10,
            existing_reassembly_gate=1e-10,shear_differences_kn_m=[141.283499,437.863419],instantaneous_work_differences_kn_m=[.013640722,.015947827],b02_passed=False,values_not_rewritten=True),
        fresh_solves=fresh,primary_768_4_peak_residuals=peak,native_peak_accepted_solver_attempt=attempt,
        independent_peak_direct_equation=dict(**linear,physical_applied_denominator_kn=peak['independent']['denominator_kn'],
            f64_direct_equation_absolute_l2=linear['true_f64_relative']*max(linear['rhs_norm'],1.),
            physical_reassembled_absolute_l2_kn=peak['independent']['residuals']['original']['absolute_l2_kn'],
            denominator_ratio=max(linear['rhs_norm'],1.)/peak['independent']['denominator_kn'],
            physical_minus_direct_soil_absolute_l2_kn=linear['physical_minus_direct_soil_relative']*max(linear['rhs_norm'],1.),
            caution='Different equations and normalization: the coupled norm also contains rotations; no scalar identity or new acceptance gate.'),
        stopping_input_controls=read(controls),recovery_and_interpolation=dict(maximum_six_point_fit_error_kn_m=max(fit),maximum_raw_vs_kinematic_error_kn_m=max(raw),
            maximum_common_point_recovered_vs_direct_kn_m=max(direct),maximum_fit_condition=max(conditions),minimum_owners=min(owners),maximum_owners=max(owners),
            maximum_work_recomposition_error_kn_m=max(work_errors),raw_gradient_shear_is_different_observable=True,no_recovery_replacement=True),
        refinement=f,finite_three_level_observed_orders=orders,load_source_localization=source,source_work_growth=growth,
        diagnosis=dict(final_reassembly_sum_order_is_not_the_main_native_excess=True,stored_stress_accumulation_is_not_the_main_native_excess=True,
            stopping_accuracy_and_stored_discrete_solution_residual_remain=True,independent_small_direct_residual_uses_larger_rhs_norm=True,
            common_point_interpolation_does_not_explain_shear_growth=True,load_steps_do_not_explain_spatial_growth=True,
            corner_beam_attachment_rotation_and_localized_shear_dominate_observed_spatial_response=True,continuum_divergence_proved=False,
            production_implementation_bug_proved=False,production_fix_applied=False,production_tolerances_references_defaults_unchanged=True),
        acceptance=dict(limited_confirmation='Dry elastic MC, planar positive offset, one linear beam: observed field gates only. Existing nine-condition historical limited record retained separately.',
            overall_passed={k:False for k in ['B-02','NC-01','P6','A-01','A-02','A-03','A-04','B-01','B-03','C-01','C-02','C-03','C-04']},
            numerical_acceptance_all_fields=None,full_field_uniform_convergence=False,windows_ci='unconfirmed after push'),
        limits=['768/8,16 native fields observed at five saved frame endpoints; no reconstructed intermediate fields.',
                'The fresh original independent solver gives 1.1557208463335835e-10, distinct from the immutable historical 1.0508229745493171e-10; no complete attribution of the old execution-environment difference.',
                '6144 and full intermediate replay at 768/8,16 are PC procedures, not cloud results.'])
    write(output,result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('analysis',type=Path);p.add_argument('localization',type=Path);p.add_argument('controls',type=Path);p.add_argument('output',type=Path);a=p.parse_args();summarize(a.formal.resolve(),a.analysis.resolve(),a.localization.resolve(),a.controls.resolve(),a.output.resolve())
