"""Compact measurements without upgrading incomplete goals or undefined gates."""
import argparse
from pathlib import Path
import json
import math
from cases import read,write
from study import sha


def collect(formal,scope,mutations,convergence,work_audit,out):
    raw=read(formal/'measurements.json');rows=[];maximum={};all_attempts=[]
    for item in raw['cases']:
        folder=formal/item['case'];r=read(folder/'native/run.json');comparison=item['comparison']
        attempts=r['geostatic']['attempts']+r['attempts'];accepted=[a for a in attempts if a['accepted']]
        force=max(a['relative_force_residual'] for a in accepted)
        linear=max((s['true_relative_residual'] for a in accepted for s in a['linear_solves']),default=0.)
        independent_linear=0.
        with (folder/'independent.jsonl').open() as stream:
            for line in stream:
                independent_linear=max(independent_linear,json.loads(line).get('direct_linear_residual',0.))
        checks=dict(fields=comparison['existing_field_gates_passed'],native_reassembly=comparison['maximum_native_reassembled_residual']<1e-10,
                    independent_reassembly=comparison['maximum_independent_reassembled_residual']<1e-10,force=force<=1e-8,linear=linear<=1e-10,independent_direct_linear=independent_linear<1e-10)
        rows.append(dict(case=item['case'],offset_m=item['offset_m'],beam=item['beam'],nodes=item['nodes'],tet10=item['tet10'],plates=item['plates'],load_steps=item['load_steps'],
                         accepted_points=comparison['accepted_points'],existing_gates=checks,existing_bounded_gates_passed=all(checks.values()),
                         maximum_force_residual=force,maximum_linear_residual=linear,maximum_independent_direct_linear_residual=independent_linear,
                         maximum_native_reassembled_residual=comparison['maximum_native_reassembled_residual'],maximum_independent_reassembled_residual=comparison['maximum_independent_reassembled_residual'],
                         maximum_by_field=comparison['maximum_by_field'],new_native_pid=item['native_process']['pid'],new_independent_pid=item['independent_process']['pid'],
                         plain_same_binary_run_exact=item['all_run_equal_to_unobserved'],input_sha256=item['input_sha256'],solver_sha256=item['solver_sha256']))
        for field,v in comparison['maximum_by_field'].items():
            old=maximum.setdefault(field,dict(maximum_absolute=0.,maximum_pointwise_relative_l2=0.,components_compared=0,gate=v['gate']))
            old['maximum_absolute']=max(old['maximum_absolute'],v['maximum_absolute']);old['maximum_pointwise_relative_l2']=max(old['maximum_pointwise_relative_l2'],v['relative_l2']);old['components_compared']+=v['components_compared']
    restarts=read(formal/'restarts.json');refusals=read(scope/'refusals.json');tamper=read(mutations/'mutations.json')
    result=dict(passed=False,b02_complete=False,p6_complete=False,nc01_complete=False,base='e940bad4015652a6cc713d0bdcde35f062c8f546',
        execution_date_jst='2026-11-12',folder_date_is_request_label=True,
        fresh_observed_native_paths=len(rows),fresh_plain_native_control_paths=len(rows),fresh_input_only_independent_paths=len(rows),accepted_pairs=sum(x['accepted_points'] for x in rows),
        cases=rows,maximum_by_field=maximum,
        checkpoint_restart=dict(real_files=len(restarts),later_accepted_states=sum(x['later_accepted_states'] for x in restarts),
                                all_later_states_exact=all(x['all_later_states_exact'] for x in restarts),
                                whole_run_exact_except_elapsed=all(x['run_exact_except_elapsed'] for x in restarts),all_distinct_pids=all(x['producer_pid']!=x['resume_pid'] for x in restarts),records=restarts),
        mutations=dict(authenticated_rehashed=sum(x['outer_hash_recomputed'] for x in tamper['records']),
                       checksum_only=sum(x['name'].startswith('unrehashed-') for x in tamper['records']),
                       typed=sum(x['name'].startswith('fake-') for x in tamper['records']),all_refused=tamper['refusal_contract_passed'],records=tamper['records']),
        scope=dict(confirmed_cli_refusals=sum(x['refused'] for x in refusals['records']),numerical_incomplete_cases=refusals['numerical_incomplete_cases'],
                   refused_output_file_count=sum(x['result_file_exists'] for x in refusals['records'] if x['refused']),records=refusals['records']),
        scope_judgments={
            'dry_elastic_MC_planar_offset_plate':'限定受入の維持・再確認（現行範囲内、ケース別の既存ゲートのみ）',
            'dry_elastic_MC_planar_offset_plate_linear_beam':'限定受入の維持・再確認（現行範囲内、ケース別の既存ゲートのみ）',
            **{k:'未達／未対応（現行拒否を保持）' for k in ('JOINT','hydrostatic_initial','hydrostatic_staged','steady_Darcy','transient_Richards','plastic_MC_history','MCC','IDH','anchors','soil_BIRTH','soil_DEATH','soil_MATSET')},
            'new_MC_plastic_response':'未達。低強度試行は数値未完了。専用拒否成功と区別する'},
        convergence=read(convergence),load_work_audit=read(work_audit),fresh_solve_sources_unchanged=read(formal/'source-binding-before.json')==read(formal/'source-binding-after.json'),
        numerical_acceptance_all_fields=None,
        note='No general all-field/reaction-resultant/work/refinement tolerance was introduced. Unsupported combinations have no invented zero field differences. Earlier failures and historical NC unexecuted records are preserved.')
    sensitivity={}
    for series,names in [('space',['beam-positive','mesh-1','mesh-2']),('load',['beam-positive','load-8','load-16'])]:
        sensitivity[series]={}
        for kind in ('native','independent'):
            tables=[{(r['stage_index'],r['progress']):r for r in result['load_work_audit']['cases'][name][kind] if r['progress'] in (.25,.5,.75,1.)} for name in names]
            assert tables[0].keys()==tables[1].keys()==tables[2].keys()
            fields={}
            for field in ('cumulative_soil_work_kn_m','cumulative_structure_work_kn_m','cumulative_total_work_kn_m'):
                changes=[math.sqrt(sum((tables[i][point][field]-tables[i+1][point][field])**2 for point in tables[0])) for i in range(2)]
                fields[field]=dict(coarse_to_medium=changes[0],medium_to_fine=changes[1],observed_order=math.log2(changes[0]/changes[1]) if all(x>0 for x in changes) else None)
            sensitivity[series][kind]=dict(common_accepted_points=len(tables[0]),fields=fields)
    result['work_refinement']=dict(series=sensitivity,accepted=None,note='Discrete endpoint cumulative work at the same 20 progress points; no continuous-path or new convergence certificate.')
    write(out/'measurements.json',result)
    write(out/'scope-refusals.json',refusals);write(out/'checkpoint-mutations.json',tamper);write(out/'convergence.json',read(convergence));write(out/'restart-audit.json',restarts);write(out/'work-audit.json',read(work_audit))
    write(out/'work-refinement.json',result['work_refinement'])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--formal',type=Path,required=True);p.add_argument('--scope',type=Path,required=True);p.add_argument('--mutations',type=Path,required=True);p.add_argument('--convergence',type=Path,required=True);p.add_argument('--work-audit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=collect(a.formal,a.scope,a.mutations,a.convergence,a.work_audit,a.output);print(r['accepted_pairs'],r['checkpoint_restart']['real_files'])
