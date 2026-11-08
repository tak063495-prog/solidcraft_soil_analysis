"""Independent saved-force/linear audits, separated from input-only 3D fields.

No Rust or constitutive solve is invoked. Numerical/physical acceptance remains
false. Independent reference fields are generated only by the separate driver.
"""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import probe_low_preload_nc as low
from probe_nc_midpoint import inputs
from measure_low_preload_nc_followup import fields, difference
from verify_high_order import internal, prepare


def write(p,value): low.write(p,value)
def read(p): return low.read(p)
def events(p): return [json.loads(line) for line in p.read_text().splitlines()]
def linear_audit(requests,responses):
    if not requests or len(requests)!=len(responses):
        raise ValueError('missing independent linear systems')
    spec=importlib.util.spec_from_file_location('prior_nc_linear_audit',ROOT/'interim/2026-11-08-nc-compensated-linear-replay/audit_linear.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.audit(requests,responses)

def force_history(mesh,ground,case,run,trace,retain=False,responses=None):
    q,offsets=prepare(mesh); state=run['last_committed']; origin=run['frames'][2]['state']
    active=np.array(state['active_elements']);continuation=run['continuation']
    rows=[];previous=None;pending=None;matrices=[];candidate_errors=[];response_index=-1
    checked_states=0;repartitions=[]
    if not trace or not any(e['kind']=='iterate' for e in trace):raise ValueError('missing Newton trace')
    for e in trace:
        kind=e['kind']
        if kind in ('iterate','line_search') and 'error' not in e:
            checked_states+=1
            if e['joint_state']!=origin['joints'] or any(p['active'] for p in e['joint_report']['points']):
                raise ValueError('inactive contact history or activity changed')
            for ei in np.flatnonzero(~active):
                if e['states'][offsets[ei]:offsets[ei+1]]!=origin['material_states'][offsets[ei]:offsets[ei+1]]:
                    raise ValueError('removed material history changed')
        if kind=='iterate':
            if e['committed']!=state:raise ValueError('trial changed immutable committed state')
            if pending is not None and retain:
                if e['displacement_m']!=previous['displacement_m'] or e['partitions']!=pending['partitions']:
                    raise ValueError('repartition failed to retain evaluated current iterate')
            pending=None;previous=e
            force,_=internal(mesh,ground,e['states'],active,q,offsets)
            joint=np.array(e['joint_internal']);
            if np.count_nonzero(joint):raise ValueError('inactive contact has nonzero force')
            full=force.ravel()+joint;target=np.array(e['target'])
            expected=np.array(continuation['end_load'])+(1-e['progress'])*np.array(continuation['holding_load'])
            if not np.array_equal(target,expected):raise ValueError('affine release load changed')
            dofs=e['dof_rows']; n=max(r for terms in dofs for r,c in terms)+1
            transform=np.zeros((len(dofs),n))
            for g,terms in enumerate(dofs):
                for r,c in terms:transform[g,r]=c
            reduced=transform.T@(target-full)
            internal_error=float(np.linalg.norm(full-np.array(e['internal'])))
            reduced_error=float(np.linalg.norm(reduced-np.array(e['reduced_residual'])))
            if max(internal_error,reduced_error)>1e-10:raise ValueError('saved force reassembly differs')
            active_indices=np.concatenate([np.arange(offsets[i],offsets[i+1]) for i in np.flatnonzero(active)])
            pc=[e['states'][i]['preconsolidation_pressure_kpa'] for i in active_indices]
            rows.append(dict(iteration=e['iteration'],relative_force_residual=e['relative_force_residual'],
                independent_stored_stress_force_residual=float(np.linalg.norm(reduced)/max(np.linalg.norm(target),1.)),
                independent_internal_difference_kn=internal_error,reduced_difference_kn=reduced_error,
                maximum_material_partition=max(e['partitions'][i] for i in active_indices),
                active_pc_range_kpa=[min(pc),max(pc)],
                softening_points=sum(e['states'][i]['preconsolidation_pressure_kpa']<state['material_states'][i]['preconsolidation_pressure_kpa'] for i in active_indices)))
        elif kind=='tangent':
            response_index+=1
            a=np.array(e['matrix_columns']).T;s=np.linalg.svd(a,compute_uv=False)
            matrices.append(dict(iteration=e['iteration'],minimum_singular_value_kn_m=float(s[-1]),condition_number_2=float(s[0]/s[-1])))
        if kind in ('line_search','repartition') and responses is not None and 'displacement_m' in e:
            response=responses[response_index]
            if response['accepted']:
                x=response['candidate_x'];u=np.asarray(previous['displacement_m']).ravel();expected=u.copy()
                for g,terms in enumerate(previous['dof_rows']):
                    expected[g]+=e['alpha']*sum(c*x[r] for r,c in terms)
                candidate_errors.append(float(np.max(np.abs(expected-np.asarray(e['displacement_m']).ravel()))))
        if kind=='repartition':pending=e;repartitions.append(e['iteration'])
    counts=Counter(e['kind'] for e in trace)
    linear=[e['event']['linear_solve'] for e in trace if e['kind']=='progress' and e['event']['kind']=='linear_solved']
    return dict(passed=False,history_converged=False,production_adopted=False,
        iterate_count=len(rows),line_search_candidate_count=counts['line_search'],linear_solves=len(linear),
        repartition_iterations=repartitions,all_repartitions_retain_current_iterate=retain,
        maximum_successful_original_linear_residual=max(x['true_relative_residual'] for x in linear),
        final_relative_force_residual=rows[-1]['relative_force_residual'],
        minimum_relative_force_residual=min(x['relative_force_residual'] for x in rows),
        maximum_contact_active_points=0,checked_material_and_contact_fields=checked_states,
        removed_material_and_contact_history_frozen=True,
        maximum_candidate_displacement_reconstruction_error_m=max(candidate_errors,default=None),
        independent_plastic_3d_solution=False,rows=rows,spectra=matrices)


def measure(root,out):
    out.mkdir(parents=True,exist_ok=False)
    exp=root/'experiment';folder=exp/'prior/evidence/prior/nc-resumed'
    mesh,ground,case=[read(folder/'input'/(n+'.json')) for n in low.NAMES]
    if (mesh,ground,case)!=inputs('midpoint'):raise ValueError('5 kPa original profile mismatch')
    old=read(folder/'job/recovery-checkpoint-0.json')['run']
    saved_sha=low.sha(folder/'job/recovery-checkpoint-0.json')
    if saved_sha!='f0b53377010227aead58171b1f5923a6639f9c1c9b1e3132bc7ea47c6b8ecd6d':raise ValueError('changed sealed checkpoint')
    controls=read(exp/'controls.json')
    if not all(controls.values()):raise ValueError('control reproduction differs')
    for label,request,response in [('gate','gate-matrices.json','gate-linear-results.json'),('trial','trial-matrices.json','trial-linear-results.json')]:
        report=linear_audit(read(exp/request),read(exp/response));write(out/(label+'-linear-audit.json'),report)
        if report!=read(exp/(label+'-linear-audit.json')):raise ValueError('saved independent linear audit differs')
    control_linear=linear_audit(read(exp/'prior/candidate-matrices.json'),read(exp/'prior/candidate-linear-results.json'))
    if control_linear!=read(exp/'prior/candidate-linear-audit.json'):raise ValueError('prior control linear audit differs')
    write(out/'control-linear-audit.json',control_linear)
    trace_reports={}
    for label in ('control','retain'):
        report=force_history(mesh,ground,case,old,events(exp/(label+'-trial/trace.jsonl')),label=='retain',
                             read(exp/'trial-linear-results.json') if label=='retain' else None)
        trace_reports[label]=report;write(out/(label+'-trace-audit.json'),report)
    fresh=read(root/'fresh-head-midpoint/observed-run.json');resumed=read(root/'fresh-head-resume/observed-run.json')
    state_equal=fresh['last_committed']==old['last_committed']
    if not state_equal:raise ValueError('fresh head state differs from sealed old state')
    current=low.inspect_run(mesh,ground,case,fresh)
    sealed=read(root/'fresh-head-resume/job/inputs/checkpoint.json')['run']
    if sealed['last_committed']!=fresh['last_committed']:raise ValueError('imported checkpoint state differs')
    audited={**fresh['last_committed'],'applied_nodal_load_kn':
        (np.array(sealed['continuation']['end_load'])+np.array(fresh['last_committed']['excavation_release_load_kn'])).tolist()}
    if audited['applied_nodal_load_kn']!=fresh['last_committed']['applied_nodal_load_kn']:
        raise ValueError('restart lambda-zero load differs')
    q,offsets=prepare(mesh);system=low.System(mesh,case['joint_system'])
    current_lambda_zero=low.force_check(mesh,ground,case,audited,q,offsets,system)
    low.compare_prefix(fresh,resumed)
    restart=read(root/'fresh-head-resume/experiment.json')['restart_observation']
    native={};reference={};sample_audits=[]
    paths=root/'partial-paths'
    for n in (128,256):
        expected=inputs('midpoint');expected[2]['stages'][2]['initial_increment']=1/n
        actual=tuple(read(paths/('native-'+str(n))/(key+'.json')) for key in low.NAMES)
        if actual!=expected:raise ValueError('partial native input mismatch')
        for progress in (.875,.9375):
            run=read(paths/('native-'+str(n))/(str(progress)+'.json'))['run'];native[n,progress]=run['last_committed']
            sample_audits.append(dict(profile=n,progress=progress,**low.inspect_run(*actual,run)))
        ref=read(paths/('reference-'+str(n))/'reference.json')
        if not ref['input_only'] or ref['uses_native_fields_or_fractions'] or ref['passed']:raise ValueError('unqualified reference')
        for progress in (.875,.9375):reference[n,progress]=read(paths/('reference-'+str(n))/('excavate_upper-progress-'+str(progress)+'.json'))
    comparisons=[]
    for progress in (.875,.9375):
        active=native[256,progress]['active_elements']
        comparisons.append(dict(progress=progress,
            native_128_to_256=difference(fields(native[128,progress]),fields(native[256,progress]),active),
            independent_128_to_256=difference(fields(reference[128,progress],True),fields(reference[256,progress],True),active),
            native_256_vs_independent_256=difference(fields(native[256,progress]),fields(reference[256,progress],True),active)))
    write(out/'partial-field-comparisons.json',dict(passed=False,history_converged=False,
        independent_input_only=True,full_path=False,comparisons=comparisons,native_sample_audits=sample_audits,
        references=[{k:v for k,v in read(paths/('reference-'+str(n))/'reference.json').items() if k!='last_accepted'} |
                    dict(last_accepted_progress=read(paths/('reference-'+str(n))/'reference.json')['last_accepted']['progress']) for n in (128,256)]))
    report=dict(schema='soil-fem-nc-strict-removal-followup/1',passed=False,history_converged=False,production_adopted=False,
        base_commit='6429e8366ad7b02ab546cc9a9f7ff321b693962b',
        highest_accepted_progress=current['stage_progress'],accepted_increments=current['accepted_increments'],
        fresh_head_run=True,fresh_head_committed_state_exact_to_old_sealed_state=state_equal,
        fresh_failure_reason=fresh['failure_reason'],fresh_process=read(root/'fresh-head-midpoint/process-observation.json'),
        history_audit=current,current_restart_lambda_zero=current_lambda_zero,
        restart_observation=restart,fresh_and_resumed_physical_run_equal_except_elapsed=True,
        fresh_restart_process=read(root/'fresh-head-resume/process-observation.json'),
        experiment=read(exp/'bounded-result.json'),
        experiment_scope='one minimum increment, not public checkpoint resume and not a new committed increment',
        experimental_retain_failure=read(exp/'retain-trial/diagnostic.json')['outcome'],
        saved_checkpoint_sha256=saved_sha,
        full_removal_reached=False,fully_released_and_frozen=False,
        independent_plastic_3d_full_path=False,both_implementations_full_path_step_sensitivity=False,
        full_path_restart_verified=False,partial_field_comparison='partial-field-comparisons.json',
        p3b_rows_complete=False,p4_rows_complete=False,structure_water_joint_complete=False,p6_complete=False,
        windows_ci_numerical_acceptance=False,priority2_new_solves=False,
        retained_wet_joint=dict(passed=False,history_converged=False,
            production_512_max_contact_increment_m=2.5e-6,
            production_512_vs_independent_4096_64_plastic_slip_percent=5.083166083259321,
            independent_2048_64_to_4096_64_plastic_slip_percent=0.9535071429812697,
            independent_512_64_to_4096_64_plastic_slip_percent=5.083092768505683,
            new_measurement=False),
        missing=['全上部撤去1.0と全解放・凍結','撤去中の全経路材料/接触履歴',
                 '入力のみの独立塑性3D全場の全経路','両実装の全経路刻み感度','完了までの途中再開一致'])
    write(out/'measurements.json',report)
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--evidence',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();r=measure(a.evidence.resolve(),a.output.resolve());print(json.dumps({k:r[k] for k in ('passed','history_converged','highest_accepted_progress','accepted_increments')},indent=2))
