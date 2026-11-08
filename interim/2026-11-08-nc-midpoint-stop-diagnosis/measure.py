"""Measure saved fields and recorded trials. No full-path acceptance is granted."""
import argparse
import json
import math
import shutil
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from probe_low_preload_nc import read, write, sha, inspect_run
from probe_nc_midpoint import inputs
from verify_high_order import internal, prepare, physical
import fem_reference as fem
from joint_reference import System


def transform(rows):
    n=1+max(i for row in rows for i,c in row)
    t=np.zeros((len(rows),n))
    for g,row in enumerate(rows):
        for i,c in row:t[g,i]=c
    return t


def spectrum(k,rhs):
    s=np.linalg.svd(k,compute_uv=False)
    sym=np.linalg.eigvalsh((k+k.T)/2)
    x=np.linalg.solve(k,rhs)
    nr=np.linalg.norm(rhs)
    rr=lambda a:float(np.linalg.norm(k@a-rhs)/max(nr,1e-30))
    history=[rr(x)]
    for _ in range(4):
        candidate=x+np.linalg.solve(k,rhs-k@x)
        history.append(rr(candidate));x=candidate
    extended=np.asarray(k,dtype=np.longdouble)@np.asarray(x,dtype=np.longdouble)-np.asarray(rhs,dtype=np.longdouble)
    return dict(singular_values_kn_m=s.tolist(),condition_number_2=float(s[0]/s[-1]),
                symmetric_part_eigenvalues_kn_m=sym.tolist(),
                numpy_true_linear_residuals=history,
                long_double_dot_residual=float(np.linalg.norm(extended)/max(nr,1e-30)),
                rhs_norm_kn=float(nr),nonsymmetry_relative=float(np.linalg.norm(k-k.T)/np.linalg.norm(k)))


def material_rows(states,old,indices,m):
    out=[]
    for pi in indices:
        s=states[pi];o=old[pi];stress=np.array(s['effective_stress_compression_positive_kpa'])
        p=stress[:3].mean();dev=stress.copy();dev[:3]-=p
        q=math.sqrt(1.5*(np.dot(dev[:3],dev[:3])+2*np.dot(dev[3:],dev[3:])))
        pc=s['preconsolidation_pressure_kpa'];M=m[pi]['critical_state_m'];f=(q*q+M*M*p*(p-pc))/(M*M*pc*pc)
        dep=np.array(s['plastic_strain_compression_positive'])-np.array(o['plastic_strain_compression_positive'])
        out.append(dict(point=int(pi),p_kpa=float(p),q_kpa=q,pc_kpa=pc,void_ratio=s['void_ratio'],
                        normalized_yield=float(f),q_over_Mp=q/(M*p),flow_volumetric_sign_term_kpa=float(2*p-pc),
                        delta_pc_kpa=pc-o['preconsolidation_pressure_kpa'],delta_plastic_volumetric_strain=float(dep[:3].sum())))
    return out


def measure(prior,traces,probe,out,reuse_material=None):
    out.mkdir(parents=True,exist_ok=False)
    folder=prior/'nc-resumed';mesh,ground,case=[read(folder/'input'/(n+'.json')) for n in ('mesh','ground','stages')]
    assert (mesh,ground,case)==inputs('midpoint'),'physical/numerical inputs changed'
    cp=read(folder/'job/recovery-checkpoint-0.json');run=cp['run'];oldrun=read(folder/'job/recovery-checkpoint-1.json')['run']
    assert run['last_committed']['stage_progress']==0.9514000415802002
    q,offset=prepare(mesh,run['integration_points']);cells=fem.cells(mesh)
    active=np.array(run['last_committed']['active_elements']);indices=np.concatenate([np.arange(offset[ei],offset[ei+1]) for ei in np.flatnonzero(active)])
    params={m['id']:m['mcc'] for m in ground['materials']}
    point_params=[dict(params[mesh['material_ids'][ei]],integration='midpoint') for ei in range(len(q)) for _ in q[ei]['weight']]
    audit=inspect_run(mesh,ground,case,read(folder/'observed-run.json'))
    requests=[];records=[];points=[];matrices=[];candidate_counts={};raw_events={}
    system=System(mesh,case['joint_system']);origin=run['frames'][2]['state'];c=run['continuation']
    for name in ('last-accepted','terminal-rejected'):
        events=[json.loads(line) for line in (traces/name/'trace.jsonl').read_text().splitlines()]
        raw_events[name]=events
        outcome=read(traces/name/'diagnostic.json')
        if name=='last-accepted':assert outcome['outcome']['state']==run['last_committed'],'accepted state replay mismatch'
        else:assert outcome['outcome']['error']=='numerical solve failed: true linear residual 1.4875587579753504e-10','terminal mismatch'
        assert outcome['saved_checkpoint_unchanged'] and not outcome['checkpoint_resume']
        old=(oldrun if name=='last-accepted' else run)['last_committed']
        candidate_counts[name]=dict(iterates=0,line_search_candidates=0,active_contact_points=0,contact_label_switches=0,contact_history_changes=0,removed_material_changes=0)
        for e in events:
            if e['kind'] in ('iterate','line_search'):
                counts=candidate_counts[name];counts['iterates' if e['kind']=='iterate' else 'line_search_candidates']+=1
                if 'error' in e:raise ValueError('candidate material evaluation failed')
                assert e['joint_state']==origin['joints'],'inactive contact history changed'
                assert all(not p['active'] for p in e['joint_report']['points']),'unexpected active contact'
                for ei in np.flatnonzero(~active):
                    assert e['states'][offset[ei]:offset[ei+1]]==origin['material_states'][offset[ei]:offset[ei+1]],'removed GP changed'
            if e['kind']=='iterate':
                t=transform(e['dof_rows']);bulk,y=internal(mesh,ground,e['states'],active,q,offset)
                joint=np.array(e['joint_internal']);assert np.count_nonzero(joint)==0
                full=bulk.ravel()+joint
                target=np.array(e['target']);rr=t.T@(target-full);scale=max(np.linalg.norm(target),1.)
                assert np.allclose(full,e['internal'],rtol=1e-12,atol=1e-11),'independent internal assembly mismatch'
                assert np.linalg.norm(rr-np.array(e['reduced_residual']))<1e-10,'independent residual mismatch'
                end=np.array(c['end_load']);holding=(1-e['progress'])*np.array(c['holding_load'])
                assert np.array_equal(target,end+holding),'recorded release path mismatch'
                rows=material_rows(e['states'],old['material_states'],indices,point_params)
                points.append(dict(attempt=name,iteration=e['iteration'],points=rows))
                start=len(requests)
                for pi in indices:requests.append(dict(parameters=point_params[pi],state=old['material_states'][pi],strain=e['strain_increment'][pi],partitions=e['partitions'][pi]))
                records.append(dict(attempt=name,iteration=e['iteration'],event=e,request_start=start,request_end=len(requests),transform=t))
                matrices.append(dict(attempt=name,iteration=e['iteration'],relative_force_residual=e['relative_force_residual'],
                    independent_relative_force_residual=float(np.linalg.norm(rr)/scale),
                    independent_internal_difference_kn=float(np.linalg.norm(full-np.array(e['internal']))),
                    reduced_material_internal_norm_kn=float(np.linalg.norm(t.T@bulk.ravel())),
                    reduced_contact_internal_norm_kn=float(np.linalg.norm(t.T@joint)),
                    reduced_physical_load_norm_kn=float(np.linalg.norm(t.T@end)),
                    reduced_remaining_holding_load_norm_kn=float(np.linalg.norm(t.T@holding)),
                    reduced_total_external_load_norm_kn=float(np.linalg.norm(t.T@target)),
                    relative_force_scale_kn=float(scale),residual_norm_kn=float(np.linalg.norm(rr)),
                    reduced_residual_kn=rr.tolist(),
                    yield_max=max(p['normalized_yield'] for p in rows),
                    yield_branch_band_count=sum(abs(p['normalized_yield'])<=1e-11 for p in rows),
                    softening_points=sum(p['delta_pc_kpa']<0 for p in rows),hardening_points=sum(p['delta_pc_kpa']>0 for p in rows),
                    dilative_side_points=int(sum(p['flow_volumetric_sign_term_kpa']<0 for p in rows)),
                    partitions_min=int(min(e['partitions'][i] for i in indices)),partitions_max=int(max(e['partitions'][i] for i in indices))))
            elif e['kind']=='tangent':
                assert matrices[-1]['attempt']==name and matrices[-1]['iteration']==e['iteration']
                matrices[-1]['native_fd_matrix']=spectrum(np.array(e['matrix_columns']).T,np.array(e['reduced_residual']))
    assert [candidate_counts[n]['iterates'] for n in ('last-accepted','terminal-rejected')]==[8,12],'recorded iterate coverage differs'
    assert [candidate_counts[n]['line_search_candidates'] for n in ('last-accepted','terminal-rejected')]==[14,34],'recorded candidate coverage differs'
    write(out/'material-requests.json',requests)
    if reuse_material is None:
        if probe is None:raise ValueError('material probe or verified saved material results required')
        started=time.monotonic()
        with (out/'material-probe.log').open('x') as log:
            result=subprocess.run([str(probe),str(out/'material-requests.json'),str(out/'material-results.json')],stdout=log,stderr=subprocess.STDOUT,timeout=120)
        write(out/'material-process.json',dict(exit_code=result.returncode,elapsed_seconds=time.monotonic()-started,solver_sha256=sha(probe)))
        assert result.returncode==0
    else:
        assert read(reuse_material/'material-requests.json')==requests,'saved material requests differ'
        if probe is not None:
            assert read(reuse_material/'material-process.json')['solver_sha256']==sha(probe),'saved material binary differs'
        shutil.copyfile(reuse_material/'material-results.json',out/'material-results.json')
        write(out/'material-process.json',dict(executed=False,source_process=read(reuse_material/'material-process.json'),source_results_sha256=sha(reuse_material/'material-results.json')))
    maps=read(out/'material-results.json');assert len(maps)==len(requests)
    for rec,row in zip(records,matrices,strict=True):
        ks=[np.zeros((len(mesh['nodes'])*3,)*2) for _ in range(2)];r=rec['event'];subset=maps[rec['request_start']:rec['request_end']]
        branches=0;diffs=[];position=0
        for ei in np.flatnonzero(active):
            nodes=cells[ei];dofs=np.array([[3*n+j for j in range(3)] for n in nodes]).ravel()
            for qi,(b,w) in enumerate(zip(q[ei]['b'],q[ei]['weight'],strict=True)):
                m=subset[position];pi=int(offset[ei]+qi);position+=1
                assert 'error' not in m,m.get('error')
                assert m['implicit_primal_exact'] and m['primal']==r['states'][pi],'material primal replay mismatch'
                fd=m['finite_differences'][-1];branches+=sum(a!=b for a,b in fd['plus_minus_plastic_substeps']);diffs.append(fd['relative_to_implicit'])
                for k,d in zip(ks,(m['implicit_tangent'],fd['tangent']),strict=True):k[np.ix_(dofs,dofs)]+=w*b.T@np.asarray(d)@b
        t=rec['transform'];ki,kf=[t.T@k@t for k in ks]
        row['implicit_matrix']=spectrum(ki,np.array(r['reduced_residual']))
        row['fd_to_implicit_relative_matrix_difference']=float(np.linalg.norm(kf-ki)/np.linalg.norm(ki))
        row['maximum_local_fd_to_implicit_difference']=max(diffs)
        row['fd_perturbation_branch_crossing_columns']=branches
        row['material_exact_replay_points']=len(subset)
        if 'native_fd_matrix' in row:
            ev=next(e for e in raw_events[rec['attempt']] if e['kind']=='tangent' and e['iteration']==rec['iteration'])
            row['independent_fd_matrix_assembly_relative_difference']=float(np.linalg.norm(kf-np.array(ev['matrix_columns']).T)/np.linalg.norm(kf))
    write(out/'material-history.json',points)
    write(out/'iterate-measurements.json',matrices)
    write(out/'contact-measurements.json',candidate_counts)
    # Rebase the retained domain at this checkpoint without advancing release.
    # Its lambda=0 load is physical load plus the CURRENT remaining holding load.
    state=run['last_committed'];end=np.array(c['end_load']);remain=(1-state['stage_progress'])*np.array(c['holding_load'])
    bulk,_=internal(mesh,ground,state['material_states'],active,q,offset)
    basis,ce=fem.staged_constraint_space(mesh,case,state)
    rr=np.concatenate([b.T@(bulk[n]-np.reshape(end+remain,(-1,3))[n]) for n,b in basis.items()])
    rebase=dict(relative_force_residual=float(np.linalg.norm(rr)/max(np.linalg.norm(end+remain),1.)),
                applied_load_exact=np.array_equal(np.array(state['applied_nodal_load_kn']),end+remain),
                remaining_holding_fraction=1-state['stage_progress'],constraint_error_m=ce,
                independent_stored_stress_reassembly=True,input_only_plastic_solution=False)
    report=dict(schema='soil-fem-nc-midpoint-stop-diagnostic/1',passed=False,history_converged=False,
        production_adopted=False,full_removal_reached=False,checkpoint_resume_performed=False,
        independent_plastic_3d_full_path=False,independent_reference=False,
        source_commit='09ff926d82e574561393be2ac877cf9f5627d750',source_parent='c4de9fdafa4d7147f4f2bd8c8c877ca969b269b6',
        checkpoint_progress=state['stage_progress'],accepted_increments=audit['accepted_increments'],
        existing_tolerances=dict(force=1e-8,true_linear=1e-10),
        original_stage_lambda_zero=audit['zero_release_check'],checkpoint_rebased_lambda_zero=rebase,
        checkpoint_stored_force=audit['last_committed_force_check'],
        contact=candidate_counts,material_map_requests=len(requests),
        last_accepted_state_exact=True,terminal_rejection_exact=True,
        original_failure=read(folder/'observed-run.json')['failure_reason'],
        saved_field_audit=audit,windows_ci='unverified',
        missing_acceptance=['full removal','post-release freeze','independent plastic 3D full path','step sensitivity of both implementations','restart through full removal'])
    report['qualified_numerical_improvement_found']=False
    report['cause_assessment']={
        'observed_proximal_cause':'MCC softening and global tangent ill-conditioning; dense LU true-residual gate fails before nonlinear force convergence',
        'contact_switching_observed':False,
        'repartition_resets_to_committed_iterate':sum(e['kind']=='repartition' for e in raw_events['terminal-rejected']),
        'terminal_fd_vs_implicit_relative_difference':matrices[-1]['fd_to_implicit_relative_matrix_difference'],
        'singularity_or_physical_collapse_proven':False,
        'offline_numpy_linear_pass_is_not_nonlinear_or_path_acceptance':True,
        'decision':'stop after diagnosis, no continuation or production adoption',
        'remaining_uncertainty':'Relative roles of softening, finite first-yield branch differences, material repartition restarts and floating-point cancellation are not uniquely isolated; no demonstrated method reaches full removal under unchanged gates'}
    write(out/'measurements.json',report)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prior',type=Path,required=True);p.add_argument('--traces',type=Path,required=True);p.add_argument('--material-probe',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--reuse-material',type=Path);a=p.parse_args()
    r=measure(a.prior.resolve(),a.traces.resolve(),a.material_probe.resolve() if a.material_probe else None,a.output.resolve(),a.reuse_material.resolve() if a.reuse_material else None);print('Diagnostic complete; passed=false')
