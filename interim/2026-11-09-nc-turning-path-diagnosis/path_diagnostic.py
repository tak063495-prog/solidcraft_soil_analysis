"""Bounded displacement-controlled equilibrium diagnostic using the native MCC map.

Not a production continuation, independent solution, checkpoint, or acceptance.
The removed elements/contact stay absent. Load progress may turn backwards.
"""
import argparse
import os
import json
import sys
import subprocess
import time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from verify_high_order import prepare

def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def run(evidence,out,step_size,count):
    out.mkdir(parents=True,exist_ok=False)
    inp=evidence/'prior/fresh-head-midpoint/input'
    mesh=read(inp/'mesh.json');ground=read(inp/'ground.json');case=read(inp/'stages.json')
    q,offsets=prepare(mesh)
    fields=read(evidence/'field-records.json')
    seed=next(r for r in fields if r['label']=='terminal-rejected:0')['event']
    cp=seed['committed']; active=cp['active_elements']
    assert not any(p['active'] for p in seed['joint_report']['points'])
    assert not np.any(seed['joint_internal'])
    nd=len(seed['reduced_residual']);R=np.zeros((3*len(mesh['nodes']),nd))
    for g,row in enumerate(seed['dof_rows']):
        for i,c in row:R[g,i]=c
    rb={};globalb={};pointmap=[]
    for ei in np.flatnonzero(active):
        eq=q[ei];gd=(eq['nodes'][:,None]*3+np.arange(3)).ravel()
        rb[int(ei)]=eq['b']@R[gd];globalb[int(ei)]=eq['b']
        for qi in range(len(eq['weight'])):pointmap.append((int(ei),qi,int(offsets[ei]+qi)))
    parameters={m['id']:{**m['mcc'],'integration':case['settings']['mcc_integration']} for m in ground['materials']}
    checkpoint=read(evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')
    assert checkpoint['run']['last_committed']==cp
    holding=np.array(checkpoint['run']['excavation_audits'][0]['initial_holding_load_kn']).ravel()
    load0=np.array(cp['applied_nodal_load_kn']).ravel();lam0=cp['stage_progress']
    L=-R.T@holding
    k0=np.array(read(evidence/'assembled-matrices.json')['last-accepted:7:ad'])
    _,_,vt=np.linalg.svd(k0);plane=vt[-1]
    rate=np.linalg.solve(k0,L)
    if plane@rate<0:plane=-plane
    write(out/'controls.json',dict(passed=False,production_adopted=False,
        seed_progress=lam0,step_size_m=step_size,maximum_steps=count,force_tolerance=1e-8,
        linear_tolerance=1e-10,geometric_constraint_tolerance_m=1e-12,maximum_newton_iterations=20,
        fixed_material_partitions=1,maximum_material_l1_strain=.001,minimum_corrector_iterations=1,plane=plane.tolist(),
        residual_constraint_row_scale_kn_m=1000.,note='diagnostic path allows decreasing removal progress'))
    old=cp['material_states'];uold=np.array(cp['displacement_m']);lam=lam0
    samples=[];events=[];call=0;processes=[];started=time.monotonic()
    binary=evidence/('target/release/examples/material_diagnostic'+('.exe' if os.name=='nt' else ''))
    def evaluate(du,states,label):
        nonlocal call
        req=[]
        for ei,qi,pi in pointmap:
            de=-rb[ei][qi]@du
            if np.abs(de).sum()>.001:raise ValueError('fixed material partition bound exceeded')
            req.append(dict(label=label,point=pi,parameters=parameters[mesh['material_ids'][ei]],
                state=states[pi],strain=de.tolist(),partitions=1,widths=[]))
        ip=out/f'call-{call:04d}-request.json';op=out/f'call-{call:04d}-result.json';write(ip,req)
        args=[str(binary),str(ip),str(op)];t=time.monotonic()
        with (out/f'call-{call:04d}.log').open('x') as f:
            p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,timeout=30)
        processes.append(dict(call=call,command=args,exit_code=p.returncode,elapsed_seconds=time.monotonic()-t))
        p.check_returncode();result=read(op)['points'];call+=1
        if any('error' in x for x in result):raise ValueError('native material error: '+str([x for x in result if 'error' in x][:1]))
        if not all(x['ad_primal_exact'] and x['substep_replay_exact'] for x in result):raise ValueError('primal/diagnostic mismatch')
        state=list(states);internal=np.zeros(nd);k=np.zeros((nd,nd))
        for (ei,qi,pi),p in zip(pointmap,result,strict=True):
            state[pi]=p['updated_state'];b=rb[ei][qi];w=q[ei]['weight'][qi]
            internal-=b.T@np.array(state[pi]['effective_stress_compression_positive_kpa'])*w
            k+=b.T@np.array(p['implicit_tangent'])@b*w
        return internal,k,state,result
    error=None
    try:
        f,k,s,_=evaluate(np.zeros(nd),old,'seed-zero')
        check=np.linalg.norm(R.T@load0-f)/max(np.linalg.norm(load0),1.)
        write(out/'seed-check.json',dict(relative_force_residual=float(check),
            state_exact=s==old,maximum_native_internal_difference_kn=float(np.max(np.abs(R.T@np.array(seed['internal']).ravel()-f)))))
        predictor_du=rate*(step_size/(plane@rate));predictor_dlambda=step_size/(plane@rate)
        for step in range(count):
            du=predictor_du.copy();trial_lam=lam+predictor_dlambda
            for it in range(21):
                f,k,states,points=evaluate(du,old,f'step-{step}:iteration-{it}')
                target=load0-(trial_lam-lam0)*holding
                residual=R.T@target-f;scale=max(np.linalg.norm(target),1.)
                rel=np.linalg.norm(residual)/scale;constraint=step_size-plane@du
                sig=np.linalg.svd(k,compute_uv=False)
                event=dict(step=step,iteration=it,material_call=call-1,progress=float(trial_lam),relative_force_residual=float(rel),
                    constraint_residual_m=float(constraint),sigma_min_kn_m=float(sig[-1]),
                    symmetric_eigen_min_kn_m=float(np.linalg.eigvalsh((k+k.T)/2)[0]),
                    displacement_increment_norm_m=float(np.linalg.norm(du)))
                events.append(event)
                if it>=1 and rel<1e-8 and abs(constraint)<1e-12:
                    uold+=(R@du).reshape(uold.shape)
                    samples.append(dict(**event,states=states,displacement_m=uold.tolist(),
                        previous_progress=float(lam),progress_increment=float(trial_lam-lam),
                        removed_material_unchanged=all(states[i]==cp['material_states'][i] for ei in range(len(q)) if not active[ei] for i in range(int(offsets[ei]),int(offsets[ei+1])))))
                    predictor_du=du.copy();predictor_dlambda=trial_lam-lam
                    old=states;lam=trial_lam
                    write(out/'samples.json',samples)
                    break
                if it==20:raise ValueError('diagnostic corrector iteration limit')
                aug=np.zeros((nd+1,nd+1));aug[:nd,:nd]=k;aug[:nd,nd]=-L;aug[nd,:nd]=1000.*plane
                rhs=np.r_[residual,1000.*constraint];dx=np.linalg.solve(aug,rhs)
                linear=np.linalg.norm(aug@dx-rhs)/max(np.linalg.norm(rhs),1e-30)
                write(out/f'linear-{call-1:04d}.json',dict(matrix_columns=aug.T.tolist(),rhs=rhs.tolist(),solution=dx.tolist(),true_relative_residual=float(linear)))
                event['augmented_true_relative_linear_residual']=float(linear)
                event['augmented_condition_number']=float(np.linalg.cond(aug))
                if not linear<1e-10:raise ValueError('diagnostic augmented linear gate failed')
                du+=dx[:nd];trial_lam+=dx[nd]
            else:raise ValueError('corrector did not complete')
    except (ValueError,subprocess.SubprocessError) as ex:error=str(ex)
    write(out/'events.json',events);write(out/'processes.json',processes)
    report=dict(passed=False,history_converged=False,production_adopted=False,independent_reference=False,
        seed_progress=lam0,equilibrium_samples=len(samples),requested_samples=count,stop_reason=error,
        progress_values=[x['progress'] for x in samples],
        progress_increments=[x['progress_increment'] for x in samples],
        maximum_progress=max([lam0]+[x['progress'] for x in samples]),
        elapsed_seconds=time.monotonic()-started,
        limitations=['displacement-controlled diagnostic, not the prescribed monotone release solve',
            'uses native constitutive map and saved history; not an independent solution',
            'not a sealed Run or production accepted increments',
            'no failure/stability verdict; no full-release acceptance'])
    write(out/'report.json',report);return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--step-size',type=float,default=1e-5);p.add_argument('--count',type=int,default=12);a=p.parse_args()
    assert 0<a.step_size<=1e-4 and 0<a.count<=64
    print(json.dumps(run(a.evidence.resolve(),a.output.resolve(),a.step_size,a.count),indent=2))
