"""Input-only 3D elastic soil + frictional JOINT reference, independent Newton.

Uses GL125 soil integration, SVD displacement constraints, and a separately
implemented subincremental Coulomb projection. Native accepted load fractions
define the compared load path; no native displacement, traction, tangent or
material history is used to generate a reference solution.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from joint_reference import System
from verify_high_order import physical
from verify_checkpoint import invoke


def contact(system,old_jump,old_traction,target,count):
    """Independent two-tangent Coulomb disk projection along a straight jump."""
    if target.ndim==2:target=target[:,None,:]
    n=system.normal[:,None,:];kn=system.kn[:,None];kt=system.kt[:,None];mu=system.mu[:,None]
    delta=(target-old_jump[:,None,:])/count;dt=delta-(delta*n).sum(axis=2)[:,:,None]*n
    tau=old_traction[:,None,:]-(old_traction*system.normal).sum(axis=1)[:,None,None]*n
    tau=np.broadcast_to(tau,target.shape).copy();slip=np.zeros_like(tau);cumulative=np.zeros(target.shape[:2])
    gn=(old_jump*system.normal).sum(axis=1)[:,None];dg=(delta*n).sum(axis=2)
    for step in range(count):
        next_g=gn+dg;was=gn<=0;closed=next_g<=0
        fraction=np.ones_like(next_g)
        opening=was & ~closed;closing=~was & closed
        fraction[~was & ~closed]=0
        np.divide(-gn,dg,out=fraction,where=opening & (dg!=0))
        np.divide(next_g,dg,out=fraction,where=closing & (dg!=0))
        tau=np.where((~was)[:,:,None],0,tau)
        trial=tau+kt[:,:,None]*dt*fraction[:,:,None]
        pressure=kn*np.maximum(-next_g,0);limit=mu*pressure
        norm=np.linalg.norm(trial,axis=2);scale=np.minimum(1,np.divide(limit,norm,out=np.ones_like(norm),where=norm>0))
        projected=trial*scale[:,:,None];flow=(trial-projected)/kt[:,:,None]
        slip+=flow;cumulative+=np.linalg.norm(flow,axis=2)
        tau=np.where(closed[:,:,None],projected,0);gn=next_g
    traction=tau+kn[:,:,None]*np.minimum(gn,0)[:,:,None]*n
    return traction,slip,cumulative


def reference(mesh,ground,case,native,count,neutral_predictor=False,diagnostics=None):
    xyz=np.asarray(mesh['nodes']);nd=xyz.size;active=np.ones(len(mesh['tetrahedra']),bool)
    assert nd<=3000 and not case.get('beam_system') and not case.get('anchor_system')
    assert all(s['kind']=='boundary' for s in case['stages'])
    system=System(mesh,case['joint_system']);c=system.c;w=system.weight
    k=np.zeros((nd,nd));seed=np.zeros(nd);operators=[]
    materials={m['id']:m for m in ground['materials']};regions={r['id']:r['body'] for r in mesh['regions']}
    geo=ground['initialization']['geostatic'];assert geo['initial_boundary_factor']==0
    for ei in range(len(active)):
        m=materials[mesh['material_ids'][ei]]['mohr_coulomb'];d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
        state=geo['seed']['states'][regions[mesh['region_ids'][ei]]];s0=np.asarray(state['effective_stress_compression_positive_kpa'])
        qp=ref.volume(mesh,ei,independent=True);ids=(3*qp['nodes'][:,None]+np.arange(3)).reshape(-1)
        k[np.ix_(ids,ids)]+=np.einsum('qai,ab,qbj,q->ij',qp['b'],d,qp['b'],qp['weight'],optimize=True)
        seed[ids]-=np.einsum('qai,a,q->i',qp['b'],s0,qp['weight'])
        operators.append((ref.volume(mesh,ei),d,s0))
    basis,_=ref.staged_constraint_space(mesh,case,{'active_elements':active,'boundary_factor':0,'displacement_m':np.zeros_like(xyz)})
    nf=sum(b.shape[1] for b in basis.values());t=np.zeros((nd,nf));start=0
    for node,b in basis.items():t[3*node:3*node+3,start:start+b.shape[1]]=b;start+=b.shape[1]
    kr=t.T@k@t;cr=np.einsum('qag,gi->qai',c,t)
    base,_=physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],True,joint_system=system)
    base=base.reshape(-1);holding=seed-base
    old_jump=np.zeros((len(w),3));old_q=np.zeros_like(old_jump);plastic=np.zeros_like(old_jump);acc=np.zeros(len(w));u=np.zeros(nd)
    def offset(factor):
        equations={}
        for condition in case['boundary']['constraints']:
            direction=np.asarray(condition['direction']);gradient=np.asarray(condition.get('gradient',[0.,0.,0.]))
            for node in mesh['node_sets'][condition['node_set']]:
                a,b=equations.setdefault(node,([],[]));a.append(direction);b.append(factor*(condition['value_m']+gradient@(xyz[node]-condition.get('reference_point_m',[0.,0.,0.]))))
        result=np.zeros_like(xyz)
        for node,(a,b) in equations.items():result[node]=np.linalg.lstsq(a,b,rcond=None)[0]
        return result.reshape(-1)
    maximum_residual=0.
    singular_steps=0;incompatible_steps=0;cold_starts=0
    def increment(load,factor):
        nonlocal u,old_jump,old_q,plastic,acc,maximum_residual,singular_steps,incompatible_steps,cold_starts
        fixed=offset(factor);free=t.T@(u-fixed);scale=max(np.linalg.norm(load),1.)
        if neutral_predictor and not np.any(old_jump) and not np.any(old_q):
            # Independent bilateral elastic starting solution, before the
            # actual unilateral Coulomb equations and support SVD are checked.
            # This is only an iterate: it contributes no committed force/history.
            nn=np.einsum('qa,qb->qab',system.normal,system.normal)
            elastic=system.kn[:,None,None]*nn+system.kt[:,None,None]*(np.eye(3)-nn)
            kj=k+np.einsum('qag,qab,qbh,q->gh',c,elastic,c,w,optimize=True)
            free=np.linalg.solve(t.T@kj@t,t.T@(load-seed-kj@fixed));cold_starts+=1
        for iteration in range(60):
            u=fixed+t@free;jump=np.einsum('qag,g->qa',c,u)
            q,slip,cum=contact(system,old_jump,old_q,jump,count);q=q[:,0];force=k@u+seed+np.einsum('qag,qa,q->g',c,q,w)
            residual=t.T@(force-load);relative=np.linalg.norm(residual)/scale
            if relative<1e-10:
                if neutral_predictor:
                    states=[]
                    for op,d,s0 in operators:
                        stresses=s0-np.einsum('ab,qbi,i->qa',d,op['b'],u[(3*op['nodes'][:,None]+np.arange(3)).reshape(-1)])
                        states.extend({'effective_stress_compression_positive_kpa':s.tolist()} for s in stresses)
                    support=system.rigid_support(case,{'active_elements':active,'displacement_m':u.reshape(-1,3),
                        'material_states':states,'joints':{'points':[{'traction_kpa':v.tolist()} for v in q]},'boundary_factor':factor})
                    assert support['unrestrained_modes']==0,'independent solution has a rigid mechanism'
                maximum_residual=max(maximum_residual,float(relative));old_jump=jump;old_q=q;plastic+=slip[:,0];acc+=cum[:,0];return
            h=2e-9;perturb=np.repeat(jump[:,None,:],7,axis=1)
            for a in range(3):perturb[:,1+2*a,a]+=h;perturb[:,2+2*a,a]-=h
            pq,_,_=contact(system,old_jump,old_q,perturb,count)
            d=np.stack([(pq[:,1+2*a]-pq[:,2+2*a])/(2*h) for a in range(3)],axis=2)
            jac=kr+np.einsum('qai,qab,qbj,q->ij',cr,d,cr,w,optimize=True)
            if neutral_predictor:
                correction,_,rank,_=np.linalg.lstsq(jac,-residual,rcond=1e-11)
                singular_steps+=rank<jac.shape[0]
                incompatible_steps+=np.linalg.norm(jac@correction+residual)/max(np.linalg.norm(residual),1e-30)>1e-8
            else:correction=np.linalg.solve(jac,-residual)
            accepted=False
            for power in range(16):
                candidate=free+correction/(2**power);uc=fixed+t@candidate;jc=np.einsum('qag,g->qa',c,uc)
                qc,_,_=contact(system,old_jump,old_q,jc,count)
                rc=t.T@(k@uc+seed+np.einsum('qag,qa,q->g',c,qc[:,0],w)-load)
                if np.linalg.norm(rc)<np.linalg.norm(residual):free=candidate;accepted=True;break
            if not accepted:raise AssertionError(('independent line search',factor,iteration,relative))
        raise AssertionError(('independent Newton iteration limit',factor))
    frames=[]
    for attempt in native['geostatic']['attempts']:
        if attempt['accepted']:increment(base+(1-attempt['progress_to'])*holding,0.)
    frames.append({'u':u.copy(),'traction':old_q.copy(),'plastic':plastic.copy(),'accumulated':acc.copy()})
    origin_factor=0.
    for stage in case['stages']:
        target=stage['target_factor']
        for attempt in native['attempts']:
            if attempt['stage']==stage['name'] and attempt['accepted']:
                factor=origin_factor+(target-origin_factor)*attempt['progress_to']
                load=base+factor*ref.conditions_load(mesh,case['boundary'],active).reshape(-1)
                increment(load,factor)
        origin_factor=target
        frames.append({'u':u.copy(),'traction':old_q.copy(),'plastic':plastic.copy(),'accumulated':acc.copy()})
    if diagnostics is not None:diagnostics.update({'singular_least_squares_steps':int(singular_steps),'incompatible_linear_predictors':int(incompatible_steps),'bilateral_elastic_starting_iterates':cold_starts})
    return frames,maximum_residual


def verify(run,output,steps,solver=None):
    mp=run/'ground/mesh.json';gp=run/'ground/ground.json';rp=run/'stage-results.json'
    mesh=json.loads(mp.read_text());ground=json.loads(gp.read_text());native=json.loads(rp.read_text());case=native['case']
    assert native['complete'] and native['schema']=='soil-fem-stage-result/9'
    assert all(p['equivalent_plastic_strain']==0 for f in native['frames'] for p in f['state']['material_states']),'elastic soil control required'
    coarse,_=reference(mesh,ground,case,native,steps);fine,residual=reference(mesh,ground,case,native,2*steps)
    rows=[]
    for f,a,b in zip(native['frames'],coarse,fine):
        actual=np.asarray(f['state']['displacement_m']).reshape(-1);scale=max(np.linalg.norm(b['u']),1e-12)
        error=np.linalg.norm(actual-b['u'])/scale;refinement=np.linalg.norm(a['u']-b['u'])/scale
        traction=np.asarray([s['traction_kpa'] for s in f['state']['joints']['points']]);te=np.linalg.norm(traction-b['traction'])/max(np.linalg.norm(b['traction']),1.)
        assert error<.005 and refinement<.001 and te<.01,(f['name'],error,refinement,te)
        rows.append({'frame':f['name'],'relative_3d_displacement_error':float(error),'reference_refinement_error':float(refinement),'relative_traction_error':float(te)})
    report={'schema':'soil-fem-independent-joint-coupling/1','reference':'GL125 elastic soil; SVD constraints; independent subincremental 3D Coulomb disk projection and Newton; input-only states, native accepted load fractions only',
        'reference_subincrements':[steps,2*steps],'frames':rows,'maximum_independent_force_residual':residual,'solver_sha256':native['provenance']['solver_binary_sha256'],
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mp,gp,rp,Path(__file__),Path(ref.__file__),Path(__file__).with_name('joint_reference.py'))}}
    if solver is not None:
        assert hashlib.sha256(solver.read_bytes()).hexdigest()==native['provenance']['solver_binary_sha256']
        refined=output.with_name(output.stem+'-refined-native');refined.mkdir(exist_ok=False);(refined/'ground').mkdir()
        rmp=refined/'ground/mesh.json';rgp=refined/'ground/ground.json';rcp=refined/'stages.json';rrp=refined/'stage-results.json'
        shutil.copyfile(mp,rmp);shutil.copyfile(gp,rgp);configuration=json.loads(json.dumps(case))
        maximum=configuration['joint_system']['maximum_increment_m'];configuration['joint_system']['maximum_increment_m']=maximum/10
        rcp.write_text(json.dumps(configuration,indent=2),encoding='utf-8')
        invoke([str(solver),'solve-stages',str(rmp),str(rgp),str(rcp),str(rrp)],refined/'solve.log')
        updated=json.loads(rrp.read_text());assert updated['complete']
        assert [(a['stage'],a['progress_to']) for a in updated['attempts'] if a['accepted']]==[(a['stage'],a['progress_to']) for a in native['attempts'] if a['accepted']], 'local-integration comparison requires the same global accepted increments'
        invoke([sys.executable,str(Path(__file__).with_name('export_stages.py')),str(rmp),str(rrp),str(refined/'results')],refined/'export.log')
        invoke([sys.executable,str(Path(__file__).with_name('verify_stage_package.py')),str(rmp),str(rgp),str(rrp),str(refined/'results'),'--solver',str(solver)],refined/'audit.log')
        refined_rows=[]
        for frame,b in zip(updated['frames'],fine):
            u=np.asarray(frame['state']['displacement_m']).reshape(-1);error=float(np.linalg.norm(u-b['u'])/max(np.linalg.norm(b['u']),1e-12))
            q=np.asarray([s['traction_kpa'] for s in frame['state']['joints']['points']]);te=float(np.linalg.norm(q-b['traction'])/max(np.linalg.norm(b['traction']),1.))
            assert error<1e-4 and te<.002,(frame['name'],error,te)
            refined_rows.append({'frame':frame['name'],'relative_3d_displacement_error':error,'relative_traction_error':te})
        before=max(r['relative_3d_displacement_error'] for r in rows);after=max(r['relative_3d_displacement_error'] for r in refined_rows)
        assert after<before/2,(before,after)
        report['native_contact_refinement']={'maximum_increment_m':[maximum,maximum/10],'same_global_accepted_path':True,
            'maximum_relative_3d_displacement_error_before_after':[before,after],'frames':refined_rows}
        report['sha256'].update({str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (rmp,rgp,rcp,rrp)})
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--steps',type=int,default=100);p.add_argument('--solver',type=Path)
    a=p.parse_args();assert a.steps>=20;verify(a.run.resolve(),a.output.resolve(),a.steps,a.solver.resolve() if a.solver else None)
