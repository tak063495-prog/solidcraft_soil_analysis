"""Independent elastic soil + tension-only anchors by active-set enumeration.

Uses independent host inversion, GL125 soil integration, SVD constraints and
input-only installation/jacking/lock-off. Rust Newton/active states/history are
used only for comparison. Small reference problems only (at most 16 segments).
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from anchor_reference import System
from verify_high_order import physical


def direct(k,seed,load,t,system,state):
    locked=[];fixed=np.zeros(len(seed));k0=t.T@k@t
    for d,s in zip(system.data,state['anchors']):
        if not s['active']:continue
        for j,(c,stiffness,reference) in enumerate(zip(d['c'],d['k'],s['reference_extensions_m'])):
            if j==0 and s['force_control_kn'] is not None:fixed+=c*s['force_control_kn']
            else:locked.append((c,c@t,stiffness,reference))
    assert len(locked)<=16
    base=t.T@(load-seed-fixed);checked=0
    for active in itertools.product((False,True),repeat=len(locked)):
        checked+=1;a=k0.copy();rhs=base.copy()
        for enabled,(_,c,kj,reference) in zip(active,locked):
            if enabled:a+=kj*np.outer(c,c);rhs+=kj*reference*c
        q=np.linalg.solve(a,rhs);u=t@q
        raw=np.array([kj*(c@u-reference) for c,_,kj,reference in locked])
        if any((v<-1e-7 if enabled else v>1e-7) for enabled,v in zip(active,raw)):continue
        updated,internal,_,reports=system.evaluate(state,u.reshape(-1,3))
        residual=t.T@(k@u+seed+internal.reshape(-1)-load)
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(load),1.));assert relative<1e-9
        return u.reshape(-1,3),updated,internal,reports,relative,checked
    raise AssertionError('no admissible tension/slack active set')


def verify(run,output):
    if output.exists():raise ValueError('reference output exists')
    mp=run/'ground/mesh.json';gp=run/'ground/ground.json';rp=run/'stage-results.json'
    mesh=json.loads(mp.read_text());ground=ref.normalize_ground(json.loads(gp.read_text()));native=json.loads(rp.read_text());case=native['case']
    assert native['complete'] and native['schema']=='soil-fem-stage-result/8' and not case.get('beam_system')
    assert all(s['kind'] in ('install_anchors','remove_anchors','jack_anchors','boundary') for s in case['stages'])
    assert ground['initialization']['method']=='equilibrate' and not case.get('initially_inactive_bodies')
    geo=ground['initialization']['geostatic'];assert geo['seed']['kind']=='body_states' and geo['initial_boundary_factor']==0
    for c in case['boundary']['constraints']:assert c['value_m']==0 and not np.asarray(c['gradient']).any()
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;active=np.ones(len(mesh['tetrahedra']),bool);k=np.zeros((ndof,ndof));seed=np.zeros(ndof);operators=[]
    assert ndof<=6500
    materials={m['id']:m for m in ground['materials']};regions={r['id']:r['body'] for r in mesh['regions']}
    for ei in range(len(active)):
        m=materials[mesh['material_ids'][ei]]['mohr_coulomb'];d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
        s=geo['seed']['states'][regions[mesh['region_ids'][ei]]];s0=np.asarray(s['effective_stress_compression_positive_kpa'])
        q=ref.volume(mesh,ei,independent=True);dofs=(3*q['nodes'][:,None]+np.arange(3)).reshape(-1)
        k[np.ix_(dofs,dofs)]+=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        total=np.broadcast_to(s0,(len(q['weight']),6)).copy();h=ground['initialization']['water_level_m']
        if h is not None:total[:,:3]+=np.maximum(h-q['position'][:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
        seed[dofs]-=np.einsum('qai,qa,q->i',q['b'],total,q['weight'],optimize=True)
        operators.append((ref.volume(mesh,ei),d,s0,s['void_ratio']))
    basis,error=ref.staged_constraint_space(mesh,case,{'active_elements':active,'boundary_factor':0,'displacement_m':np.zeros_like(xyz)})
    assert error==0;count=sum(b.shape[1] for b in basis.values());t=np.zeros((ndof,count));start=0
    for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
    system=System(mesh,case['anchor_system']);state=system.initial_state();u=np.zeros_like(xyz);factor=0.;frames=[]
    for fi,frame in enumerate(native['frames']):
        selected=[]
        if fi:
            stage=case['stages'][fi-1]
            if stage['kind']=='boundary':factor=stage['target_factor']
            elif stage['kind']=='jack_anchors':
                for i,d in enumerate(system.data):
                    if d['definition']['id'] in stage['anchors']:selected.append(i);state['anchors'][i]['force_control_kn']=stage['target_force_kn']
            else:state=system.change(state,u,stage)
        f,_=physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False))
        f+=factor*ref.conditions_load(mesh,case['boundary'],active)
        u,state,ai,ar,residual,enumerated=direct(k,seed,f.reshape(-1),t,system,state)
        for i in selected:
            s=state['anchors'][i];d=system.data[i];force=s['force_control_kn']
            extension=float(d['c'][0]@u.reshape(-1))
            s['reference_extensions_m'][0]=(max(s['reference_extensions_m'][0],extension) if force==0 else extension-force/d['k'][0])
            s['last_lock_off_force_kn']=force;s['force_control_kn']=None
        state,ai,_,ar=system.evaluate(state,u)
        actual=frame['state'];assert [s['active'] for s in actual['anchors']['anchors']]==[s['active'] for s in state['anchors']]
        du=float(np.linalg.norm(np.asarray(actual['displacement_m'])-u)/max(np.linalg.norm(u),1e-12));assert du<1e-7,du
        force_error=0.;reference_error=0.
        for a,b in zip(actual['anchors']['anchors'],state['anchors']):
            force_error=max(force_error,float(np.max(np.abs(np.asarray(a['segment_forces_kn'])-b['segment_forces_kn']))))
            reference_error=max(reference_error,float(np.max(np.abs(np.asarray(a['reference_extensions_m'])-b['reference_extensions_m']))))
        assert force_error<1e-5 and reference_error<1e-9,(force_error,reference_error)
        stress=[];void=[]
        for qp,d,s0,e0 in operators:
            strain=np.einsum('qai,i->qa',qp['b'],u[qp['nodes']].reshape(-1));stress.extend(s0-strain@d.T);void.extend((1+e0)*np.exp(strain[:,:3].sum(axis=1))-1)
        points=actual['material_states'];assert all(p['equivalent_plastic_strain']==0 for p in points),'linear reference entered MC plasticity'
        stress=np.asarray(stress);ds=float(np.linalg.norm(np.asarray([p['effective_stress_compression_positive_kpa'] for p in points])-stress)/np.linalg.norm(stress))
        de=float(np.max(np.abs(np.asarray([p['void_ratio'] for p in points])-void)));assert ds<1e-7 and de<1e-8
        balance=(k@u.reshape(-1)+seed-f.reshape(-1)).reshape(-1,3)+ai;support=np.zeros_like(xyz)
        for n,b in basis.items():support[n]=balance[n]-b@(b.T@balance[n])
        assert np.allclose(support.sum(axis=0),frame['support_reaction_kn'],atol=1e-5,rtol=1e-9)
        moment=np.cross(xyz,f+support).sum(axis=0);assert np.linalg.norm(moment)<1e-6
        frames.append({'frame':frame['name'],'soil_free_dofs':count,'enumerated_active_sets':enumerated,'relative_3d_displacement_error':du,
            'relative_independent_equilibrium_residual':residual,'maximum_tension_error_kn':force_error,'maximum_installation_lock_reference_error_m':reference_error,
            'relative_all_point_stress_error':ds,'maximum_void_ratio_error':de,'global_moment_balance_kn_m':moment.tolist()})
    report={'schema':'soil-fem-independent-anchor-coupling/1','reference':'input-only GL125 elastic soil + independently inverted host interpolation; exhaustive tension/slack active-set enumeration, SVD, NumPy solve; installation and force-control/lock-off references from independent displacements',
        'frames':frames,'solver_sha256':native['provenance']['solver_binary_sha256'],
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mp,gp,rp,Path(__file__),Path(ref.__file__),Path(__file__).with_name('anchor_reference.py'))}}
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);a=p.parse_args();verify(a.run.resolve(),a.output.resolve())
