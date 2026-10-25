"""Independent 3D linear-history reference for BIRTH/MATSET/death/replacement.

All stresses and displacements originate in inputs. Each new material keeps its
own strain origin; changing E does not recompute old stress from cumulative u.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from verify_high_order import physical


def verify(folder,output):
    if output.exists():raise ValueError('reference output exists')
    mp=folder/'ground/mesh.json';gp=folder/'ground/ground.json';rp=folder/'stage-results.json'
    mesh=json.loads(mp.read_text());ground=ref.normalize_ground(json.loads(gp.read_text()));result=json.loads(rp.read_text());case=result['case']
    assert result['complete'] and result['schema'] in ('soil-fem-stage-result/5','soil-fem-stage-result/6')
    assert ground['initialization']['water_level_m'] is None
    cfg=ground['initialization']['geostatic'];assert cfg['seed']['kind']=='body_states'
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;assert ndof<=6500
    case={**case,'boundary':case.get('boundary') or {'constraints':[],'surface_loads':[],'nodal_loads':[]}}
    for c in case['boundary']['constraints']:assert c['value_m']==0 and not np.asarray(c['gradient']).any()
    material={m['id']:m for m in ground['materials']};assert all('mohr_coulomb' in m for m in material.values())
    regions={r['id']:r['body'] for r in mesh['regions']};body=[regions[r] for r in mesh['region_ids']]
    independent=[ref.volume(mesh,i,independent=True) for i in range(len(body))]
    actual_rule=[ref.volume(mesh,i) for i in range(len(body))]
    rules=[independent,actual_rule];stress=[[],[]];void=[[],[]]
    seeds=cfg['seed']['states']
    for ri,rule in enumerate(rules):
        for ei,q in enumerate(rule):
            s=seeds[body[ei]];stress[ri].append(np.tile(s['effective_stress_compression_positive_kpa'],(len(q['weight']),1)).astype(float))
            void[ri].append(np.full(len(q['weight']),s['void_ratio']))
    ids=list(mesh['material_ids']);active=np.asarray([b not in case.get('initially_inactive_bodies',[]) for b in body]);u=np.zeros_like(xyz)
    pressure=ground['initialization']['surface_surcharge_kpa'];factor=cfg['initial_boundary_factor'];reports=[];display_reference=None;explicit_display=False
    for index,frame in enumerate(result['frames']):
        if index:
            stage=case['stages'][index-1];kind=stage['kind'];selection=[i for i,b in enumerate(body) if b in stage.get('bodies',[])]
            if kind=='pressure':pressure=stage['target_pressure_kpa']
            elif kind=='boundary':factor=stage['target_factor']
            elif kind=='excavate':active[selection]=False;pressure=stage['remaining_surface_pressure_kpa']
            elif kind=='reset_display':display_reference=u.copy();explicit_display=True
            elif kind in ('birth','material_set','reset_state'):
                if kind=='birth':
                    existing=np.unique(ref.cells(mesh)[active]);new=np.setdiff1d(np.unique(ref.cells(mesh)[selection]),existing)
                    if stage.get('placement')=='cad_new_nodes':
                        u[new]=0.
                    if explicit_display:display_reference[new]=u[new]
                    assert not active[selection].any();active[selection]=True;replacement=stage['states']
                elif kind=='material_set':
                    for ei in selection:ids[ei]=stage['material_id']
                    replacement=stage['history'].get('states')
                else:replacement=stage['states']
                if replacement is not None:
                    for ri,rule in enumerate(rules):
                        for ei in selection:
                            s=replacement[body[ei]];stress[ri][ei][:]=s['effective_stress_compression_positive_kpa'];void[ri][ei][:]=s['void_ratio']
            else:raise ValueError('unsupported reference operation')
        k=np.zeros((ndof,ndof));old_internal=np.zeros(ndof)
        for ei in np.flatnonzero(active):
            q=independent[ei];p=material[ids[ei]]['mohr_coulomb'];d=ref.elastic_tensor(p['young_kpa'],p['poisson']);dofs=(q['nodes'][:,None]*3+np.arange(3)).reshape(-1)
            k[np.ix_(dofs,dofs)]+=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
            old_internal[dofs]-=np.einsum('qai,qa,q->i',q['b'],stress[0][ei],q['weight'],optimize=True)
        support_state={'active_elements':active.tolist(),'boundary_factor':factor,'displacement_m':np.zeros_like(xyz)}
        basis,error=ref.staged_constraint_space(mesh,case,support_state);assert error==0
        count=sum(b.shape[1] for b in basis.values());t=np.zeros((ndof,count));start=0
        for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
        f,_=physical(mesh,ground,active,pressure,case.get('pressure_on_original_surface_only',False),ids)
        if factor:f+=factor*ref.conditions_load(mesh,case['boundary'],active)
        noop=index and case['stages'][index-1]['kind']=='reset_display'
        du=np.zeros_like(u) if noop else (t@np.linalg.solve(t.T@k@t,t.T@(f.reshape(-1)-old_internal))).reshape(-1,3)
        for ei in ([] if noop else np.flatnonzero(active)):
            p=material[ids[ei]]['mohr_coulomb'];d=ref.elastic_tensor(p['young_kpa'],p['poisson'])
            for ri,rule in enumerate(rules):
                q=rule[ei];eps=np.einsum('qai,i->qa',q['b'],du[q['nodes']].reshape(-1));stress[ri][ei]-=eps@d.T
                void[ri][ei]=(1+void[ri][ei])*np.exp(eps[:,:3].sum(axis=1))-1
        u+=du;state=frame['state'];assert state['active_elements']==active.tolist() and state['element_material_ids']==ids
        if index==0:display_reference=u.copy()
        native_reference=np.asarray(state.get('display_reference_displacement_m') or result['reference_displacement_m'])
        display_absolute_error=float(np.linalg.norm((np.asarray(state['displacement_m'])-native_reference)-(u-display_reference)))
        display_norm=float(np.linalg.norm(u-display_reference));display_error=display_absolute_error/display_norm if display_norm>1e-9 else None
        assert display_absolute_error<=1e-12+1e-7*display_norm,(frame['name'],display_absolute_error,display_norm)
        assert all(s['equivalent_plastic_strain']==0 for s in state['material_states']),'linear reference requires an elastic branch'
        actual=np.asarray(state['displacement_m']);displacement_error=float(np.linalg.norm(actual-u)/max(np.linalg.norm(u),1e-12))
        predicted=np.concatenate(stress[1]);actual_stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in state['material_states']])
        stress_error=float(np.linalg.norm(actual_stress-predicted)/max(np.linalg.norm(predicted),1e-12))
        void_error=float(np.max(np.abs(np.asarray([s['void_ratio'] for s in state['material_states']])-np.concatenate(void[1]))))
        reaction=(k@du.reshape(-1)+old_internal-f.reshape(-1)).reshape(-1,3)
        assert np.allclose(reaction.sum(axis=0),frame['support_reaction_kn'],atol=1e-5,rtol=1e-9)
        assert displacement_error<1e-7 and stress_error<1e-7 and void_error<1e-8,(frame['name'],displacement_error,stress_error,void_error)
        reports.append({'stage':frame['name'],'active_elements':int(active.sum()),'free_dofs':count,'relative_3d_displacement_error':displacement_error,'relative_display_displacement_error':display_error,'absolute_display_displacement_error_m':display_absolute_error,'relative_all_point_stress_error':stress_error,'maximum_void_ratio_error':void_error})
    report={'schema':'soil-fem-independent-construction/1','reference':'NumPy 3D elastic stress history, independent GL/Duffy integration, nodal SVD; no Rust states in reference construction',
        'frames':reports,'solver_sha256':result['provenance']['solver_binary_sha256'],
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mp,gp,rp,Path(__file__),Path(ref.__file__))}}
    with output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);a=p.parse_args();verify(a.run,a.output)
