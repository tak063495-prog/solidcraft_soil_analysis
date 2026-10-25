"""Independent full-3D linear system for an elastic MC geostatic slope.

Uses input seed stress and loads, NumPy/SVD, and different quadrature. Never
uses Rust stresses/displacements to construct the reference solution.
"""
import argparse
import hashlib
import json
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from verify_high_order import physical


def verify(run,output):
    if output.exists():raise ValueError('reference output exists')
    mp=run/'ground/mesh.json';gp=run/'ground/ground.json';rp=run/'stage-results.json'
    mesh=json.loads(mp.read_text());ground=ref.normalize_ground(json.loads(gp.read_text()));result=json.loads(rp.read_text());case=result['case']
    case={**case,'boundary':case.get('boundary') or {'constraints':[],'nodal_loads':[],'surface_loads':[]}}
    assert result['complete'] and result['schema']=='soil-fem-stage-result/4'
    assert ground['initialization']['geostatic']['seed']['kind']=='body_states'
    for c in case.get('boundary',{}).get('constraints',[]):assert c['value_m']==0 and not np.asarray(c['gradient']).any()
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size
    if ndof>6500:raise ValueError('independent dense reference limited to 6500 physical DOFs')
    material={m['id']:m for m in ground['materials']};regions={r['id']:r for r in mesh['regions']};seeds=ground['initialization']['geostatic']['seed']['states']
    operators=[]
    for ei in range(len(mesh['tetrahedra'])):
        p=material[mesh['material_ids'][ei]]['mohr_coulomb'];d=ref.elastic_tensor(p['young_kpa'],p['poisson']);q=ref.volume(mesh,ei,independent=True)
        dofs=(q['nodes'][:,None]*3+np.arange(3)).reshape(-1)
        ke=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        seed=np.array(seeds[regions[mesh['region_ids'][ei]]['body']]['effective_stress_compression_positive_kpa'])
        total=np.broadcast_to(seed,(len(q['weight']),6)).copy();h=ground['initialization']['water_level_m']
        if h is not None:total[:,:3]+=np.maximum(h-q['position'][:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
        fe=-np.einsum('qai,qa,q->i',q['b'],total,q['weight'],optimize=True)
        operators.append((dofs,ke,fe,len(ref.volume(mesh,ei)['weight'])))
    reports=[];references=[]
    for frame in result['frames']:
        state=frame['state'];assert all(s['equivalent_plastic_strain']==0 for s in state['material_states']),'this linear reference is valid only for elastic MC branches'
        active=np.asarray(state['active_elements']);used=np.unique(np.asarray(mesh['tetrahedra10'] if mesh.get('tetrahedra10') else mesh['tetrahedra'])[active])
        k=np.zeros((ndof,ndof));fseed=np.zeros(ndof);point_ids=[];offset=0
        for ei,(dofs,ke,fe,nq) in enumerate(operators):
            if active[ei]:
                k[np.ix_(dofs,dofs)]+=ke;fseed[dofs]+=fe;point_ids.extend(range(offset,offset+nq))
            offset+=nq
        basis,error=ref.staged_constraint_space(mesh,case,{**state,'displacement_m':np.zeros_like(xyz)});assert error==0
        count=sum(b.shape[1] for b in basis.values());t=np.zeros((ndof,count));start=0
        for n,b in basis.items():t[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
        f,_=physical(mesh,ground,active,state['surface_pressure_kpa'],case.get('pressure_on_original_surface_only',False))
        f+=state['boundary_factor']*ref.conditions_load(mesh,case.get('boundary',{}),active)
        rhs=t.T@(f.reshape(-1)-fseed);u=t@np.linalg.solve(t.T@k@t,rhs);references.append(u.reshape(-1,3))
        actual=np.asarray(state['displacement_m']);scale=max(np.linalg.norm(u.reshape(-1,3)[used]),1e-12);du=float(np.linalg.norm(actual[used]-u.reshape(-1,3)[used])/scale)
        assert du<1e-7,du
        reaction=(k@u+fseed-f.reshape(-1)).reshape(-1,3);assert np.allclose(reaction.sum(axis=0),frame['support_reaction_kn'],atol=1e-5,rtol=1e-9)
        stress=[];void=[]
        for ei in np.flatnonzero(active):
            q=ref.volume(mesh,ei);p=material[mesh['material_ids'][ei]]['mohr_coulomb'];d=ref.elastic_tensor(p['young_kpa'],p['poisson'])
            s0=seeds[regions[mesh['region_ids'][ei]]['body']];eps=np.einsum('qai,i->qa',q['b'],u.reshape(-1,3)[q['nodes']].reshape(-1));stress.extend(np.asarray(s0['effective_stress_compression_positive_kpa'])-eps@d.T)
            void.extend((1+s0['void_ratio'])*np.exp(eps[:,:3].sum(axis=1))-1)
        actual_points=[state['material_states'][i] for i in point_ids]
        actual_stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in actual_points]);stress=np.asarray(stress)
        ds=float(np.linalg.norm(actual_stress-stress)/np.linalg.norm(stress));de=float(np.max(np.abs(np.asarray([s['void_ratio'] for s in actual_points])-void)))
        assert ds<1e-7 and de<1e-8,(ds,de)
        reports.append({'frame':frame['name'],'free_dofs':count,'active_elements':int(active.sum()),'relative_3d_displacement_error':du,'relative_all_point_stress_error':ds,'maximum_void_ratio_error':de})
    reference=np.asarray(result['reference_displacement_m']);assert np.linalg.norm(reference-references[0])/np.linalg.norm(references[0])<1e-7
    report={'schema':'soil-fem-independent-geostatic-slope/1','reference':'full 3D NumPy dense solve, nodal SVD, Gauss-Legendre/Duffy volume and face quadrature; input seed/loads only',
        'frames':reports,'solver_sha256':result['provenance']['solver_binary_sha256'],
        'artifacts_sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mp,gp,rp,Path(__file__),Path(ref.__file__))}}
    output.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);a=p.parse_args();verify(a.run,a.output)
