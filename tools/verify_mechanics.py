"""Independent high-order elasticity: different quadrature, SVD constraints,
NumPy dense global solve. Also validates raw stress/strain and reactions."""
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


def verify(mesh_path,case_path,result_path,report_path,reference_order=5):
    if report_path.exists():raise ValueError('verification report exists')
    mesh=json.loads(mesh_path.read_text());case=json.loads(case_path.read_text());result=json.loads(result_path.read_text())
    assert result['schema']=='soil-fem-mechanics-result/1'
    assert result['solver_binary_sha256']==hashlib.sha256((Path(__file__).resolve().parents[1]/'target/release/soil-fem-3d.exe').read_bytes()).hexdigest()
    ref.normalize_conditions(case['boundary']);assert case==result['case']
    xyz=np.array(mesh['nodes']);ndof=xyz.size
    if ndof>6500:raise ValueError('independent dense mechanics reference is limited to 6500 physical DOFs')
    materials={m['id']:m for m in case['materials']};k=np.zeros((ndof,ndof));f=np.zeros(ndof)
    for ei in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,ei,independent=True,independent_order=reference_order);m=materials[mesh['material_ids'][ei]]
        d=ref.elastic_tensor(m['young_kpa'],m['poisson']);dofs=(3*q['nodes'][:,None]+np.arange(3)).reshape(-1)
        ke=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        k[np.ix_(dofs,dofs)]+=ke
        body=(q['weight']@q['shape'])[:,None]*m['body_force_kn_m3'];f[dofs]+=body.reshape(-1)
    for load in case['boundary'].get('nodal_loads',[]):
        nodes=np.array(mesh['node_sets'][load['node_set']]);f.reshape(-1,3)[nodes]+=np.asarray(load['total_force_kn'])/len(nodes)
    owners=ref.topology(mesh)
    for load in case['boundary'].get('surface_loads',[]):
        for face in mesh['face_sets'][load['face_set']]:
            adjacent=owners[tuple(sorted(face))];assert len(adjacent)==1
            q=ref.face(mesh,adjacent[0],face);dx=q['position']-load.get('reference_point_m',[0,0,0]);area=q['area_vector']
            if load['kind']=='pressure':
                pressure=load['pressure_kpa']+dx@np.asarray(load.get('gradient_kpa_per_m',[0,0,0]));force=-pressure[:,None]*area
            else:
                traction=load['traction_kpa']+dx@np.asarray(load.get('gradient_kpa_per_m',np.zeros((3,3)))).T
                force=traction*np.linalg.norm(area,axis=1)[:,None]
            f.reshape(-1,3)[q['nodes']]+=q['shape'].T@force
    eq=[[] for _ in xyz];values=[[] for _ in xyz]
    for c in case['boundary']['constraints']:
        for node in mesh['node_sets'][c['node_set']]:
            eq[node].append(c['direction'])
            values[node].append(c['value_m']+np.asarray(c.get('gradient',[0,0,0]))@(xyz[node]-c.get('reference_point_m',[0,0,0])))
    nulls=[];offset=np.zeros(ndof);count=0
    for node in range(len(xyz)):
        if eq[node]:
            a=np.array(eq[node]);b=np.array(values[node]);u0=np.linalg.lstsq(a,b,rcond=1e-12)[0]
            assert np.linalg.norm(a@u0-b)<1e-10
            _,s,vh=np.linalg.svd(a,full_matrices=True);rank=int((s>1e-10*max(s.max(),1)).sum());basis=vh[rank:].T
            offset[3*node:3*node+3]=u0
        else:basis=np.eye(3)
        nulls.append((count,basis));count+=basis.shape[1]
    t=np.zeros((ndof,count))
    for node,(start,basis) in enumerate(nulls):t[3*node:3*node+3,start:start+basis.shape[1]]=basis
    reduced=t.T@k@t;rhs=t.T@(f-k@offset)
    u=offset+(t@np.linalg.solve(reduced,rhs) if count else 0)
    actual=np.array(result['displacement_m']).reshape(-1);scale=max(np.linalg.norm(u),1e-12)
    difference=np.linalg.norm(actual-u)/scale
    assert difference<1e-6,difference
    assert np.allclose(f.reshape(-1,3),result['applied_nodal_load_kn'],rtol=1e-10,atol=1e-9)
    reaction=k@u-f;actual_reaction=np.array(result['reaction_kn']).reshape(-1)
    reaction_error=np.linalg.norm(reaction-actual_reaction)/max(np.linalg.norm(reaction),1.)
    assert reaction_error<1e-7,reaction_error
    stress=[];strain=[]
    for ei in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,ei);m=materials[mesh['material_ids'][ei]];d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
        eps=np.einsum('qai,i->qa',q['b'],u.reshape(-1,3)[q['nodes']].reshape(-1));strain.extend(eps);stress.extend(eps@d.T)
    stress=np.array(stress);strain=np.array(strain)
    stress_error=np.linalg.norm(np.array(result['stress_tension_positive_kpa'])-stress)/max(np.linalg.norm(stress),1.)
    strain_error=np.linalg.norm(np.array(result['strain_engineering'])-strain)/max(np.linalg.norm(strain),1e-12)
    assert stress_error<1e-6 and strain_error<1e-6,(stress_error,strain_error)
    report={'schema':'soil-fem-independent-mechanics/1','reference':f'independent NumPy basis + Gauss-Legendre/Duffy {reference_order**3}-point volume integration + nodal SVD + dense solve',
            'free_dofs':count,'relative_displacement_difference':float(difference),'relative_reaction_difference':float(reaction_error),
            'relative_stress_difference':float(stress_error),'relative_strain_difference':float(strain_error),
            'force_balance_error_kn':(f+reaction).reshape(-1,3).sum(axis=0).tolist(),
            'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (mesh_path,case_path,result_path)}}
    report_path.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
    return u,stress


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('case',type=Path);p.add_argument('result',type=Path);p.add_argument('report',type=Path)
    p.add_argument('--reference-order',type=int,choices=range(5,10),default=5)
    a=p.parse_args();verify(a.mesh,a.case,a.result,a.report,a.reference_order)
