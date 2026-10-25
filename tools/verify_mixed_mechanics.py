"""Independent dense u/P1 calculation, raw-field comparison and pressure rank audit.

Uses NumPy basis, GL/Duffy quadrature and per-node SVD constraints, never the
native displacement, matrix, pressure or reaction to construct the solution.
The discrete pressure-rank check is evidence for this mesh only.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def constraint_map(mesh,conditions):
    xyz=np.asarray(mesh['nodes']);equations=[[] for _ in xyz];values=[[] for _ in xyz]
    for c in conditions['constraints']:
        for node in mesh['node_sets'][c['node_set']]:
            equations[node].append(c['direction'])
            values[node].append(c['value_m']+np.asarray(c.get('gradient',[0,0,0]))@(xyz[node]-c.get('reference_point_m',[0,0,0])))
    offset=np.zeros_like(xyz);bases=[];count=0
    for node in range(len(xyz)):
        if equations[node]:
            a=np.asarray(equations[node]);b=np.asarray(values[node]);offset[node]=np.linalg.lstsq(a,b,rcond=1e-12)[0]
            if np.linalg.norm(a@offset[node]-b)>1e-10:raise ValueError('inconsistent reference constraints')
            _,s,vh=np.linalg.svd(a,full_matrices=True);rank=int((s>1e-10*max(s.max(),1.)).sum());basis=vh[rank:].T
        else:basis=np.eye(3)
        bases.append((count,basis));count+=basis.shape[1]
    t=np.zeros((xyz.size,count))
    for node,(start,basis) in enumerate(bases):t[3*node:3*node+3,start:start+basis.shape[1]]=basis
    return t,offset.reshape(-1)


def solve_reference(mesh,case,order=5):
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;materials={m['id']:m for m in case['materials']}
    keys=sorted({(node,mat) for nodes,mat in zip(mesh['tetrahedra'],mesh['material_ids']) for node in nodes})
    index={key:i for i,key in enumerate(keys)};npres=len(keys);total=ndof+npres
    if total>2600:raise ValueError('independent mixed dense reference is limited to 2600 physical plus pressure DOFs')
    k=np.zeros((total,total));h1=np.zeros((ndof,ndof));mass=np.zeros((npres,npres))
    force=ref.conditions_load(mesh,case['boundary'],np.ones(len(mesh['tetrahedra']),dtype=bool)).reshape(-1)
    for element in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,element,independent=True,independent_order=order);mat=materials[mesh['material_ids'][element]]
        young=mat['young_kpa'];nu=mat['poisson'];shear=young/(2*(1+nu));bulk=young/(3*(1-2*nu))
        dev=np.zeros((6,6));dev[:3,:3]=-2*shear/3
        for a in range(3):dev[a,a]+=2*shear;dev[a+3,a+3]=shear
        dofs=(3*q['nodes'][:,None]+np.arange(3)).reshape(-1)
        pi=np.array([index[(node,mat['id'])] for node in mesh['tetrahedra'][element]])
        l=np.column_stack((1-q['natural'].sum(axis=1),q['natural']))
        div=q['b'][:,:3].sum(axis=1)
        a=np.einsum('qai,ab,qbj,q->ij',q['b'],dev,q['b'],q['weight'],optimize=True)
        b=-np.einsum('qa,qi,q->ai',l,div,q['weight'],optimize=True)
        m=np.einsum('qa,qb,q->ab',l,l,q['weight'],optimize=True)
        k[np.ix_(dofs,dofs)]+=a;k[np.ix_(ndof+pi,dofs)]+=b;k[np.ix_(dofs,ndof+pi)]+=b.T
        k[np.ix_(ndof+pi,ndof+pi)]-=m/bulk;mass[np.ix_(pi,pi)]+=m
        gradient=np.stack([q['b'][:,a,a::3] for a in range(3)],axis=2)
        scalar=np.einsum('qia,qja,q->ij',gradient,gradient,q['weight'],optimize=True)
        for a in range(3):h1[np.ix_(dofs[a::3],dofs[a::3])]+=scalar
        force[dofs]+=((q['weight']@q['shape'])[:,None]*mat['body_force_kn_m3']).reshape(-1)
    t,u0=constraint_map(mesh,case['boundary']);count=t.shape[1]
    lift=np.r_[u0,np.zeros(npres)];transform=np.zeros((total,count+npres));transform[:ndof,:count]=t;transform[ndof:,count:]=np.eye(npres)
    rhs=transform.T@(np.r_[force,np.zeros(npres)]-k@lift);reduced=transform.T@k@transform
    coordinates=np.linalg.solve(reduced,rhs);state=lift+transform@coordinates
    residual=float(np.linalg.norm(reduced@coordinates-rhs)/max(np.linalg.norm(rhs),1e-30))
    if residual>1e-10:raise ValueError('independent mixed true residual exceeded tolerance')
    coupling=t.T@k[:ndof,ndof:];sv=np.linalg.svd(coupling,compute_uv=False)
    rank=int((sv>max(sv.max(initial=0.)*1e-10,1e-14)).sum())
    beta=None
    if count:
        schur=coupling.T@np.linalg.solve(t.T@h1@t,coupling);lower=np.linalg.cholesky(mass)
        normalized=np.linalg.solve(lower,np.linalg.solve(lower,schur).T).T
        eigen=np.linalg.eigvalsh((normalized+normalized.T)/2);positive=eigen[eigen>1e-10*max(eigen.max(),1.)]
        if len(positive):beta=float(np.sqrt(positive.min()))
    imbalance=k[:ndof]@state-force;reaction=np.zeros(ndof)
    for node in range(len(xyz)):
        basis=t[3*node:3*node+3];block=imbalance[3*node:3*node+3]
        reaction[3*node:3*node+3]=block-basis@(basis.T@block)
    return dict(displacement=state[:ndof].reshape(-1,3),pressure=state[ndof:],pressure_keys=keys,
                reaction=reaction.reshape(-1,3),load=force.reshape(-1,3),true_relative_residual=residual,
                pressure_rank=rank,pressure_nullity=npres-rank,pressure_inf_sup_positive_beta=beta,free_displacement_dofs=count)


def verify(mesh_path,case_path,result_path,solver,report_path,order):
    if report_path.exists():raise ValueError('refusing to overwrite a verification report')
    mesh=read(mesh_path);case=read(case_path);native=read(result_path)
    if native['schema']!='soil-fem-mixed-mechanics-result/1':raise ValueError('wrong native formulation')
    if native['solver_binary_sha256']!=sha(solver):raise ValueError('native solver binary does not match')
    ref.normalize_conditions(case['boundary'])
    if case!=native['case']:raise ValueError('native input case differs')
    expected=solve_reference(mesh,case,order);actual_pressure={(p['node_index'],p['material_id']):p['compression_positive_kpa'] for p in native['pressure_dofs']}
    if set(actual_pressure)!=set(expected['pressure_keys']):raise ValueError('pressure DOF identity differs')
    pressure=np.array([actual_pressure[key] for key in expected['pressure_keys']])
    def metric(a,b,floor):
        a=np.asarray(a);b=np.asarray(b)
        return dict(relative_l2=float(np.linalg.norm(a-b)/max(np.linalg.norm(b),floor)),maximum_absolute=float(np.max(np.abs(a-b))))
    metrics={key:metric(native[field],expected[key],1.) for key,field in (('reaction','reaction_kn'),('load','applied_nodal_load_kn'))}
    metrics['displacement']=metric(native['displacement_m'],expected['displacement'],1e-12)
    metrics['pressure']=metric(pressure,expected['pressure'],1.)
    index={key:i for i,key in enumerate(expected['pressure_keys'])};stress=[];strain=[];gp_pressure=[];materials={m['id']:m for m in case['materials']}
    for element in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,element);mat=materials[mesh['material_ids'][element]];g=mat['young_kpa']/(2*(1+mat['poisson']))
        eps=np.einsum('qai,i->qa',q['b'],expected['displacement'][q['nodes']].reshape(-1))
        l=np.column_stack((1-q['natural'].sum(axis=1),q['natural']));pc=l@expected['pressure'][[index[(node,mat['id'])] for node in mesh['tetrahedra'][element]]]
        sig=eps*np.array([2*g,2*g,2*g,g,g,g]);sig[:,:3]-=(2*g*eps[:,:3].sum(axis=1)/3+pc)[:,None]
        stress.extend(sig);strain.extend(eps);gp_pressure.extend(pc)
    metrics.update(stress=metric(native['stress_tension_positive_kpa'],stress,1.),strain=metric(native['strain_engineering'],strain,1e-12),
                   gp_pressure=metric(native['pressure_compression_positive_kpa'],gp_pressure,1.))
    if not all(m['relative_l2']<1e-6 for m in metrics.values()):raise ValueError(('independent field comparison failed',metrics))
    files=[mesh_path,case_path,result_path,solver,Path(__file__),Path(ref.__file__)]
    report=dict(passed=True,goal_complete=False,p6_complete=False,
        scope='Independent linear mixed 3D field comparison for this mesh; pressure rank is not a general stability proof',
        quadrature_order=order,metrics=metrics,reference_true_relative_residual=expected['true_relative_residual'],
        pressure_dofs=len(pressure),free_displacement_dofs=expected['free_displacement_dofs'],
        pressure_coupling_rank=expected['pressure_rank'],pressure_coupling_nullity=expected['pressure_nullity'],
        pressure_inf_sup_positive_beta=expected['pressure_inf_sup_positive_beta'],sha256={str(p):sha(p) for p in files})
    with report_path.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('case',type=Path);p.add_argument('result',type=Path)
    p.add_argument('solver',type=Path);p.add_argument('report',type=Path);p.add_argument('--reference-order',type=int,choices=range(5,10),default=5);a=p.parse_args()
    verify(*(v.resolve() for v in (a.mesh,a.case,a.result,a.solver,a.report)),a.reference_order)
