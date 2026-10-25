"""Assess displacement-element volumetric locking against a smooth exact 3D field.

Manufactured solution: u = curl(0, 0, f(x) f(y) f(z)), f(t)=t^2(1-t)^2.
It vanishes on the unit-cube boundary and div(u)=0. With constant shear
modulus G, b=-G Laplacian(u), independent of Poisson's ratio. Equivalent nodal
body forces are integrated independently, never generated from a stiffness
matrix or the exact nodal displacements. This assesses existing elements;
it does not introduce or certify a mixed incompressible formulation.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import fem_reference as ref

G=1000.


def fields(xyz):
    x,y,z=np.asarray(xyz).T
    def f(t):return t*t*(1-t)*(1-t)
    def d(t):return 2*t-6*t*t+4*t**3
    def dd(t):return 2-12*t+12*t*t
    def ddd(t):return -12+24*t
    u=np.column_stack((f(x)*d(y)*f(z),-d(x)*f(y)*f(z),np.zeros_like(x)))
    strain=np.column_stack((d(x)*d(y)*f(z),-d(x)*d(y)*f(z),np.zeros_like(x),
        (f(x)*dd(y)-dd(x)*f(y))*f(z),-d(x)*f(y)*d(z),f(x)*d(y)*d(z)))
    body=-G*np.column_stack((dd(x)*d(y)*f(z)+f(x)*ddd(y)*f(z)+f(x)*d(y)*dd(z),
        -ddd(x)*f(y)*f(z)-d(x)*dd(y)*f(z)-d(x)*f(y)*dd(z),np.zeros_like(x)))
    return u,strain,body


def mesh(n,quadratic):
    xyz=[[x/n,y/n,z/n] for z in range(n+1) for y in range(n+1) for x in range(n+1)]
    tet=[]
    for z,y,x in itertools.product(range(n),repeat=3):
        def node(a,b,c):return x+a+(n+1)*(y+b+(n+1)*(z+c))
        c=[node(*p) for p in ((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))]
        tet.extend([[c[i] for i in t] for t in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6))])
    result=dict(schema='soil-fem-mesh/3' if quadratic else 'soil-fem-mesh/2',units='m-kN-kPa',axis='z_up',
                crs='LOCAL',origin_m=[0.]*3,nodes=xyz,tetrahedra=tet,material_ids=[1]*len(tet),node_sets={},
                region_ids=[1]*len(tet),regions=[dict(id=1,body='ManufacturedCube',material_id=1,
                    bounds_m=[[0.]*3,[1.]*3],volume_m3=1.,
                    source_sha256=hashlib.sha256(b'analytic unit cube; six tetrahedra per subcube').hexdigest())])
    if quadratic:
        edges={};cells=[]
        for t in tet:
            cell=t.copy()
            for a,b in ref.EDGES:
                key=tuple(sorted((t[a],t[b])))
                if key not in edges:
                    edges[key]=len(xyz);xyz.append(((np.asarray(xyz[key[0]])+xyz[key[1]])/2).tolist())
                cell.append(edges[key])
            cells.append(cell)
        result.update(tetrahedra10=cells,quadrature_order=3)
    result['node_sets']={'boundary':[i for i,p in enumerate(xyz) if any(t in (0.,1.) for t in p)]}
    for i in range(len(xyz)):result['node_sets'][f'n{i}']=[i]
    return result


def body_load(mesh,order):
    out=np.zeros((len(mesh['nodes']),3))
    for i in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,i,independent=True,independent_order=order)
        out[q['nodes']]+=np.einsum('q,qi,qa->ia',q['weight'],q['shape'],fields(q['position'])[2])
    return out


def error(mesh,u,nu,exact_u_norm2,exact_strain_norm2):
    # Global exact squared norms are evaluated separately. The remaining cross
    # terms are degree <=13; GL8/Duffy integrates them exactly on straight cells.
    unorm=0.;ucross=0.;enorm=0.;ecross=0.;divnorm=0.
    shear_metric=np.array([2.,2.,2.,1.,1.,1.])
    for i in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,i,independent=True,independent_order=8)
        uh=q['shape']@u[q['nodes']];eh=np.einsum('qai,i->qa',q['b'],u[q['nodes']].reshape(-1))
        ue,ee,_=fields(q['position']);w=q['weight']
        unorm+=float(w@np.square(uh).sum(axis=1));ucross+=float(w@np.einsum('qa,qa->q',uh,ue))
        enorm+=float(w@(np.square(eh)*shear_metric).sum(axis=1))
        ecross+=float(w@(eh*ee*shear_metric).sum(axis=1))
        divnorm+=float(w@np.square(eh[:,:3].sum(axis=1)))
    strain_error2=max(0.,enorm-2*ecross+exact_strain_norm2)
    return dict(displacement_relative_l2=float(np.sqrt(max(0.,unorm-2*ucross+exact_u_norm2)/exact_u_norm2)),
        strain_relative_l2=float(np.sqrt(strain_error2/exact_strain_norm2)),
        displacement_norm_ratio=float(np.sqrt(unorm/exact_u_norm2)),divergence_l2=float(np.sqrt(divnorm)),
        energy_error_relative=float(np.sqrt((G*strain_error2+2*G*nu/(1-2*nu)*divnorm)/(G*exact_strain_norm2))))


def run(out,subdivisions,poissons,solver=None,command='solve-mechanics'):
    root=Path(__file__).resolve().parents[1];solver=(solver or root/'target/release/soil-fem-3d.exe').resolve()
    files=[Path(__file__),Path(ref.__file__),solver,root/'tools/run_logged.py']
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    fingerprints={str(p):digest(p) for p in files}
    out.mkdir(exist_ok=False)
    def write(path,value):
        with path.open('x',encoding='utf-8') as f:json.dump(value,f)
    # Differentiate the analytic displacement independently before using b.
    sample=np.array([[.23,.37,.61],[.59,.41,.73]]);h=1e-3;laplace=np.zeros_like(sample)
    for axis in range(3):
        delta=np.eye(3)[axis]*h
        laplace+=(-fields(sample+2*delta)[0]+16*fields(sample+delta)[0]-30*fields(sample)[0]
                  +16*fields(sample-delta)[0]-fields(sample-2*delta)[0])/(12*h*h)
    np.testing.assert_allclose(-G*laplace,fields(sample)[2],rtol=1e-7,atol=1e-8)
    r,w=np.polynomial.legendre.leggauss(8);r=(r+1)/2;w=w/2
    points=np.array(list(itertools.product(r,repeat=3)));weights=np.array([a*b*c for a,b,c in itertools.product(w,repeat=3)])
    exact_u,exact_e,_=fields(points)
    exact_u_norm2=float(weights@np.square(exact_u).sum(axis=1))
    exact_strain_norm2=float(weights@(np.square(exact_e)*[2.,2.,2.,1.,1.,1.]).sum(axis=1))
    rows=[]
    for quadratic in ((True,) if command=='solve-mixed-mechanics' else (False,True)):
        for n in subdivisions:
            m=mesh(n,quadratic);name=f'tet{10 if quadratic else 4}-n{n}';mp=out/(name+'-mesh.json');write(mp,m)
            subprocess.run([str(solver),'validate',str(mp)],cwd=root,check=True)
            force=body_load(m,8)
            if n==subdivisions[0]:
                alternative=body_load(m,7)
                np.testing.assert_allclose(force,alternative,rtol=1e-10,atol=1e-13)
            for nu in poissons:
                tag=name+'-nu'+str(nu).replace('.','_');cp=out/(tag+'-case.json');rp=out/(tag+'-result.json')
                case=dict(schema='soil-fem-mechanics/1',units='m-kN-kPa',
                    materials=[dict(id=1,young_kpa=2*G*(1+nu),poisson=nu,body_force_kn_m3=[0.]*3)],
                    boundary=dict(constraints=[dict(node_set='boundary',direction=np.eye(3)[a].tolist(),value_m=0.) for a in range(3)],
                        nodal_loads=[dict(node_set=f'n{i}',total_force_kn=f.tolist()) for i,f in enumerate(force)]),
                    relative_linear_tolerance=1e-10,maximum_free_dofs=100000,
                    linear_solver=dict(method='sparse_gmres_ilu0',options=dict(restart=100,maximum_iterations=5000,pivot_tolerance=1e-14)))
                write(cp,case)
                done=subprocess.run([sys.executable,str(root/'tools/run_logged.py'),str(out/(tag+'.log')),
                    str(solver),command,str(mp),str(cp),str(rp)],cwd=root)
                row=dict(element='Tet10' if quadratic else 'Tet4',subdivisions=n,poisson=nu,nodes=len(m['nodes']),elements=len(m['tetrahedra']),exit_code=done.returncode)
                if done.returncode==0:
                    native=json.loads(rp.read_text());u=np.asarray(native['displacement_m'])
                    assert native['true_relative_free_force_residual']<1e-10
                    assert np.max(np.abs(u[m['node_sets']['boundary']]))<1e-13
                    row.update(error(m,u,nu,exact_u_norm2,exact_strain_norm2))
                    row.update(true_relative_free_force_residual=native['true_relative_free_force_residual'],linear=native['linear'])
                rows.append(row);write(out/(tag+'-metrics.json'),row)
                print(json.dumps(row),flush=True)
    assert fingerprints=={str(p):digest(p) for p in files},'solver/source changed during assessment'
    record=dict(schema='soil-fem-incompressibility-assessment/1',assessment_complete=True,
        all_cases_solved=all(row['exit_code']==0 for row in rows),goal_complete=False,p6_complete=False,
        scope='Straight unit cube; smooth divergence-free manufactured field; '+command,
        analytic_body_force_checked_by_independent_finite_difference=True,gl7_gl8_body_load_agreement=True,
        exact_displacement_squared_norm=exact_u_norm2,exact_strain_squared_norm=exact_strain_norm2,
        input_and_source_sha256=fingerprints,cases=rows,
        limitations='Element assessment only; finite bulk modulus; no plasticity, CAD integration or general locking-free claim',
        evidence_sha256={p.name:digest(p) for p in out.iterdir() if p.is_file()})
    write(out/'assessment.json',record)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--subdivisions',type=int,nargs='+',default=[2,4,6])
    p.add_argument('--poissons',type=float,nargs='+',default=[.3,.49,.4999]);p.add_argument('--solver',type=Path)
    p.add_argument('--command',choices=['solve-mechanics','solve-mixed-mechanics'],default='solve-mechanics');a=p.parse_args()
    if len(set(a.subdivisions))!=len(a.subdivisions) or min(a.subdivisions)<2 or any(not(-1<v<.5) for v in a.poissons):
        p.error('distinct subdivisions >=2 and -1<nu<0.5 required')
    run(a.output.resolve(),a.subdivisions,a.poissons,a.solver,a.command)
