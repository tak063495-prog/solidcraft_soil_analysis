"""Independent Q4 degenerate director shell; NumPy/SciPy only.

5 DOFs/node: global u[3] and tangent director increments beta[2].
MITC-style covariant transverse shear tying; full membrane integration.
This is a reference experiment, not the production curved Tri6 implementation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent


def catalog():
    return json.loads((HERE / 'benchmarks.json').read_text(encoding='utf-8'))


def shape(r, s):
    signs = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    rr, ss = signs.T
    return ((1 + rr*r)*(1 + ss*s)/4,
            np.column_stack((rr*(1 + ss*s), ss*(1 + rr*r)))/4)


def mesh(case, n):
    """Regular quarter roof/hemisphere or eighth cylinder, analytic directors."""
    if not isinstance(n, int) or n < 1:
        raise ValueError('n must be a positive integer')
    g = case['geometry']
    xyz, directors, bases, fixed = [], [], [], set()
    for j in range(n+1):
        for i in range(n+1):
            k = j*(n+1)+i
            if case['id'] == 'hemisphere_hole':
                phi = math.pi/2 * i/n
                theta = math.pi/2 - math.radians(90-g['hole_polar_angle_deg'])*j/n
                d = np.array([math.sin(theta)*math.cos(phi),
                              math.sin(theta)*math.sin(phi), math.cos(theta)])
                b1 = np.array([-math.sin(phi), math.cos(phi), 0.])
                b2 = np.cross(d, b1)
                x = g['radius']*d
                if i in (0, n):
                    fixed.update((5*k+(1 if i == 0 else 0), 5*k+3))
                if i == 0 and j == 0:
                    fixed.add(5*k+2)  # eliminate vertical rigid translation only
            else:
                extent = g['half_angle_deg'] if case['id'] == 'scordelis_lo' else 90.
                theta = math.radians(extent)*i/n
                d = np.array([math.sin(theta), 0., math.cos(theta)])
                b1 = np.array([math.cos(theta), 0., -math.sin(theta)])
                b2 = np.array([0., 1., 0.])
                x = g['radius']*d + [0., g['length']/2*j/n, 0.]
                if i == 0:
                    fixed.update((5*k, 5*k+3))
                if i == n and case['id'] == 'pinched_cylinder':
                    fixed.update((5*k+2, 5*k+3))
                if j == 0:
                    fixed.update((5*k+1, 5*k+4))
                if j == n:
                    fixed.update((5*k, 5*k+2, 5*k+3))
            xyz.append(x)
            directors.append(d)
            bases.append(np.column_stack((b1, b2)))
    quads = []
    for j in range(n):
        for i in range(n):
            a = j*(n+1)+i
            quads.append([a, a+1, a+n+2, a+n+1])
    return np.array(xyz), np.array(directors), np.array(bases), np.array(quads), sorted(fixed)


def kinematics(x, d, bases, r, s, z):
    N, dN = shape(r, s)
    J = np.column_stack(((dN.T @ (x+z*d)).T, N@d))
    deriv = np.zeros((3, 3, 20))  # displacement component, natural derivative, DOF
    for a in range(4):
        for c in range(3):
            deriv[c, :2, 5*a+c] = dN[a]
        deriv[:, :2, 5*a+3:5*a+5] = z*bases[a, :, None, :]*dN[a, None, :, None]
        deriv[:, 2, 5*a+3:5*a+5] = N[a]*bases[a]
    return N, J, deriv


def covariant_shear(x, d, bases, r, s, alpha):
    _, J, deriv = kinematics(x, d, bases, r, s, 0.)
    return J[:, 2]@deriv[:, alpha, :] + J[:, alpha]@deriv[:, 2, :]


def element(x, d, bases, E, nu, thickness):
    """2x2 in surface, 2 through thickness; plane stress, shear kappa=5/6.

    X=sum N(Xa+z da), U=sum N(ua+z Ba beta_a).
    Membrane strains are projected from sym(dU/dX) onto local tangents.
    Covariant gamma_rz is tied at (0,+/-1), gamma_sz at (+/-1,0).
    No drilling DOF, penalty, hourglass stabilization, or fitted coefficients.
    """
    C = E/(1-nu*nu)*np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])
    shear_modulus = E/(2*(1+nu))*5/6
    kr = [covariant_shear(x,d,bases,0.,s,0) for s in (-1.,1.)]
    ks = [covariant_shear(x,d,bases,r,0.,1) for r in (-1.,1.)]
    ke, load = np.zeros((20,20)), np.zeros(20)
    for r in (-1/math.sqrt(3), 1/math.sqrt(3)):
        for s in (-1/math.sqrt(3), 1/math.sqrt(3)):
            N, J0, _ = kinematics(x,d,bases,r,s,0.)
            normal = np.cross(J0[:,0],J0[:,1])
            area = np.linalg.norm(normal)
            e3 = normal/area
            e1 = J0[:,0]/np.linalg.norm(J0[:,0])
            e2 = np.cross(e3,e1)
            # Map the assumed covariant strain tensor back to Cartesian.
            cov = np.zeros((3,3,20))
            cov[0,2] = cov[2,0] = ((1-s)*kr[0]+(1+s)*kr[1])/4
            cov[1,2] = cov[2,1] = ((1-r)*ks[0]+(1+r)*ks[1])/4
            inv0 = np.linalg.inv(J0)
            eps = np.einsum('ai,abk,bj->ijk',inv0,cov,inv0)
            Bs = np.array([2*np.einsum('i,ijk,j->k',e1,eps,e3),
                           2*np.einsum('i,ijk,j->k',e2,eps,e3)])
            ke += thickness*abs(np.linalg.det(J0))*shear_modulus*(Bs.T@Bs)
            for z in (-thickness/math.sqrt(12), thickness/math.sqrt(12)):
                _, J, deriv = kinematics(x,d,bases,r,s,z)
                grad = np.einsum('iak,aj->ijk',deriv,np.linalg.inv(J))
                B = np.array([np.einsum('i,ijk,j->k',e1,grad,e1),
                              np.einsum('i,ijk,j->k',e2,grad,e2),
                              np.einsum('i,ijk,j->k',e1,grad,e2)+
                              np.einsum('i,ijk,j->k',e2,grad,e1)])
                ke += thickness/2*abs(np.linalg.det(J))*(B.T@C@B)
            load[2::5] -= N*area  # unit global -z surface traction
    return (ke+ke.T)/2, load


def assemble(case, n):
    x,d,b,q,fixed = mesh(case,n)
    m = case['material']
    rows,cols,data = [],[],[]
    f = np.zeros(5*len(x))
    for nodes in q:
        dofs = (5*nodes[:,None]+np.arange(5)).ravel()
        ke, fe = element(x[nodes],d[nodes],b[nodes],m['E'],m['nu'],case['geometry']['thickness'])
        rows.extend(np.repeat(dofs,20)); cols.extend(np.tile(dofs,20)); data.extend(ke.ravel())
        if case['id'] == 'scordelis_lo':
            f[dofs] += case['load']['surface_traction_magnitude']*fe
    if case['id'] == 'pinched_cylinder':
        f[2] = -case['load']['point_force_full']/4
    if case['id'] == 'hemisphere_hole':
        f[0] = case['load']['point_force_full']/2
        f[5*n+1] = -case['load']['point_force_full']/2
    K = sparse.coo_matrix((data,(rows,cols)),shape=(len(f),len(f))).tocsr()
    return K,f,np.array(fixed),x,d,b


def rigid_modes(x,d,b):
    modes = np.zeros((5*len(x),6))
    for a in range(len(x)):
        modes[5*a:5*a+3,:3] = np.eye(3)
        for k in range(3):
            omega = np.eye(3)[k]
            modes[5*a:5*a+3,3+k] = np.cross(omega,x[a])
            modes[5*a+3:5*a+5,3+k] = b[a].T@np.cross(omega,d[a])
    return modes


def eigen_audit(case, n=2):
    """Free patch and supported patch, dimensionless scaling, full dense spectra.

    n=2 is deliberately small: this checks all eigenvalues, not only a sparse
    solver's selected modes. Threshold is numerical rank, not B-01 acceptance.
    """
    K,_,fixed,x,d,b = assemble(case,n)
    scale = np.tile([case['geometry']['radius']]*3+[1.,1.],len(x))
    A = K.toarray()*scale[:,None]*scale[None,:]
    A /= np.linalg.norm(A,2)
    eig = np.linalg.eigvalsh(A)
    tol = 1e-9
    R = rigid_modes(x,d,b)/scale[:,None]
    R /= np.linalg.norm(R,axis=0)
    residuals = np.linalg.norm(A@R,axis=0)
    free = np.setdiff1d(np.arange(len(A)),fixed)
    supported = np.linalg.eigvalsh(A[np.ix_(free,free)])
    passed = (np.count_nonzero(abs(eig)<tol)==6 and np.min(eig)>-tol and
              max(residuals)<tol and min(supported)>tol)
    return dict(n=n, numerical_rank_tolerance=tol,
                free_eigenvalues=eig.tolist(), supported_eigenvalues=supported.tolist(),
                rigid_residuals=residuals.tolist(),
                rigid_modes=6, detected_zero_modes=int(sum(abs(eig)<tol)),
                extra_zero_modes=max(0,int(sum(abs(eig)<tol))-6),
                passed=bool(passed))


def solve(case,n):
    start=time.perf_counter()
    K,f,fixed,x,_,_ = assemble(case,n)
    free=np.setdiff1d(np.arange(len(f)),fixed)
    u=np.zeros(len(f))
    u[free]=spsolve(K[free][:,free],f[free])
    if not np.isfinite(u).all():
        raise RuntimeError('nonfinite solution')
    residual=float(np.linalg.norm((K@u-f)[free])/np.linalg.norm(f[free]))
    if residual>1e-7:
        raise RuntimeError(f'equation residual too large: {residual}')
    probe=5*n+2 if case['id']=='scordelis_lo' else (2 if case['id']=='pinched_cylinder' else 0)
    value=float(u[probe])
    target=case['reference']['signed_value']
    row=dict(n=n,nodes=len(x),dofs=len(f),free_dofs=len(free),h_over_radius=case['h_span_over_radius']/n,
             displacement=value,reference=target,relative_error=abs((value-target)/target),
             equation_residual=residual,strain_energy=float(u@(K@u)/2),
             external_work=float(u@f),elapsed_seconds=time.perf_counter()-start)
    return row,x,u.reshape(-1,5)[:,:3]


def convergence(rows):
    """Reference-error orders and three-grid difference orders are diagnostics."""
    for i,row in enumerate(rows):
        row['reference_error_order']=None
        row['difference_order']=None
        row['richardson_limit']=None
        if i:
            prev=rows[i-1]
            if row['relative_error']>0 and prev['relative_error']>0:
                row['reference_error_order']=math.log(prev['relative_error']/row['relative_error'])/math.log(row['n']/prev['n'])
        if i>=2:
            a,b,c=rows[i-2:i+1]
            ratio=b['n']/a['n']
            if not math.isclose(c['n']/b['n'],ratio):
                continue
            d1=b['displacement']-a['displacement']
            d2=c['displacement']-b['displacement']
            if d1*d2>0 and abs(d2)<abs(d1):
                p=math.log(abs(d1/d2))/math.log(ratio)
                row['difference_order']=p
                row['richardson_limit']=c['displacement']+d2/(ratio**p-1)
    return rows


def write_nodes(path,x,u):
    with path.open('w',newline='',encoding='utf-8') as out:
        w=csv.writer(out);w.writerow(['node_id','x','y','z','ux','uy','uz'])
        for k,(point,disp) in enumerate(zip(x,u)):
            w.writerow([k,*point,*disp])


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--levels',type=int,nargs='+',default=[4,8,16])
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--allow-large',action='store_true',help='PC only, explicitly allow n>32')
    args=p.parse_args(argv)
    if args.levels!=sorted(set(args.levels)) or min(args.levels)<2:
        p.error('levels must be unique increasing integers >=2')
    if max(args.levels)>32 and not args.allow_large:
        p.error('n>32 requires --allow-large (user PC)')
    if args.out.exists() and any(args.out.iterdir()):
        p.error('output directory must be new or empty; preserve previous evidence')
    args.out.mkdir(parents=True,exist_ok=True)
    results=dict(schema_version=1,production_solver_called=False,B_01_passed=False,
                 benchmark_catalog_sha256=hashlib.sha256((HERE/'benchmarks.json').read_bytes()).hexdigest(),
                 implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),
                 benchmarks={})
    for case in catalog()['benchmarks']:
        rows=[]
        audits=[eigen_audit(case,n) for n in (2,4)]
        if not all(audit['passed'] for audit in audits):
            raise RuntimeError(f"{case['id']}: eigen audit failed")
        for n in args.levels:
            row,x,u=solve(case,n);rows.append(row)
            write_nodes(args.out/f"{case['id']}-n{n}.csv",x,u)
            print(case['id'],json.dumps(row),flush=True)
        results['benchmarks'][case['id']]=dict(rows=convergence(rows),eigen_audits=audits,
                                             benchmark_acceptance='未定義')
    (args.out/'results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
