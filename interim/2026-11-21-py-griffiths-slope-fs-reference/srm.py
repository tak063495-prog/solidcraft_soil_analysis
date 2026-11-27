"""Independent small-strain plane-strain CST FEM / MC strength reduction.

NumPy/SciPy only. Constitutive backward Euler active-set return in principal
stress, all six MC planes; no smoothing, hardening or production solver calls.
"""
import itertools
import time
import warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, MatrixRankWarning


class Unresolved(RuntimeError):
    pass


def elastic(E=100000., nu=.3):
    g=E/(2*(1+nu)); la=E*nu/((1+nu)*(1-2*nu))
    d=np.zeros((4,4)); d[:3,:3]=la; d[np.arange(3),np.arange(3)]+=2*g; d[3,3]=g
    return d


class MohrCoulomb:
    def __init__(self, factor=1., psi=0., E=100000., nu=.3, cohesion=10., phi=20.):
        if factor<=0 or not 0<=psi<=phi or not 0<=phi<90 or cohesion<0:
            raise ValueError('invalid material')
        self.D=elastic(E,nu); self.invD=np.linalg.inv(self.D)
        # psi is held fixed under reduction; default zero matches fixture.
        friction=np.arctan(np.tan(np.deg2rad(phi))/factor)
        self.k=2*cohesion/factor*np.cos(friction)
        self.A=[]; flow=[]
        for i,j in itertools.permutations(range(3),2):
            a=np.zeros(3); a[i]=1+np.sin(friction); a[j]=-(1-np.sin(friction)); self.A.append(a)
            b=np.zeros(3); b[i]=1+np.sin(np.deg2rad(psi)); b[j]=-(1-np.sin(np.deg2rad(psi))); flow.append(b)
        self.A=np.array(self.A); flow=np.array(flow); self.active=[]
        dp=self.D[:3,:3]
        for size in [1,2,3]:
            for ids in itertools.combinations(range(6),size):
                a=self.A[list(ids)]; g=flow[list(ids)]; m=a@dp@g.T
                if np.linalg.cond(m)<1e12:
                    self.active.append((list(ids),np.linalg.inv(m),dp@g.T))

    def update(self, trial):
        trial=np.atleast_2d(trial)
        tensor=np.zeros((len(trial),3,3)); tensor[:,0,0]=trial[:,0]; tensor[:,1,1]=trial[:,1]
        tensor[:,2,2]=trial[:,2]; tensor[:,0,1]=tensor[:,1,0]=trial[:,3]
        val,vec=np.linalg.eigh(tensor)
        corrected=val.copy(); f=val@self.A.T-self.k
        yielded=np.max(f,axis=1)>1e-8
        remaining=yielded.copy()
        for ids,inv,dpgt in self.active:
            rows=np.flatnonzero(remaining)
            if not len(rows): break
            multipliers=(f[rows][:,ids])@inv.T
            sig=val[rows]-multipliers@dpgt.T
            valid=(np.min(multipliers,axis=1)>=-1e-12)&(np.max(sig@self.A.T-self.k,axis=1)<=1e-7)
            corrected[rows[valid]]=sig[valid]; remaining[rows[valid]]=False
        if np.any(remaining):
            raise Unresolved('constitutive active set unresolved (including zero-dilation apex)')
        tensor=np.einsum('nik,nk,njk->nij',vec,corrected,vec)
        result=tensor[:,[0,1,2,0],[0,1,2,1]]
        # Eigen reconstruction roundoff must not create plasticity on an elastic branch.
        result[~yielded]=trial[~yielded]
        increment=(trial-result)@self.invD.T
        increment[~yielded]=0.
        return result,increment

    def tangent(self, trial):
        # Numerical consistent tangent of the local map with engineering shear.
        tangent=np.empty((len(trial),4,4)); eps=1e-8
        for j in range(4):
            delta=self.D[:,j]*eps
            tangent[:,:,j]=(self.update(trial+delta)[0]-self.update(trial-delta)[0])/(2*eps)
        return tangent


def mesh(n):
    """Mapped structured triangles on (0,0),(32,0),(12,10),(0,10)."""
    if n<1: raise ValueError('n must be positive')
    nx=3*n; ny=n
    xy=np.array([(i/nx*(32-2*(10*j/ny)),10*j/ny) for j in range(ny+1) for i in range(nx+1)])
    triangles=[]
    for j in range(ny):
        for i in range(nx):
            a=j*(nx+1)+i; b=a+1; d=a+nx+1; c=d+1
            triangles.extend([[a,b,c],[a,c,d]])
    conn=np.array(triangles); pts=xy[conn]
    twice=np.linalg.det(np.stack((pts[:,1]-pts[:,0],pts[:,2]-pts[:,0]),axis=2))
    area=twice/2
    grad=np.stack((pts[:,[1,2,0],1]-pts[:,[2,0,1],1],pts[:,[2,0,1],0]-pts[:,[1,2,0],0]),axis=2)/twice[:,None,None]
    B=np.zeros((len(conn),4,6))
    B[:,0,::2]=grad[:,:,0]; B[:,1,1::2]=grad[:,:,1]
    B[:,3,::2]=grad[:,:,1]; B[:,3,1::2]=grad[:,:,0]
    dofs=(2*conn[:,:,None]+np.arange(2)).reshape(-1,6)
    fixed=np.zeros((len(xy),2),bool)
    fixed[:,0]=(np.isclose(xy[:,0],0)|np.isclose(xy[:,1],0)); fixed[:,1]=np.isclose(xy[:,1],0)
    load=np.zeros(2*len(xy)); np.add.at(load,dofs[:,1::2],np.broadcast_to(-20*area[:,None]/3,(len(conn),3)))
    h=float(np.max(np.linalg.norm(pts-np.roll(pts,1,axis=1),axis=2)))
    return dict(xy=xy,conn=conn,area=area,B=B,dofs=dofs,free=np.flatnonzero(~fixed.ravel()),load=load,h=h)


def trial_solve(m, factor, steps=8, iterations=35, seconds=20., psi=0.):
    start=time.monotonic(); mat=MohrCoulomb(factor,psi)
    B,dofs,area,load,free=[m[k] for k in ['B','dofs','area','load','free']]
    nd=len(load); ne=len(area); u=np.zeros(nd); stress=np.zeros((ne,4)); plastic=np.zeros((ne,4))
    history=[]; status='iteration_limit'; reason=''; last_residual=None
    row=np.broadcast_to(dofs[:,:,None],(ne,6,6)).ravel(); col=np.broadcast_to(dofs[:,None,:],(ne,6,6)).ravel()
    def evaluate(du,old,target):
        strain=np.einsum('eij,ej->ei',B,du[dofs]); tr=old+strain@mat.D.T
        sig,dp=mat.update(tr)
        internal=np.zeros(nd); np.add.at(internal,dofs,(area[:,None]*np.einsum('eij,ei->ej',B,sig)))
        residual=target*load-internal
        norm=float(np.linalg.norm(residual[free])/max(np.linalg.norm(target*load[free]),1.))
        return norm,residual,tr,sig,dp
    try:
        for step in range(1,steps+1):
            du=np.zeros(nd); target=step/steps
            for iteration in range(iterations):
                if time.monotonic()-start>seconds: raise TimeoutError('trial wall time budget')
                norm,r,tr,sig,dp=evaluate(du,stress,target); last_residual=norm
                if time.monotonic()-start>seconds: raise TimeoutError('trial budget exceeded during residual evaluation')
                history.append(dict(step=step,iteration=iteration,relative_residual=norm))
                if norm<1e-8:
                    u+=du; stress=sig; plastic+=dp; break
                tangent=mat.tangent(tr)
                ke=area[:,None,None]*np.einsum('eik,eij,ejl->ekl',B,tangent,B)
                K=coo_matrix((ke.ravel(),(row,col)),shape=(nd,nd)).tocsr()
                with warnings.catch_warnings():
                    warnings.simplefilter('error',MatrixRankWarning)
                    correction=spsolve(K[free][:,free],r[free])
                if not np.all(np.isfinite(correction)): raise Unresolved('nonfinite Newton update')
                improved=False
                for alpha in [1.,.5,.25,.125,.0625,.03125]:
                    candidate=du.copy(); candidate[free]+=alpha*correction
                    try: newnorm=evaluate(candidate,stress,target)[0]
                    except Unresolved: continue
                    if newnorm<norm:
                        du=candidate; improved=True; break
                if not improved: raise Unresolved('line search stalled')
            else:
                raise Unresolved('Newton iteration limit')
        status='equilibrated'; reason='all gravity increments equilibrated'
    except TimeoutError as exc:
        status='timeout'; reason=str(exc)
    except (Unresolved,MatrixRankWarning,np.linalg.LinAlgError) as exc:
        status='iteration_limit' if 'iteration limit' in str(exc) else 'unresolved'; reason=str(exc)
    return dict(factor=float(factor),status=status,reason=reason,seconds=time.monotonic()-start,
                relative_residual=last_residual,completed_steps=sum(1 for x in history if x['relative_residual']<1e-8),
                history=history,displacement=u.reshape(-1,2).tolist(),stress=stress.tolist(),plastic_strain=plastic.tolist(),
                fields_scope='last_committed_gravity_step',failure_certified=False)


def bisect_srm(n, low=.8, high=2., width=.05, max_trials=8, seconds=20., steps=8, psi=0.):
    """Bisect computational equilibrium boundary, NEVER infer physical failure.

    Unknown endpoints remain unknown. A finite certified Fs upper bound needs
    an independent failure certificate, absent from this exploratory FEM.
    """
    if not 0<low<high or width<=0 or max_trials<2: raise ValueError('bad bracket/budget')
    m=mesh(n); trials=[]
    def run(f):
        r=trial_solve(m,f,steps=steps,seconds=seconds,psi=psi); trials.append(r); return r
    left=run(low); right=run(high)
    eq=low if left['status']=='equilibrated' else None
    upper=None
    if eq is not None and right['status']!='equilibrated':
        upper=high
        while upper-eq>width and len(trials)<max_trials:
            mid=(eq+upper)/2
            r=run(mid)
            if r['status']=='equilibrated': eq=mid
            else: upper=mid
    elif right['status']=='equilibrated': eq=high
    stop=('no_equilibrated_lower' if eq is None else 'no_unresolved_upper' if upper is None else
          'width_reached' if upper-eq<=width else 'trial_limit')
    return dict(n=n,h_max_m=m['h'],nodes=len(m['xy']),elements=len(m['conn']),psi_deg=psi,
                equilibrium_lower=eq,unresolved_upper=upper,
                computational_transition=[eq,upper] if upper is not None else None,
                certified_fs_interval=None,accepted=None,stop_reason=stop,trial_cap_reached=len(trials)>=max_trials,
                interpretation='computational endpoints only; unresolved is not failure; no Fs midpoint',
                mesh=dict(xy=m['xy'].tolist(),triangles=m['conn'].tolist()),trials=trials)
