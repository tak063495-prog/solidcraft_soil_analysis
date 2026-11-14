"""Input-only diagnostic for the existing two exact force-recovery gates."""
import numpy as np


def recover(system, residual, free=None):
    nd=system.c.shape[2]
    residual=np.asarray(residual)
    if residual.shape!=(nd,) or not np.isfinite(residual).all():
        raise ValueError('finite complete nodal vector required')
    free=np.eye(nd) if free is None else free
    pairs=sorted({(a,b) for p in system.points for a,b in zip(p['nodes_a'],p['nodes_b']) if a!=b})
    pm=np.zeros((nd,3*len(pairs)))
    for j,(a,b) in enumerate(pairs):
        pm[3*a:3*a+3,3*j:3*j+3]=-np.eye(3)
        pm[3*b:3*b+3,3*j:3*j+3]=np.eye(3)
    multipliers=np.linalg.lstsq(free.T@pm,free.T@residual,rcond=1e-12)[0]
    target=pm@multipliers
    first=float(np.max(np.abs(free.T@(target-residual)),initial=0))
    if first>1e-10:raise ValueError('unreproducible free nodal residual')
    matrix=np.einsum('qag,q->gqa',system.c,np.sqrt(system.weight)).reshape(nd,-1)
    x=np.linalg.lstsq(matrix,target,rcond=1e-12)[0]
    second=float(np.max(np.abs(matrix@x-target),initial=0))
    if second>1e-10:raise ValueError('unreproducible physical nodal transfer')
    q=x.reshape(-1,3)/np.sqrt(system.weight)[:,None]
    qn=np.einsum('qa,qa->q',q,system.normal)
    shear=q-qn[:,None]*system.normal
    if np.any(qn>0) or np.any(np.linalg.norm(shear,axis=1)>system.mu*(-qn)):
        raise ValueError('no-tension/Coulomb cone')
    return q,dict(free_nodal_maximum_kn=first,physical_nodal_maximum_kn=second)
