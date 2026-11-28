"""Bound the same 10 nm ball after a fixed inverse action in the free coordinates."""
import numpy as np


def check(red, refined, elastic, f, err, force_guard, kn, kt, mu, radius, archive):
    free=red.free
    indices=np.full(len(free),-1)
    indices[free]=np.arange(free.sum())
    K=red.g.k[np.ix_(free,free)].copy()
    scale=float(np.linalg.norm(K))
    for sl,_,ids,_,mass,_,_ in red.blocks:
        gram=mass.sum(axis=0).reshape(15,15)
        for d,k in enumerate((kt,kt,kn)):
            for i,a in enumerate(ids):
                for j,b in enumerate(ids):
                    ia,ib=indices[3*a+d],indices[3*b+d]
                    assert np.all(ia>=0) and np.all(ib>=0)
                    K[np.ix_(ia,ib)]+=(1 if i==j else -1)*k*gram
        # Absolute sum bounds Gram integration and nodal assembly roundoff.
        scale+=4*max(kn,kt)*np.linalg.norm(np.abs(mass).sum(0))
    A=np.linalg.solve(K,np.eye(len(K)))
    norm_A=float(np.linalg.norm(A)*(1+1e-12))
    # Frobenius norms bound the induced Euclidean norms. 1e-9 exceeds the
    # Bernstein restriction, moment assembly, summation and dense product errors.
    assembly=1e-9*norm_A*scale+1e-12
    base_error=float(np.linalg.norm(np.eye(len(K))-A@K)*(1+1e-12)+assembly)
    perturb=0.
    for (sl,_,ids,_,_,_,_),(_,w,_,_,area,gain2) in zip(red.blocks,refined,strict=True):
        sel=~elastic[sl]
        W=w.reshape(-1,15,15)[sel]
        if not len(W):
            continue
        # D(x) maps free DOFs to a 3-component jump. Its rows have disjoint
        # supports. A D_t(x)^T is a convex combination of its Bernstein
        # coefficients, including both signed sides of the same child.
        ampl2=np.zeros(W.shape[:2])
        for d in (0,1):
            V=(A[:,indices[3*ids[1]+d]]-A[:,indices[3*ids[0]+d]]).T
            actions=W@V
            ampl2+=np.sum(actions**2,axis=-1)
        right=np.sqrt(gain2[sel])
        left=np.sqrt(ampl2.max(axis=1))+1e-10*norm_A*right
        # The disk loses at most kt in its tangential block and acquires a
        # normal column of norm <= mu*kn. The normal force is still elastic.
        perturb+=float(area[sel]@(left*right))*np.hypot(kt,mu*kn)
    q=base_error+perturb*(1+1e-9)
    fv=f.ravel()[free]
    ev=err.ravel()[free]
    applied=A@fv
    spread=abs(A)@ev
    np.savez_compressed(archive,free_dofs=np.flatnonzero(free),center_kn=fv,center_radius_kn=ev,
                        after_inverse_m=applied,after_inverse_radius_m=spread)
    eta=float(np.linalg.norm(abs(applied)+spread)*(1+1e-10)+
              norm_A*(force_guard+1e-12*np.linalg.norm(fv)))
    ok=bool(q<1 and eta+q*radius<radius)
    return dict(reference_inverse_frobenius_upper_m_kn=norm_A,elastic_inverse_defect=base_error,
                assembly_guard=assembly,projected_tangent_defect=float(perturb*(1+1e-9)),
                contraction_upper=float(q),center_after_inverse_upper_m=eta,
                ball_image_upper_m=eta+q*radius,included=ok,
                root_radius_upper_m=eta/(1-q)*(1+1e-12) if ok else None)
