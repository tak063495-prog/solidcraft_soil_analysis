"""Integrate consecutive nominal disk projections at common material positions."""
import inspect
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")
sys.path.insert(0, str(P14))
import face_chain as chain

loc = chain.loc


def fixed_inverse(red, kn, kt):
    free = red.free
    ix = np.full(free.size,-1)
    ix[free] = np.arange(free.sum())
    K = red.g.k[np.ix_(free,free)].copy()
    for _,_,ids,_,mass,_,_ in red.blocks:
        gram = mass.sum(0).reshape(15,15)
        for d,k in enumerate((kt,kt,kn)):
            for i,a in enumerate(ids):
                for j,b in enumerate(ids):
                    ia,ib = ix[3*a+d],ix[3*b+d]
                    assert min(ia.min(),ib.min())>=0
                    K[np.ix_(ia,ib)] += (1 if i==j else -1)*k*gram
    return np.linalg.solve(K,np.eye(len(K))),ix


def integrate(red, refined, center_box, displacements, kn, kt, mu, f, pf, pe,
              guard, old, total, out):
    A,ix = fixed_inverse(red,kn,kt)
    norm_A = old["reference_inverse_frobenius_upper_m_kn"]
    assert norm_A == np.linalg.norm(A)*(1+1e-12)
    basis = np.array([math.factorial(4)/math.prod(math.factorial(int(i)) for i in alpha)/3**4
                      for alpha in chain.face.base.prev.ALPHA])
    mats = loc.subdivision()
    # Two further four-way subdivisions for nominal integration only.
    restriction = (mats[:,None]@mats[None]).reshape(16,15,15)
    force = f.copy()
    nodal_width = np.zeros_like(f)
    direct_error = 0.
    counts = 0
    for (sl,_,ids,_,_,mean,_),(_,weights,_,sign,areas,gain2) in zip(red.blocks,refined,strict=True):
        lo,hi = [b[sl] for b in center_box]
        old_correction = mean.T@((lo+hi)/2)
        terms = np.zeros((15,2))
        spreads = np.zeros((15,2))
        jumps = [(u.reshape(-1,3)[ids[1]]-u.reshape(-1,3)[ids[0]])*[1,1,sign]
                 for u in displacements]
        # A leaf whose nominal enclosure is exactly zero contributes exactly zero.
        active = np.flatnonzero(np.any((lo!=0)|(hi!=0),axis=1))
        W = weights.reshape(-1,15,15)
        V = [(A[:,ix[3*ids[1]+d]]-A[:,ix[3*ids[0]+d]]).T for d in (0,1)]
        for start in range(0,len(active),24):
            take = active[start:start+24]
            cells = (restriction[None]@W[take,None]).reshape(-1,15,15)
            sizes = np.repeat(areas[take]/16,16)
            h = np.zeros((len(cells),2))
            enc = h.copy(),h.copy()
            for jump,u in zip(jumps,displacements,strict=True):
                coeff = cells@jump
                point = basis@coeff
                v = kt*point[:,:2]+h
                r = -mu*kn*point[:,2]
                pad = 1e-11*abs(u).max()+1e-90
                assert (coeff[:,:,2].max(1)+pad<0).all()
                h -= v*np.maximum(0,1-r/np.maximum(np.linalg.norm(v,axis=1),1e-90))[:,None]
                enc,_ = chain.project((coeff.min(1)-pad,coeff.max(1)+pad),enc,kt,mu*kn)
            width = np.maximum(abs(h-enc[0]),abs(enc[1]-h))*(1+1e-12)+1e-85
            terms += cells.mean(1).T@(sizes[:,None]*h)
            spreads += np.max(abs(cells),axis=1).T@(sizes[:,None]*width)
            # Each A D_t(x)^T is a convex combination of its Bernstein
            # coefficients. Use its Frobenius norm before nodal intervalization.
            amps = np.sqrt(sum(np.sum((W[take]@v)**2,axis=2) for v in V).max(1))
            amps += 1e-10*norm_A*np.sqrt(gain2[take])
            direct_error += (sizes*np.repeat(amps,16))@np.linalg.norm(width,axis=1)
            counts += len(cells)
        for nodes,side in zip(ids,(-1,1),strict=True):
            force[nodes,:2] += side*(terms-old_correction)
            nodal_width[nodes,:2] += spreads
    def test(past, past_width):
        fv = (force+past).ravel()[red.free]
        ev = past_width.ravel()[red.free]
        action = A@fv
        spread = abs(A)@ev
        roundoff = norm_A*(guard+1e-12*np.linalg.norm(fv))
        residual = np.linalg.norm(abs(action)+spread)
        eta = (residual+direct_error)*(1+1e-9)+roundoff
        q = total["contraction_upper"]
        image = eta+q*chain.R
        ok = bool(q<1 and image<chain.R)
        nodal = np.linalg.norm(abs(action)+spread+abs(A)@nodal_width.ravel()[red.free])
        return dict(center_after_inverse_upper_m=float(eta),ball_image_upper_m=float(image),
                    contraction_upper=q,included=ok,root_radius_upper_m=float(eta/(1-q)*(1+1e-12)) if ok else None,
                    direct_center_remainder_m=float(direct_error*(1+1e-9)),force_guard_m=float(roundoff),
                    nodalized_same_integral_image_m=float(nodal*(1+1e-9)+roundoff+q*chain.R))
    zero = np.zeros_like(force)
    cen, all_hist = test(zero,zero),test(pf,pe)
    np.savez_compressed(out,center_kn=force,center_nodal_radius_kn=nodal_width,past_kn=pf,past_radius_kn=pe,
                        direct_error_m=direct_error,tiles=counts,center_after_inverse_m=A@force.ravel()[red.free])
    return cen,all_hist,counts


def run(sign,out):
    # Reuse p14's certified prefix, paired-history propagation and tangent check.
    # Replace only the nominal residual test, starting at the first p14 failure.
    source = inspect.getsource(chain.run)
    anchor = "        upper = (fc[1]+fp[1])"
    assert source.count(anchor)==1
    source = source.replace(anchor, '''        if idx+1 >= 831:
            row.update(nodal_center=center_only,nodal_full_history=inv)
            center_only, inv, tiles = integrate(red,refined,center_def,
                data["u"][p13["step"]-1:idx+1],kn,kt,mu,f,pf,pe,guard,center_only,inv,
                out/f"{sign}-{idx+1}-direct.npz")
            row["integration_tiles"] = tiles
''' + anchor)
    scope = dict(vars(chain),integrate=integrate)
    exec(compile(source,str(P14/"face_chain.py"),"exec"),scope)
    scope["run"](sign,out)


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
