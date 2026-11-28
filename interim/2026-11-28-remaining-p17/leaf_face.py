"""Split only the largest p16 Bernstein leaf's nominal integral at a stopped step."""
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.with_name("2026-11-28-remaining-p16")))
import paired_face as p16

p15 = p16.p15_face
baseline = p15.integrate


def integrate(red, refined, center_box, displacements, kn, kt, mu, f, pf, pe, guard, old, total, out):
    center, all_hist, count = baseline(red, refined, center_box, displacements, kn, kt, mu,
                                       f, pf, pe, guard, old, total, out)
    sign, step, _ = out.stem.split("-")
    if int(step) < {"plus":833, "minus":832}[sign]:
        return center, all_hist, count
    A, ix = p15.fixed_inverse(red,kn,kt)
    norm_a = old["reference_inverse_frobenius_upper_m_kn"]
    mats = p15.loc.subdivision()
    grid = (mats[:,None] @ mats[None]).reshape(-1,15,15)
    finer = (mats[:,None] @ grid[None]).reshape(-1,15,15)
    basis = np.array([math.factorial(4)/math.prod(math.factorial(int(i)) for i in alpha)/3**4
                      for alpha in p16.chain.face.base.prev.ALPHA])

    def terms(grid, weights, areas, jumps, amps):
        n = len(grid)
        cells = (grid[None] @ weights[:,None]).reshape(-1,15,15)
        sizes = np.repeat(areas/n,n)
        h = np.zeros((len(cells),2))
        enc = h.copy(),h.copy()
        for jump,u in zip(jumps,displacements,strict=True):
            coeff = cells @ jump
            point = basis @ coeff
            v = kt*point[:,:2]+h
            radius = -mu*kn*point[:,2]
            pad = 1e-11*abs(u).max()+1e-90
            assert (coeff[:,:,2].max(1)+pad<0).all()
            h -= v*np.maximum(0,1-radius/np.maximum(np.linalg.norm(v,axis=1),1e-90))[:,None]
            enc,_ = p16.chain.project((coeff.min(1)-pad,coeff.max(1)+pad),enc,kt,mu*kn)
        width = np.maximum(abs(h-enc[0]),abs(enc[1]-h))*(1+1e-12)+1e-85
        error = sizes*np.repeat(amps,n)*np.linalg.norm(width,axis=1)
        return error.reshape(-1,n),cells.mean(1).T@(sizes[:,None]*h),np.max(abs(cells),axis=1).T@(sizes[:,None]*width)

    scores, keys, contexts = [], [], []
    for block_index,(block,mesh) in enumerate(zip(red.blocks,refined,strict=True)):
        sl,_,ids,*_ = block
        _,weights,_,orient,areas,gain2 = mesh
        weights = weights.reshape(-1,15,15)
        jumps = [(u.reshape(-1,3)[ids[1]]-u.reshape(-1,3)[ids[0]])*[1,1,orient] for u in displacements]
        lo,hi = (v[sl] for v in center_box)
        active = np.flatnonzero(np.any((lo!=0)|(hi!=0),axis=1))
        columns = [(A[:,ix[3*ids[1]+d]]-A[:,ix[3*ids[0]+d]]).T for d in (0,1)]
        amps = np.zeros(len(areas))
        for start in range(0,len(active),24):
            take = active[start:start+24]
            amps[take] = np.sqrt(sum(np.sum((weights[take]@v)**2,axis=2) for v in columns).max(1))
            amps[take] += 1e-10*norm_a*np.sqrt(gain2[take])
            error,_,_ = terms(grid,weights[take],areas[take],jumps,amps[take])
            scores.extend(error)
            keys.extend((block_index,int(k)) for k in take)
        contexts.append((ids,weights,areas,jumps,amps))
    scores = np.asarray(scores)
    totals = scores.sum(axis=1)
    winner = int(np.argmax(totals))
    block,leaf = keys[winner]
    ids,weights,areas,jumps,amps = contexts[block]
    subset = slice(leaf,leaf+1)
    args = weights[subset],areas[subset],jumps,amps[subset]
    before,term0,width0 = terms(grid,*args)
    after,term1,width1 = terms(finer,*args)
    np.testing.assert_allclose(before[0],scores[winner],rtol=1e-13,atol=1e-40)
    with np.load(out) as data:
        saved = dict(data)
    error0 = float(saved["direct_error_m"])
    np.testing.assert_allclose(totals.sum(),error0,rtol=2e-13,atol=0.)
    error = math.fsum((error0,-before.sum(),after.sum()))*(1+1e-14)
    force,spread = saved["center_kn"].copy(),saved["center_nodal_radius_kn"].copy()
    for nodes,side in zip(ids,(-1,1),strict=True):
        force[nodes,:2] += side*(term1-term0)
        spread[nodes,:2] += width1-width0
    action = A@force.ravel()[red.free]
    q = total["contraction_upper"]
    def test(past,past_width):
        fv = (force+past).ravel()[red.free]
        ev = past_width.ravel()[red.free]
        image = abs(A@fv)+abs(A)@ev
        rounding = norm_a*(guard+1e-12*np.linalg.norm(fv))
        eta = (np.linalg.norm(image)+error)*(1+1e-9)+rounding
        bound = eta+q*p16.chain.R
        ok = bool(q<1 and bound<p16.chain.R)
        return dict(center_after_inverse_upper_m=eta,ball_image_upper_m=bound,contraction_upper=q,
                    included=ok,root_radius_upper_m=eta/(1-q)*(1+1e-12) if ok else None,
                    direct_center_remainder_m=error*(1+1e-9),force_guard_m=rounding,
                    nodalized_same_integral_image_m=np.linalg.norm(image+abs(A)@spread.ravel()[red.free])*(1+1e-9)+rounding+q*p16.chain.R)
    zero = np.zeros_like(force)
    new_count = count+len(finer)-len(grid)
    np.savez_compressed(out.with_name(f"{sign}-{step}-leaves.npz"),block_leaf=keys,
                        tile_remainder_m=scores,selected_index=winner,child_tile_remainder_m=after,
                        prior_direct_remainder_m=error0,updated_direct_remainder_m=error,
                        roots_per_block=[len(c[2]) for c in contexts],center_delta_m=action-saved["center_after_inverse_m"])
    np.savez_compressed(out,**(saved|dict(center_kn=force,center_nodal_radius_kn=spread,
                        center_after_inverse_m=action,direct_error_m=error,tiles=new_count)))
    return test(zero,zero),test(pf,pe),new_count


if __name__ == "__main__":
    p15.integrate = integrate
    p16.run(sys.argv[1],Path(sys.argv[2]),"projected")
