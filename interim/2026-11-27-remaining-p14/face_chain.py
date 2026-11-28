"""Carry signed closed-face traction defects on p13's Bernstein leaves."""
import json
import os
from pathlib import Path
import sys

import numpy as np

from plans import P12, P13, save

sys.path.insert(0, str(P13))
# local_face imports P12 and save from plans (both modules define these).
import local_face as loc

face = loc.face
R = loc.RADIUS


def add_history(box, history, kt):
    shifted = [x.copy() for x in box]
    for a, h in zip(shifted, history, strict=True):
        v = h / kt
        a[:, :2] += v
    pad = 1e-13 * np.maximum(abs(box[0]), abs(box[1])) + 1e-90
    pad[:,:2] += 1e-13*np.maximum(abs(history[0]),abs(history[1]))/kt
    return shifted[0]-pad, shifted[1]+pad


def project(box, past, kt, km):
    trial = add_history(box, past, kt)
    stick = (kt*face.norm_hi(np.maximum(-trial[0][:,:2], trial[1][:,:2]))*(1+1e-12)
             < -km*trial[1][:,2]*(1-1e-12))
    dl, dh = loc.defect(trial, kt, km)
    # When elastic, P(kt*j+h) - kt*j equals h exactly.
    lo, hi = past[0]+dl, past[1]+dh
    guard = 1e-12*np.maximum(abs(lo), abs(hi)) + 1e-80
    return (np.where(stick[:,None], past[0], lo-guard),
            np.where(stick[:,None], past[1], hi+guard)), stick


def hull_by_parent(parent, h, n):
    lo, hi = np.full((n,2), np.inf), np.full((n,2), -np.inf)
    np.minimum.at(lo, parent, h[0])
    np.maximum.at(hi, parent, h[1])
    return lo, hi


def assemble(red, u, phase, t, center, past, kn, kt, guard0):
    g = red.g
    f = (g.k@u + g.initial-g.load - ((1-t)*(g.initial-g.load) if phase == 0 else 0)).reshape(-1,3)
    err = np.zeros_like(f)
    pf = np.zeros_like(f)
    pe = np.zeros_like(f)
    scale = np.linalg.norm(abs(g.k)@abs(u)+2*(abs(g.initial)+abs(g.load)))
    for sl, _, ids, _, mass, mean, amax in red.blocks:
        jump = u.reshape(-1,3)[ids[1]]-u.reshape(-1,3)[ids[0]]
        M = mass.sum(0).reshape(15,15)
        spring = M@(jump*[kt,kt,kn])
        scale += np.linalg.norm(abs(M)@(abs(jump)*[kt,kt,kn]))
        def nodal(h):
            lo, hi = h[0][sl], h[1][sl]
            return mean.T@((lo+hi)/2), amax.T@((hi-lo)/2)
        mid, spread = nodal(center)
        pm, ps = nodal(past)
        spring[:,:2] += mid
        for nodes, side in zip(ids, (-1,1), strict=True):
            f[nodes] += side*spring
            err[nodes,:2] += spread
            pf[nodes,:2] += side*pm
            pe[nodes,:2] += ps
    guard = max(guard0, float(8192*loc.EPS*(1+scale)))
    def limits(v,e,g):
        lo = np.maximum(abs(v)-e,0.).ravel()[red.free]
        hi = (abs(v)+e).ravel()[red.free]
        return [max(0.,float(np.linalg.norm(lo)*(1-1e-10)-g)), float(np.linalg.norm(hi)*(1+1e-10)+g)]
    return f, err, pf, pe, guard, limits(f,err,guard), limits(pf,pe,0.)


def run(sign, out):
    out.mkdir(exist_ok=True, parents=True)
    orig, _, gains, law, data = face.setup(sign)
    p13 = json.loads((P13/"c02"/(sign+".json")).read_text())
    kn, kt, mu = (law[k] for k in ("normal_stiffness_kpa_m","tangential_stiffness_kpa_m","friction_coefficient"))
    km = mu*kn
    n = len(gains)
    mask = np.zeros(n,bool)
    mask[[c["child"] for c in p13["signed_children"]]] = True
    refined, parent = loc.refine(orig, mask)
    jump = np.sqrt(np.concatenate([r[-1] for r in refined]))*(1+1e-12)
    zero = np.zeros((len(parent),2))
    u0 = data["u"][p13["step"]-1]
    radius = p13["certified_radius_m"]
    center0 = loc.box(refined,u0)
    archived_lo,archived_hi = loc.defect(center0,kt,km)
    for c in p13["signed_children"]:
        sel = parent == c["child"]
        np.testing.assert_array_equal([archived_lo[sel].min(0),archived_hi[sel].max(0)],c["center_defect_kpa"])
    hnom, _ = project(center0,(zero,zero),kt,km)
    htrue, _ = project(face.expand(center0,jump*radius),(zero,zero),kt,km)
    # At fixed x the previous P_r(v)-v is 1-Lipschitz in v and
    # 1-Lipschitz in r. Keep this paired bound separate from spatial variation.
    pair = np.where(np.any((htrue[0]!=0)|(htrue[1]!=0),axis=1), (kt+km)*jump*radius*(1+1e-12), 0.)
    steps = []
    for idx in range(p13["step"],len(data["u"])):
        phase,t,width = data["clock"][idx]
        prev, u = data["u"][idx-1:idx+1]
        coarse = face.expand(loc.box(orig,u),gains*R)
        _, stick = project(coarse,hull_by_parent(parent,htrue,n),kt,km)
        newmask = mask | ~stick
        if not np.array_equal(mask,newmask):
            repeats = np.where((~mask & newmask)[parent],4**loc.DEPTH,1)
            hnom, htrue = [tuple(np.repeat(x,repeats,axis=0) for x in h) for h in (hnom,htrue)]
            pair = np.repeat(pair,repeats)
            mask = newmask
            refined, newparent = loc.refine(orig,mask)
            np.testing.assert_array_equal(newparent,np.repeat(parent,repeats))
            parent = newparent
        area = np.concatenate([r[-2] for r in refined])
        jump = np.sqrt(np.concatenate([r[-1] for r in refined]))*(1+1e-12)
        center = loc.box(refined,u)
        ball = face.expand(center,jump*R)
        old = face.expand(loc.box(refined,prev),jump*radius)
        dg = 1e-13*(abs(prev).max()+abs(u).max())+1e-90
        delta = loc.box(refined,u-prev)
        _, cmax = loc.count_bounds((delta[0]-dg,delta[1]+dg),R+radius,jump,width)
        truth_trial = add_history(ball,htrue,kt)
        closed = (old[1][:,2]<0)&(ball[1][:,2]<0)
        _, elastic = project(ball,htrue,kt,km)
        interior = loc.prior_elastic(add_history(old,htrue,kt),truth_trial,cmax,kt,km)
        eligible = closed & (elastic|interior)
        row = dict(step=idx+1, t=t, previous_radius_m=radius,
                   uncertified_stick_children=np.unique(parent[~elastic]).tolist(),
                   uncertain_leaves=int((~elastic).sum()),max_substeps=int(cmax.max()),
                   unresolved_constitutive_leaves=int((~eligible).sum()))
        if not eligible.all():
            row.update(included=False,reason="closed single final projection not enclosed")
            steps.append(row)
            break
        center_def,_ = project(center,hnom,kt,km)
        actual_at_center,_ = project(center,htrue,kt,km)
        # At identical current j,r the two projected tractions differ by at
        # most the previous paired bound, including when either branch sticks.
        pl = np.maximum(-pair[:,None],actual_at_center[0]-center_def[1]-1e-12*(abs(actual_at_center[0])+abs(center_def[1])))
        ph = np.minimum(pair[:,None],actual_at_center[1]-center_def[0]+1e-12*(abs(actual_at_center[1])+abs(center_def[0])))
        assert (pl<=ph).all()
        os.environ["P11_TEST_WEIGHT"]="1"
        red = face.Reduction(refined,area)
        a,_ = face.first_lower(red,elastic,np.zeros(len(parent),bool),kn,kt,mu)
        b,_,_ = loc.tangent.integrate(red,truth_trial,elastic,parent,mask,kn,kt,mu)
        lower = max(a,b)
        f,err,pf,pe,guard,fc,fp = assemble(red,u,phase,t,center_def,(pl,ph),kn,kt,p13["guard_kn"])
        center_only = loc.inverse.check(red,refined,elastic,f,err,guard,kn,kt,mu,R,out/f"{sign}-{idx+1}-center.npz")
        inv = loc.inverse.check(red,refined,elastic,f+pf,err+pe,guard,kn,kt,mu,R,out/f"{sign}-{idx+1}-residual.npz")
        upper = (fc[1]+fp[1])*(1+1e-12)
        bounds = [upper/lower*(1+1e-12)] if upper<lower*R else []
        if inv["included"]:
            bounds.append(inv["root_radius_upper_m"])
        row.update(center_residual_kn=fc,past_difference_kn=fp,lower_kn_m=lower,
                   lower_times_radius_kn=lower*R,center_only_inverse_test=center_only,inverse_test=inv,included=bool(bounds),
                   certified_radius_m=min(bounds) if bounds else None)
        sel=mask[parent]
        payload=dict(parent=parent[sel],past_nominal_kpa=np.stack(hnom)[:,sel],
                     past_enclosed_kpa=np.stack(htrue)[:,sel],paired_previous_upper_kpa=pair[sel],
                     paired_at_center_kpa=np.stack([pl,ph])[:,sel],center_defect_kpa=np.stack(center_def)[:,sel],
                     leaf_stick=elastic[sel])
        steps.append(row)
        if bounds:
            radius=min(bounds)
            accepted=face.expand(center,jump*radius)
            htrue,stick_new=project(accepted,htrue,kt,km)
            hnom=center_def
            # Fixed-r projection is nonexpansive in history. Altering
            # the current displacement changes the defect by <= (kt+km)*|Dj|.
            pair=np.where(stick_new,pair,pair+(kt+km)*jump*radius)*(1+1e-12)
            payload.update(out_nominal_kpa=np.stack(hnom)[:,sel],out_enclosed_kpa=np.stack(htrue)[:,sel],
                           out_paired_upper_kpa=pair[sel])
        np.savez_compressed(out/f"{sign}-{idx+1}-history.npz",**payload)
        if not bounds:
            break
    save(out/(sign+".json"),dict(start=p13["step"],start_radius_m=p13["certified_radius_m"],steps=steps,
         consecutive=p13["step"]+sum(r["included"] for r in steps),total=len(data["clock"]),
         threshold=.001,full_history_relative_upper=None))


if __name__ == "__main__":
    run(sys.argv[1],Path(sys.argv[2]))
