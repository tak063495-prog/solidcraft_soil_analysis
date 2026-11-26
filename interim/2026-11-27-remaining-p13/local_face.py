"""Integrate closed-face spring and projected branches only on the p12 stop children."""
from collections import Counter
import json
import importlib.util
import os
from pathlib import Path
import sys

import numpy as np

from plans import P12, save
import tangent
import inverse
sys.path.insert(0, str(P12))
spec=importlib.util.spec_from_file_location("p12_center",P12/"face.py")
face=importlib.util.module_from_spec(spec)
spec.loader.exec_module(face)

EPS = face.EPS
RADIUS = 1e-8
DEPTH = 4


def subdivision():
    # Blossom restriction from the parent Bernstein basis to each of its four children.
    alpha = face.base.prev.ALPHA
    unit = np.eye(3)
    a, b, c = unit
    ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
    matrices = []
    for vertices in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)):
        rows = []
        for powers in alpha:
            weights = {(0, 0, 0): 1.}
            for vertex in np.repeat(vertices, powers, axis=0):
                child = Counter()
                for key, val in weights.items():
                    for j in np.flatnonzero(vertex):
                        k = list(key)
                        k[j] += 1
                        child[tuple(k)] += val * vertex[j]
                weights = child
            rows.append([weights.get(k, 0.) for k in alpha])
        matrices.append(rows)
    matrices = np.array(matrices)
    assert np.all(matrices >= 0) and np.all(matrices.sum(-1) == 1)
    return matrices


def refine(original, mask):
    mats = subdivision()
    refined, parents = [], []
    offset = 0
    for ids, w, _, sign, areas, _ in original:
        weights, sizes = [], []
        for i, (v, area) in enumerate(zip(w.reshape(-1,15,15), areas, strict=True)):
            cells = v[None]
            if mask[offset+i]:
                for _ in range(DEPTH):
                    cells = (mats[:,None] @ cells[None]).reshape(-1,15,15)
            weights.append(cells)
            sizes.extend([area/len(cells)]*len(cells))
            parents.extend([offset+i]*len(cells))
        W = np.concatenate(weights)
        gain2 = 2 * np.sum(abs(W).max(1)**2,axis=1) * (1+1e-12)
        refined.append((ids,W.reshape(-1,15),abs(W).sum(2).max()*(1+1e-12),sign,np.array(sizes),gain2))
        offset += len(areas)
    return refined, np.array(parents)


def box(refined, u):
    lo,hi=face.enclosure(refined,u)
    guard=1e-12*np.max(abs(u))+1e-90
    return lo-guard, hi+guard


def count_bounds(delta, radius, jump, width):
    l,h=face.expand(delta,jump*radius)
    near=np.maximum(np.maximum(l,-h),0.)
    far=np.maximum(-l,h)
    return (np.maximum(1,np.ceil(np.linalg.norm(near,axis=1)*(1-1e-12)/width)).astype(int),
            np.maximum(1,np.ceil(face.norm_hi(far)*(1+1e-12)/width)).astype(int))


def prior_elastic(old, end, m, kt, knmu):
    # Both the old point and this last possible interior point are admissible.
    # Convexity of the Coulomb cone then covers every earlier material substep.
    s = (1.-1./m)[:,None]
    l = (1-s)*old[0]+s*end[0]
    h = (1-s)*old[1]+s*end[1]
    pad=1e-13*np.maximum.reduce([abs(old[0]),abs(old[1]),abs(end[0]),abs(end[1])])+1e-90
    q=kt*face.norm_hi(np.maximum(-l[:,:2],h[:,:2])+pad[:,:2])
    return (m==1) | (q*(1+1e-12) < -knmu*(h[:,2]+pad[:,2])*(1-1e-12))


def defect(end, kt, knmu):
    # P_r(v)-v = -v * max(0,1-r/|v|), with signed component bounds.
    l,h=end
    vl,vh=l[:,:2]*kt,h[:,:2]*kt
    qmax=face.norm_hi(np.maximum(-vl,vh))
    qmin=np.linalg.norm(np.maximum(np.maximum(vl,-vh),0),axis=1)*(1-1e-12)
    low=np.maximum(0,1-(-knmu*l[:,2])*(1+1e-12)/np.maximum(qmin,1e-90)-1e-12)
    high=np.maximum(0,1-(-knmu*h[:,2])*(1-1e-12)/np.maximum(qmax,1e-90)+1e-12)
    low=low[:,None]
    high=high[:,None]
    a,b=face.base.prev.scale(-vh,-vl,low,high)
    pad=1e-12*np.maximum(abs(a),abs(b))+1e-90
    return a-pad,b+pad


def run(sign, out):
    original, reds, gains, law, data=face.setup(sign)
    rows=[json.loads(s) for s in (P12/"c02"/(sign+".jsonl")).read_text().splitlines()]
    last=rows[-1]
    for r in rows[:-1]:
        assert any(v["included"] and v["elastic"]==6080 for v in r["tests"])
    i=last["step"]-1
    u,prev=data["u"][i],data["u"][i-1]
    phase,t,width=data["clock"][i]
    r0=last["previous_radius_m"]
    kn,kt,mu=[law[k] for k in ("normal_stiffness_kpa_m","tangential_stiffness_kpa_m","friction_coefficient")]
    o=face.expand(face.enclosure(original,prev),gains*r0)
    e=face.expand(face.enclosure(original,u),gains*RADIUS)
    _,stick,op=face.defects(o,e,np.zeros(len(gains)),width,False,kn,kt,mu)
    target=~stick
    assert not op.any() and target.sum()=={"plus":12,"minus":2}[sign]
    refined, parent=refine(original,target)
    area=np.concatenate([r[-2] for r in refined])
    jump=np.sqrt(np.concatenate([r[-1] for r in refined]))*(1+1e-12)
    center=box(refined,u)
    old=face.expand(box(refined,prev),jump*r0)
    ball=face.expand(center,jump*RADIUS)
    assert np.all(old[1][:,2]<0) and np.all(ball[1][:,2]<0)
    delta=box(refined,u-prev)
    # Subtraction of the nodal displacement vectors precedes the polynomial evaluation.
    dg=1e-13*(abs(u).max()+abs(prev).max())+1e-90
    delta=delta[0]-dg,delta[1]+dg
    cmin,cmax=count_bounds(delta,RADIUS+r0,jump,width)
    elastic=kt*face.norm_hi(np.maximum(-ball[0][:,:2],ball[1][:,:2]))*(1+1e-12)<-mu*kn*ball[1][:,2]*(1-1e-12)
    assert np.all(elastic[~target[parent]])
    simple=elastic | prior_elastic(old,ball,cmax,kt,kn*mu)
    assert simple.all(), (int((~simple).sum()),int(cmax.max()))
    # On these cells nominal and every enclosed past give the identical fixed-center
    # projection. Their paired difference is exactly zero, in both signed directions.
    dl,dh=defect(center,kt,kn*mu)
    os.environ["P11_TEST_WEIGHT"]="1"
    red=face.Reduction(refined,area)
    lower,margins=face.first_lower(red,elastic,np.zeros(len(area),bool),kn,kt,mu)
    direction_lower,direction_margins,integrals=tangent.integrate(red,ball,elastic,parent,target,kn,kt,mu)
    lower=max(lower,direction_lower)
    g=red.g
    force=(g.k@u+g.initial-g.load-((1-t)*(g.initial-g.load) if phase==0 else 0.)).reshape(-1,3)
    uncertainty=np.zeros_like(force)
    scale=np.linalg.norm(abs(g.k)@abs(u)+2*(abs(g.initial)+abs(g.load)))
    child_records=[]
    for sl,pos,ids,sgn,mass,mean,amax in red.blocks:
        jump_u=u.reshape(-1,3)[ids[1]]-u.reshape(-1,3)[ids[0]]
        M=mass.sum(axis=0).reshape(15,15)
        spring=M@(jump_u*[kt,kt,kn])
        scale+=np.linalg.norm(abs(M)@(abs(jump_u)*[kt,kt,kn]))
        correction=mean.T@((dh[sl]+dl[sl])/2)
        err=amax.T@((dh[sl]-dl[sl])/2)
        spring[:,:2]+=correction
        for nodes,side in zip(ids,(-1,1),strict=True):
            force[nodes]+=side*spring
            uncertainty[nodes,:2]+=err
        for k in np.unique(parent[sl][target[parent[sl]]]):
            sel=parent[sl]==k
            lo=dl[sl][sel]
            hi=dh[sl][sel]
            child_records.append(dict(child=int(k),leaves=int(sel.sum()),area=float(area[sl][sel].sum()),
                count_min=int(cmin[sl][sel].min()),count_max=int(cmax[sl][sel].max()),
                elastic_leaves=int(elastic[sl][sel].sum()),center_defect_kpa=[lo.min(0).tolist(),hi.max(0).tolist()],
                past_difference_kpa=[[0.,0.],[0.,0.]],
                signed_nodal_center_kn=(mean[sel].T@((lo+hi)/2)).tolist(),
                nodal_radius_kn=(amax[sel].T@((hi-lo)/2)).tolist(),
                signed_nodal_past_kn=np.zeros((15,2)).tolist(),
                integrated_sliding_lower_gram=(mass[sel & ~elastic[sl]].sum(0).reshape(15,15)).tolist()))
    guard=max(last["assembly_guard_kn"],float(8192*EPS*(1+scale)))
    upper=float(np.linalg.norm((abs(force)+uncertainty).ravel()[red.free])*(1+1e-10)+guard)
    inner=float(max(0,np.linalg.norm(np.maximum(abs(force)-uncertainty,0.).ravel()[red.free])*(1-1e-10)-guard))
    direct=bool(upper<lower*RADIUS)
    inv=inverse.check(red,refined,elastic,force,uncertainty,guard,kn,kt,mu,RADIUS,out/(sign+"-residual.npz"))
    bounds=([upper/lower*(1+1e-12)] if direct else [])+([inv["root_radius_upper_m"]] if inv["included"] else [])
    ok=bool(bounds)
    result=dict(step=i+1,phase=int(phase),t=t,radius_m=RADIUS,previous_radius_m=r0,
        original_children=6080,target_children=int(target.sum()),depth=DEPTH,uncertain_refined=int((~elastic).sum()),
        signed_children=child_records,isotropic_schur_margins=margins,directional_schur_margins=direction_margins,
        child_directional_tangent_lower_integrals_kn_m=integrals,lower_kn_m=lower,
        center_residual_interval_kn=[inner,upper],paired_past_upper_kn=0.,guard_kn=guard,
        lower_times_radius_kn=lower*RADIUS,unscaled_included=direct,inverse_test=inv,included=ok,certified_radius_m=min(bounds) if bounds else None,
        consecutive=i+int(ok),total_increments=len(data["clock"]),threshold=.001,
        reequilibrated_full_history_relative_upper=None)
    save(out/(sign+".json"),result)


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
