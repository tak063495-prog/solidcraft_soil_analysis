"""Score archived and matched reference histories using disjoint physical bands."""
import json
from pathlib import Path
import numpy as np
from freeze import HERE, P7, save, digest

NATIVE = HERE.parent/"2026-11-25-remaining-p3/c02/p4"


def masks(position):
    x,y=position[:,:2].T
    distances=np.array([x,1-x,y,1-y]).T
    near=distances.argmin(axis=1)
    d=distances.min(axis=1)
    mask={"all":np.ones(len(x),dtype=bool)}
    for side,name in enumerate(("west","east","south","north")):
        mask[name]=(near==side)&(d<.1)
    mask["interior"]=d>=.1
    for lo,hi in ((0,.025),(.025,.05),(.05,.1)):
        mask[f"d[{lo},{hi})"]=(d>=lo)&(d<hi)
    return mask


def score(a,b,weight=1.):
    delta=np.linalg.norm((a-b)*weight)
    denom=np.linalg.norm(b*weight)
    return dict(numerator=float(delta), denominator=float(denom),
                relative=float(delta/denom) if denom else None)


def arrays(file):
    with np.load(file) as z:
        return {k:z[k] for k in z.files}


def integral(slip,w,sel):
    return np.einsum("tqc,q->tc",slip[:,sel],w[sel])


def compare(n,na,r,ra):
    for k in ("position","parent","weight"):
        np.testing.assert_array_equal(n[k],r[k])
    rm, nm, am = masks(r["position"]), masks(na["position"]), masks(ra["position"])
    rows={}
    histories={}
    for key, mask in rm.items():
        field=score(n["slip"][:,mask],r["slip"][:,mask],np.sqrt(r["weight"][mask])[None,:,None])
        ni=integral(na["slip"],na["weight"],nm[key])
        ri=integral(ra["slip"],ra["weight"],am[key])
        nc=integral(n["slip"],n["weight"],mask)
        rc=integral(r["slip"],r["weight"],mask)
        rows[key]=dict(common_points=int(mask.sum()),common_area=r["weight"][mask].sum(),
            native_active_points=int(nm[key].sum()),reference_active_points=int(am[key].sum()),
            native_sampled_area=na["weight"][nm[key]].sum(),reference_sampled_area=ra["weight"][am[key]].sum(),
            slip=field,own_active_signed_integral=score(ni,ri),common_signed_integral=score(nc,rc))
        histories[key]=np.stack((ni,ri,nc,rc))
    for j,original in enumerate((n,r)):
        np.testing.assert_allclose(histories["all"][j], original["signed_area_integral"],atol=1e-17,rtol=1e-12)
    parts=["west","east","south","north","interior"]
    np.testing.assert_allclose(sum(histories[k] for k in parts),histories["all"],atol=1e-17)
    return rows,histories


def summarize_residual(folder):
    v=arrays(folder/"residual.npz")
    out={}
    for band, mask in masks(v["centroid"]).items():
        rows={}
        for label,choice in (("all",np.ones_like(v["mixed"])),("mixed_at_anchor",v["mixed"]),
                             ("switch_since_last_anchor",v["interval_switch"])):
            chosen=choice & mask[None,:]
            err=np.where(chosen[:,:,None],v["signed"],0.)
            child=np.where(chosen[:,:,None],v["child"],0.)
            # Time-L2 of the sum of per-leaf magnitudes exposes cancellation.
            numerator=float(np.linalg.norm(err.sum(axis=1)))
            denominator=float(np.linalg.norm(child.sum(axis=1)))
            rows[label]=dict(leaf_anchors=int(chosen.sum()),unique_leaves=int(chosen.any(axis=0).sum()),
                total_leaf_anchors=int(mask.sum()*len(choice)),
                signed_resultant=dict(numerator=numerator,denominator=denominator,
                    relative=numerator/denominator if denominator else None),
                l2_time_sum_local_norms=float(np.linalg.norm(np.linalg.norm(err,axis=2).sum(axis=1))),
                max_leaf_residual=float(np.linalg.norm(err,axis=2).max()))
        out[band]=rows
    return out


def run(native_active,root):
    n=arrays(NATIVE/"principal.npz")
    na=arrays(native_active)
    report={}
    for name in ("plus","minus"):
        folder=root/name
        r,ra=(arrays(folder/f) for f in ("principal.npz","slip.npz"))
        rows,hist=compare(n,na,r,ra)
        np.savez_compressed(folder/"band-integrals.npz",**hist)
        report[name]=dict(bands=rows,residual=summarize_residual(folder),
            whole_slip_within=rows["all"]["slip"]["relative"]<=.001,
            whole_area_within=rows["all"]["own_active_signed_integral"]["relative"]<=.001)
    old=arrays(P7/"c02/principal.npz")
    pos=arrays(root/"plus/principal.npz")
    report["plus_replay_max_abs"]={k:float(np.max(np.abs(pos[k]-old[k]))) for k in ("slip","signed_area_integral")}
    report["source_hashes"]={str(p.name):digest(p) for p in (NATIVE/"principal.npz",native_active,NATIVE/"increments.npz")}
    save(root/"comparison.json",report)


if __name__=="__main__":
    import sys
    run(Path(sys.argv[1]),Path(sys.argv[2]))
