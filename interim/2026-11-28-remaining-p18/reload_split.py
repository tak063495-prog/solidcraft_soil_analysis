"""Resolve signed equilibrium shares along accepted reload half-intervals."""
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P16=HERE.with_name("2026-11-28-remaining-p16")
sys.path.insert(0,str(P16))
import attribution as prior

FIELDS=("u","jump","reaction")
POINTS=(1,25,21)


def run(raw,out):
    op,L=prior.operators()
    by_route={}
    grid=np.arange(1,257)/256
    all_times=[]
    for route,folder in (("outer",prior.P15),("middle",prior.P14)):
        levels=[]
        for n in prior.LEVELS:
            with np.load(raw/f"{route}-{n}.npz") as f:
                v=dict(f)
            t=v["t"]
            archived=prior.restore(folder/"c01/alternate"/route,n,[*FIELDS,"traction"])
            v["reaction"]=np.array([(op.k@u.ravel()+op.seed+
                np.einsum("qag,qa,q->g",op.system.c,q,op.system.weight)-load.ravel()).reshape(-1,3)
                for u,q,load in zip(v["u"],v["traction"],v["load"],strict=True)])
            indices=np.searchsorted(t,grid)
            np.testing.assert_array_equal(t[indices],grid)
            for k in archived:
                np.testing.assert_array_equal(v[k][indices],archived[k][-512:-256])
                np.testing.assert_array_equal(v[k][0],archived[k][-513])
            levels.append(v)
            all_times.append(t)
        by_route[route]=levels

    common=all_times[0]
    for t in all_times[1:]:
        common=np.intersect1d(common,t)
    assert common[0]==0 and common[-1]==1 and np.isin(grid,common).all()
    saved={}
    for route,levels in by_route.items():
        indices=[np.searchsorted(v["t"],common) for v in levels]
        def delta(k):
            return -np.diff(np.stack([v[k][ix] for v,ix in zip(levels,indices,strict=True)]),axis=0)
        n=levels[0]["normals"]
        for v in levels:
            np.testing.assert_array_equal(n,v["normals"])
        np.testing.assert_allclose(np.linalg.norm(n,axis=1),1,rtol=1e-15)
        dt=delta("traction")
        normal=(dt*n).sum(-1)[...,None]*n
        components=np.stack([normal,dt-normal],axis=-1)
        rt=dict(progress=common,normals=n,traction_difference=dt)
        result={}
        # For every matched endpoint give the precise enclosing accepted half in each level.
        for level,v,ix in zip(prior.LEVELS,levels,indices,strict=True):
            rt[f"accepted_spans_{level}"]=v["spans"]
            rt[f"half_{level}"]=np.column_stack((v["t"][np.maximum(ix-1,0)],v["t"][ix]))
            # Native traction changes over the immediately preceding accepted half.
            # Cross-level growth itself is evaluated only at shared load coordinates.
            traction = v["traction"]
            inc = traction[ix]-traction[np.maximum(ix-1,0)]
            rt[f"native_last_half_tangent_change_{level}"] = inc - (inc*n).sum(-1)[...,None]*n
        gi=np.searchsorted(common,grid)
        with np.load(P16/f"c01/alternate/{route}-1792.npz") as previous:
            for k in FIELDS:
                d=delta(k).reshape(2,len(common),-1)
                operator=L[k].reshape(96,-1)
                proj=(d@operator.T).reshape(2,len(common),32,3)
                shares=np.einsum("ltpc,ltpcv->ltpv",proj,components)
                growth=shares[1]-shares[0]
                sq=(d*d).sum(-1)
                remainder=(d*(d-dt.reshape(2,len(common),96)@operator)).sum(-1)
                np.testing.assert_allclose(shares.sum((2,3))+remainder,sq,rtol=1e-10,atol=1e-40)
                np.testing.assert_allclose(growth[gi].sum(-1),previous[k+"_point_growth"][-512:-256],
                    atol=1e-35 if k!="reaction" else 1e-22,rtol=1e-9)
                rt[k+"_difference"]=d
                rt[k+"_growth"]=growth
                rt[k+"_residual"]=remainder[1]-remainder[0]
                rt[k+"_squared_norms"]=sq
                res={}
                for i in POINTS:
                    res[i]={}
                    for label,series in (("normal",growth[:,i,0]),("tangent",growth[:,i,1]),
                                         ("sum",growth[:,i].sum(-1))):
                        def interval(j):
                            return common[j-1:j+1].tolist()
                        pos=np.flatnonzero(series[1:]>0)+1
                        inc=np.flatnonzero(np.diff(series)>0)+1
                        cross=np.flatnonzero((series[1:]>0)&(series[:-1]<=0))+1
                        res[i][label]=dict(at_start=series[0],first_value=series[1],
                            first_positive=interval(pos[0]) if len(pos) else None,
                            first_increase=interval(inc[0]) if len(inc) else None,
                            positive_crossings=[interval(j) for j in cross],
                            largest_increase=interval(np.argmax(np.diff(series))+1),
                            max_increase=np.diff(series).max(),grid256_sum=series[gi].sum())
                result[k]=res
        np.savez_compressed(out/f"{route}-reload.npz",**rt)
        with np.load(HERE.with_name("2026-11-28-remaining-p17")/f"c01/alternate/{route}-return.npz") as ret:
            for k in FIELDS:
                np.testing.assert_array_equal(rt[k+"_difference"][:,-1],ret[k+"_difference"][:,0])
                np.testing.assert_allclose(rt[k+"_growth"][-1],ret[k+"_growth"][0],rtol=1e-11,atol=1e-35)
                for point in POINTS:
                    series = rt[k+"_growth"][:,point,1]
                    j = np.flatnonzero(series<=0)
                    ret_series = ret[k+"_growth"][:,point,1]
                    peak = np.argmax(np.diff(ret_series))+1
                    result[k][point]["join"] = dict(
                        reload_last_half={str(n):rt[f"half_{n}"][-1].tolist() for n in prior.LEVELS},
                        penultimate=float(series[-2]), at_reload_end=float(series[-1]),
                        at_return_start=float(ret_series[0]), at_return_first=float(ret_series[1]),
                        last_nonpositive_then_positive=common[j[-1]:j[-1]+2].tolist() if len(j) and j[-1]+1<len(common) else None,
                        return_max_span=ret["progress"][peak-1:peak+1].tolist(),
                        return_max_increase=float(ret_series[peak]-ret_series[peak-1]))
        saved[route]=result
    (out/"summary.json").write_text(json.dumps(dict(common_endpoints=len(common),reload_gate_endpoints=256,
        points=POINTS,components=["normal","tangent"],routes=saved),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]))
