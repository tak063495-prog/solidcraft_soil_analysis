"""Attribute point 7 and 9 endpoint products to the accepted recontact half steps."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np
from plans import P11, P12, save
sys.path.insert(0,str(P12))
from alternate_plan import halves, LEVELS


def run(out):
    folders={"p11":P11/"c01/alternate","p12":P12/"c01/alternate",
             "point14":out/"point14","point26":out/"point26"}
    paths={}
    grids=[]
    for route,src in folders.items():
        levels=[]
        for n in LEVELS:
            data=[(i,h,l,r,np.array([x[k][1] for k in (7,9)]),np.array([y[k][1] for k in (7,9)]))
                  for i,h,l,r,x,y in halves(src/str(n)/"trace.jsonl.gz")]
            clock={data[0][2]:data[0][4]}
            for _,_,l,r,x,y in data:
                np.testing.assert_array_equal(clock[l],x)
                clock[r]=y
            grids.append(set(clock))
            levels.append((data,clock))
        paths[route]=levels
    common=sorted(set.intersection(*grids))
    end=common[-1]
    results={}
    with gzip.open(out/"signed-halves-7-9.jsonl.gz","xt") as f:
        for route,levels in paths.items():
            A,B=-np.diff([clock[end] for _,clock in levels],axis=0)
            weights=[B/2,(A-B)/2,-A/2]
            initial_a,initial_b=-np.diff([clock[0.] for _,clock in levels],axis=0)
            initial=((initial_a*B+initial_b*A)/2).sum(-1)
            endpoint=(A*B).sum(-1)
            bins=np.zeros((len(common)-1,2))
            for n,(path,_),w in zip(LEVELS,levels,weights,strict=True):
                for i,h,l,r,x,y in path:
                    if r>end:
                        break
                    j=np.searchsorted(common,r,side="left")-1
                    assert common[j]<=l<r<=common[j+1]
                    term=((y-x)*w).sum(-1)
                    bins[j]+=term
                    print(json.dumps(dict(path=route,level=n,attempt=i,half=h,span=[l,r],
                        delta_m=(y-x).tolist(),weight_m=w.tolist(),term_m2=term.tolist()),allow_nan=False),file=f)
            np.testing.assert_allclose(initial+bins.sum(0),endpoint,rtol=1e-8,atol=1e-28)
            results[route]=dict(initial_m2=initial.tolist(),endpoint_m2=endpoint.tolist(),terms_m2=bins.tolist())
    effects={}
    orig=np.array(results["p11"]["terms_m2"])
    for route in ("point14","point26","p12"):
        delta=np.array(results[route]["terms_m2"])-orig
        worst=delta.argmin(0)
        effects[route]=[dict(point=k,span=common[j:j+2],change_m2=float(delta[j,col])) for col,(k,j) in enumerate(zip((7,9),worst,strict=True))]
    save(out/"intervals-7-9.json",dict(points=[7,9],endpoint_progress=end,common=common,routes=results,
         worst_change_from_p11=effects))


if __name__=="__main__":
    run(Path(sys.argv[1]))
