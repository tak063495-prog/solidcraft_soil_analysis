"""Score the six measured histories on their intersection of accepted load coordinates."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np
from plans import P11, P12, P13, ROUTES, save

sys.path.insert(0,str(P11))
from alternate_score import score, STAGES
sys.path.insert(0,str(P12))
from alternate_plan import halves, LEVELS


def accepted(path):
    rows = []
    clock = {}
    for i,h,l,r,x,y in halves(path):
        x,y = (np.array([v[1] for v in states]) for states in (x,y))
        if l in clock:
            np.testing.assert_array_equal(clock[l],x)
        else:
            assert not clock and l == 0.
            clock[l] = x
        clock[r] = y
        rows.append((i,h,l,r,(y-x)[[7,9]]))
    return rows,clock


def score_route(work,out,route):
    dest = out/route
    score(work/route,dest)
    with np.load(dest/"signed-32.npz") as f:
        v,progress = f["p11"],f["progress"]
    np.savez_compressed(dest/"signed-32.npz",progress=progress,**{route:v})
    generic = json.loads((dest/"contributions.json").read_text())
    save(dest/"contributions.json",{s:r["p11"] for s,r in generic.items()})


def summary(rows):
    return dict(count=len(rows),sum_m2=float(rows.sum()),by_point_m2=rows.sum(0).tolist(),
                last_progress=len(rows)/256,endpoint_total_m2=float(rows[-1].sum()),endpoint_m2=rows[-1].tolist())


def compare(out):
    folders = dict(p11=P11/"c01/alternate",p12=P12/"c01/alternate",
                   point14=P13/"c01/alternate/point14",point26=P13/"c01/alternate/point26")
    fixed = {}
    folders.update({name:out/name for name in ROUTES})
    for name,src in folders.items():
        source = folders["p12"] if name == "p11" else src
        with np.load(source/"signed-32.npz") as f:
            fixed[name] = f[name]
            np.testing.assert_array_equal(f["progress"],np.tile(np.arange(1,257)/256,7)[:len(f[name])])
    paths = {name:[accepted(src/str(n)/"trace.jsonl.gz") for n in LEVELS] for name,src in folders.items()}
    count = min(map(len,fixed.values()))
    stages = {}
    completed_stages = {}
    for j,stage in enumerate(STAGES):
        stages_data = {}
        completed_stages[stage] = {}
        for name,v in fixed.items():
            start = j*256
            x = v[start:min(count,start+256)]
            if len(x):
                stages_data[name] = summary(x)
            full = v[start:start+256]
            completed_stages[stage][name] = summary(full) if len(full)==256 else None
        if stages_data:
            stages[stage] = stages_data
    # Also retain the denser shared recontact end loads (exclude half-only events).
    end_grid = sorted(set.intersection(*(set(r for _,h,_,r,_ in rows if h==2)
                                         for levels in paths.values() for rows,_ in levels)))
    with gzip.open(out/"common-recontact.jsonl.gz","xt") as stream:
        for t in end_grid:
            item = {"t":t}
            for name,levels in paths.items():
                a,b = -np.diff([clock[t] for _,clock in levels],axis=0)
                item[name] = (a*b).sum(1).tolist()
            print(json.dumps(item),file=stream)
    shared_effects = attribute(paths,end_grid[-1],out)
    full_paths = {name:levels for name,levels in paths.items() if all(1. in clock for _,clock in levels)}
    full_dir = out/"at-recontact-end"
    full_dir.mkdir(exist_ok=True)
    full_effects = attribute(full_paths,1.,full_dir) if full_paths else {}
    save(out/"comparison.json",dict(common_coordinates=count,common_recontact_ends=len(end_grid),
         matched_endpoint=end_grid[-1],stages=stages,p11_differences=shared_effects,
         full_stage_comparison=completed_stages,full_endpoint_p11_differences=full_effects))


def attribute(paths,endpoint,out):
    grid = sorted(t for t in set.intersection(*(set(clock) for levels in paths.values() for _,clock in levels))
                  if t<=endpoint)
    assert grid[0]==0 and grid[-1]==endpoint
    terms = {}
    with gzip.open(out/"signed-halves-7-9.jsonl.gz","xt") as stream:
        for name,levels in paths.items():
            a,b = -np.diff([c[endpoint][[7,9]] for _,c in levels],axis=0)
            a0,b0 = -np.diff([c[0.][[7,9]] for _,c in levels],axis=0)
            initial = ((a0*b+b0*a)/2).sum(1)
            weights = [b/2,(a-b)/2,-a/2]
            bins = np.zeros((len(grid)-1,2))
            for n,(rows,_),w in zip(LEVELS,levels,weights,strict=True):
                for i,h,l,r,delta in rows:
                    if r>endpoint:
                        break
                    k = np.searchsorted(grid,r)-1
                    assert grid[k]<=l<r<=grid[k+1]
                    v = (delta*w).sum(1)
                    bins[k] += v
                    print(json.dumps(dict(path=name,level=n,attempt=i,half=h,span=[l,r],
                                delta_m=delta.tolist(),weight_m=w.tolist(),term_m2=v.tolist())),file=stream)
            final = (a*b).sum(1)
            np.testing.assert_allclose(initial+bins.sum(0),final,rtol=1e-8,atol=1e-28)
            terms[name] = dict(initial_m2=initial.tolist(),endpoint_m2=final.tolist(),terms_m2=bins.tolist())
    effects = {}
    for name in paths:
        d = np.array(terms[name]["terms_m2"])-terms["p11"]["terms_m2"]
        item = dict(endpoint_change_m2=(np.array(terms[name]["endpoint_m2"])-terms["p11"]["endpoint_m2"]).tolist())
        for idx,(pt,l,r) in enumerate([(7,.9296875,.9307782966271231),(9,.92578125,.9286986368319801)]):
            k = d[:,idx].argmin()
            item[str(pt)] = dict(p13_span=[l,r],change_m2=float(d[grid.index(l):grid.index(r),idx].sum()),
                                 worst_span=grid[k:k+2],worst_m2=float(d[k,idx]))
        effects[name] = item
    save(out/"intervals-7-9.json",dict(points=[7,9],common=grid,endpoint_progress=endpoint,routes=terms))
    return effects


if __name__ == "__main__":
    work,out=map(Path,sys.argv[1:])
    for route in ROUTES:
        score_route(work,out,route)
    compare(out)
