"""All 32 signed point contributions at actual endpoints of all three runs."""
import json
from pathlib import Path
import sys

import numpy as np
from alt_plan import collect, read, P10

def run(out):
    plan=json.loads((out/"plan.json").read_text())
    grids=[]
    original=[]
    for n in (256,512,1024):
        rows=read(out/str(n)/"trace.jsonl.gz")
        candidates={(r["p11_span"][0],r["p11_span"][2]):r["all"][-1] for r in rows if "p11_span" in r}
        grid={r["to"]:candidates[r["from"],r["to"]] for r in rows if r.get("stage")=="recontact" and r["accepted"]}
        if max(grid)>=max(plan["knots"]):
            assert set(plan["knots"])<=grid.keys()
        grids.append(grid)
        original.append({r:after for _,h,_,r,_,after in collect(out,n) if h==2})
    shared=sorted(set.intersection(*(set(g) for g in grids)))
    matches=[]
    for t in shared:
        record=dict(t=t)
        for name,paths in (("p11",grids),("p10",original)):
            if not all(t in g for g in paths):
                continue
            v=np.array([[p[1] for p in g[t]] for g in paths])
            a,b=-np.diff(v,axis=0)
            record[name]=dict(dot=(a*b).sum(-1).tolist(),status=[[p[0] for p in g[t]] for g in paths])
        matches.append(record)
    with (out/"matched-32.jsonl").open("w") as stream:
        for row in matches:
            print(json.dumps(row,allow_nan=False),file=stream)
if __name__=="__main__":
    run(Path(sys.argv[1]))
