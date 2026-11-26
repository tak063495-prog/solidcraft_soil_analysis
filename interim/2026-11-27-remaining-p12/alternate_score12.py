"""Score native fields with the unchanged p11 routine, compare to p11 at matched loads."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

from alternate_plan import BASE, HERE, LEVELS, halves
sys.path.insert(0,str(BASE.parents[1]))
from alternate_score import score, save, STAGES


def run(raw,out):
    score(raw,out)  # identical 22-field decision and reversible XOR against archived p10 bytes
    with np.load(BASE/"signed-32.npz") as old, np.load(out/"signed-32.npz") as new:
        progress=new["progress"].copy()
        n=len(progress)
        np.testing.assert_array_equal(progress,old["progress"][:n])
        before=old["p11"][:n].copy()
        after=new["p11"].copy()
    np.savez_compressed(out/"signed-32.npz",progress=progress,p11=before,p12=after)
    report={}
    for j,s in enumerate(STAGES):
        sl=slice(j*256,(j+1)*256)
        if not len(progress[sl]):
            break
        report[s]={}
        for label,v in (("p11",before),("p12",after)):
            x=v[sl]
            report[s][label]=dict(coordinates=len(x),sum_m2=float(x.sum()),by_point_m2=x.sum(0).tolist(),
                                  endpoint=x[-1].tolist(),endpoint_total=float(x[-1].sum()),
                                  endpoint_negative=int((x[-1]<0).sum()),
                                  endpoint_other_30=float(np.delete(x[-1],[14,26]).sum()))
    save(out/"contributions.json",report)
    grids={}
    for label,base in (("p11",BASE),("p12",out)):
        grids[label]=[{r:y for _,h,_,r,_,y in halves(base/str(n)/"trace.jsonl.gz") if h==2} for n in LEVELS]
    plan=json.loads((out/"plan.json").read_text())
    shared={k:set.intersection(*(set(path) for path in v)) for k,v in grids.items()}
    assert all(t in shared["p12"] for t in plan["knots"] if t<=min(max(p) for p in grids["p12"]))
    with gzip.open(out/"matched-32.jsonl.gz","xt") as stream:
        for t in sorted(set.union(*shared.values())):
            row=dict(t=t)
            for label,paths in grids.items():
                if t in shared[label]:
                    a,b=-np.diff([[x[1] for x in p[t]] for p in paths],axis=0)
                    row[label]=(a*b).sum(-1).tolist()
            print(json.dumps(row,allow_nan=False),file=stream)


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]))
