"""Compare four paths at the same accepted loads using the archived p11 decomposition."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np
from plans import P11, P12, save
from signed_halves import run as attribute_halves

sys.path.insert(0,str(P11))
from alternate_score import score, STAGES
sys.path.insert(0,str(P12))
from alternate_plan import halves, LEVELS


def run(work, out):
    ref=P12/"c01/alternate"
    with np.load(ref/"signed-32.npz") as f:
        progress=f["progress"]
        dots={name:f[name] for name in ("p11","p12")}
    clocks={}
    gates={}
    for point in (14,26):
        label=f"point{point}"
        dest=out/label
        score(work/label,dest)
        with np.load(dest/"signed-32.npz") as f:
            v=dots[label]=f["p11"]
            np.testing.assert_array_equal(f["progress"],progress[:len(v)])
        # Replace p11's generic output keys with the name of the measured route.
        np.savez_compressed(dest/"signed-32.npz", progress=progress[:len(v)], **{label:v})
        legacy=json.loads((dest/"contributions.json").read_text())
        save(dest/"contributions.json",{s:x["p11"] for s,x in legacy.items()})
        gates[label]=json.loads((dest/"score.json").read_text())["passed_fields"]
        levels=[]
        for n in LEVELS:
            levels.append({r:np.array([x[1] for x in y]) for _,half,_,r,_,y in halves(dest/str(n)/"trace.jsonl.gz") if half==2})
        common=set.intersection(*(set(g) for g in levels))
        knots=json.loads((dest/"plan.json").read_text())["knots"]
        assert all(k in common for k in knots if k<=min(max(g) for g in levels))
        clock={}
        for t in common:
            a,b=-np.diff([g[t] for g in levels],axis=0)
            clock[t]=(a*b).sum(1)
        clocks[label]=clock
    with gzip.open(ref/"matched-32.jsonl.gz","rt") as f:
        old={r["t"]:r for line in f if "p11" in (r:=json.loads(line)) and "p12" in r}
    shared=sorted(set.intersection(set(old),*(set(g) for g in clocks.values())))
    with gzip.open(out/"common-recontact.jsonl.gz","xt") as f:
        for t in shared:
            row=dict(old[t], **{label:g[t].tolist() for label,g in clocks.items()})
            print(json.dumps(row, allow_nan=False),file=f)
    count=min(map(len,dots.values()))
    summary={}
    for i,s in enumerate(STAGES):
        sl=slice(i*256,min(count,(i+1)*256))
        if sl.start>=sl.stop:
            break
        summary[s]={}
        for label,all_values in dots.items():
            x=all_values[sl]
            summary[s][label]=dict(sum_m2=float(x.sum()),endpoint_m2=float(x[-1].sum()),
                endpoint=x[-1].tolist(),by_point=x.sum(0).tolist(),count=len(x),negative=int((x[-1]<0).sum()))
    e={k:np.array(r["endpoint"]) for k,r in summary["recontact"].items()}
    attribution={str(j):dict(p14_alone=float(e["point14"][j]-e["p11"][j]),
        p26_alone=float(e["point26"][j]-e["p11"][j]),
        interaction=float((e["p12"]-e["point14"]-e["point26"]+e["p11"])[j]))
        for j in (7,9,14,26)}
    save(out/"comparison.json",dict(common_coordinates=count,common_recontact=len(shared),
        passed_fields=gates,stages=summary,endpoint_effects_m2=attribution))
    attribute_halves(out)


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]))
