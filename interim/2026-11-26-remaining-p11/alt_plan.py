"""Attribute the pre-intervention endpoint to accepted half steps and freeze target loads."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P10=HERE.parent/"2026-11-26-remaining-p10/c01/alternate"

def read(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(s) for s in f if s.startswith("{")]

def collect(folder,n):
    rows=read(folder/f"observed-{n}.jsonl.gz")
    spans={(r["p11_span"][0],r["p11_span"][2]):r for r in rows if "p11_span" in r}
    def accepted(records):
        return [r for r in records if r.get("stage")=="recontact"]
    attempts=accepted(rows)
    prior=accepted(read(P10/str(n)/"trace.jsonl.gz"))
    for x,y in zip(attempts,prior,strict=True):
        for k in ("from","to","accepted","reason"):
            assert x[k]==y[k]
    halves=[]
    for j,r in enumerate(attempts,1):
        if r["accepted"]:
            s=spans[r["from"],r["to"]]
            for h in (0,1):
                halves.append((j,h+1,*s["p11_span"][h:h+2],*s["all"][h:h+2]))
    return halves

def run(folder):
    archive=json.loads((P10/"endpoint-32.json").read_text())
    plan=json.loads((P10/"plan.json").read_text())
    knots=set(plan["knots"])
    transitions=[]
    clocks=[]
    all_halves={}
    attribution=[]
    for level,n in enumerate((256,512,1024)):
        halves=all_halves[n]=collect(folder,n)
        clock={halves[0][2]:halves[0][-2]}
        for _,_,l,r,x,y in halves:
            np.testing.assert_array_equal([v[1] for v in clock[l]],[v[1] for v in x])
            clock[r]=y
        clocks.append(clock)
        for point in (17,23):
            v=np.array(archive[point]["p10"]["slip"])
            a,b=v[0]-v[1],v[1]-v[2]
            np.testing.assert_array_equal(halves[-1][-1][point][1],v[level])
            # Symmetric, exact telescoping decomposition of a.b across levels.
            coefficient=(b/2,(a-b)/2,-a/2)[level]
            initial=float(coefficient@np.array(halves[0][-2][point][1]))
            contrib=[]
            for attempt,h,l,r,x,y in halves:
                sx,ux=x[point]
                sy,uy=y[point]
                increment=np.array(uy)-ux
                contrib.append(float(coefficient@increment))
                if sx!=sy:
                    transitions.append(dict(n=n,point=point,attempt=attempt,half=h,span=[l,r],
                        transition=[sx,sy],delta=increment.tolist(),endpoint_term=contrib[-1]))
                    if n==512 and "open" in sx.lower() and sy=="slip":
                        knots.update((l,r))
            attribution.append(dict(n=n,point=point,initial=initial,half_sum=float(sum(contrib)),endpoint=float(coefficient@v[level])))
    assert {r["point"] for r in transitions if r["n"]==512 and "open" in r["transition"][0].lower()}=={17,23}
    # At common accepted half endpoints both adjacent differences are evaluated
    # at precisely the same load. The endpoint attribution telescopes over these
    # bins, even when the three native subdivisions between loads differ.
    for point in (17,23):
        changes=[r["span"] for r in transitions if r["point"]==point and r["transition"]==["open","slip"]]
        assert len(changes)==3
        knots.update((min(l for l,_ in changes),max(r for _,r in changes)))
    shared=sorted(set.intersection(*(set(c) for c in clocks)))
    assert shared[0]==0 and shared[-1]==1
    intervals={}
    for point in (17,23):
        v=np.array(archive[point]["p10"]["slip"])
        final_a,final_b=-np.diff(v,axis=0)
        records=[]
        for t in shared:
            levels=[c[t][point] for c in clocks]
            a,b=-np.diff([r[1] for r in levels],axis=0)
            records.append(dict(t=t,status=[r[0] for r in levels],dot=float(a@b),
                endpoint_share=float((a@final_b+b@final_a)/2)))
        spans=[dict(left=a,right=b,delta_endpoint=b["endpoint_share"]-a["endpoint_share"])
               for a,b in zip(records,records[1:])]
        localized=[]
        for group in sorted(spans,key=lambda s:s["delta_endpoint"])[:3]:
            lo,hi=group["left"]["t"],group["right"]["t"]
            group["accepted_halves"]={}
            for n in (256,512,1024):
                group["accepted_halves"][n]=[dict(attempt=j,half=h,l=l,r=r,labels=[x[point][0],y[point][0]],
                    increment=(np.array(y[point][1])-x[point][1]).tolist())
                    for j,h,l,r,x,y in all_halves[n] if lo<=l<r<=hi]
            localized.append(group)
        intervals[point]=dict(worst=localized,initial=records[0],end=records[-1],
            net=float(sum(r["delta_endpoint"] for r in spans)))
    result=dict(signed_attribution=intervals,knots=sorted(knots),p10_knots=plan["knots"],transitions=transitions,attribution=attribution)
    (folder/"plan.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
if __name__=="__main__":
    run(Path(sys.argv[1]))
