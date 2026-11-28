"""Attribute outer and supplemental evaluations to each pretrial parent."""
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")


def load_trace(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(s) for s in f if s.startswith("{")]


def cost(path):
    trace = load_trace(path)
    rows = [r for r in trace if r.get("stage")=="recontact" and "accepted" in r]
    extras = [r for r in trace if r.get("stage")==3 and "charge" in r]
    charges = defaultdict(Counter)
    for i,e in enumerate(extras,1):
        # Outer trials are also charged, first, in actual execution order.
        trial = e["spent"]-i
        assert 1<=trial<=len(rows)
        charges[trial][e["charge"]] += 1
    bins, parents = defaultdict(Counter),defaultdict(Counter)
    last = 0.
    for i,r in enumerate(rows,1):
        assert last==r["from"]
        c=Counter(accepted=int(r["accepted"]),rejected=int(not r["accepted"]),**charges[i])
        if r["accepted"]:
            last=r["to"]
        else:
            assert r["parent_unchanged"]
        c["total"]=sum(c.values())
        for dest in (bins[math.floor(r["from"]*1024)],parents[r["from"]]):
            dest.update(c)
    assert sum(c["total"] for c in bins.values())==len(rows)+len(extras)<=1787
    return dict(last=last,outer=len(rows),extra=len(extras),rejected=sum(not r["accepted"] for r in rows),
                total=len(rows)+len(extras),last_accepted=rows[max(i for i,r in enumerate(rows) if r["accepted"])]),bins,parents


def run(out):
    dest = out/"c01/alternate"
    sources = {"p14_outer":P14/"c01/alternate/outer","middle":P14/"c01/alternate/middle"}
    if (dest/"outer/1024/trace.jsonl.gz").exists():
        sources["p15_outer"]=dest/"outer"
    data={s:cost(p/"1024/trace.jsonl.gz") for s,p in sources.items()}
    bands={}
    for k in sorted(set.union(*(set(b) for _,b,_ in data.values()))):
        values={s:x[1][k] for s,x in data.items()}
        diff = values["p14_outer"].copy()
        diff.subtract(values["middle"])
        if any(diff.values()) or len({tuple(sorted(v.items())) for v in values.values()})>1:
            bands[k]=dict(span=[k/1024,(k+1)/1024],counts=values,p14_outer_minus_middle=dict(diff))
    for name,(_,_,parents) in data.items():
        with gzip.open(dest/(name+"-parents.jsonl.gz"),"wt") as f:
            for t,c in parents.items():
                print(json.dumps(dict(parent=t,**c)),file=f)
    (dest/"cost.json").write_text(json.dumps(dict(bin_width=1/1024,totals={s:v[0] for s,v in data.items()},
        differing_bands=bands),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
