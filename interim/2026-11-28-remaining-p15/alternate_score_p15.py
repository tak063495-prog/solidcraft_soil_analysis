"""Reapply the inherited 22-field decision and compare at identical recontact loads."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P14=HERE.with_name("2026-11-27-remaining-p14")
sys.path.insert(0,str(P14))
from alternate import accepted,score_route
from plans import save
from trial_cost import run as cost


def run(work,out):
    dest=out/"c01/alternate"
    score_route(work,dest,"outer")
    cost(out)
    dirs={"p14_middle":P14/"c01/alternate/middle","p14_outer":P14/"c01/alternate/outer","p15_outer":dest/"outer"}
    summaries={}
    clocks={}
    ends={}
    for route,folder in dirs.items():
        s=json.loads((folder/"score.json").read_text())
        summaries[route]={k:v for k,v in s.items() if k!="fields"}
        summaries[route]["failed_fields"]=[k for k,v in s["fields"].items() if not v["full"]["passes"]]
        clocks[route]=[]
        ends[route]=[]
        for n in (256,512,1024):
            halves,clock=accepted(folder/str(n)/"trace.jsonl.gz")
            clocks[route].append(clock)
            ends[route].append({r for _,h,_,r,_ in halves if h==2})
    def matched(routes,file):
        grid=sorted(set.intersection(*(e for route in routes for e in ends[route])))
        with gzip.open(dest/file,"wt") as f:
            for t in grid:
                data={}
                for route in routes:
                    a,b=-np.diff([c[t] for c in clocks[route]],axis=0)
                    data[route]=(a*b).sum(1).tolist()
                print(json.dumps(dict(t=t,points_m2=data)),file=f)
        return dict(count=len(grid),last=t,endpoint_m2={s:sum(v) for s,v in data.items()})
    save(dest/"comparison.json",dict(routes=summaries,all_shared=matched(list(dirs),"common-all.jsonl.gz"),
         p15_vs_middle=matched(["p14_middle","p15_outer"],"common-middle.jsonl.gz")))


if __name__=="__main__":
    run(*map(Path,sys.argv[1:]))
