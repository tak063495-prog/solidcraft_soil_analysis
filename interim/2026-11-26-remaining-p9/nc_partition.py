"""Supplement the prespecified search near the observed material partition change."""
import gzip
import json
import os
from pathlib import Path
import sys

import numpy as np

from nc_scan import P8,p8


def run(work,out):
    summary=json.loads((out/"summary.json").read_text())
    near=summary["nearby_parents"]
    with gzip.open(out/"point-inputs.jsonl.gz","rt") as f:
        inputs={int(r["trial"]):r["input"] for r in map(json.loads,f) if r["trial"].isdigit() and
                int(r["trial"]) in near and r["evaluation"]==1 and r["input"]["point"]==272}
    arc={r["attempt"]:r for r in map(json.loads,(p8.P7/"nc/continuation.jsonl").read_text().splitlines())}
    start=max(int(path.name.split("-")[1].split(".")[0]) for path in (P8/"nc").glob("parent-*.json.gz")
              if int(path.name.split("-")[1].split(".")[0])<=min(near))
    old=p8.prior.old
    parent=old.read(P8/f"nc/parent-{start}.json.gz")
    state,lam,t=parent["active_material_states"],parent["load"],np.array(parent["secant"])
    probe=work/"probe/target/release/examples"
    os.environ["NC_WORKER"]=str(probe/"p5_homogeneous")
    m=old.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    watched=old.nc_extended.Extended(work/"inputs",P8/"capture.py")
    cap=work/"partition-capture"
    os.environ["P8_CAPTURE"]=str(cap)
    report=[]
    calls=0
    with gzip.open(out/"partition-observations.jsonl.gz","xt") as f:
        for a in range(start,max(near)+1):
            row=arc[a]
            assert abs(lam-row["parent_lambda"])<1e-13
            if a in near:
                # At the first global predictor the strain is linear in ds.
                # Above this width the existing rule changes count from 1 to 4.
                cut=row["ds_m"]*m.width/np.abs(inputs[a]["strain"]).sum()
                widths=[cut*x for x in (.99,.995,.9975,.999,.99999,1.00001)]
                if a==summary["first"]["group"]:
                    low,high=[b for b in summary["brackets"] if b["field"]=="negative272" and b["group"]==a]
                    widths+=np.linspace(low["upper"],high["lower"],9)[1:-1].tolist()
                samples=[]
                for width in widths:
                    cap.write_text("")
                    result,it,reason=old.nc_branch.correct(watched,state,lam,t,width,1e-12)
                    calls+=1
                    signs=[]
                    for j,call in enumerate(map(json.loads,cap.read_text().splitlines()),1):
                        neg={}
                        for rec in call["points"]:
                            obs=p8.point(probe,rec)
                            neg[rec["point"]]={s["step"] for s in obs["steps"] if s["direction"]<0}
                            print(json.dumps(dict(group=a,width=width,evaluation=j,input=rec,observation=obs)),file=f)
                        signs.append(bool(neg[351]&neg[272]))
                    samples.append(dict(width=width,accepted=result is not None,reason=reason,both=any(signs)))
                report.append(dict(parent=a,load=lam,partition_width=cut,samples=samples))
            if a<max(near):
                result,_,reason=old.nc_branch.correct(m,state,lam,t,row["ds_m"],1e-12)
                calls+=1
                assert (result is not None)==row["accepted"] and reason==row.get("reason")
                if result is not None:
                    dz,state,_=result
                    lam+=dz[-1]
                    t=dz/np.linalg.norm(dz)
    total=calls+summary["total_charged"]
    assert total<=2000
    (out/"partition.json").write_text(json.dumps(dict(parents=report,additional_calls=calls,
        total_charged_conservatively=total),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
