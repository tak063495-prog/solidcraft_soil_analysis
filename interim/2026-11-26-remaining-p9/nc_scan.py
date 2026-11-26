"""Trace both material inputs on every global corrector call and bracket sign changes."""
import gzip
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P8 = HERE.parent / "2026-11-26-remaining-p8"
sys.path.insert(0,str(P8))
import nc as p8


def run(work,out):
    out.mkdir(parents=True,exist_ok=True)
    probe = work/"probe/target/release/examples"
    os.environ["NC_WORKER"] = str(probe/"p5_homogeneous")
    old = p8.prior.old
    m = old.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    observed = old.nc_extended.Extended(work/"inputs",P8/"capture.py")
    archive = [json.loads(l) for l in (p8.P7/"nc/continuation.jsonl").read_text().splitlines()]
    start = max(int(path.name.split("-")[1].split(".")[0]) for path in (P8/"nc").glob("parent-*.json.gz")
                if int(path.name.split("-")[1].split(".")[0]) <= 460)
    parent = old.read(P8/f"nc/parent-{start}.json.gz")
    state,lam,t = parent["active_material_states"],parent["load"],np.array(parent["secant"])
    u = np.ravel(parent["displacement_m"])
    assert len(m.free)==348 and len(state)==1296
    snap = {}
    traced = []
    cache = {}
    cap = work/"p9capture.jsonl"
    os.environ["P8_CAPTURE"] = str(cap)
    count = 0
    def trial(key,state,lam,t,width,group,kind):
        nonlocal count
        cap.write_text("")
        count += 1
        assert 505 + start-106 + count <= 2000
        result,it,reason = old.nc_branch.correct(observed,state,lam,t,width,1e-12)
        flags={351:False,272:False}
        simultaneous = []
        starts={351:0,272:0}
        for j,line in enumerate(cap.read_text().splitlines(),1):
            call=json.loads(line)
            local={}
            for rec in call["points"]:
                i=rec["point"]
                assert rec["state"]==state[i]
                sample=p8.point(probe,rec)
                starts[i]+=len(sample["steps"])
                local[i]={s["step"] for s in sample["steps"] if s["direction"]<0}
                flags[i] |= bool(local[i])
                print(json.dumps(dict(trial=key,evaluation=j,input=rec,worker_exit=call["exit"],
                                      observation=sample),allow_nan=False),file=detail)
            simultaneous.extend((j,s) for s in sorted(local[351]&local[272]))
        row=dict(trial=key,group=group,parent=lam,width=width,kind=kind,accepted=result is not None,
                 reason=reason,iterations=it,negative351=flags[351],negative272=flags[272],
                 same_evaluation_and_substep=simultaneous,starts=starts)
        print(json.dumps(row,allow_nan=False),file=log,flush=True)
        cache[group,width]=row
        return result,row

    with gzip.open(out/"point-inputs.jsonl.gz","xt") as detail, (out/"trials.jsonl").open("x") as log:
        group=start
        ds=None
        for row in archive:
            a=row["attempt"]
            if a<start: continue
            if ds is None: ds=row["ds_m"]
            assert ds==row["ds_m"] and abs(lam-row["parent_lambda"])<1e-13
            if a>=460:
                if a==group:
                    snap[group]=(state,lam,t.copy())
                result,record=trial(str(a),state,lam,t,ds,group,"path")
                traced.append(record)
                reason=record["reason"]
            else:
                result,it,reason=old.nc_branch.correct(m,state,lam,t,ds,1e-12)
                count+=1
            assert (result is not None)==row["accepted"] and reason==row.get("reason")
            if result is None:
                ds/=2
            else:
                dz,state,diag=result
                lam+=dz[-1]
                u[m.free]+=dz[:-1]
                t=dz/np.linalg.norm(dz)
                ds=min(.04,ds*(1.3 if len(row["iterations"])<=8 else 1.))
                group=a+1
        assert len(traced)==178
        onset=next(r for r in traced if r["negative272"])
        parents=sorted(snap)
        k=parents.index(onset["group"])
        nearby=parents[max(k-2,0):k+3]
        width=onset["width"]
        brackets=[]
        for group in nearby:
            def sample(w):
                if (group,w) not in cache:
                    trial(f"{group}@{w.hex()}",*snap[group],w,group,"width_probe")
                return cache[group,w]
            grid=[sample(width*j/16) for j in range(8,21)]
            for key in ("negative351","negative272","same_evaluation_and_substep"):
                for left,right in zip(grid,grid[1:]):
                    if bool(left[key])==bool(right[key]): continue
                    lo,hi=left,right
                    for _ in range(10):
                        mid=sample((lo["width"]+hi["width"])/2)
                        if bool(mid[key])==bool(lo[key]): lo=mid
                        else: hi=mid
                    brackets.append(dict(group=group,field=key,lower=lo["width"],upper=hi["width"],
                                         lower_sign=bool(lo[key]),upper_sign=bool(hi[key])))
        result=dict(first=onset,previous=traced[traced.index(onset)-1] if onset!=traced[0] else None,
            path_counts={k:sum(bool(r[k]) for r in traced) for k in ("negative351","negative272","same_evaluation_and_substep")},
            nearby_parents=nearby,grid_base=width,grid_multipliers=[j/16 for j in range(8,21)],
            brackets=brackets,prior_to_start=505+start-106,global_calls=count,
            total_charged=505+start-106+count,lambda_end=lam,tangent_end=diag["lambda_tangent"])
        (out/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
