"""Repeat each stopped global trial once and observe the real material calls."""
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P16=HERE.with_name("2026-11-28-remaining-p16")
sys.path.insert(0,str(P16))
import parent_nc as prior


def run(work,out):
    probe=work/"probe/target/release/examples"
    model=prior.old.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    assert len(model.free)==348 and model.b.shape[:2]==(48,27)
    ref=json.loads((P16/"nc/parents.json").read_text())
    total=ref["cumulative"]
    records={}
    ledger=out/"trials.jsonl"
    assert not ledger.exists()
    def note(entry):
        with ledger.open("a") as stream:
            print(json.dumps(entry,allow_nan=False),file=stream)
    for path,previous in ref["paths"].items():
        end=prior.old.read(prior.P13/f"nc/{path}-end.json.gz")
        state,lam,t,u=end["active_material_states"],end["load"],np.array(end["secant"]),np.array(end["displacement_m"])
        parent=state,lam,t,u,end["next_width"]
        key=prior.digest(parent)
        assert key==end["state_hash"]==previous["parent_hash"]
        w=previous["width_m"]
        flat=np.zeros_like(model.base)
        flat[model.free]=t[:-1]
        strain=-np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)[351]*w
        trace=out/(path+"-material.jsonl")
        directions=out/(path+"-directions.jsonl")
        assert not (trace.exists() or directions.exists())
        os.environ.update(P16_STRAIN=json.dumps(strain.tolist()),P16_LU=str(probe/"local_lu"),
            P16_PYTHON=sys.executable,P16_COMPENSATE=str(P16/"compensate.py"),P16_RECORD=str(directions),
            P17_NC_LOG=str(trace))
        total+=1
        note(dict(event="start",path=path,total=total,parent_hash=key,width_m=w))
        sol,it,reason=prior.old.nc_branch.correct(model,state,lam,t,w,1e-12)
        note(dict(event="verdict",path=path,total=total,accepted=sol is not None,iterations=it,
                  reason=reason,parent_hash_after=prior.digest(parent)))
        assert prior.digest(parent)==key
        assert (sol is None and it==previous["iterations"] and reason==previous["reason"])
        assert directions.read_bytes()==(P16/f"nc/{path}-global-direction.jsonl").read_bytes()
        rows=[json.loads(s) for s in trace.read_text().splitlines()]
        split=max(i for i,r in enumerate(rows) if "evaluation" in r)
        reached=rows[split:]
        pending=max(i for i,r in enumerate(reached) if "iteration" in r)
        iteration=reached[pending]
        candidates=[r for r in reached[pending+1:] if "alpha" in r]
        assert len(candidates)==30
        counts={}
        for r in candidates:
            cond=dict(finite=r["finite"],log_pc=abs(r["log_pc"])<100,
                multiplier=r["multiplier"]>=0,pressure=r["pressure"]>0)
            declines=r["candidate_norm"] is not None and r["candidate_norm"]<iteration["norm"]
            assert not(all(cond.values()) and declines)
            for field,ok in cond.items():
                counts[field]=counts.get(field,0)+int(ok)
            counts["residual_decrease"]=counts.get("residual_decrease",0)+int(declines)
            counts["admissible"]=counts.get("admissible",0)+int(all(cond.values()))
        records[path]=dict(parent=previous["parent"],parent_hash=key,end_hash=prior.digest(parent),
            lambda_before=lam,lambda_after=lam,width_m=w,accepted=False,iterations=it,reason=reason,
            material_calls=sum("evaluation" in r for r in rows),failed_call=reached[0],
            failed_iteration=iteration["iteration"],residual=iteration["norm"],candidate_count=len(candidates),
            condition_pass_counts=counts,first=candidates[0],last=candidates[-1],
            compensated_checks=len(directions.read_text().splitlines()),additional_global_trials=1,cumulative=total)
        (out/"summary.json").write_text(json.dumps(dict(paths=records,rtol=5e-9,atol=1e-14,
            start_count=ref["cumulative"],cumulative=total),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
