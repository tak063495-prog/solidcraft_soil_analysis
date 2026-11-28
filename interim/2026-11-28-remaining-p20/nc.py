"""Continue the four immutable p19 parents round-robin under the original shared budget."""
import gzip
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P19 = HERE.with_name("2026-11-28-remaining-p19")
sys.path.insert(0,str(P19))
from nc_next import P16, parent_at, prior


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def run(work, out):
    probe = work/"probe/target/release/examples"
    model = prior.old.nc_extended.Extended(work/"inputs", probe/"p5_homogeneous")
    assert len(model.free)==348 and model.b.shape[:2]==(48,27)
    os.environ.update(P16_LU=str(probe/"local_lu"), P16_PYTHON=sys.executable,
                      P16_COMPENSATE=str(P16/"compensate.py"))
    old = json.loads((P19/"nc/summary.json").read_text())
    total = old["cumulative"]
    assert total==1406
    states, stats = {}, {}
    for name, ref in old["paths"].items():
        parent = parent_at(P19/f"nc/{name}-next.json.gz")
        assert parent[1]==ref["end_lambda"] and prior.digest(parent)==ref["end_hash"]
        states[name] = parent
        stats[name] = dict(start_lambda=parent[1],end_lambda=parent[1],start_hash=ref["end_hash"],
            end_hash=ref["end_hash"],arc_m=0.,chord_m=0.,accepted=0,attempts=0,next_trial_width=parent[-1],
            first_positive_tangent=None,first_ascending_interval=None,stop=None,reached_one=False)

    with (out/"global.jsonl").open("x",buffering=1) as ledger, gzip.open(out/"point351.jsonl.gz","xt") as archive:
        def record(**fields):
            print(json.dumps(fields,allow_nan=False),file=ledger)
        while total<2000 and any(s["stop"] is None for s in stats.values()):
            for name,s in stats.items():
                if s["stop"] or total>=2000:
                    continue
                parent = states[name]
                width = s["next_trial_width"]
                total += 1
                s["attempts"] += 1
                key = prior.digest(parent)
                record(event="start",total=total,path=name,parent_hash=key,width_m=width,arc_m=s["arc_m"],
                       lambda_before=parent[1],predictor_secant_lambda=parent[2][-1])
                flat=np.zeros_like(model.base)
                flat[model.free]=parent[2][:-1]
                slope=-np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)[351]
                trace, check = (work/f"p20-{total}-{part}.jsonl" for part in ("trace","check"))
                os.environ.update(P17_NC_LOG=str(trace),P16_RECORD=str(check),P16_STRAIN=json.dumps((slope*width).tolist()))
                sol,it,reason=prior.old.nc_branch.correct(model,*parent[:3],width,1e-12)
                assert prior.digest(parent)==key
                candidate=prior.advance(model,parent,width,sol,it) if sol is not None else parent
                if sol is not None:
                    s["accepted"]+=1
                    s["arc_m"]+=width
                    s["chord_m"]+=float(np.linalg.norm(sol[0]))
                    for field, test in [("first_positive_tangent",sol[2]["lambda_tangent"]>0),
                                        ("first_ascending_interval",candidate[1]>parent[1])]:
                        if test and s[field] is None:
                            s[field]=dict(total=total,arc_m=s["arc_m"],parent=parent[1],load=candidate[1])
                    s.update(end_hash=prior.digest(candidate),end_lambda=candidate[1],next_trial_width=candidate[-1])
                    states[name]=candidate
                    if candidate[1]>=1:
                        s.update(stop="reached lambda=1",reached_one=True)
                else:
                    s["next_trial_width"]=width/2
                    if width/2<model.settings["minimum_increment"]:
                        s["stop"]="minimum increment"
                directions=[]
                if trace.exists():
                    # Save every reached 351 evaluation, including rejected trials.
                    for line in trace.read_text().splitlines():
                        row=json.loads(line)
                        print(json.dumps(dict(total=total,trace=row),separators=(",",":")),file=archive)
                if check.exists():
                    for line in check.read_text().splitlines():
                        audit=json.loads(line)
                        directions.append(dict(denominator=audit["denominator"],direction=audit["refined"][-1],
                            schur=audit["schur"],matched=audit["passes"],allowance=audit["allowance"]))
                    assert all(d["matched"] for d in directions) or sol is None
                record(event="verdict",total=total,path=name,accepted=sol is not None,reason=reason,
                    iterations=it,arc_m=s["arc_m"],end_hash=s["end_hash"],lambda_after=candidate[1],
                    secant_lambda_after=candidate[2][-1],equilibrium_tangent=sol[2] if sol is not None else None,delta_lambda=candidate[1]-parent[1],
                    point351_live=directions)
                for p in (trace,check):
                    p.unlink(missing_ok=True)
    for name,s in stats.items():
        if not s["stop"]:
            s["stop"]="cumulative global trial limit 2000"
        p=states[name]
        with gzip.open(out/f"{name}-end.json.gz","xt") as dest:
            json.dump(dict(active_material_states=p[0],load=p[1],secant=p[2].tolist(),
                displacement_m=p[3].tolist(),next_width=p[-1],state_hash=prior.digest(p)),dest,allow_nan=False)
    dump(out/"summary.json",dict(start=1406,cumulative=total,additional=total-1406,paths=stats,
        minimum_increment=model.settings["minimum_increment"],rtol=5e-9,atol=1e-14))


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
