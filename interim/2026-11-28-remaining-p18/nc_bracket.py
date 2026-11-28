"""Probe widths on frozen parents; classify every parent only by the original global corrector."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P17 = HERE.with_name("2026-11-28-remaining-p17")
sys.path.insert(0, str(P17))
from nc_parent import prior, P16


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def run(work, out):
    probe = work/"probe/target/release/examples"
    model = prior.old.nc_extended.Extended(work/"inputs", probe/"p5_homogeneous")
    assert len(model.free) == 348 and model.b.shape[:2] == (48, 27)
    base = json.loads((P17/"nc/summary.json").read_text())
    assert base["cumulative"] == 1368
    p13 = [json.loads(x) for x in (prior.P13/"nc/global.jsonl").read_text().splitlines()]
    ledger = out/"global.jsonl"
    total = base["cumulative"]
    if ledger.exists():
        total = max(json.loads(s)["total"] for s in ledger.read_text().splitlines())
    def record(data):
        with ledger.open("a") as stream:
            print(json.dumps(dict(total=total, **data),allow_nan=False),file=stream)
    summary = {}
    os.environ.update(P16_LU=str(probe/"local_lu"), P16_PYTHON=sys.executable,
                      P16_COMPENSATE=str(P16/"compensate.py"))
    for path in ("ordinary","predicted"):
        end = prior.old.read(prior.P13/f"nc/{path}-end.json.gz")
        state,lam,t,u = end["active_material_states"],end["load"],np.array(end["secant"]),np.array(end["displacement_m"])
        parent = state,lam,t,u,end["next_width"]
        key = prior.digest(parent)
        assert key == base["paths"][path]["parent_hash"]
        pred = next(r["frozen_prediction"] for r in reversed(p13) if r["path"]==path and "frozen_prediction" in r)
        br, = pred["unresolved_brackets"]
        left,mid,right = br["left"]["width"],br["failed_midpoint"],br["right"]["width"]
        flat = np.zeros_like(model.base)
        flat[model.free] = t[:-1]
        slope = -np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)[351]
        samples = {}

        def trial(width):
            nonlocal total
            # Re-read completed trials when resuming. A started trial remains charged.
            filename = f"{path}-{width.hex()}"
            result = out/(filename+".json")
            if result.exists():
                row = json.loads(result.read_text())
                samples[width] = row
                return row
            total += 1
            record(dict(event="start",path=path,parent=key,width=width))
            trace = work/f"nc-{total}.jsonl"
            os.environ.update(P16_STRAIN=json.dumps((slope*width).tolist()),
                              P16_RECORD=str(work/f"nc-{total}-p16.jsonl"),P17_NC_LOG=str(trace))
            sol,it,why = prior.old.nc_branch.correct(model,state,lam,t,width,1e-12)
            assert prior.digest(parent) == key
            accepted = sol is not None
            candidate = prior.advance(model,parent,width,sol,it) if accepted else parent
            row = dict(path=path,total=total,width_m=width,accepted=accepted,iterations=it,reason=why,
                       parent_hash=key,end_hash=prior.digest(candidate),lambda_before=lam,lambda_after=candidate[1])
            record(dict(event="verdict",**{k:v for k,v in row.items() if k!="total"}))
            if accepted and not (out/f"{path}-accepted.json.gz").exists():
                with gzip.open(out/f"{path}-accepted.json.gz","xt") as f:
                    json.dump(dict(active_material_states=candidate[0],load=candidate[1],secant=candidate[2].tolist(),
                              displacement_m=candidate[3].tolist(),next_width=candidate[-1],state_hash=row["end_hash"],width_m=width),f)
            calls = []
            with trace.open() as stream:
                for line in stream:
                    r = json.loads(line)
                    if "evaluation" in r:
                        calls.append(dict(strain=r["strain"],count=r["count"]))
                    elif r.get("iteration")==0 and "matrix" not in calls[-1] and len(calls)<=2:
                        call = calls[-1]
                        matrix = dict(matrix=r["jac_rows"],rhs=r["r"],native=r["dx"])
                        check = work/f"nc-{total}-{len(calls)}-check.jsonl"
                        proc = subprocess.run([sys.executable,P16/"compensate.py"],input=json.dumps(matrix),text=True,
                            capture_output=True,env=dict(os.environ,P16_RECORD=str(check)))
                        audit = json.loads(check.read_text())
                        assert audit["passes"] == (proc.returncode==0)
                        call.update(matrix=audit["matrix"],rhs=audit["rhs"],native_multiplier=r["dx"][-1],
                                    multiplier=audit["refined"][-1],schur=audit["schur"],denominator=audit["denominator"],
                                    allowance=audit["allowance"],matched=audit["passes"])
            row["calls"] = calls[:2]
            row["material_evaluations_at_351"] = len(calls)
            dump(result,row)
            samples[width] = row
            return row

        for width in (left,mid,right):
            trial(width)
        def signed(width,call):
            calls = trial(width)["calls"]
            if len(calls)>call and calls[call].get("matched"):
                return np.sign(calls[call]["multiplier"])
            return None

        def bracket(call):
            valid = [(w,signed(w,call)) for w in sorted(samples) if signed(w,call) is not None]
            return next(((a,b) for (a,sa),(b,sb) in zip(valid,valid[1:]) if sa>0>sb),None)

        # Resolve the corrected side from actual equilibria, never by scaling a corrected strain.
        # If needed extend to the lower width already sampled by p13's point search.
        for j in range(1,12):
            if bracket(1):
                break
            trial(left-(mid-left)*2**j)
        for call in (0,1):
            for _ in range(5):
                pair = bracket(call)
                if not pair:
                    break
                a,b = pair
                h = (a+b)/2
                if h in samples or signed(h,call) is None:
                    break
        with gzip.open(prior.P13/f"nc/{path}-local.jsonl.gz","rt") as stream:
            archived = [r for line in stream if (r:=json.loads(line)).get("point")==351
                        and r["parent"]==pred["parent"] and r.get("partition")==1 and r.get("first")]
        prior_low = max((r for r in archived if r["reason"] is None and r["width"]<=min(samples)),
                        key=lambda r:r["width"])
        prior_high = min((r for r in archived if r["reason"] and r["width"]>=max(samples)),
                         key=lambda r:r["width"])
        summary[path] = dict(input_widths=[left,mid,right],parent_hash=key,
            predictor=bracket(0),corrected=bracket(1),
            previous_local_success=prior_low,previous_local_failure=prior_high,
            tested_widths=sorted(samples),accepted_widths=sorted(w for w,v in samples.items() if v["accepted"]))
        dump(out/"summary.json",dict(paths=summary,point=351,rtol=5e-9,atol=1e-14,start=1368,cumulative=total,
                                     additional=total-1368))
        print(path,summary[path],flush=True)


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
