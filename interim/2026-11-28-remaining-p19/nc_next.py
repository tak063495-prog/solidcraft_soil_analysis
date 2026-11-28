"""Try one successor from each named and serialized p18 equilibrium."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P18 = HERE.with_name("2026-11-28-remaining-p18")
sys.path.insert(0, str(P18))
from nc_bracket import prior, P16


def parent_at(path):
    with gzip.open(path, "rt") as f:
        row = json.load(f)
    parent = (row["active_material_states"], row["load"], np.array(row["secant"]),
              np.array(row["displacement_m"]), row["next_width"])
    assert prior.digest(parent) == row["state_hash"]
    return parent


def run(work, out):
    bin_dir = work/"probe/target/release/examples"
    model = prior.old.nc_extended.Extended(work/"inputs", bin_dir/"p5_homogeneous")
    assert len(model.free) == 348 and model.b.shape[:2] == (48,27)
    os.environ.update(P16_LU=str(bin_dir/"local_lu"),P16_PYTHON=sys.executable,
                      P16_COMPENSATE=str(P16/"compensate.py"))
    total = json.loads((P18/"nc/summary.json").read_text())["cumulative"]
    assert total == 1396
    targets = dict(ordinary=0.5481291446371566,predicted=0.5481291312021739)
    results = {}
    log_path = out/"global.jsonl"
    history = [json.loads(s) for s in log_path.read_text().splitlines()] if log_path.exists() else []
    total += sum(r["event"]=="start" for r in history)
    with log_path.open("a",buffering=1) as ledger:
        def trial(path, purpose, parent, width):
            nonlocal total
            total += 1
            assert total <= 2000
            key = prior.digest(parent)
            identity = dict(total=total,path=path,purpose=purpose,width_m=width,parent_hash=key,
                            lambda_before=parent[1])
            def record(**fields):
                print(json.dumps(dict(**identity,**fields),allow_nan=False),file=ledger)
            record(event="start")
            flat = np.zeros_like(model.base)
            flat[model.free] = parent[2][:-1]
            strain = -np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)[351]*width
            trace = work/f"material-{total}.jsonl"
            os.environ.update(P17_NC_LOG=str(trace),P16_RECORD=str(out/f"{total}-live.jsonl"),
                              P16_STRAIN=json.dumps(strain.tolist()))
            sol, it, why = prior.old.nc_branch.correct(model,*parent[:3],width,1e-12)
            assert prior.digest(parent) == key
            candidate = prior.advance(model,parent,width,sol,it) if sol is not None else parent
            record(event="verdict",accepted=sol is not None,reason=why,iterations=it,
                   lambda_after=candidate[1],end_hash=prior.digest(candidate))
            directions = []
            if trace.exists():
                with trace.open("rb") as src, gzip.open(out/f"{total}-point351.jsonl.gz","wb") as dest:
                    shutil.copyfileobj(src,dest)
                call = 0
                for line in trace.read_text().splitlines():
                    r = json.loads(line)
                    if "evaluation" in r:
                        call += 1
                    if r.get("iteration") != 0:
                        continue
                    payload = dict(matrix=r["jac_rows"],rhs=r["r"],native=r["dx"])
                    check = out/f"{total}-checks.jsonl"
                    proc = subprocess.run([sys.executable,P16/"compensate.py"],input=json.dumps(payload),
                        text=True,capture_output=True,env=dict(os.environ,P16_RECORD=str(check)))
                    audit = json.loads(check.read_text().splitlines()[-1])
                    assert audit["passes"] == (proc.returncode==0)
                    directions.append(dict(call=call,residual=r["norm"],direction=audit["refined"][-1],
                        denominator=audit["denominator"],schur=audit["schur"],matched=audit["passes"]))
            record(event="diagnostics",directions=directions)
            return candidate,sol is not None

        for path, target in targets.items():
            archived = [json.loads(f.read_text()) for f in (P18/"nc").glob(f"{path}-*.json")]
            ref, = [r for r in archived if r["accepted"] and r["lambda_after"]==target]
            stored = parent_at(P18/f"nc/{path}-accepted.json.gz")
            for source in ("named","serialized"):
                key = path+"-"+source
                completed = out/(key+"-next.json.gz")
                if completed.exists():
                    candidate = parent_at(completed)
                    rows = [r for r in history if r["path"]==path and r["purpose"]==source and r["event"]=="verdict"]
                    assert rows[-1]["accepted"] and rows[-1]["end_hash"]==prior.digest(candidate)
                    results[key] = dict(start_lambda=rows[0]["lambda_before"],end_lambda=candidate[1],
                        start_hash=rows[0]["parent_hash"],end_hash=prior.digest(candidate),attempts=len(rows),
                        accepted=True,stop="requested next parent accepted",reached_one=bool(candidate[1]>=1))
                    continue
                if source=="named":
                    ancestor = parent_at(prior.P13/f"nc/{path}-end.json.gz")
                    assert prior.digest(ancestor) == ref["parent_hash"]
                    parent,ok = trial(path,"reconstruct_named",ancestor,ref["width_m"])
                    assert ok and prior.digest(parent)==ref["end_hash"] and parent[1]==target
                else:
                    parent = stored
                    if prior.digest(parent)==ref["end_hash"]:
                        continue
                start = total
                width = parent[-1]
                candidate,ok = parent,False
                while width >= model.settings["minimum_increment"]:
                    candidate,ok = trial(path,source,parent,width)
                    if ok:
                        break
                    width /= 2
                key = path+"-"+source
                if ok:
                    with gzip.open(out/(key+"-next.json.gz"),"xt") as f:
                        json.dump(dict(active_material_states=candidate[0],load=candidate[1],
                            secant=candidate[2].tolist(),displacement_m=candidate[3].tolist(),
                            next_width=candidate[-1],state_hash=prior.digest(candidate)),f,allow_nan=False)
                results[key] = dict(start_lambda=parent[1],end_lambda=candidate[1],
                    start_hash=prior.digest(parent),end_hash=prior.digest(candidate),attempts=total-start,
                    accepted=ok,stop="requested next parent accepted" if ok else "minimum increment",
                    reached_one=bool(candidate[1]>=1))
                (out/"summary.json").write_text(json.dumps(dict(start=1396,cumulative=total,
                    additional=total-1396,paths=results,rtol=5e-9,atol=1e-14),indent=2)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
