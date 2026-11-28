"""Continue from the P5 parent, then replay the first adjacent different-point reject/accept pair."""
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P6 = HERE.parent/"2026-11-26-remaining-p6"
P5 = HERE.parent/"2026-11-26-remaining-p5"
sys.path.insert(0,str(P5))
import nc_trace as old
sys.path.insert(0,str(P6))
import score_nc


def save(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def trace(binary, record, stream, evaluation):
    proc = subprocess.run([binary], text=True, capture_output=True, check=True,
                          env=dict(os.environ,P5_POINT_TRACE="1"),
                          input=json.dumps(dict(record,mode="homogeneous")))
    point = json.loads(proc.stdout)
    assert bool(point["reason"]) == bool(record["worker_exit"])
    summaries = []
    step = None
    for item in map(json.loads,proc.stderr.splitlines()):
        if "substep" in item:
            step = item["substep"]+1
            summaries.append(dict(step=step))
        if "predictor" in item:
            summaries[-1].update(f0=item["f0"],f1=item["f1"],loading=item["loading"])
        if "first_yield" in item:
            summaries[-1].update(crossing_fraction=item["fraction"])
        if item.get("iteration")==0:
            J,r = np.array(item["jac_rows"]),np.array(item["r"])
            slope = J[7,7] - J[7,:7] @ np.linalg.solve(J[:7,:7],J[:7,7])
            rhs = -r[7] + J[7,:7] @ np.linalg.solve(J[:7,:7],r[:7])
            np.testing.assert_allclose(rhs/slope,item["dx"][7],rtol=5e-9,atol=1e-14)
            summaries[-1].update(elastic_predictor_yield=r[7], multiplier_dx=item["dx"][7],
                                 schur_tangent=slope,reduced_rhs=rhs)
        print(json.dumps(dict(evaluation=evaluation,step=step,record=item),allow_nan=False),file=stream)
    return dict(evaluation=evaluation,count=record["count"],reason=point["reason"],substeps=summaries)


def run(work,out):
    out.mkdir(parents=True)
    archive = old.P4/"nc/extended-004/trials.jsonl"
    rows = list(map(json.loads,archive.read_text().splitlines()))
    prior = old.read(P5/"nc/parent.json.gz")
    state,lam = prior["active_material_states"],prior["load"]
    t = np.array(prior["secant"])
    u = np.array(prior["displacement_m"]).ravel()
    binary = work/"probe/target/release/examples"
    os.environ["NC_WORKER"] = str(binary/"p5_homogeneous")
    m = old.nc_extended.Extended(work/"inputs",binary/"p5_homogeneous")
    assert len(state)==1296 and len(m.free)==348 and t[-1]<0
    # Re-use the actual restart width and the P4 restart update after each decision.
    # Two extra global calls are reserved for a bitwise replay of the selected pair.
    ds = rows[105]["ds_m"]
    inherited = 400 + 105
    prefix = []
    pair = None
    with (out/"continuation.jsonl").open("x",buffering=1) as log:
        for j in range(2000-inherited-2):
            result,it,reason = old.nc_branch.correct(m,state,lam,t,ds,1e-12)
            trial = dict(attempt=106+j,parent_lambda=lam,ds_m=ds,iterations=it,accepted=result is not None)
            if result is None:
                trial["reason"] = reason
            else:
                dz,new,diag = result
                trial.update(diag,load=lam+dz[-1],delta_lambda=dz[-1])
                assert dz[-1]<0 and diag["lambda_tangent"]<0 and lam+dz[-1]>0
            print(json.dumps(trial,allow_nan=False),file=log)
            if result is not None and prefix and not prefix[-1]["accepted"]:
                match = re.search(r"active point (\d+)",prefix[-1]["reason"])
                if (match and int(match[1])!=351 and "line search failed" in prefix[-1]["reason"]
                    and lam!=prior["load"] and prefix[-1]["parent_lambda"]==lam):
                    point = int(match[1])
                    pair = [prefix[-1],trial]
                    prefix.append(trial)
                    break
            prefix.append(trial)
            if result is None:
                ds /= 2
                assert ds>=m.settings["minimum_increment"],reason
            else:
                state,lam = new,trial["load"]
                u[m.free] += dz[:-1]
                t = dz/np.linalg.norm(dz)
                ds = min(.04,ds*(1.3 if len(it)<=8 else 1.))
    if pair is None:
        raise RuntimeError("no adjacent pair within the original attempt limit")
    with gzip.open(out/"parent.json.gz","wt") as f:
        json.dump(dict(load=lam,secant=t.tolist(),displacement_m=u.tolist(),active_material_states=state),f)
    m = old.nc_extended.Extended(work/"inputs",HERE/"capture_point.py")
    report = []
    for row in pair:
        assert row["parent_lambda"]==lam
        attempt = row["attempt"]
        cap = out/f"inputs-{attempt}.jsonl"
        os.environ.update(P7_CAPTURE=str(cap),P7_POINT=str(point))
        result,it,reason = old.nc_branch.correct(m,state,lam,t,row["ds_m"],1e-12)
        if it!=row["iterations"] or reason!=row.get("reason") or (result is not None)!=row["accepted"]:
            save(out/"mismatch.json",dict(expected=row,iterations=it,reason=reason,accepted=result is not None))
            raise AssertionError(f"archive mismatch at attempt {attempt}")
        evaluations = []
        with gzip.open(out/f"substeps-{attempt}.jsonl.gz","wt") as stream:
            for k,line in enumerate(cap.read_text().splitlines(),1):
                record = json.loads(line)
                assert record["state"]==state[point] and record["point"]==point
                evaluations.append(trace(binary/"p5_point",record,stream,k))
        report.append(dict(archive=row,point=point,evaluations=evaluations))
        if result is not None:
            dz,state,diag = result
            u[m.free] += dz[:-1]
            lam += dz[-1]
            assert lam==row["load"] and all(diag[k]==row[k] for k in diag)
            t = dz/np.linalg.norm(dz)
    save(out/"summary.json",report)
    original = {r["attempt"]:r for r in rows}
    historical = dict(compared=len(prefix),
        same_width=sum(r["ds_m"]==original[r["attempt"]]["ds_m"] for r in prefix),
        same_decision=sum(r["accepted"]==original[r["attempt"]]["accepted"] for r in prefix),
        same_reason=sum(r.get("reason")==original[r["attempt"]].get("reason") for r in prefix),
        bitwise_equal_iterations=sum(r["iterations"]==original[r["attempt"]]["iterations"] for r in prefix),
        max_parent_lambda_difference=max(abs(r["parent_lambda"]-original[r["attempt"]]["parent_lambda"]) for r in prefix),
        selected_archived=[original[r["attempt"]] for r in pair],selected_new=pair)
    save(out/"historical-comparison.json",historical)
    save(out/"replay.json",dict(point=point,inherited=inherited,continuation=len(prefix),
        global_trials=len(prefix)+2,total_including_inherited=inherited+len(prefix)+2,
        pair_replayed_exactly=True,unit480_trials=0,original_cap=2000,mesh_cap=16520,
        local_replays=sum(len(r["evaluations"]) for r in report)))
    score_nc.run(out)


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
