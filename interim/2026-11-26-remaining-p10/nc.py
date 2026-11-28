"""First-predictor Schur complement vs width on the physical continuation branch."""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-26-remaining-p9"))
from nc_scan import P8, p8

POINTS = (351, 272)


def run(work, out):
    old = p8.prior.old
    probe = work/"probe/target/release/examples"
    model = old.nc_extended.Extended(work/"inputs", probe/"p5_homogeneous")
    archive = {r["attempt"]:r for r in map(json.loads, (p8.P7/"nc/continuation.jsonl").read_text().splitlines())}
    parent = old.read(P8/"nc/parent-632.json.gz")
    state, lam = parent["active_material_states"], parent["load"]
    t = np.array(parent["secant"])
    assert len(state)==1296 and len(model.free)==348 and len(model.b)==48
    count = 0
    summaries = []
    with (out/"global.jsonl").open("x", buffering=1) as global_log, (out/"local.jsonl").open("x", buffering=1) as local:
        def trial(width, label):
            nonlocal count
            assert 1153+count+1 <= 2000
            count += 1
            result, iterations, reason = old.nc_branch.correct(model,state,lam,t,width,1e-12)
            print(json.dumps(dict(parent=index, load=lam, control=label, width=width, accepted=result is not None,
                                  reason=reason, iterations=iterations), allow_nan=False),file=global_log)
            return result, iterations, reason

        index = 632
        width = archive[index]["ds_m"]
        while len(summaries)<6:
            if index in archive:
                assert abs(lam - archive[index]["parent_lambda"])<1e-13 and width==archive[index]["ds_m"]
            if index>=634:
                flat = np.zeros_like(model.base)
                flat[model.free] = t[:-1]
                slope = -np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)
                roots = []
                cuts = {}
                for i in POINTS:
                    cut = model.width/np.abs(slope[i]).sum()
                    cuts[str(i)] = cut
                    cache = {}
                    def sample(w):
                        if w not in cache:
                            strain = slope[i]*w
                            needed = max(1, int(np.ceil(np.abs(strain).sum()/model.width)))
                            partition = 1 if needed==1 else 2*needed
                            rec = dict(point=i,law=model.laws[model.ids[i]],state=state[i],strain=strain.tolist(),count=partition)
                            obs = p8.point(probe, rec)
                            step = next((s for s in obs["steps"] if s["step"]==1), None)
                            cache[w] = dict(width=w,partition=partition,first=step)
                            print(json.dumps(dict(parent=index,load=lam,point=i,**cache[w],reason=obs["result"]["reason"]),
                                             allow_nan=False),file=local)
                        return cache[w]
                    grid = np.r_[np.linspace(.05,.99,36),.995,.999,.99999,1-1e-9,1+1e-9,1.00001,1.001,1.005,
                                 np.linspace(1.01,1.99,36)]
                    vals = [sample(cut*f) for f in grid]
                    for a,b in zip(vals,vals[1:]):
                        if a["partition"]!=b["partition"] or not a["first"] or not b["first"]:
                            continue
                        def tangent(v):
                            return v["first"]["denominator"]
                        if tangent(a)*tangent(b)>=0:
                            continue
                        for _ in range(14):
                            mid=sample((a["width"]+b["width"])/2)
                            assert mid["first"] and a["partition"]==mid["partition"]==b["partition"]
                            if tangent(mid)*tangent(a)>0:
                                a=mid
                            else:
                                b=mid
                        roots.append(dict(point=i, left=a, right=b))
                proposed=min(width,.9*min([*cuts.values(),*[r["left"]["width"] for r in roots]]))
                row=dict(parent=index,load=lam,partition_boundaries=cuts,sign_brackets=roots,ordinary=width,predicted=proposed)
                print(json.dumps(dict(frozen_prediction=row),allow_nan=False),file=global_log)
                def search(w,label):
                    initial=count
                    while True:
                        outcome,it,why=trial(w,label)
                        if outcome is not None:
                            return outcome,it,w,count-initial
                        w/=2
                        assert w>=model.settings["minimum_increment"]
                _,_,wp,npred=search(proposed,"predicted")
                result,it,accepted,no=search(width,"ordinary")
                row.update(predicted_accepted_width=wp,predicted_attempts=npred,ordinary_accepted_width=accepted,ordinary_attempts=no)
                summaries.append(row)
            else:
                result,it,why=trial(width,"replay")
                assert (result is not None)==archive[index]["accepted"] and why==archive[index].get("reason")
                if result is None:
                    width/=2
                    index+=1
                    continue
                accepted=width
                no=1
            dz,state,diag=result
            assert dz[-1]<0 and diag["lambda_tangent"]<0
            lam+=dz[-1]
            t=dz/np.linalg.norm(dz)
            width=min(.04,accepted*(1.3 if len(it)<=8 else 1.))
            index+=no
    (out/"summary.json").write_text(json.dumps(dict(parents=summaries,load_end=lam,tangent_end=diag["lambda_tangent"],
        new_global_calls=count,total_including_p9=1153+count,mesh_budget=16520,default_budget=2000,
        local_calls=sum(1 for _ in (out/"local.jsonl").open()),passed=False),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
