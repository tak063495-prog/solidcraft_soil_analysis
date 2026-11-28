"""One original lambda-unknown corrector attempt at each stopped midpoint."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P13 = HERE.with_name("2026-11-27-remaining-p13")
P15 = HERE.with_name("2026-11-28-remaining-p15")
sys.path.insert(0,str(P13))
from continue_nc import advance,digest,old


def run(work,out):
    probe = work/"probe/target/release/examples"
    model = old.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    assert len(model.free)==348 and model.b.shape[:2]==(48,27)
    src = [json.loads(s) for s in (P13/"nc/global.jsonl").read_text().splitlines()]
    p15 = json.loads((P15/"nc/schur.json").read_text())
    total = p15["new_global_trials"] + json.loads((P13/"nc/comparison.json").read_text())["cumulative"]
    assert total==1364
    result = {}
    for path in ("ordinary","predicted"):
        end = old.read(P13/f"nc/{path}-end.json.gz")
        state,lam,t,u = end["active_material_states"],end["load"],np.array(end["secant"]),np.array(end["displacement_m"])
        parent = state,lam,t,u,end["next_width"]
        key = digest(parent)
        last = next(r for r in reversed(src) if r["path"]==path and "frozen_prediction" in r)
        assert last["parent_hash"] == end["state_hash"] == key
        (bracket,) = last["frozen_prediction"]["unresolved_brackets"]
        w = bracket["failed_midpoint"]
        assert bracket["point"] == 351 and w == p15["paths"][path]["width_m"]
        flat = np.zeros_like(model.base)
        flat[model.free] = t[:-1]
        slope = -np.einsum("eqai,ei->eqa",model.b,flat[model.dofs]).reshape(-1,6)
        strain = slope[351]*w
        assert abs(strain).sum() <= model.width
        os.environ.update(P16_STRAIN=json.dumps(strain.tolist()),P16_LU=str(probe/"local_lu"),
                          P16_PYTHON=sys.executable,P16_COMPENSATE=str(HERE/"compensate.py"))
        log = out/(path+"-local-direction.jsonl")
        os.environ["P16_RECORD"] = str(log)
        point = dict(mode="homogeneous",point=351,count=1,state=state[351],law=model.laws[model.ids[351]],
                     strain=strain.tolist())
        with (out/(path+"-local-trace.jsonl")).open("w") as stream:
            proc = subprocess.run([probe/"p5_point"],env=dict(os.environ,P5_POINT_TRACE="1"),
                input=json.dumps(point),text=True,stdout=subprocess.PIPE,stderr=stream,check=True)
        rows = [json.loads(s) for s in log.read_text().splitlines()]
        assert len(rows)==1 and rows[0]["passes"]
        assert rows[0]["refined"][-1] == p15["paths"][path]["lu"]
        assert rows[0]["schur"] == p15["paths"][path]["variants"]["both"]["direction"]
        total += 1
        glog = out/(path+"-global-direction.jsonl")
        os.environ["P16_RECORD"] = str(glog)
        sol,it,reason = old.nc_branch.correct(model,state,lam,t,w,1e-12)
        assert digest(parent)==key
        updated = advance(model,parent,w,sol,it) if sol is not None else parent
        if sol is not None:
            with gzip.open(out/(path+"-accepted.json.gz"),"xt") as f:
                json.dump(dict(active_material_states=updated[0],load=updated[1],secant=updated[2].tolist(),
                          displacement_m=updated[3].tolist(),next_width=updated[4],state_hash=digest(updated)),f)
        result[path] = dict(parent=last["frozen_prediction"]["parent"],parent_hash=key,width_m=w,point=351,
             local=json.loads(proc.stdout),global_direction_calls=len(glog.read_text().splitlines()) if glog.exists() else 0,
             iterations=it,accepted=sol is not None,reason=reason,parent_updated=sol is not None,
             end_hash=digest(updated),lambda_before=lam,lambda_after=updated[1],additional_global_trials=1,cumulative=total)
        (out/"parents.json").write_text(json.dumps(dict(paths=result,rtol=5e-9,atol=1e-14,
               start_count=1364,cumulative=total),indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
