"""Restart each branch by verified replay, then continue over the preregistered arc."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P11 = HERE.parent / "2026-11-26-remaining-p11"
sys.path.insert(0, str(P11))
from nc_paths import p8, P8
from nc_prediction import predict

old = p8.prior.old


def digest(parent):
    state, lam, t, u, _ = parent
    data = json.dumps([state,lam,t.tolist(),u.tolist()], separators=(",",":"), allow_nan=False)
    return hashlib.sha256(data.encode()).hexdigest()


def advance(model, parent, w, result, iterations):
    _, lam, _, u, _ = parent
    dz, state, _ = result
    u = u.copy()
    if len(u) == len(dz)-1:
        u += dz[:-1]
    else:
        u[model.free] += dz[:-1]
    next_width = min(.04, w * (1.3 if len(iterations)<=8 else 1.))
    return state, lam+dz[-1], dz/np.linalg.norm(dz), u, next_width


def run(work, out):
    arc = json.loads((HERE/"design.json").read_text())["nc"]["common_arc_m"]
    probe = work/"probe/target/release/examples"
    model = old.nc_extended.Extended(work/"inputs", probe/"p5_homogeneous")
    assert len(model.free)==348 and len(model.b)==48
    base = old.read(P8/"nc/parent-632.json.gz")
    log11 = [json.loads(s) for s in (P11/"nc/global.jsonl").read_text().splitlines()]
    ref = json.loads((P11/"nc/comparison.json").read_text())
    initial = (base["active_material_states"],base["load"],np.array(base["secant"]),
               np.ravel(base["displacement_m"]),log11[0]["width"])
    assert len(initial[0]) == 1296
    attempt_logs=sorted((HERE/"nc").glob("*-stop/global.jsonl"))
    initial_attempts=sum("accepted" in json.loads(s) for path in attempt_logs for s in path.read_text().splitlines())
    cumulative=ref["total_including_p10"]+initial_attempts
    paths = {}
    with (out/"global.jsonl").open("x", buffering=1) as log, gzip.open(out/"local.jsonl.gz","xt") as local:
        def record(**row):
            print(json.dumps(row,allow_nan=False),file=log)
        def solve(path, parent, width, replay=False):
            nonlocal cumulative
            cumulative += 1
            assert cumulative <= min(2000,16520)
            state, lam, direction, _, _ = parent
            before = digest(parent)
            result, it, reason = old.nc_branch.correct(model,state,lam,direction,width,1e-12)
            assert before == digest(parent)
            record(path=path, replay=replay, parent_hash=before, width=width, load=lam,
                   accepted=result is not None, iterations=it, reason=reason, cumulative=cumulative)
            return result,it
        def restore(path,parent):
            for row in (r for r in log11 if r["path"]==path):
                if "accepted_parent" in row:
                    assert digest(parent)==row["state_hash"]
                if "accepted" not in row:
                    continue
                assert digest(parent)==row["hash"]
                result,it=solve(path,parent,row["width"],True)
                assert row["accepted"]==(result is not None)
                if result is not None:
                    parent=advance(model,parent,row["width"],result,it)
            return parent
        common=restore("replay",initial)
        assert digest(common)==ref["initial_hash"]
        for path in ("ordinary","predicted"):
            parent=start=restore(path,common)
            assert start[1]==ref["paths"][path]["load"]
            stage_start=cumulative
            length=chord=0.
            accepted=0
            cuts=[]
            stop=None
            unresolved_parents=0
            while length<arc:
                state,lam,t,_,ordinary=parent
                row=predict(model,probe,state,lam,t,ordinary,accepted+1,local)
                cuts.append(row["partition_boundaries"])
                record(path=path,parent_hash=digest(parent),arc_before=length,frozen_prediction=row)
                unresolved_parents+=int(bool(row["unresolved_brackets"]))
                if not row["unique_width"]:
                    stop="unresolved p10 sign bracket can determine the minimum"
                    break
                w=min(ordinary if path=="ordinary" else row["predicted"],arc-length)
                while True:
                    result,it=solve(path,parent,w)
                    if result is not None:
                        break
                    w/=2
                    assert w>=model.settings["minimum_increment"]
                parent=advance(model,parent,w,result,it)
                accepted+=1
                length+=w
                chord+=np.linalg.norm(result[0])
                record(path=path,accepted_parent=accepted,arc=length,chord=chord,
                       state_hash=digest(parent),load=parent[1],next_width=parent[-1])
            paths[path]=dict(stop=stop,unresolved_parents=unresolved_parents,start_hash=digest(start),end_hash=digest(parent),start_load=start[1],
                             end_load=parent[1],delta_lambda=parent[1]-start[1],arc=length,chord=chord,
                             displacement_change_m=float(np.linalg.norm(parent[3]-start[3])),
                             boundary_ranges_m={k:[min(c[k] for c in cuts),max(c[k] for c in cuts)] for k in cuts[0]},
                             first_boundaries_m=cuts[0],last_boundaries_m=cuts[-1],
                             accepted=accepted,attempts=cumulative-stage_start,rejected=cumulative-stage_start-accepted)
            with gzip.open(out/(path+"-end.json.gz"),"xt") as f:
                state,lam,t,u,width=parent
                json.dump(dict(active_material_states=state,load=lam,secant=t.tolist(),displacement_m=u.tolist(),
                               next_width=width,state_hash=digest(parent)),f,allow_nan=False)
    (out/"comparison.json").write_text(json.dumps(dict(common_arc=arc,paths=paths,
            initial_stop_attempts=initial_attempts,cumulative_including_replay=cumulative),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]))
