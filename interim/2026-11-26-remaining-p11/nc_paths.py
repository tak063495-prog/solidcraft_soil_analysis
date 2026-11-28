"""First-predictor Schur complement vs width on the physical continuation branch."""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-26-remaining-p9"))
from nc_scan import P8, p8

POINTS = (351, 272)


def predict(model, probe, state, lam, t, width, index, local):
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
    return row

def run(work, out):
    import hashlib
    import gzip
    old = p8.prior.old
    probe = work / "probe/target/release/examples"
    model = old.nc_extended.Extended(work/"inputs", probe/"p5_homogeneous")
    assert len(model.free)==348 and len(model.b)==48
    archive = {r["attempt"]:r for r in map(json.loads,(p8.P7/"nc/continuation.jsonl").read_text().splitlines())}
    p = old.read(P8/"nc/parent-632.json.gz")
    state,lam,t = p["active_material_states"],p["load"],np.array(p["secant"])
    u=np.ravel(p["displacement_m"])
    assert len(state)==1296
    width=archive[632]["ds_m"]
    count=0
    arc_limit=0.004
    def digest(state,lam,t,u):
        payload=json.dumps([state,lam,t.tolist(),u.tolist()],separators=(",",":"),allow_nan=False)
        return hashlib.sha256(payload.encode()).hexdigest()
    def attempt(path,index,state,lam,t,u,w):
        nonlocal count
        assert 1182+count+1 <= min(2000,16520)
        count+=1
        parent=digest(state,lam,t,u)
        result,it,reason=old.nc_branch.correct(model,state,lam,t,w,1e-12)
        assert parent==digest(state,lam,t,u)
        print(json.dumps(dict(path=path,parent=index,hash=parent,load=lam,width=w,
            accepted=result is not None,iterations=it,reason=reason),allow_nan=False),file=log)
        return result,it,reason
    def advance(w,result,it,u,lam):
        dz,state,diag=result
        # Only active free displacement entries participate in the continuation.
        u=u.copy()
        if len(u)==len(dz)-1:
            u+=dz[:-1]
        else:
            u[model.free]+=dz[:-1]
        return state,lam+dz[-1],dz/np.linalg.norm(dz),u,min(.04,w*(1.3 if len(it)<=8 else 1.))
    summary={}
    with (out/"global.jsonl").open("x",buffering=1) as log, gzip.open(out/"local.jsonl.gz","xt") as local:
        for i in range(632,642):
            result,it,reason=attempt("replay",i,state,lam,t,u,width)
            if i in archive:
                assert archive[i]["accepted"]==(result is not None)
            if result is None:
                width/=2
            else:
                state,lam,t,u,width=advance(width,result,it,u,lam)
        ref=json.loads((HERE.parent/"2026-11-26-remaining-p10/nc/summary.json").read_text())
        assert abs(lam-ref["load_end"])<1e-13
        initial=(state,lam,t,u,width)
        for path in ("ordinary","predicted"):
            state,lam,t,u,width=initial
            start=count
            length=chord=0.
            accepted=0
            while length<arc_limit:
                index=accepted+1
                if path=="predicted":
                    row=predict(model,probe,state,lam,t,width,index,local)
                    print(json.dumps(dict(frozen_prediction=row,path=path,parent_hash=digest(state,lam,t,u)),allow_nan=False),file=log)
                    width=row["predicted"]
                w=min(width,arc_limit-length)
                while True:
                    result,it,_=attempt(path,index,state,lam,t,u,w)
                    if result is not None:
                        break
                    w/=2
                    assert w>=model.settings["minimum_increment"]
                chord+=np.linalg.norm(result[0])
                state,lam,t,u,width=advance(w,result,it,u,lam)
                length+=w
                accepted+=1
                print(json.dumps(dict(path=path,accepted_parent=index,arc=length,chord=chord,
                    load=lam,displacement_change=float(np.linalg.norm(u-initial[3])),
                    state_hash=digest(state,lam,t,u)),allow_nan=False),file=log)
            summary[path]=dict(arc=length,chord=chord,accepted=accepted,rejected=count-start-accepted,
                attempts=count-start,load=lam,delta_lambda=lam-initial[1],
                displacement_change=float(np.linalg.norm(u-initial[3])),secant_lambda=t[-1])
    (out/"comparison.json").write_text(json.dumps(dict(start_load=initial[1],common_arc=arc_limit,
        initial_hash=digest(*initial[:4]),paths=summary,total_including_p10=1182+count),indent=2)+"\n")

if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
