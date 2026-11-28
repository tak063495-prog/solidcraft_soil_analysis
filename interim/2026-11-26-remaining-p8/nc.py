"""Follow the archived branch and inspect only the preselected parent/width series."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
from collections import Counter

import numpy as np
from freeze import HERE, P7, save, digest
sys.path.insert(0,str(P7))
import nc_pair as prior


def point(probe, record):
    proc = subprocess.run([probe/"p5_point"], input=json.dumps(dict(record,mode="homogeneous")),
        env=dict(os.environ,P5_POINT_TRACE="1"),text=True,capture_output=True,check=True)
    steps = []
    sub = 0
    for line in proc.stderr.splitlines():
        j = json.loads(line)
        if "substep" in j:
            sub = j["substep"]+1
        if j.get("iteration") != 0:
            continue
        assert j["x"][7] == 0
        J,r = np.array(j["jac_rows"]), np.array(j["r"])
        v = -np.linalg.solve(J[:7,:7],J[:7,7])
        stress, hard = J[7,:6] @ v[:6], J[7,6]*v[6]
        denominator = J[7,7] + stress + hard
        rhs = -r[7] + J[7,:7] @ np.linalg.solve(J[:7,:7],r[:7])
        np.testing.assert_allclose(rhs/denominator,j["dx"][7],rtol=5e-9,atol=1e-14)
        steps.append(dict(step=sub,stress=stress,hardening=hard,direct=J[7,7],
            denominator=denominator,rhs=rhs,direction=j["dx"][7]))
    return dict(point=record["point"],result=json.loads(proc.stdout),steps=steps,partition=record["count"])


def run(work,out):
    selection=json.loads((HERE/"nc/selection.json").read_text())
    archive=P7/"nc/continuation.jsonl"
    assert digest(archive)==selection["source_sha256"]
    selected={a:k for k,s in enumerate(selection["intervals"]) for g in s["selected"]
              for a in range(g["first"],g["last"]+1)}
    firsts = {g["first"] for s in selection["intervals"] for g in s["selected"]}
    probe=work/"probe/target/release/examples"
    os.environ["NC_WORKER"]=str(probe/"p5_homogeneous")
    old=prior.old
    m=old.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    watched=old.nc_extended.Extended(work/"inputs",HERE/"capture.py")
    p=old.read(prior.P5/"nc/parent.json.gz")
    state,lam,t=p["active_material_states"],p["load"],np.array(p["secant"])
    u=np.ravel(p["displacement_m"])
    assert len(state)==1296 and len(m.free)==348
    totals=[Counter() for _ in selection["intervals"]]
    out.mkdir(exist_ok=True,parents=True)
    ds=None
    with gzip.open(out/"materials.jsonl.gz","xt") as details, (out/"global.jsonl").open("x") as trace:
        for row in map(json.loads,archive.read_text().splitlines()):
            a=row["attempt"]
            if ds is None:
                ds=row["ds_m"]
            assert ds==row["ds_m"] and abs(lam-row["parent_lambda"])<1e-13
            if a in firsts:
                with gzip.open(out/f"parent-{a}.json.gz","xt") as f:
                    json.dump(dict(load=lam,secant=t.tolist(),displacement_m=u.tolist(),
                                   active_material_states=state),f)
            cap=work/f"captured-{a}.jsonl"
            os.environ["P8_CAPTURE"]=str(cap)
            result,it,reason=old.nc_branch.correct(watched if a in selected else m,state,lam,t,ds,1e-12)
            ok=result is not None
            evidence=dict(attempt=a,ds=ds,parent=lam,accepted=ok,reason=reason,
                          archive_iterations_equal=it==row["iterations"],iterations=it)
            print(json.dumps(evidence),file=trace,flush=True)
            assert ok==row["accepted"] and reason==row.get("reason"),a
            if a in selected:
                summary=totals[selected[a]]
                summary["trials"]+=1
                summary["rejected"]+=int(not ok)
                flags={i:[] for i in (351,272)}
                for e,call in enumerate(map(json.loads,cap.read_text().splitlines()),1):
                    for rec in call["points"]:
                        i=rec["point"]
                        assert rec["state"]==state[i]
                        info=point(probe,rec)
                        steps=info["steps"]
                        neg=[s for s in steps if s["direction"]<0]
                        flags[i].extend(neg)
                        summary[f"{i}_evaluations"]+=1
                        summary[f"{i}_newton_starts"]+=len(steps)
                        summary[f"{i}_negative_starts"]+=len(neg)
                        summary[f"{i}_negative_evaluations"]+=int(bool(neg))
                        summary[f"{i}_local_failure_evaluations"]+=int(bool(info["result"]["reason"]))
                        print(json.dumps(dict(attempt=a,evaluation=e,input=rec,worker_exit=call["exit"],
                            worker_error=call["error"],observation=info),allow_nan=False),file=details)
                for i,neg in flags.items():
                    summary[f"{i}_negative_trials"]+=int(bool(neg))
                    summary[f"{i}_negative_rejected"]+=int(bool(neg) and not ok)
                    summary[f"{i}_negative_accepted"]+=int(bool(neg) and ok)
                    summary[f"{i}_stress_hardening_sign_mechanism_trials"]+=int(any(
                        s["rhs"]<0 and s["denominator"]>0 and s["stress"]<0 and s["hardening"]>-s["stress"]
                        for s in neg))
            if ok:
                dz,state,diag=result
                assert dz[-1]<0 and diag["lambda_tangent"]<0
                lam+=dz[-1]
                u[m.free]+=dz[:-1]
                t=dz/np.linalg.norm(dz)
                ds=min(.04,ds*(1.3 if len(it)<=8 else 1.))
            else:
                ds/=2
                assert ds>=m.settings["minimum_increment"]
    save(out/"frequencies.json",dict(intervals=totals,load=lam,lambda_tangent=diag["lambda_tangent"],
        attempts=a-105,inherited=505,total=505+a-105,mesh_limit=16520,default_limit=2000,comparison480=0))


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
