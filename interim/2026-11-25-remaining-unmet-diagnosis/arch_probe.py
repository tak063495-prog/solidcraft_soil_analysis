"""Observe individual fold conditions on return from the unmodified classifier."""
from pathlib import Path
import json
import sys

import numpy as np

root = Path(sys.argv[1])
sys.path[:0] = [str(root/"interim/2026-11-24-t9d-benchmarks"),
                str(root/"interim/2026-11-22-impl-t7-physical-instability")]
from arch import beam_curve
from continuation import simple_point

def observe(frame,event,result):
    if frame.f_code is simple_point.__code__ and event=="return":
        sys.setprofile(None)
        d=frame.f_locals
        row={"kind":result["kind"], "reason":result.get("reason")}
        if "aug" in d:
            ev,scale,v,tol=d["ev"],d["scale"],d["v"],d["tol"]
            row.update(force_ratio=float(np.linalg.norm(ev.residual)/ev.denominator),
                       null_residual_over_scale=float(np.linalg.norm(ev.k@v)/scale),
                       unit_error=float(abs(v@v-1)),along=float(d["along"]),
                       eligible=bool(ev.eligible),same_event=bool(ev.event==d["ea"].event and all(p.event==ev.event for p in d["probes"])),
                       smallest_two_over_scale=(d["sv"][-2:]/scale).tolist(),
                       nullity=int(sum(d["sv"]<=tol*scale)),
                       transversality=float(abs(d["w"]@ev.q)/np.linalg.norm(ev.q)),
                       curvature_over_scale=float(abs(d["curvature"])/scale),
                       bordered_rank=int(result["bordered_rank"]),required_rank=len(v)+1,tol=tol)
            curvature_ratios={}
            for factor in (0.125,0.5,1.,2.,8.):
                step=factor*d["h"]
                plus=d["evaluate"](d["parent"],d["z"]+np.r_[step*v,0.])
                minus=d["evaluate"](d["parent"],d["z"]-np.r_[step*v,0.])
                curvature_ratios[str(factor)]=float(abs(d["w"]@(plus.k-minus.k)@v)/(2*step*scale))
            row["curvature_over_scale_by_difference_step_factor"]=curvature_ratios
        print(json.dumps(row),flush=True)

for n in (12,60,96):
    sys.setprofile(observe)
    r=beam_curve(n,extrapolate=True)
    print(json.dumps({"elements":n,"trials":len(r["trials"]),"accepted":sum(bool(t["accepted"]) for t in r["trials"]),"peak":r["fold_load"]}),flush=True)
