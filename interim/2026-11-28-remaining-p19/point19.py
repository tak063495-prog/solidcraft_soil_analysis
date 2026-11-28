"""Observe the frozen p18 trial 877 with causal two-, three- and four-end extrapolants."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
P18 = HERE.with_name("2026-11-28-remaining-p18")


def run(out):
    case = json.loads((HERE.with_name("2026-11-26-remaining-p7")/"c01/original1024/input.json").read_text())
    law = case["joint_system"]["joints"][0]["material"]
    with gzip.open(P18/"c01/reload/trace.jsonl.gz","rt") as f:
        trace = [json.loads(s) for s in f if s.startswith("{")]
    # Read only accepted-parent observations available before this trial.
    states = {}
    for note in trace:
        if "p18_parent" in note:
            states[note["p18_parent"]] = note["points"][19]
            if note["attempt"] == 877:
                break
    assert note["attempt"]==877
    start, end = note["p18_parent"],note["nominal"]
    width = end-start
    status,jump,n,tau = states[start]
    n,tau,jump = map(np.asarray,(n,tau,jump))
    shear = tau - n*(n@tau)
    cap = law["friction_coefficient"]*law["normal_stiffness_kpa_m"]
    guard = 1/1024/64
    trials = []
    for count,degree in ((2,1),(3,2),(4,2)):
        times = list(states)[-count:]
        assert all(states[t][0] == status == "stick" for t in times)
        # Fit changes around the present endpoint to anchor the observed current state.
        scaled = (np.array(times[:-1])-start)/width
        design = np.column_stack([scaled**d for d in range(1,degree+1)])
        coeff = np.linalg.lstsq(design,np.array([states[t][1] for t in times[:-1]])-jump,rcond=None)[0]
        def margin(f):
            delta = np.array([f**d for d in range(1,degree+1)])@coeff
            normal = n@delta
            return np.linalg.norm(shear+law["tangential_stiffness_kpa_m"]*(delta-n*normal))+cap*(n@jump+normal)
        assert margin(0)<0<margin(1)
        root = start+width*brentq(margin,0,1,xtol=1e-14)
        trials.append(dict(endpoints=times,degree=degree,coefficients=coeff.tolist(),root=root,
                           end_distance=end-root,within_guard=min(end-root,root-start)<=guard))
    # Join the actual verdict only after the forecasts have been produced.
    rows = [r for r in trace if r.get("stage")=="reload" and "accepted" in r]
    event,point,*_ = rows[876]["audit"]["event"]
    assert point==19 and (rows[876]["from"],rows[876]["to"])==(start,end)
    for trial in trials:
        trial["signed_error"] = trial["root"]-event
    (out/"point19.json").write_text(json.dumps(dict(attempt=877,point=19,start=start,end=end,
        event=event,exclusion=guard,observation_only=True,trials=trials),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
