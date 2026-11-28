"""Scalar consistency slope of the original 8x8 local linearization."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np


def run(out):
    report=[]
    for entry in json.loads((out/"summary.json").read_text()):
        a=entry["archive"]["attempt"]
        results={}
        with gzip.open(out/f"substeps-{a}.jsonl.gz","rt") as stream:
            for line in stream:
                row=json.loads(line)
                j=row["record"]
                if j.get("iteration")!=0:
                    continue
                assert j["x"][7]==0
                J=np.array(j["jac_rows"])
                v=-np.linalg.solve(J[:7,:7],J[:7,7])
                stress=float(J[7,:6] @ v[:6])
                hardening=J[7,6]*v[6]
                slope=J[7,7]+stress+hardening
                rhs=-j["r"][7]+J[7,:7] @ np.linalg.solve(J[:7,:7],j["r"][:7])
                np.testing.assert_allclose(rhs/slope,j["dx"][7],rtol=5e-9)
                results.setdefault(row["evaluation"],[]).append(dict(step=row["step"],
                    stress_part=stress,log_pc_part=hardening,direct=J[7,7],schur=slope,
                    rhs=rhs,direction=j["dx"][7]))
        for e in entry["evaluations"]:
            steps=results[e["evaluation"]]
            crossing=[s["step"] for s in e["substeps"] if "crossing_fraction" in s]
            report.append(dict(width=entry["archive"]["ds_m"],evaluation=e["evaluation"],
                partition=e["count"],executed=len(e["substeps"]),first_yield_splits=crossing,
                negative_initial_direction=[s["step"] for s in steps if s["direction"]<0],
                steps=steps))
    (out/"consistency-slope.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
