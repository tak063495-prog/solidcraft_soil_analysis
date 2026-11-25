"""Replay only attempts 106–109 from the archived complete committed state."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P5 = HERE.parent / "2026-11-26-remaining-p5"
sys.path.insert(0, str(P5))
import nc_trace as p5


def run(work, out):
    out.mkdir(parents=True)
    binary = work/"probe/target/release/examples"
    os.environ["NC_WORKER"] = str(binary/"p5_homogeneous")
    m = p5.nc_extended.Extended(work/"inputs", HERE/"capture351.py")
    parent = p5.read(P5/"nc/parent.json.gz")
    state, lam, secant = parent["active_material_states"], parent["load"], np.array(parent["secant"])
    assert len(state) == 1296 and lam == 0.5481398948309102
    archived = [json.loads(line) for line in (p5.P4/"nc/extended-004/trials.jsonl").open()]
    report = []
    for row in archived[105:109]:
        attempt = row["attempt"]
        assert row["parent_lambda"] == lam
        cap = out/f"inputs-{attempt}.jsonl"
        os.environ["P6_CAPTURE"] = str(cap)
        result, iterations, reason = p5.nc_branch.correct(m,state,lam,secant,row["ds_m"],1e-12)
        assert iterations == row["iterations"] and reason == row.get("reason")
        assert (result is not None) == row["accepted"]
        if result is not None:
            dz, _, diag = result
            assert lam+dz[-1] == row["load"] and all(diag[k] == row[k] for k in diag)
        evaluations = []
        with gzip.open(out/f"substeps-{attempt}.jsonl.gz","wt") as stream:
            for j,line in enumerate(cap.read_text().splitlines(), 1):
                local = json.loads(line)
                assert local["state"] == state[351]
                proc = subprocess.run([binary/"p5_point"],text=True,capture_output=True,check=True,
                    env=dict(os.environ,P5_POINT_TRACE="1"),input=json.dumps(dict(local,mode="homogeneous")))
                point = json.loads(proc.stdout)
                assert bool(point["reason"]) == bool(local["worker_exit"])
                step = None
                summaries = []
                for item in map(json.loads,proc.stderr.splitlines()):
                    if "substep" in item:
                        step = item["substep"]+1
                        summary = dict(step=step)
                        summaries.append(summary)
                    if "f0" in item and "predictor" in item:
                        summary.update(f0=item["f0"], f1=item["f1"], loading=item["loading"])
                    if "first_yield" in item:
                        summary.update(crossing_fraction=item["fraction"], yield_at_crossing=item["first_yield"])
                    if item.get("iteration") == 0:
                        J = np.array(item["jac_rows"])
                        r = np.array(item["r"])
                        # Eliminate the six stresses and log-hardening variable.
                        A,b,c = J[:7,:7], J[:7,7], J[7,:7]
                        S = J[7,7] - c @ np.linalg.solve(A,b)
                        rhs = -r[7] + c @ np.linalg.solve(A,r[:7])
                        dx = item["dx"][7]
                        np.testing.assert_allclose(rhs/S,dx,rtol=5e-9,atol=1e-14)
                        summary.update(elastic_predictor_yield=r[7], multiplier_dx=dx,
                                       schur_tangent=S, reduced_rhs=rhs, linear_defect=item["linear_defect"])
                    stream.write(json.dumps(dict(evaluation=j,step=step,record=item),allow_nan=False)+"\n")
                evaluations.append(dict(evaluation=j, count=local["count"], reason=point["reason"],
                                        substeps=summaries))
        report.append(dict(archive=row,evaluations=evaluations))
    assert [r["archive"]["ds_m"] for r in report] == [.04,.02,.01,.005]
    (out/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]))
