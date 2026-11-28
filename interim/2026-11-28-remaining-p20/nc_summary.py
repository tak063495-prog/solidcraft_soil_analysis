"""Tabulate every global verdict and both lambda slopes in accepted-arc order."""
from collections import Counter
import csv
import json
from pathlib import Path
import sys


def span(values):
    return [min(values),max(values)] if values else None


def run(out):
    state=json.loads((out/"summary.json").read_text())
    data=[json.loads(r) for r in (out/"global.jsonl").read_text().splitlines()]
    rows=[]
    for a,b in zip(data[::2],data[1::2],strict=True):
        assert a["total"]==b["total"]
        directions=b["point351_live"]
        d=directions[0] if directions else {}
        gl=b["iterations"][-1] if b["iterations"] else {}
        rows.append(dict(path=a["path"],trial=a["total"],accepted=b["accepted"],arc_before_m=a["arc_m"],
            attempted_width_m=a["width_m"],arc_after_m=b["arc_m"],lambda_before=a["lambda_before"],
            lambda_after=b["lambda_after"],delta_lambda=b["delta_lambda"],predictor_secant=a["predictor_secant_lambda"],
            endpoint_tangent=(b["equilibrium_tangent"] or {}).get("lambda_tangent"),relative_force=gl.get("relative_force"),
            arc_residual_m=gl.get("arc_residual_m"),point351_denominator=d.get("denominator"),
            point351_multiplier_direction=d.get("direction"),matched=d.get("matched"),reason=b["reason"]))
    with (out/"arc.csv").open("w") as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    result={}
    for name in state["paths"]:
        rr=[r for r in rows if r["path"]==name]
        ok=[r for r in rr if r["accepted"]]
        summary=dict(first=ok[0],last=ok[-1],last_trial=rr[-1],max_relative_force=max(r["relative_force"] for r in ok),
            max_abs_arc_residual_m=max(abs(r["arc_residual_m"]) for r in ok),
            tangent_range=span([r["endpoint_tangent"] for r in ok]),
            accepted_delta_lambda_range=span([r["delta_lambda"] for r in ok]),
            rejection_reasons=Counter(r["reason"] for r in rr if not r["accepted"]))
        for flag in (True,False):
            cohort=[r for r in rr if r["accepted"]==flag]
            summary["accepted_local" if flag else "rejected_local"]=dict(count=len(cohort),
                measured=sum(r["matched"] is not None for r in cohort),
                matched=sum(r["matched"] is True for r in cohort),denominator_range=span([r["point351_denominator"]
                    for r in cohort if r["point351_denominator"] is not None]),multiplier_range=span([
                    r["point351_multiplier_direction"] for r in cohort if r["point351_multiplier_direction"] is not None]))
        result[name]=summary
    (out/"ranges.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
