"""Trace the public deep arch on two grids with the existing T7 solvers."""
import argparse
from collections import Counter
import json
import gzip
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent


def beam_curve(n, *, extrapolate=False):
    """Continue from T7's verified restart, releasing full Newton matrices in batches."""
    from arch_elements import Arch
    from continuation import Controls, trace, simple_point
    from references import force_unit

    model = Arch(n)
    origin = np.zeros(len(model.free)+1)
    forward = origin.copy()
    forward[-1] = 1.
    points, trials, checkpoint = [origin], [], None
    highest = 0.
    while len(points)-1 < 1500:
        block = trace(model, {}, points[-1], forward,
                      Controls(0.18, min(4,1501-len(points)), 1e-9, 1e-10, 20),
                      restart=checkpoint, extrapolate=extrapolate,
                      attempt_limit=(2065*n+11)//12-len(trials) if extrapolate else None,
                      stop=lambda pts: pts[-1]["z"][-1] < max(highest, *(p["z"][-1] for p in pts))-0.45)
        new = [p["z"] for p in block["accepted"][1:]]
        points.extend(new)
        highest = max(highest, *(p[-1] for p in points[-5:]))
        trials.extend(block["trials"])
        reason = block["termination_reason"]
        checkpoint = block["restart"]
        del block
        if reason != "step budget" or not new:
            break
    curve = [dict(load=z[-1]*force_unit, downward_crown=model({},z).details["crown"]) for z in points]
    i = int(np.argmax([r["load"] for r in curve]))
    fold = (simple_point(model, {}, points[i-1], points[i+1], 1e-9)
            if 0 < i < len(points)-1 else {"kind":"D"})
    return dict(curve=curve,trials=trials,fold=fold,termination=reason,
                fold_load=fold["z"][-1]*force_unit if fold["kind"]=="simple fold" else None)


def summarize(case, peak_key, band):
    curve = [{"load": row["load"], "downward_crown": row["downward_crown"]}
             for row in case["curve"]]
    loads = [row["load"] for row in curve]
    i = int(np.argmax(loads))
    reversal = 0 < i < len(curve)-1 and loads[i-1] < loads[i] > loads[i+1]
    peak = case[peak_key]
    accepted = [r for r in case["trials"] if r["accepted"]]
    return {"fold": case["fold"]["kind"], "fold_load": peak,
            "termination": case["termination"], "accepted": len(accepted),
            "rejected": len(case["trials"])-len(accepted),
            "rejection_reasons": dict(Counter(r["reason"] for r in case["trials"] if not r["accepted"])),
            "max_accepted_relative_force": max(r["relative_force"] for r in accepted),
            "reversal": bool(reversal), "bracket": curve[max(0,i-1):i+2],
            "figure7_peak_and_reversal_passed": bool(reversal and peak is not None
                                                    and band[0] <= peak <= band[1]),
            "curve": curve}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--extrapolate", action="store_true", help="seed beam corrections from accepted points on the current arc plane")
    parser.add_argument("--root", type=Path, default=HERE.parents[1])
    args = parser.parse_args()
    sys.path.insert(0, str(args.root/"interim/2026-11-22-impl-t7-physical-instability"))
    from references import arch

    band = json.loads((HERE/"criteria.json").read_text())["arch"]["first_peak_load_reading"]
    results, histories = {}, {}
    columns = ["accepted", "reason", "lambda_parent", "delta_lambda", "ds",
               "iterations", "relative_force", "true_linear_residual"]
    for n in (12, 60, 96):
        print("Tracing", n, flush=True)
        results[str(n)], histories[str(n)] = {}, {}
        for method, calculate, key in [("beam", beam_curve, "fold_load"),
                                       ("shooting", arch, "estimated_peak_load")]:
            case = calculate(n, extrapolate=args.extrapolate) if method == "beam" else calculate(n)
            results[str(n)][method] = summarize(case, key, band)
            histories[str(n)][method] = [[r[k] for k in columns] for r in case["trials"]]
    encoded = json.dumps(results, indent=2, allow_nan=False)
    trace_bytes = json.dumps({"columns":columns,"cases":histories},
                             separators=(",",":"),allow_nan=False).encode()
    with args.output.with_name(args.output.stem+"-trials.json.gz").open("xb") as f:
        f.write(gzip.compress(trace_bytes,mtime=0))
    with args.output.open("x") as f:
        f.write(encoded+"\n")
