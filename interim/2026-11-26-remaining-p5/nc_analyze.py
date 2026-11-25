"""Audit recorded local line searches and the normal/subnormal transition."""
import json
import math
from pathlib import Path
import sys

import numpy as np


def run(out):
    load = lambda p: json.loads((out/p).read_bytes())
    trace = load("line-search.json")["trace"]
    substep = None
    search = []
    for row in trace:
        if "substep" in row:
            substep = row["substep"]
        elif "iteration" in row:
            search.append(dict(substep=substep, **row, candidates=[]))
        elif "alpha" in row:
            search[-1]["candidates"].append(row)
    final = search[-1]
    candidates = final["candidates"]
    assert [r["alpha"] for r in candidates] == [2.**-k for k in range(30)]
    def failures(candidate):
        result = []
        if not candidate["finite"]:
            result.append("finite")
        if abs(candidate["log_pc"]) >= 100:
            result.append("log_pc")
        if candidate["multiplier"] < 0:
            result.append("multiplier")
        if candidate["pressure"] <= 0:
            result.append("pressure")
        if candidate["candidate_norm"] is None or candidate["candidate_norm"] >= final["norm"]:
            result.append("residual_decrease")
        return result
    rejected = [failures(r) for r in candidates]
    assert all(rejected)
    history = load("pc-history.json")
    cross = next(i for i,r in enumerate(history) if r["scaled_square_subnormal"])
    assert cross > 0 and history[-1]["scaled_square_subnormal"]
    first_yield = {}
    for mode in ("units448","comparison480"):
        cases = [load(f"yield-{i}-{mode}.json") for i in range(1,18)]
        events = [next((r for r in c["trace"] if "first_yield" in r), None) for c in cases]
        entry = [r for r in events if r is not None]
        square = cases[0]["result"]["pc_squared"]
        first_yield[mode] = dict(calls_with_yield_bracket=len(entry), initial_event_abs_range=(
            [min(abs(r["first_yield"]) for r in entry),max(abs(r["first_yield"]) for r in entry)]
            if entry else None), first_event=entry[0] if entry else None, pc_squared=square,
            pc_square_ulp=math.ulp(square),relative_square_ulp=math.ulp(square)/square,
            completed_local_calls=sum(c["result"]["reason"] is None for c in cases))
    counts = load("attempts.json")
    retries = load("matched-retries.json")
    spent = (counts["extended_prefix"] + counts["homogeneous_point"]
             + counts["auxiliary_including_fold"] + counts["auxiliary_frozen_local_calls"]
             + len(retries) + sum(not r["accepted"] for r in retries))
    assert spent <= counts["budget"]
    report = dict(total_trials_including_local_calls=spent,
        component_order=["xx","yy","zz","xy","yz","zx","log_hardening","yield"],
        final_substep_zero_based=final["substep"],final_iteration_zero_based=final["iteration"],
        x=final["x"],residual=final["r"],newton_direction=final["dx"],residual_norm=final["norm"],
        linear_defect=final["linear_defect"],residual_largest_component=int(np.argmax(np.abs(final["r"]))),
        alpha_full=candidates[0],alpha_min=candidates[-1], failed_conditions=[
            dict(alpha=r["alpha"],failed=f) for r,f in zip(candidates,rejected)],
        per_condition={key:sum(key in f for f in rejected) for key in
                       ("finite","log_pc","multiplier","pressure","residual_decrease")},
        pc_first_subnormal_index_after_fold=cross,pc_first_pair=history[cross-1:cross+1],
        pc_end=history[-1],local_unit_diagnostics=first_yield)
    (out/"analysis.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]))
