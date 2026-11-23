"""Observe native open/reload trials without changing the production defaults."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import subprocess
import sys

import ijson

from c01 import HERE, PHASES, ROUTES, observation


def run(steps, route, solver, folder):
    folder.mkdir(parents=True, exist_ok=False)
    case_file = folder / "case.json"
    locate, substeps = ROUTES[route]
    probe = HERE.parent / "2026-11-22-impl-c01-c03-nonlinear/probe_c01.py"
    subprocess.run([sys.executable, probe, "--steps", str(steps), "--event-proposal", "--local-substeps",
                    str(substeps), "--native-input", "--output", case_file], check=True)
    case = json.loads(case_file.read_bytes())
    control = case["joint_system"]["global_step_control"]["path_error"]
    control["event_comparison"] = "observe" if locate else "bisect_rejections"
    case_file.write_text(json.dumps(case, allow_nan=False), encoding="utf-8")
    inputs = HERE.parents[1] / "examples"
    raw = folder / "run.json"
    subprocess.run([solver, inputs/"joint-water-column-tet10-mesh.json", inputs/"ground-joint-water-mc.json",
                    case_file, "reload", "1.0", raw], check=True, timeout=21600)
    counts, last, rejected = Counter(), {}, {}
    trials = {phase: [] for phase in PHASES}
    for prefix in ("geostatic.attempts.item", "attempts.item"):
        with raw.open("rb") as stream:
            for row in ijson.items(stream, prefix, use_float=True):
                stage = "geostatic" if prefix.startswith("geostatic") else row["stage"]
                counts[stage] += 1
                if counts[stage] > 2000 or last.get(stage, 0) != row["progress_from"]:
                    raise ValueError("stage budget or accepted history differs")
                if row["accepted"]:
                    last[stage] = row["progress_to"]
                else:
                    if not row["committed_state_unchanged_on_rejection"]:
                        raise ValueError("rejection changed committed history")
                    rejected[stage] = [row["progress_from"], row["progress_to"]]
                if stage in trials:
                    audit = row.get("joint_path_error")
                    if audit:
                        audit.setdefault("first_gap_zero", None)
                    trials[stage].append(observation(row, case["settings"]["relative_force_tolerance"]))
    with raw.open("rb") as stream:
        stop = next(ijson.items(stream, "failure_reason"))
    return dict(steps=steps, route=route, stop=stop, counts=dict(counts), last_accepted=last,
                last_rejected=rejected, trials=trials, passed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, choices=(512,1024), required=True)
    parser.add_argument("--route", choices=ROUTES, required=True)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    value = run(args.steps, args.route, args.solver.resolve(), args.work)
    with gzip.open(args.output, "xt", encoding="utf-8") as stream:
        json.dump(value, stream, allow_nan=False, separators=(",", ":"))
