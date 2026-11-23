"""Opt-in 12 Tet10 independent path comparison with the original stage budget."""
import argparse
from collections import Counter
import gzip
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROUTES = {"event4": (True, 4), "bisect4": (False, 4), "event1": (True, 1)}
PHASES = ("open", "reload")


def proposal_observer(original, locate):
    """Keep the accepted-step rule identical. Only the rejected-step locator varies."""
    def propose(audit, gaps, progress, steps, minimum):
        h = audit["actual_increment"]
        candidate = original(audit, gaps, progress, steps, minimum)
        roots = []
        for half, (start, end) in enumerate(zip(gaps[:-1], gaps[1:])):
            for gp in np.flatnonzero((start > 0) != (end > 0)):
                fraction = (half + start[gp] / (start[gp] - end[gp])) / 2
                roots.append((float(progress - h + fraction * h), int(gp)))
        audit["first_gap_zero"] = min(roots, default=None)
        return h / 2 if audit["reject"] and not locate else candidate
    return propose


def observation(row, force_tolerance=1e-10):
    start, end = row["progress_from"], row["progress_to"]
    result = dict(start=start, end=end, accepted=row["accepted"], interval="unknown")
    audit = row.get("joint_path_error")
    if audit is None:
        if row["accepted"]:
            raise ValueError("accepted trial lacks a path audit")
        return result
    statuses = np.array([audit[k] for k in
                         ("old_status", "midpoint_status", "fine_status", "coarse_status")])
    root = audit["first_gap_zero"]
    switched = bool((statuses != statuses[0]).any())
    result.update(interval="crossing" if root or switched else "smooth", gap_zero=root,
                  state_switch=switched, proposed_increment=audit["proposed_increment"])
    for name, key, budget in (("slip", "slip_errors_m", "slip_budget_m"),
                              ("dissipation", "dissipation_errors_kpa_m", "dissipation_budget_kpa_m")):
        array = np.abs(audit[key])
        index = np.unravel_index(array.argmax(), array.shape)
        value = float(array[index])
        if row["accepted"] and value > audit[budget]:
            raise ValueError("accepted trial exceeded path budget")
        result[name] = dict(value=value, gp=int(index[0]))
    if row["accepted"] and max(row["relative_force_residual"], audit["midpoint_relative_force_residual"],
                               audit["fine_relative_force_residual"]) > force_tolerance:
        raise ValueError("accepted equilibrium exceeds force tolerance")
    return result


def run(steps, route):
    path = HERE.parent / "2026-12-12-c01-switch-event-localization/study.py"
    spec = importlib.util.spec_from_file_location("c01_window_study", path)
    study = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(study)
    inputs = HERE.parents[1] / "examples"
    mesh = study.read(inputs / "joint-water-column-tet10-mesh.json")
    ground = study.fem.normalize_ground(study.read(inputs / "ground-joint-water-mc.json"))
    stages = study.fem.normalize_stage_case(study.read(inputs / "stages-joint-water-cycle.json"))
    case = study.prior.configuration(stages, steps, 5e-5)
    case["settings"]["minimum_increment"] = min(case["settings"]["minimum_increment"], 1/(steps*65536))
    locate, substeps = ROUTES[route]
    control = dict(study.PATH_ERROR, event_proposal=True, minimum_contact_substeps=substeps)
    case["joint_system"]["global_step_control"] = dict(path_error=control)
    reference = study.reference  # private module instance for this run
    reference.event_proposal = proposal_observer(reference.event_proposal, locate)
    trials = {phase: [] for phase in PHASES}
    counts = Counter()
    last, rejected = {}, {}

    def record(row):
        stage = row["stage"]
        counts[stage] += 1
        if counts[stage] > 2000 or row["progress_from"] != last.get(stage, 0):
            raise ValueError("stage budget or accepted history differs")
        if row["accepted"]:
            last[stage] = row["progress_to"]
        else:
            rejected[stage] = [row["progress_from"], row["progress_to"]]
        if stage in trials:
            trials[stage].append(observation(row))

    try:
        reference._reference(mesh, ground, case, steps, 1, 7, False, True, lambda *_: None,
                             maximum_m=5e-5, path_error=control, attempt_sink=record)
        stop = "completed"
    except AssertionError as exc:
        if str(exc) != "independent stage attempt budget":
            raise
        stop = str(exc)
    return dict(steps=steps, route=route, stop=stop, counts=dict(counts), last_accepted=last,
                last_rejected=rejected, trials=trials, passed=False)


def summarize(rows):
    result = {}
    for label in ("crossing", "smooth", "unknown"):
        group = [r for r in rows if r["interval"] == label]
        accepted = [r for r in group if r["accepted"]]
        lengths = [r["end"] - r["start"] for r in accepted]
        summary = dict(accepted=len(accepted), rejected=len(group)-len(accepted),
                       accepted_length_sum=sum(lengths), accepted_length_min=min(lengths, default=None),
                       accepted_length_max=max(lengths, default=None))
        for quantity in ("slip", "dissipation"):
            worst = max((r for r in group if quantity in r),
                        key=lambda r: r[quantity]["value"], default=None)
            summary[quantity] = (dict(**worst[quantity], start=worst["start"], end=worst["end"],
                                      accepted=worst["accepted"]) if worst else None)
        result[label] = summary
    return result


def compare(runs):
    by_route = {r["route"]: r for r in runs}
    if set(by_route) != set(ROUTES) or len(runs) != 3 or len({r["steps"] for r in runs}) != 1:
        raise ValueError("three one-axis routes at one N required")
    result = {}
    for axis, other in (("gap_locator", "bisect4"), ("local_substeps", "event1")):
        pair = [by_route[other], by_route["event4"]]
        for stage in PHASES:
            traces = [run["trials"][stage] for run in pair]
            common = set.intersection(*({r["end"] for r in rows if r["accepted"]} for rows in traces))
            if not common:
                result[f"{axis}/{stage}"] = dict(common_end=None, routes=None)
                continue
            end = max(common)
            rows_by_route = {}
            for route, rows in zip((other, "event4"), traces, strict=True):
                # Include overshooting rejected trials actually spent before landing on the anchor.
                stop = next((i+1 for i,r in enumerate(rows) if r["accepted"] and r["end"] == end), 0)
                rows_by_route[route] = summarize(rows[:stop])
            result[f"{axis}/{stage}"] = dict(common_end=end, routes=rows_by_route)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, required=True, choices=(512, 1024))
    parser.add_argument("--route", required=True, choices=ROUTES)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    value = run(args.steps, args.route)
    with gzip.open(args.output, "xt", encoding="utf-8") as stream:
        json.dump(value, stream, allow_nan=False, separators=(",", ":"))
