"""Seal the full factorial C-01 matrix, including every terminated path."""

import argparse
import csv
import json
import math
from pathlib import Path

from run_case import STAGES, context, seal, sha

MODES = ("base", "w", "s", "a", "ws", "wa", "sa", "wsa")
SIDES = ("native", "independent")


def write(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def inspect(folder, side, n, mode, expected):
    flags = set() if mode == "base" else set(mode)
    case = json.loads((folder / "stages.json").read_bytes())
    if "s" in flags:
        if (
            case["joint_system"]["global_step_control"]["path_error"].pop("semismooth_contact", False)
            is not True
        ):
            raise ValueError("active set option missing")
    if case != expected:
        raise ValueError(f"physical input or error control differs in {folder.name}")
    result = json.loads((folder / "result.json").read_bytes())
    if (result.get("seed_queries", 0) > 0) != ("w" in flags) or (result.get("arc_reviews", 0) > 0) != (
        "a" in flags
    ):
        raise ValueError(f"applied method differs in {folder.name}")
    rows = result["stages"]
    if tuple(r["stage"] for r in rows) != STAGES:
        raise ValueError("stage inventory differs")
    for row in rows:
        attempts = row["attempts"]
        if (
            not 0 <= row["accepted"] <= attempts <= 2000
            or row["accepted"] + sum(row["rejections"].values()) != attempts
        ):
            raise ValueError("invalid stage budget")
        if not math.isfinite(row["max_force"]) or row["max_force"] > (
            case["settings"]["relative_force_tolerance"] if side == "native" else 1e-10
        ):
            raise ValueError("original equilibrium tolerance exceeded")
    complete = all(row["last_accepted"] == 1 for row in rows)
    if complete != result["complete"]:
        raise ValueError("completion differs from observed endpoints")
    stopped = next((r for r in rows if r["last_accepted"] != 1), None)
    if stopped and (not result["stop"] or not stopped["attempts"]):
        raise ValueError("missing stop provenance")
    per_stage = {r["stage"]: r for r in rows}
    return dict(
        solver=side,
        N=n,
        mode=mode,
        complete=complete,
        total_attempts=sum(r["attempts"] for r in rows),
        reload_attempts=per_stage["reload"]["attempts"],
        reload_accepted=per_stage["reload"]["last_accepted"],
        stop_stage=stopped["stage"] if stopped else None,
        stop_progress=stopped["last_accepted"] if stopped else None,
        result_sha256=sha(folder / "result.json"),
        case_sha256=sha(folder / "stages.json"),
        details=result,
    )


def generate(root, out):
    out.mkdir()
    matrix = []
    seeds = {}
    for side in SIDES:
        for n in (256, 512):
            folder = root / f"{side}-{n}-base"
            ident = seal(folder, side, n)
            if ident != json.loads((folder / "identity.json").read_bytes()):
                raise ValueError("unbound seed")
            result = json.loads((folder / "result.json").read_bytes())
            seeds[folder.name] = dict(**ident, attempts=sum(r["attempts"] for r in result["stages"]))
        for n in (512, 1024):
            expected = context(n)[3]
            for mode in MODES:
                matrix.append(inspect(root / f"{side}-{n}-{mode}", side, n, mode, expected))
    with (out / "matrix.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "solver",
                "N",
                "mode",
                "complete",
                "total_attempts",
                "reload_attempts",
                "reload_accepted",
                "stop_stage",
                "stop_progress",
                *STAGES,
            ]
        )
        for r in matrix:
            writer.writerow(
                [
                    r[k]
                    for k in (
                        "solver",
                        "N",
                        "mode",
                        "complete",
                        "total_attempts",
                        "reload_attempts",
                        "reload_accepted",
                        "stop_stage",
                        "stop_progress",
                    )
                ]
                + [s["attempts"] for s in r["details"]["stages"]]
            )
    write(out / "routes.json", {"coarse_seeds": seeds, "routes": matrix})
    qualifying = [
        mode
        for mode in MODES
        if all(row["complete"] for row in matrix if row["N"] == 1024 and row["mode"] == mode)
    ]
    write(
        out / "validation.json",
        dict(
            measured_cases=len(matrix),
            full_factorial_measured=True,
            both_solvers_complete_at_1024=bool(qualifying),
            jointly_completed_modes=qualifying,
            stage_attempt_budget=2000,
            mesh_elements=12,
            scopes=[
                dict(
                    id="paired_method_trials",
                    status="達成",
                    scope="同一の12 Tet10で、対照と全7選択をN512/N1024・両実装で計測。",
                ),
                dict(
                    id="spatial_refinement",
                    status="未達",
                    scope="Nは荷重進捗の分割数。異なる空間メッシュ間の補間と要素数比例の予算は未検証。",
                ),
                dict(
                    id="full_history_convergence",
                    status="未達",
                    scope="同じ受理時刻における全場の刻み収束と元のC-01全体。",
                ),
            ],
            c01_complete=False,
            p6_complete=False,
            counts_sha256=sha(out / "matrix.csv"),
            routes_sha256=sha(out / "routes.json"),
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    generate(args.root, args.output)
