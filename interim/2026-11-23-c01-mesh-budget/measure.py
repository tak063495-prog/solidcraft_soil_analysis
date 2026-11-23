"""Seal the actual accepted C-01 histories and compare completed load grids."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "2026-11-23-c01-method-improvements"
sys.path.insert(0, str(PREVIOUS))
from run_case import STAGES, context, seal

sys.path.insert(0, str(HERE.parent / "2026-12-19-c04-transient-water-cause-diagnosis"))
from c04cause_metrics import cosine, metric


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def write(path, obj):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(obj, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def on_grid(history, n):
    keys = [(stage, i / n) for stage in range(len(STAGES)) for i in range(1, n + 1)]
    if any(key not in history for key in keys):
        raise ValueError("missing exact accepted coordinate")
    return np.array([history[k] for k in keys])


def convergence(histories):
    pairs = {}
    for n in (256, 512):
        pairs[f"{n}-{2 * n}"] = metric(on_grid(histories[n], n), on_grid(histories[2 * n], n))
    a, b, c = (on_grid(histories[n], 256) for n in (256, 512, 1024))
    ab, bc = a - b, b - c
    coarse, fine = map(float, (np.linalg.norm(ab), np.linalg.norm(bc)))
    return dict(
        adjacent=pairs,
        common_times=len(a),
        coarse_norm=coarse,
        fine_norm=fine,
        cosine=cosine(ab, bc),
        observed_order=float(np.log2(coarse / fine)) if min(coarse, fine) > 0 else None,
    )


def inspect(folder, out, side, n, mode, expected, old):
    cfg = json.loads((folder / "stages.json").read_bytes())
    policy = cfg["settings"].pop("attempt_budget_policy", None)
    if policy != "mesh_proportional":
        raise ValueError("missing budget policy")
    active = cfg["joint_system"]["global_step_control"]["path_error"].pop("semismooth_contact", False)
    if active != (mode == "ws") or cfg != expected:
        raise ValueError("physical input or control differs")
    result = json.loads((folder / "result.json").read_bytes())
    if (result.get("seed_queries", 0) > 0) != active or result["arc_reviews"] != 0:
        raise ValueError("method selection differs")
    rows = result["stages"]
    if [r["stage"] for r in rows] != list(STAGES):
        raise ValueError("stage order differs")
    data = (folder / "dissipation.jsonl").read_bytes()
    hist = {}
    counts = [0] * len(rows)
    last = [None] * len(rows)
    for line in data.splitlines():
        stage, time, field = json.loads(line)
        if not (
            0 <= stage < len(rows)
            and (last[stage] or 0) < time <= 1
            and len(field) == 32
            and np.isfinite(field).all()
            and min(field) >= 0
        ):
            raise ValueError("invalid accepted field or coordinate")
        last[stage] = time
        counts[stage] += 1
        hist[stage, time] = field
    for i, r in enumerate(rows):
        if (
            not 0 <= r["attempts"] <= 2065
            or r["accepted"] != counts[i]
            or r["last_accepted"] != last[i]
            or sum(r["rejections"].values()) + counts[i] != r["attempts"]
        ):
            raise ValueError("trial counts or accepted endpoints differ from history")
        if (
            not 0
            <= r["max_force"]
            <= (cfg["settings"]["relative_force_tolerance"] if side == "native" else 1e-10)
        ):
            raise ValueError("force tolerance exceeded")
    complete = all(t == 1.0 for t in last)
    stopped = next((r for r in rows if r["last_accepted"] != 1.0), None)
    if complete != result["complete"] or bool(result["stop"]) == complete:
        raise ValueError("inconsistent termination")
    if stopped and (
        stopped["attempts"] != 2065 or any(r["attempts"] for r in rows[rows.index(stopped) + 1 :])
    ):
        raise ValueError("unexpected termination; inspect source run")
    archived = folder.name + ".jsonl.gz"
    with (out / archived).open("xb") as stream:
        with gzip.GzipFile(fileobj=stream, filename="", mode="wb", mtime=0) as z:
            z.write(data)
    previous = old.get((side, n, mode))
    unchanged = rows == previous["details"]["stages"] if previous else None
    return hist, dict(
        solver=side,
        N=n,
        mode=mode,
        complete=complete,
        stop=result["stop"],
        stopped_stage=stopped["stage"] if stopped else None,
        total=sum(r["attempts"] for r in rows),
        stages=rows,
        same_stage_records_as_fixed=unchanged,
        archive=archived,
        accepted_history_sha256=digest(data),
        case_sha256=digest((folder / "stages.json").read_bytes()),
        result_sha256=digest((folder / "result.json").read_bytes()),
        seed_queries=result["seed_queries"],
        arc_reviews=result["arc_reviews"],
        independent_active_tangents=result.get("active_tangents"),
    )


def measure(root, seeds, out):
    out.mkdir()
    old_list = json.loads((PREVIOUS / "results/routes.json").read_bytes())["routes"]
    old = {(r["solver"], r["N"], r["mode"]): r for r in old_list}
    routes = {}
    orders = {}
    for side in ("native", "independent"):
        for mode in ("base", "ws"):
            histories = {}
            for n in (512, 1024, 256):
                if n == 256 and not routes[f"{side}-1024-{mode}"]["complete"]:
                    continue
                name = f"{side}-{n}-{mode}"
                hist, route = inspect(root / name, out, side, n, mode, context(n)[3], old)
                if mode == "ws":
                    source = seeds / f"{side}-{n // 2}-base"
                    receipt = seal(source, side, n // 2)
                    if receipt != json.loads((source / "identity.json").read_bytes()):
                        raise ValueError("unbound coarse seed")
                    route["coarse_seed"] = receipt
                histories[n] = hist
                routes[name] = route
            orders[f"{side}-{mode}"] = convergence(histories) if 256 in histories else None
    with (out / "stages.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "solver",
                "N",
                "mode",
                "stage",
                "attempts",
                "accepted",
                "rejected",
                "last_accepted",
                "outcome",
                "stop",
            ]
        )
        for route in routes.values():
            for row in route["stages"]:
                stopped = row["stage"] == route["stopped_stage"]
                writer.writerow(
                    [
                        route["solver"],
                        route["N"],
                        route["mode"],
                        row["stage"],
                        row["attempts"],
                        row["accepted"],
                        row["attempts"] - row["accepted"],
                        row["last_accepted"],
                        "complete"
                        if row["last_accepted"] == 1
                        else "stopped"
                        if stopped
                        else "not_started",
                        route["stop"] if stopped else "",
                    ]
                )
    write(out / "routes.json", routes)
    write(out / "dissipation.json", orders)
    completed = [name for name, r in routes.items() if r["N"] == 1024 and r["complete"]]
    regression = all(r["same_stage_records_as_fixed"] for r in routes.values() if r["N"] == 512)
    scopes = [
        dict(
            id="n512_regression",
            status="達成" if regression else "未達",
            scope="12 Tet10・N512の両実装・両方式で従来予算の全工程記録と一致。",
        ),
        dict(id="measured", status="達成", scope="12 Tet10・湿潤7工程のevent4とWSを本体・独立で実測。"),
        dict(
            id="other_meshes", status="未達", scope="他の空間メッシュでの比例予算の十分性・計算費用。"
        ),
    ]
    for name, r in routes.items():
        if r["N"] == 1024:
            scopes.append(
                dict(
                    id=name,
                    status="達成" if r["complete"] else "未達",
                    scope="12 Tet10・N1024・この方式の7工程完走。",
                )
            )
    for key, result in orders.items():
        scopes.append(
            dict(
                id=key + "-dissipation-measured",
                status="達成" if result else "未達",
                scope="同一受理座標の全履歴散逸を3水準で比較し、有限刻みの自己差と次数を計測。",
            )
        )
        if result is not None:
            scopes.append(
                dict(
                    id=key + "-finite-dissipation-decrease",
                    status="達成" if result["fine_norm"] < result["coarse_norm"] else "未達",
                    scope="12 Tet10・この湿潤経路の3水準で、同一1792座標の散逸自己差が縮小。",
                )
            )
    scopes.append(
        dict(
            id="general_convergence",
            status="未達",
            scope="漸近域、連続解誤差、全場の刻み収束および元C-01全体。",
        )
    )
    write(
        out / "validation.json",
        dict(
            scopes=scopes,
            elements=12,
            stage_limit=2065,
            default_limit=2000,
            completed_n1024=completed,
            completed_n1024_count=len(completed),
            regression_n512_all_stage_records_match=regression,
            c01_complete=False,
            p6_complete=False,
            routes_sha256=digest((out / "routes.json").read_bytes()),
            stages_sha256=digest((out / "stages.csv").read_bytes()),
            dissipation_sha256=digest((out / "dissipation.json").read_bytes()),
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seeds", type=Path, required=True)
    args = parser.parse_args()
    measure(args.root, args.seeds, args.output)
