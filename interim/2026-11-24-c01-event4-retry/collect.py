"""Validate bounded C-01 route records and assemble full-field evidence."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil

import numpy as np

from c01_retry_metrics import History, align, compare, sha256
from c01_retry_fields import HERE, STAGES
from run_case import context, seal

PREVIOUS = HERE.parent / "2026-11-23-c01-mesh-budget" / "results"


def write(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def inspect(folder, side, n, mode, seeds, previous):
    cfg = json.loads((folder / "stages.json").read_bytes())
    study, mesh, _, expected, ctrl = context(n)
    expected["settings"]["attempt_budget_policy"] = "mesh_proportional"
    flags = {
        "ws": "semismooth_contact",
        "local": "contact_history_bisection",
        "retry": "error_scaled_retry",
    }
    if mode in flags:
        ctrl[flags[mode]] = True
    if mode in ("local-ws", "local-active", "friction-events","resolved-slip"):
        ctrl.update(semismooth_contact=True, contact_history_bisection=True)
    if mode == "resolved-slip":
        ctrl["resolve_slip_path"] = True
    if mode in ("friction-events","resolved-slip"):
        ctrl["localize_friction_events"] = True
    warm = mode in ("ws", "local-ws")
    if cfg != expected or study.reference.stage_attempt_budget(mesh, cfg["settings"]) != 2065:
        raise ValueError("case input or attempt budget changed")
    route = json.loads((folder / "result.json").read_bytes())
    if (route["seed_queries"] > 0) != warm or route["arc_reviews"] != 0:
        raise ValueError("unrequested method")
    if side == "independent" and (route["active_tangents"] > 0) != bool(ctrl.get("semismooth_contact")):
        raise ValueError("independent tangent does not match the requested method")
    if warm:
        seed_dir = seeds / f"{side}-{n // 2}-base"
        receipt = seal(seed_dir, side, n // 2)
        if receipt != json.loads((seed_dir / "identity.json").read_bytes()):
            raise ValueError("coarse seed identity differs")
        route["coarse_seed"] = receipt
    rows = route["stages"]
    if [s["stage"] for s in rows] != list(STAGES):
        raise ValueError("stage order differs")
    stopped = next((i for i, s in enumerate(rows) if s["last_accepted"] != 1), len(rows))
    if route["complete"] != (stopped == len(rows)) or bool(route["stop"]) == route["complete"]:
        raise ValueError("inconsistent completion")
    for i, s in enumerate(rows):
        if not (
            0 <= s["accepted"] <= s["attempts"] <= 2065
            and s["attempts"] == s["accepted"] + sum(s["rejections"].values())
            and 0
            <= s["max_force"]
            <= (cfg["settings"]["relative_force_tolerance"] if side == "native" else 1e-10)
        ):
            raise ValueError("stage count or force gate")
        if (
            i == stopped
            and (s["attempts"] != 2065 or "budget" not in route["stop"])
            or i > stopped
            and s["attempts"]
        ):
            raise ValueError("unexpected stop")
    history = folder / "dissipation.jsonl"
    history_hash = None
    if history.exists():
        counts = [0] * len(rows)
        last = [None] * len(rows)
        with history.open("rb") as stream:
            for line in stream:
                i, t, d = json.loads(line)
                if not (
                    0 <= i < len(rows)
                    and (last[i] or 0) < t <= 1
                    and len(d) == 32
                    and all(0 <= x < float("inf") for x in d)
                ):
                    raise ValueError("invalid accepted dissipation history")
                counts[i] += 1
                last[i] = t
        if counts != [s["accepted"] for s in rows] or last != [s["last_accepted"] for s in rows]:
            raise ValueError("accepted trace differs from stage records")
        history_hash = sha256(history)
    elif (mode, n) != ("retry", 1024):
        raise ValueError("missing accepted history")
    case_digest = sha256(folder / "stages.json")
    result_digest = sha256(folder / "result.json")
    replay = None
    if previous is not None:
        replay = dict(
            case_bytes=case_digest == previous["case_sha256"],
            result_bytes=result_digest == previous["result_sha256"],
            coarse_seed=route.get("coarse_seed") == previous.get("coarse_seed"),
            stages=rows == previous["stages"],
            termination=(route["complete"], route["stop"]) == (previous["complete"], previous["stop"]),
            accepted_dissipation=history_hash == previous["accepted_history_sha256"],
        )
    route.update(
        solver=side,
        N=n,
        mode=mode,
        completed_stages=stopped,
        attempts=sum(s["attempts"] for s in rows),
        accepted=sum(s["accepted"] for s in rows),
        case_sha256=case_digest,
        result_sha256=result_digest,
        accepted_history_sha256=history_hash,
        previous_replay=replay,
    )
    return route


def check_dissipation_anchors(series, history):
    grid = series.grid("plastic_dissipation", series.n)
    count = 0
    with history.open("rb") as src:
        for line in src:
            stage, progress, values = json.loads(line)
            if (float(progress) * series.n).is_integer():
                phase, i = divmod(count, series.n)
                if (stage, progress) != (phase, (i + 1) / series.n) or grid[count].tolist() != values:
                    raise ValueError("field stream differs from committed dissipation at anchor")
                count += 1
    if count != series.count:
        raise ValueError("committed dissipation anchors missing")
    return count


def measure(root, seeds, out):
    out.mkdir()
    old = json.loads((PREVIOUS / "routes.json").read_bytes())
    routes, histories, manifests = {}, {}, {}
    for mode in ("base", "ws", "local", "local-ws", "local-active", "retry"):
        for side in ("native", "independent"):
            for n in (
                (256, 512, 1024) if mode in ("ws", "local", "local-ws", "local-active") else (512, 1024)
            ):
                key = f"{side}-{n}-{mode}"
                folder = root / key
                r = inspect(folder, side, n, mode, seeds, old.get(key))
                routes[key] = r
                if mode in ("ws", "local", "local-ws", "local-active"):
                    series = History(root / (key + "-fields"), side, n, r["stages"])
                    r["confirmed_dissipation_anchors"] = check_dissipation_anchors(
                        series, folder / "dissipation.jsonl"
                    )
                    histories[side, n, mode] = series
                    manifests[key] = series.meta
                if mode in ("local", "local-ws", "local-active"):
                    with (
                        (out / (key + ".jsonl.gz")).open("xb") as blob,
                        gzip.GzipFile(fileobj=blob, filename="", mode="wb", mtime=0) as archive,
                        (folder / "dissipation.jsonl").open("rb") as src,
                    ):
                        shutil.copyfileobj(src, archive)
    write(out / "routes.json", routes)
    write(out / "streams.json", manifests)
    with (out / "stages.csv").open("x", encoding="utf-8", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(
            [
                "solver",
                "N",
                "method",
                "stage",
                "attempts",
                "accepted",
                "rejected",
                "last_accepted",
                "max_force",
                "result",
            ]
        )
        for r in routes.values():
            for s in r["stages"]:
                outcome = (
                    "complete"
                    if s["last_accepted"] == 1
                    else r["stop"]
                    if s["attempts"]
                    else "not_started"
                )
                w.writerow(
                    [
                        r["solver"],
                        r["N"],
                        r["mode"],
                        s["stage"],
                        s["attempts"],
                        s["accepted"],
                        s["attempts"] - s["accepted"],
                        s["last_accepted"],
                        s["max_force"],
                        outcome,
                    ]
                )
    result = compare(histories)
    for key, entry in result["refinements"].items():
        side, mode = key.split("-", 1)
        series = [histories[side, n, mode] for n in (256, 512, 1024)]
        a, b, c = align(series, "plastic_slip", 256)
        path = out / (key + "-slip-differences.f64.gz")
        with path.open("xb") as blob, gzip.GzipFile(fileobj=blob, filename="", mode="wb", mtime=0) as z:
            for diff in (a - b, b - c):
                z.write(np.ascontiguousarray(diff, dtype="<f8").tobytes())
        entry["slip_differences"] = dict(
            file=path.name,
            dtype="<f8",
            shape=[2, *a.shape],
            order=["256-512", "512-1024"],
            common_grid=256,
        )
    write(out / "fields.json", result)
    replay = all(
        all(r["previous_replay"].values()) for r in routes.values() if r["previous_replay"] is not None
    )
    scopes = [
        dict(
            id="control_replay",
            status="達成" if replay else "未達",
            scope="比較対象10経路の全工程記録・停止・全受理散逸の既報との一致。",
        )
    ]
    for key, r in routes.items():
        if r["N"] == 1024:
            scopes.append(
                dict(
                    id=key,
                    status="達成" if r["complete"] else "未達",
                    scope="12 Tet10 湿潤7工程の完走。",
                )
            )
    for key, r in result["refinements"].items():
        scopes.append(
            dict(
                id=key + "-measured",
                status="達成" if r["complete"] else "未達",
                scope="3水準共通の全7工程で連続20場と分類2場を照合。",
            )
        )
        continuous = [v["common"]["trend"] for v in r["fields"].values() if "trend" in v["common"]]
        categorical = [v["common"] for v in r["fields"].values() if "coarse_mismatches" in v["common"]]
        aligned = all(t in ("identical", "decreasing_aligned") for t in continuous)
        stable = all(v["coarse_mismatches"] == v["fine_mismatches"] == 0 for v in categorical)
        scopes.append(
            dict(
                id=key + "-finite-trends",
                status="達成" if r["complete"] and aligned and stable else "未達",
                scope="同じ有限3水準で全連続場が一致または正の次数・差方向、かつ両分類の不一致がゼロ。",
            )
        )
    scopes.extend(
        [
            dict(
                id="full_convergence",
                status="未達",
                scope="全場での刻み収束、漸近域および連続解との差。",
            ),
            dict(
                id="other_meshes", status="未達", scope="他の格子・荷重経路での比例予算と手法の妥当性。"
            ),
        ]
    )
    files = {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()}
    write(
        out / "validation.json",
        dict(
            scopes=scopes,
            default_limit=2000,
            stage_limit=2065,
            elements=12,
            c01_complete=False,
            regression=replay,
            sha256=files,
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("runs", type=Path)
    p.add_argument("seeds", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    measure(a.runs, a.seeds, a.output)
