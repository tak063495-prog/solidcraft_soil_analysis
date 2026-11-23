"""Compare the sealed producer and independent contact histories on common GP coordinates."""

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "2026-12-19-c04-transient-water-cause-diagnosis"
sys.path.insert(0, str(HERE.parent / "2026-11-22-impl-t8-acceptance-criteria"))
sys.path.insert(0, str(SOURCE))
from sealed import archive
from c04cause_metrics import metric
from c04cause_observe import last_update, reference


def accepted(records, stages):
    expected = [(name, i / 32) for name in stages for i in range(1, 33)]
    actual = [(r["row"]["stage"], r["row"]["progress"]) for r in records]
    if actual != expected:
        raise ValueError("incomplete or misaligned accepted stage/progress coordinates")


def increments(values):
    return np.diff(values, axis=0, prepend=np.zeros_like(values[:1]))


def observables(records, system, partitions):
    def field(key):
        values = np.asarray([r["fields"][key] for r in records])
        if not np.isfinite(values).all() or values.shape[:2] != partitions.shape:
            raise ValueError(f"incomplete or nonfinite GP field {key}")
        return values

    jump, traction = field("jump"), field("traction")
    signed_normal = np.einsum("tgi,gi->tg", traction, system.normal)
    normal_traction = signed_normal[:, :, None] * system.normal
    tangential = traction - normal_traction
    trial = np.empty_like(signed_normal)
    # Replay the last constitutive substep for each system using its saved partition count.
    if partitions.min() < 1:
        raise ValueError("invalid contact partition count")
    old_jump = np.concatenate((np.zeros_like(jump[:1]), jump[:-1]))
    old_traction = np.concatenate((np.zeros_like(traction[:1]), traction[:-1]))
    for count in np.unique(partitions):
        times, gps = np.nonzero(partitions == count)
        part = SimpleNamespace(
            **{k: getattr(system, k)[gps] for k in ("normal", "kn", "kt", "mu")}
        )
        trial[times, gps] = last_update(
            part,
            old_jump[times, gps],
            old_traction[times, gps],
            jump[times, gps],
            int(count),
        )["trial_margin_kpa"]
    work = field("plastic_dissipation")
    area_work = work @ system.weight
    return dict(
        gap_m=field("gap"),
        effective_pressure_kpa=field("effective_contact_pressure"),
        normal_a_kn=-normal_traction * system.weight[:, None],
        normal_b_kn=normal_traction * system.weight[:, None],
        tangent_a_kn=-tangential * system.weight[:, None],
        tangent_b_kn=tangential * system.weight[:, None],
        trial_margin_kpa=trial,
        returned_margin_kpa=field("friction_limit")
        - np.linalg.norm(tangential, axis=2),
        plastic_slip_increment_m=increments(field("plastic_slip")),
        interval_kpa_m=increments(work),
        cumulative_kpa_m=work,
        interval_area_kn_m=increments(area_work),
        cumulative_area_kn_m=area_work,
    )


def wet_history():
    with archive(SOURCE, "archive-index.json") as z:

        def read(name):
            return json.loads(z.read("data/global-final/" + name))

        def lines(name):
            with z.open("data/global-final/" + name) as stream:
                return [json.loads(line) for line in stream]

        case = read("native-w32/case.json")
        if case != read("adaptive-w32/case.json"):
            raise ValueError("native and independent input differ")
        model = reference.Model(read("mesh.json"), read("ground.json"), case)
        system = model.system
        native, other = [
            lines(n + "/fields.jsonl") for n in ("native-w32", "adaptive-w32")
        ]
        counts = np.asarray(
            [
                r["partitions"][len(model.seeds) :]
                for r in lines("native-w32/raw/constitutive.jsonl")
            ]
        )
    stages = ["geostatic", *(s["name"] for s in case["stages"])]
    accepted(native, stages)
    accepted(other, stages)
    if any(
        (a["row"]["level_m"], a["row"]["factor"])
        != (b["row"]["level_m"], b["row"]["factor"])
        for a, b in zip(native, other, strict=True)
    ):
        raise ValueError("accepted water level or load factor differs")
    a = observables(native, system, counts)
    b = observables(
        other, system, np.asarray([r["row"]["contact_partitions"] for r in other])
    )
    pairs = np.array(
        [
            [r["fields"]["contact_status"] for r in records]
            for records in (native, other)
        ]
    )
    labels, count = np.unique(pairs.reshape(2, -1).T, axis=0, return_counts=True)
    return dict(
        stages=stages,
        accepted=len(native),
        contact_gp=len(system.weight),
        counts={f"{i}/{j}": int(n) for (i, j), n in zip(labels, count, strict=True)},
        history={k: metric(a[k], b[k]) for k in a},
        by_stage={
            s: metric(
                a["cumulative_kpa_m"][32 * i : 32 * (i + 1)],
                b["cumulative_kpa_m"][32 * i : 32 * (i + 1)],
            )
            for i, s in enumerate(stages)
        },
        relative_width=None,
        passed=None,
    )
