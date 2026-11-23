"""Integrate signed MC plastic work on the accepted material subincrements of each route."""

import argparse
import heapq
import json
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-23-ac-scope-decisions"))
from classifications import Model, STAGES, coordinate, material_state, mc, read, stream

sys.path.insert(0, str(HERE.parent / "2026-11-23-c03-mcc-plastic-work"))
from independent_work import components


def totals(stresses, increments, weights):
    terms = components(stresses, increments)
    return [math.fsum(column) for column in (terms * np.array(weights)[:, None]).T]


def independent(folder, model, resolve_edges):
    assert read(folder / "result.json")["complete"]
    states = model.seeds.copy()
    previous_u = np.zeros(model.nd)
    anchors = stream(folder / "fields.jsonl")
    extras = stream(folder / "substeps.jsonl")
    for sample in heapq.merge(anchors, extras, key=lambda x: coordinate(x["row"])):
        row, fields = sample["row"], sample["fields"]
        u = np.ravel(fields["u"])
        strains = -model.strains(u - previous_u)
        previous_u = u
        stresses, increments, weights = [], [], []
        for i, (old, mat, de, count) in enumerate(
            zip(states, model.mats, strains, row["material_partitions"], strict=True)
        ):
            state = old
            for _ in range(count):
                previous = state["plastic_strain_compression_positive"]
                state = mc.step(mat, state, de / count, resolve_edges=resolve_edges)
                stresses.append(state["effective_stress_compression_positive_kpa"])
                increments.append(np.subtract(state["plastic_strain_compression_positive"], previous))
                weights.append(model.weights[i])
            if state != material_state(fields, i):
                raise ValueError(("independent substep replay", row["stage"], row["progress"], i))
            states[i] = state
        yield dict(
            stage=row["stage"],
            progress=row["progress"],
            substeps=len(weights),
            work_kn_m=totals(stresses, increments, weights),
        )


def native(matrix, route, model, probe):
    folder = matrix / (route + "-trace")
    raw = read(matrix / (route + "-native/raw/run.json"))
    assert raw["complete"]
    weights = [p["weight_m3"] for p in raw["integration_points"]]
    states = model.seeds.copy()
    with TemporaryDirectory() as tmp:
        request, response = (Path(tmp) / n for n in ("request.json", "response.json"))
        for frame, step in zip(
            stream(folder / "accepted.jsonl"), stream(folder / "constitutive.jsonl"), strict=True
        ):
            stage, progress = step["stage"], step["progress"]
            assert frame["stage_index"] == stage and (
                frame["stage_progress"] == progress
                or (stage, progress, frame["stage_progress"]) == (0, 1, 0)
            )
            requests = [
                dict(
                    route=route,
                    increment=stage,
                    gp=i,
                    parameters=mat,
                    old=old,
                    strain=de,
                    count=step["partitions"][i],
                    save_substeps=True,
                )
                for i, (mat, old, de) in enumerate(
                    zip(model.mats, states, step["strain_increment"], strict=True)
                )
            ]
            request.write_text(json.dumps(requests), encoding="utf-8")
            subprocess.run([probe, request, response], check=True)
            points = read(response)
            response.unlink()
            stresses, increments, volumes = [], [], []
            states = frame["material_states"]
            for i, (record, accepted) in enumerate(zip(points, states, strict=True)):
                if record["gp"] != i or record["state"] != accepted:
                    raise ValueError(("native substep replay", stage, progress, i))
                for sub in record["substeps"]:
                    stresses.append(sub["stress_kpa"])
                    increments.append(sub["plastic_increment"])
                    volumes.append(weights[i])
            yield dict(
                stage=STAGES[stage],
                progress=progress,
                substeps=len(volumes),
                work_kn_m=totals(stresses, increments, volumes),
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("matrix", type=Path)
    p.add_argument("route", choices=("m0c1", "m1c0", "m1c1"))
    p.add_argument("output", type=Path)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--native-probe", type=Path)
    mode.add_argument("--independent", choices=("original", "resolved"))
    a = p.parse_args()
    folder = a.matrix / (a.route + ("-native" if a.native_probe else "-independent"))
    model = Model(
        read(a.matrix / "mesh.json"), read(a.matrix / "ground.json"), read(folder / "stages.json")
    )
    rows = (
        native(a.matrix, a.route, model, a.native_probe.resolve())
        if a.native_probe
        else independent(folder, model, a.independent == "resolved")
    )
    with a.output.open("x", encoding="utf-8") as f:
        for record in rows:
            print(json.dumps(record, allow_nan=False), file=f, flush=True)
