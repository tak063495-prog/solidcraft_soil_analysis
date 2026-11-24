"""Extract native FEM inputs and work from accepted steps; verify native replay exactly."""

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from c03_point import POINTS, tensor, work


def read(path):
    return json.loads(path.read_bytes())


def extract(root, route, probe, tmp):
    folder = root / route
    mesh, ground, run = [
        read(p) for p in (root / "mesh.json", root / "ground.json", folder / "run.json")
    ]
    if not run["complete"] or len(run["integration_points"]) != 324:
        raise ValueError("incomplete 12 Tet10 run")
    attempts = run["attempts"] + run["geostatic"]["attempts"]
    if len(attempts) != 224 or not all(a["accepted"] for a in attempts):
        raise ValueError("expected 32 accepted increments per stage without retries")
    regions = {r["id"]: r["body"] for r in mesh["regions"]}
    materials = {m["id"]: m["mohr_coulomb"] for m in ground["materials"]}
    seeds = ground["initialization"]["geostatic"]["seed"]["states"]
    elements = [run["integration_points"][gp]["element_index"] for gp in POINTS]
    models = [materials[mesh["material_ids"][e]] for e in elements]
    if any(m != models[0] for m in models):
        raise ValueError("expected matching material properties at selected points")
    old = [seeds[regions[mesh["region_ids"][e]]] for e in elements]
    initial = [s["effective_stress_compression_positive_kpa"] for s in old]
    rows, requests, accepted = [], [], []
    with (
        (folder / "accepted.jsonl").open() as frames,
        (folder / "constitutive.jsonl").open() as strains,
    ):
        for index, (f, d) in enumerate(zip(frames, strains, strict=True)):
            state, step = json.loads(f), json.loads(d)
            stage, progress = step["stage"], step["progress"]
            if (
                (stage, progress) != div_coordinate(index)
                or state["stage_index"] != stage
                or (
                    state["stage_progress"] != progress
                    and (stage, progress, state["stage_progress"]) != (0, 1, 0)
                )
            ):
                raise ValueError("accepted coordinate mismatch")
            strain = [step["strain_increment"][gp] for gp in POINTS]
            for j, gp in enumerate(POINTS):
                requests.append(
                    dict(
                        route=route,
                        increment=index,
                        gp=gp,
                        parameters=models[j],
                        old=old[j],
                        strain=strain[j],
                        count=step["partitions"][gp],
                        save_substeps=True,
                    )
                )
                accepted.append(state["material_states"][gp])
                old[j] = accepted[-1]
            rows.append(
                dict(
                    stage=stage,
                    progress=progress,
                    strain=strain,
                    fem_label=[],
                    fem_substeps=[],
                    fem_work_increment_kpa=[],
                )
            )
    if len(rows) != 224:
        raise ValueError("incomplete native trace")
    request, response = tmp / "request.json", tmp / "response.json"
    request.write_text(json.dumps(requests), encoding="utf8")
    subprocess.run([probe, request, response], check=True)
    records = read(response)
    response.unlink()
    max_work_roundoff = 0.0
    for req, rec, end in zip(requests, records, accepted, strict=True):
        if (
            any(req[k] != rec[k] for k in ("increment", "gp", "route", "count"))
            or rec["state"] != end
        ):
            raise ValueError("native replay differs from accepted FEM state")
        subs = rec["substeps"]
        if len(subs) != req["count"] or subs[-1]["label"] != end["last_return"]:
            raise ValueError("native partition mismatch")
        w = work(
            tensor([s["stress_kpa"] for s in subs]),
            tensor([s["plastic_increment"] for s in subs], engineering=True),
        )
        # Rust integrates its own principal-stress work at each native substep.
        # Cross-check the Voigt/tensor reconstruction against that state variable.
        native_total = end["plastic_work_kpa"]
        replay_total = req["old"]["plastic_work_kpa"] + math.fsum(w[:, 0])
        roundoff = abs(native_total - replay_total)
        max_work_roundoff = max(max_work_roundoff, roundoff)
        if roundoff > 64 * len(subs) * math.ulp(max(1.0, abs(native_total))):
            raise ValueError("reconstructed work differs from native plastic work")
        row = rows[rec["increment"]]
        row["fem_substeps"].append(len(subs))
        row["fem_work_increment_kpa"].append(w.sum(axis=0).tolist())
        row["fem_label"].append(end["last_return"])
    return dict(
        max_native_work_roundoff_kpa=max_work_roundoff,
        material=models[0],
        points=POINTS,
        initial_stress_kpa=initial,
        frames=rows,
        trace_sha256={
            p: hash_file(folder / p) for p in ("accepted.jsonl", "constitutive.jsonl")
        },
    )


def hash_file(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def div_coordinate(index):
    stage, i = divmod(index, 32)
    return stage, (i + 1) / 32


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("root", "probe", "output"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    with TemporaryDirectory() as directory:
        data = {
            r: extract(a.root, r, a.probe.resolve(), Path(directory))
            for r in ("m0c1", "m1c0", "m1c1")
        }
    with gzip.GzipFile(str(a.output), "wb", mtime=0) as f:
        f.write(json.dumps(data, allow_nan=False, separators=(",", ":")).encode())
