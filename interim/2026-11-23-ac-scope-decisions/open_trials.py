"""Re-evaluate saved trial coordinates without advancing the accepted history."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "2026-11-23-impl-root-causes-a-nc"
sys.path.insert(0, str(PREVIOUS))
import branch as b


def inspect(op, parent, z, expected, model, engine):
    ev = op(parent, z)
    if engine == "native":
        old = [parent["joints"]["points"][i] for i in op.order]
        before, q, rest = (
            np.array(v)
            for v in (
                [p["jump_m"] for p in old],
                [p["traction_kpa"] for p in old],
                [p["birth"]["elastic_offset_m"] for p in old],
            )
        )
        new = [ev.candidate["joints"]["points"][i] for i in op.order]
        target = np.array([p["jump_m"] for p in new])
        stored = [p["status"] for p in new]
    else:
        old = parent["joint"]
        before, q, rest = (old[k] for k in ("jump", "traction", "rest"))
        target = ev.candidate["joint"]["jump"]
        stored = [("open", "stick", "slip")[k] for k in ev.candidate["joint"]["status"]]
    s = model.system
    # This calculation uses the single-subincrement branch of the existing
    # contact law. Refuse a replay whose last trial would require another one.
    assert all(len(v["path"]) == 1 for v in ev.details["branches"])
    n = s.normal
    normal = lambda a: np.einsum("ij,ij->i", a, n)
    g0, g = normal(before - rest), normal(target - rest)
    dj = target - before
    dt = dj - normal(dj)[:, None] * n
    closing = (g0 > 0) & (g <= 0)
    fraction = np.ones_like(g)
    fraction[closing] = g[closing] / (g[closing] - g0[closing])
    shear = q - normal(q)[:, None] * n
    shear[g0 > 0] = 0.0
    trial = shear + s.kt[:, None] * dt * fraction[:, None]
    excess = np.linalg.norm(trial, axis=1) - s.mu * s.kn * np.maximum(-g, 0)
    computed = np.where(g > 0, "open", np.where(excess > 0, "slip", "stick"))
    if list(computed) != list(ev.event):
        raise ValueError("replayed physical labels differ from the continuation arrival check")
    changes = []
    for i, (previous, new) in enumerate(zip(expected, ev.event, strict=True)):
        if previous != new:
            p = s.points[i]
            changes.append(
                dict(
                    gp=i,
                    joint=p["joint_id"],
                    face=p["face_index"],
                    owners=[p["owner_a"], p["owner_b"]],
                    position_m=p["position_m"],
                    previous=previous,
                    candidate=new,
                    gap_m=g[i],
                    friction_excess_kpa=excess[i],
                    stored_label=stored[i],
                )
            )
    return dict(
        progress=z[-1],
        relative_force=float(np.linalg.norm(ev.residual) / ev.denominator),
        events=list(ev.event),
        stored_label_mismatches=[i for i, v in enumerate(stored) if v != ev.event[i]],
        changed=changes,
        min_abs_gap_m=float(np.min(abs(g))),
        min_abs_closed_friction_excess_kpa=float(np.min(abs(excess[g <= 0]))),
    )


def measure(folder, bridge, engine):
    model = b.Model(*(b.read(folder / (k + ".json")) for k in ("mesh", "ground", "case")))
    origin, history = [b.unpack(folder / (k + ".gz")) for k in ("origins", engine)]
    path = origin["path"]
    path["holding"] = np.array(path["holding"])
    run = b.read(folder / "native.json")
    run["last_committed"] = origin["native"]
    op = (
        b.Native(bridge, folder, model, run, 1) if engine == "native" else b.independent(model, path, 1)
    )
    tol = model.case["settings"]["relative_force_tolerance"] if engine == "native" else 1e-10
    rows = []
    try:
        for t in history["trials"]:
            if t["lambda_parent"] < 0.02:
                continue
            (idx,) = [i for i, p in enumerate(history["points"]) if p["z"][-1] == t["lambda_parent"]]
            point = history["points"][idx]
            parent = b.checkpoint(point["state"], engine)
            assert b.identity(parent) == t["parent"]
            z = np.array(point["z"])
            secant = z - np.array(history["points"][idx - 1]["z"])
            probes = {
                "tangent": z + t["ds"] * np.array(t["direction"]),
                "secant": z + t["ds"] * secant / np.linalg.norm(secant),
                "corrected": np.array(t["z"]),
            }
            sampled = {
                name: inspect(op, parent, coords, point["event"], model, engine)
                for name, coords in probes.items()
            }
            assert (not sampled["corrected"]["changed"]) == t["accepted"]
            assert sampled["corrected"]["relative_force"] <= tol and b.identity(parent) == t["parent"]
            rows.append(dict(parent_progress=z[-1], ds=t["ds"], accepted=t["accepted"], **sampled))
    finally:
        if engine == "native":
            op.close()
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bridge", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    archive = PREVIOUS / "observations.zip"
    assert (
        hashlib.sha256(archive.read_bytes()).hexdigest()
        == b.read(PREVIOUS / "validation.json")["archive_sha256"]
    )
    data = {}
    with TemporaryDirectory() as tmp, zipfile.ZipFile(archive) as z:
        for size in ("18", "144"):
            folder = Path(tmp) / size
            folder.mkdir()
            for name in (
                "mesh.json",
                "ground.json",
                "case.json",
                "native.json",
                "origins.gz",
                "native.gz",
                "independent.gz",
            ):
                (folder / name).write_bytes(z.read(f"a01/{size}/{name}"))
            data[size] = {
                e: measure(folder, args.bridge.resolve(), e) for e in ("native", "independent")
            }
    with args.output.open("x") as out:
        json.dump(data, out, default=b.serial, allow_nan=False, indent=2)
        out.write("\n")
