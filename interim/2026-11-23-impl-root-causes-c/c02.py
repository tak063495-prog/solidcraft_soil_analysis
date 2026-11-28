"""Compare time and mesh independently on the same physical points and accepted times."""
import argparse
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-22-impl-t8-acceptance-criteria"))
import two_mesh

diag = two_mesh.diag
STAGES = ("geostatic", "shear", "open", "recontact", "reverse", "reload", "return")


def common_indices(rows, steps):
    times = [(row["stage"], row["progress"]) for row in rows]
    if len(times) != len(set(times)):
        raise ValueError("duplicate accepted time")
    index = {key: i for i, key in enumerate(times)}
    keys = [(stage, i / steps) for stage in STAGES for i in range(1, steps + 1)]
    if any(key not in index for key in keys):
        raise ValueError("missing exact accepted time")
    indices = [index[key] for key in keys]
    if indices != sorted(indices):
        raise ValueError("accepted history out of order")
    return keys, indices


def pair_metrics(left, right, weights, keys, points, roots):
    measures = two_mesh.metrics(left, right, weights, keys, points)
    _, gp = np.unravel_index(np.abs(left-right).argmax(), left.shape)[:2]
    measures["maximum_at"].update(parent_face=int(roots[gp]), point_index=int(gp))
    area = diag.area_shape(weights, left)
    surface_left = (left*area).sum(axis=1)
    surface_right = (right*area).sum(axis=1)
    # Retain signed components at the time of the largest *integrated* difference.
    delta = surface_left - surface_right
    i = int(np.unravel_index(np.abs(delta).argmax(), delta.shape)[0])
    measures["signed_integral_at_maximum"] = dict(time=keys[i], left=surface_left[i].tolist(),
                                                  right=surface_right[i].tolist(), delta=delta[i].tolist())
    return measures


def measure(source):
    legacy = diag.legacy()
    ground, stages = [diag.read(source / f"{name}.json") for name in ("ground", "stages")]
    meshes = [diag.read(source / f"mesh-h{level}.json") for level in (0, 1)]
    systems = [legacy.old.Water(mesh, ground, stages).contact for mesh in meshes]
    points = np.array([p["position_m"] for system in systems for p in system.points])
    roots = np.array([p["face_index"] // 4**h for h, system in enumerate(systems) for p in system.points])
    weights = np.array([p["weight_m2"]/2 for system in systems for p in system.points])
    plans = [two_mesh.ContactPlan(mesh, systems[h], h, points, roots, "cubic") for h, mesh in enumerate(meshes)]
    if len(points) != 160 or any(not np.array_equal(plan.owners // 4**h, roots)
                                 for h, plan in enumerate(plans)):
        raise ValueError("common point ancestry changed")
    pairs = [("time", (0,64), (0,128)), ("time", (1,64), (1,128)),
             ("mesh", (0,64), (1,64)), ("mesh", (0,128), (1,128))]
    results = []
    for solver in ("native", "independent"):
        mapped = {}
        for h, plan in enumerate(plans):
            for n in (64, 128):
                rows, fields = diag.load_contact(legacy, source / f"h{h}-g{n}" / solver)
                keys, ids = common_indices(rows, 64)
                mapped[h,n] = {name: plan.history(fields[name][ids]) for name in diag.FIELDS}
        for axis, left, right in pairs:
            results.append(dict(axis=axis, solver=solver, left=left, right=right, accepted=len(keys),
                                fields={name: pair_metrics(mapped[left][name], mapped[right][name],
                                                            weights, keys, points, roots) for name in diag.FIELDS}))
    return dict(method="cubic", steps_of_common_times=64, points=len(points), area_m2=float(sum(weights)),
                parent_faces=np.unique(roots).tolist(), comparisons=results,
                continuum_reference=None, spatial_order=None, gci=None, passed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    diag.write(args.output, measure(args.source))
