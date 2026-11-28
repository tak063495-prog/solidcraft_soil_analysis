"""Measure the two permitted meshes at their union of physical contact points."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-12-17-c02-spatial-cause-diagnosis"))
import c02diag_driver as diag
from c02diag_mapping import ContactPlan
from c02diag_metrics import difference, time_identity


def metrics(coarse, fine, weights, keys, points):
    delta = np.abs(coarse - fine)
    index = np.unravel_index(delta.argmax(), delta.shape)
    measure = diag.area_shape(weights, coarse)
    return {
        "point": difference(coarse, fine),
        "area_weighted": difference(coarse, fine, measure),
        "surface_integral": difference(
            (coarse * measure).sum(axis=1), (fine * measure).sum(axis=1)
        ),
        "maximum_at": {
            "accepted": keys[index[0]],
            "position_m": points[index[1]].tolist(),
            "component": int(index[2]) if len(index) > 2 else None,
        },
    }


def measure(source):
    legacy = diag.legacy()
    meshes = [diag.read(source / f"mesh-h{h}.json") for h in (0, 1)]
    ground, stages = [diag.read(source / f"{name}.json") for name in ("ground", "stages")]
    systems = [legacy.old.Water(mesh, ground, stages).contact for mesh in meshes]
    points = np.array([p["position_m"] for system in systems for p in system.points])
    roots = np.array([p["face_index"] // 4**h for h, system in enumerate(systems)
                      for p in system.points])
    weights = np.array([p["weight_m2"] / 2 for system in systems for p in system.points])
    plans = {method: [ContactPlan(mesh, systems[h], h, points, roots, method)
                      for h, mesh in enumerate(meshes)] for method in diag.METHODS}
    results, time_exchange, constraints = {}, {}, {}
    normal = np.array(systems[0].normal[0])
    mu = systems[0].mu[0]
    if any(not (np.allclose(system.normal, normal, rtol=0, atol=1e-15) and
                np.all(system.mu == mu)) for system in systems):
        raise ValueError("two-mesh law observation requires a common planar normal and friction")
    for solver in ("native", "independent"):
        previous = {}
        for steps in (32, 64, 128):
            histories = [diag.load_contact(legacy, source / f"h{h}-g{steps}" / solver)
                         for h in (0, 1)]
            (keys, coarse), (fine_keys, fine) = [
                diag.sampled(rows, fields, steps) for rows, fields in histories
            ]
            if keys != fine_keys:
                raise ValueError("mesh histories do not share the accepted coordinates")
            for method, (left, right) in plans.items():
                label = f"{solver}/{steps}/{method}"
                fields, mapped = {}, {}
                gap_labels = [plan.history(src["gap"]) > 0 for plan, src in
                              zip(plans["nearest"], (coarse, fine), strict=True)]
                for field in diag.FIELDS:
                    a, b = left.history(coarse[field]), right.history(fine[field])
                    mapped[field] = (a, b)
                    fields[field] = metrics(a, b, weights, keys, points)
                    if field == "plastic_dissipation":
                        fields[field]["negative_reconstructions"] = [
                            int(np.count_nonzero(values < 0)) for values in (a, b)
                        ]
                results[label] = fields
                constraints[label] = []
                for index in range(2):
                    traction = mapped["traction"][index]
                    pressure = -np.einsum("...i,i->...", traction, normal)
                    shear = np.linalg.norm(traction + pressure[..., None]*normal, axis=-1)
                    gap = mapped["gap"][index]
                    constraints[label].append(dict(
                        friction_excess_max_kpa=float(np.max(shear-mu*np.maximum(pressure, 0), initial=0)),
                        traction_in_reconstructed_open_max_kpa=float(np.max(
                            np.linalg.norm(traction[gap>0], axis=-1), initial=0))))
                if method in previous:
                    prior_keys, prior_fields, prior_labels = previous[method]
                    lookup = {key:i for i,key in enumerate(keys)}
                    ids = [lookup[(stage, 2*i)] for stage,i in prior_keys]
                    time_exchange[label] = {
                        k:time_identity(*prior_fields[k], *(v[ids] for v in mapped[k]),
                                        [*prior_labels, *(v[ids] for v in gap_labels)])
                        for k in diag.FIELDS
                    }
                previous[method] = keys, mapped, gap_labels
    return dict(mesh_tet10=[12, 96], steps=[32, 64, 128], common_points=len(points),
                area_m2=float(weights.sum()), comparisons=results, time_mesh_exchange=time_exchange,
                reconstruction_law=constraints, observed_order=None,
                continuum_error=None, passed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    value = measure(args.source)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
