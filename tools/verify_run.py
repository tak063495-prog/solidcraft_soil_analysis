"""Independent checks of a known-volume self-weight bootstrap run."""
import argparse
import json
from pathlib import Path
import meshio
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("run", type=Path)
parser.add_argument("--volume-m3", type=float, required=True)
parser.add_argument("--expected-bounds-m", type=float, nargs=6, required=True,
                    metavar=("XMIN", "YMIN", "ZMIN", "XMAX", "YMAX", "ZMAX"))
args = parser.parse_args()
mesh = json.loads((args.run / "mesh.json").read_text(encoding="utf-8"))
case = json.loads((args.run / "case.json").read_text(encoding="utf-8"))
solution = json.loads((args.run / "solution.json").read_text(encoding="utf-8"))
all_points = np.asarray(mesh["nodes"])
bounds = np.concatenate([all_points.min(axis=0), all_points.max(axis=0)])
assert np.allclose(bounds, args.expected_bounds_m, rtol=1e-12, atol=1e-10), bounds
points = np.asarray(mesh["nodes"])[np.asarray(mesh["tetrahedra"])]
jacobians = (points[:, 1:] - points[:, :1]).transpose(0, 2, 1)
volumes = np.linalg.det(jacobians) / 6
assert np.all(volumes > 0)
volume = float(volumes.sum())
assert np.isclose(volume, args.volume_m3, rtol=1e-12)
weight = args.volume_m3 * case["unit_weight_kn_m3"]
applied = np.asarray(solution["applied_force_kn"])
reaction = np.asarray(solution["support_reaction_kn"])
assert np.isclose(applied[2], -weight, rtol=1e-12)
assert np.linalg.norm(applied + reaction) / weight < 1e-8
assert solution["true_relative_residual"] <= case["relative_tolerance"]
vtk = meshio.read(args.run / "result.vtu")
assert len(vtk.points) == len(mesh["nodes"])
assert np.allclose(vtk.point_data["displacement_m"], solution["displacement_m"])
assert np.array_equal(vtk.cells_dict["tetra"], mesh["tetrahedra"])
assert np.all(np.isfinite(vtk.cell_data_dict["stress_tension_positive_kpa"]["tetra"]))
report = {"volume_m3": volume, "weight_kn": weight,
          "bounds_m": bounds.tolist(),
          "relative_force_balance_error": float(np.linalg.norm(applied + reaction) / weight),
          "true_relative_residual": solution["true_relative_residual"],
          "maximum_displacement_m": float(np.linalg.norm(vtk.point_data["displacement_m"], axis=1).max()),
          "nodes": len(mesh["nodes"]), "tetrahedra": len(mesh["tetrahedra"]),
          "vtu_roundtrip": "passed"}
(args.run / "verification.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
