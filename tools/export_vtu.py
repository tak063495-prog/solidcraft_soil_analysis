import argparse
import json
from pathlib import Path
import meshio
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("mesh", type=Path)
parser.add_argument("solution", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
if args.output.exists():
    parser.error("output already exists")
mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
result = json.loads(args.solution.read_text(encoding="utf-8"))
u = np.asarray(result["displacement_m"])
s = np.asarray(result["stress_tension_positive_kpa"])
meshio.write_points_cells(
    args.output, mesh["nodes"], [("tetra", mesh["tetrahedra"])],
    point_data={"displacement_m": u, "displacement_magnitude_m": np.linalg.norm(u, axis=1)},
    cell_data={"material_id": [np.asarray(mesh["material_ids"])],
               "stress_tension_positive_kpa": [s], "stress_compression_positive_kpa": [-s]},
)
