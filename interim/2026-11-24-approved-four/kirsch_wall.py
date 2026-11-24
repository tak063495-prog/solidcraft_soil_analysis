"""Run the staged Kirsch case on radially graded meshes and recover axis wall stress."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "2026-11-24-t9d-benchmarks"))
from kirsch import A, R, study, tensor


def extrapolate(result, nt, nr):
    """Fit the GP hoop stresses in the first two radial and angular cells of each axis."""
    points = result["integration_points"]
    # Element order is angular sector, radial cell, then six Tet4s per prism pair.
    # All recovery meshes have one axial layer. Use connectivity order to include
    # every GP in the patch, including those inside the inscribed polygon's arc.
    cell = np.array([p["element_index"] // 6 for p in points])
    xy = np.array([p["position_m"][:2] for p in points])
    radius = np.linalg.norm(xy, axis=1)
    angle = np.arctan2(xy[:, 1], xy[:, 0])
    tangent = np.column_stack((-xy[:, 1], xy[:, 0])) / radius[:, None]
    stress = np.array([tensor(s) for s in result["stress_tension_positive_kpa"]])
    hoop = -np.einsum("ij,ijk,ik->i", tangent, stress, tangent)
    volume = np.array([p["weight_m3"] for p in points])
    dr = A * ((R / A) ** (1 / nr) - 1)
    dt = np.pi / (2 * nt)
    estimates = []
    for end, theta in ((0, 0), (nt - 1, np.pi / 2)):
        take = (np.abs(cell // nr - end) < 2) & (cell % nr < 2)
        r = (radius[take] - A) / dr
        # Reflection symmetry makes the hoop stress even about the axis.
        design = np.column_stack((np.ones_like(r), r, r*r, ((angle[take] - theta) / dt)**2))
        weights = np.sqrt(volume[take])
        coefficients, _, rank, singular = np.linalg.lstsq(
            design * weights[:, None], hoop[take] * weights, rcond=None
        )
        if rank != design.shape[1]:
            raise ValueError("rank-deficient wall recovery patch")
        estimates.append(dict(wall_kpa=float(coefficients[0]), gp_count=int(take.sum()),
                              condition=float(singular[0] / singular[-1])))
    return estimates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for nt, nr in ((4, 2), (8, 4), (16, 8), (24, 16)):
        folder = args.output / f"theta{nt}-radial{nr}"
        folder.mkdir()
        report = study(args.solver.resolve(), folder, None, divisions=(nt, nr, 1))
        stages = []
        for j, stage in enumerate(report["stages"]):
            result = json.loads((folder / f"release-{j}" / "result.json").read_bytes())
            estimates = extrapolate(result, nt, nr)
            stages.append(dict(release=stage["release"], free_relative=stage["free_relative"],
                               linear_relative=stage["linear_relative"],
                               wall_gp_average_kpa=stage["wall_hoop_recovered_kpa"],
                               wall_reference_kpa=stage["wall_hoop_exact_kpa"],
                               extrapolation=estimates))
        rows.append(dict(theta_cells=nt, radial_cells=nr, axial_cells=1,
                         first_ring_width_m=A*((R/A)**(1/nr)-1),
                         retained_tets=report["retained_tets"], initial=report["initial"],
                         stages=stages))
    (args.output / "kirsch-wall.json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
