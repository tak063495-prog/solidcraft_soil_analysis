"""Re-evaluate the single sub-3% grid from the three predetermined grading families."""

import json
from pathlib import Path
import sys

import numpy as np

root = Path.cwd()
sys.path[:0] = [
    str(root / "interim" / name)
    for name in (
        "2026-11-25-approved-six",
        "2026-11-24-approved-four",
        "2026-11-24-t9d-benchmarks",
    )
]
from kirsch import A, R, study
from kirsch_refinement import score
from kirsch_wall import extrapolate

folder = Path(sys.argv[1])
folder.mkdir()
nt, nr = 22, 21
radii = A / np.linspace(1, A / R, nr + 1)
radii[[0, -1]] = [A, R]
case = study(
    root / "target/release/soil-fem-3d",
    folder,
    None,
    divisions=(nt, nr, 1),
    radii=radii,
)
read = lambda name: json.loads((folder / name).read_bytes())
final = read("release-2/result.json")
fields = score(read("release-2/mesh.json"), final, read("intact/result.json"))
wall = [e["wall_kpa"] for e in extrapolate(final, nt, nr)]
print(
    json.dumps(
        dict(
            theta=nt,
            radial=nr,
            retained_elements=case["retained_tets"],
            removed_elements=case["removed_tets"],
            radii=radii.tolist(),
            linear_solves=4,
            wall_kpa=wall,
            maximum_free_force=max(
                s["free_relative"] for s in [case["initial"], *case["stages"]]
            ),
            maximum_true_linear=max(
                s["linear_relative"] for s in [case["initial"], *case["stages"]]
            ),
            **fields,
        ),
        indent=2,
    )
)
