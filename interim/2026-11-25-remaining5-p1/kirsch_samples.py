"""Archive the native raw tensors at the twelve original scoring coordinates."""

import json
from pathlib import Path
import sys
import numpy as np

folder = Path(sys.argv[1]) / "release-2"
mesh, result = [
    json.loads((folder / f"{n}.json").read_bytes()) for n in ("mesh", "result")
]
cells = np.array(mesh["tetrahedra"])
v = np.array(mesh["nodes"])[cells]
inv = np.linalg.inv((v[:, 1:] - v[:, :1]).transpose(0, 2, 1))
stress = {
    p["element_index"]: s
    for p, s in zip(
        result["integration_points"], result["stress_tension_positive_kpa"], strict=True
    )
}
rows = []
for r in (1.1, 1.5, 2.25):
    for deg in (7, 27, 52, 79):
        angle = np.deg2rad(deg)
        x = np.array([r * np.cos(angle), r * np.sin(angle), 0.37])
        natural = np.einsum("eij,ej->ei", inv, x - v[:, 0])
        bary = np.column_stack((1 - natural.sum(axis=1), natural))
        (elem,) = np.flatnonzero((bary > 0).all(axis=1))
        rows.append(
            dict(
                x_m=x.tolist(),
                element=int(elem),
                barycentric=bary[elem].tolist(),
                raw_stress_tension_positive_kpa=stress[elem],
            )
        )
print(json.dumps(rows, indent=2))
