"""SSLS111A paired solves and fixed-displacement readouts for both eccentric layers."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-23-b01-b02-refined-cases"))
from refined import nodal_values, projected, relative, shape

MODES = ("centroid", "element_mean")


def plate(example, grid, mode):
    args = [example, grid]
    if mode == "element_mean":
        args.append("--shear-mean")
    args.append("--compare-shear")
    raw = json.loads(subprocess.check_output(args))
    layers = []
    for row in raw["layers"]:
        q, m = nodal_values(row["root_samples"])[0]
        opposite = nodal_values(row["other_interpolation_same_solution"])[0]
        exact = row["exact_nodal_field"]
        center = shape(1/3, 1/3)[0] @ nodal_values(exact)
        layers.append(dict(offset_m=row["offset_m"], qx=relative(q, -500., .001),
                           mxx=relative(m, 5000., .001), other_readout=opposite.tolist(),
                           exact_nodes_center_qx=float(center[0]),
                           exact_nodes_gp_qx=projected(exact)[:, 0].tolist()))
    return dict(elements=raw["triangles"], layers=layers,
                a3=relative(raw["tip_dz_m"], -3.90625e-5, .005),
                root_reactions=[raw["root_force_z_kn"], raw["root_moment_y_kn_m"]])


def compare(example, grid):
    old, new = (plate(example, grid, mode) for mode in MODES)
    contributions = []
    for a, b in zip(old["layers"], new["layers"], strict=True):
        assert a["offset_m"] == b["offset_m"]
        split = {}
        for i, field in enumerate(("qx", "mxx")):
            before, after = a[field]["measured"], b[field]["measured"]
            # Traverse the two axes from old/old via new solve/old readout.
            mixed = b["other_readout"][i]
            initial_error = before - a[field]["reference"]
            remaining_error = after - a[field]["reference"]
            split[field] = dict(initial_error=initial_error, remaining_error=remaining_error,
                                solve_change=mixed-before, readout_change=after-mixed,
                                readout_change_on_old_solution=a["other_readout"][i]-before,
                                removed_absolute_error_percent=100*(1-abs(remaining_error/initial_error)))
        contributions.append(dict(offset_m=a["offset_m"], **split))
    return dict(centroid=old, element_mean=new, contributions=contributions)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("--release", type=Path, default=Path("target/release"))
    args = p.parse_args()
    exe = args.release.resolve() / "examples" / ("t10_layered_cantilever" + (".exe" if os.name == "nt" else ""))
    result = {grid: compare(exe, grid) for grid in ("12x6", "18x9")}
    with args.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
