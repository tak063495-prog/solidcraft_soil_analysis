"""Run the four approved cases with the T10 references and widths."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from cylinders import HERE, cylinder_input, standard

ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "2026-11-23-impl-t10-reference-tolerance"))
from t10_compare import relative, shape

CYLINDERS = ("b01-regular16", "b01-edf364")
PLATES = ((12, 6), (18, 9))


def cylinder(case_id, cli):
    case, probe = cylinder_input(case_id)
    with tempfile.TemporaryDirectory() as tmp:
        frame, error = standard.execute(cli, Path(tmp) / "cylinder.json", case)
    if error is not None:
        raise RuntimeError(f"{case_id}: {error}")
    return dict(case_id=case_id, scope="C displacement only", nodes=len(case["nodes"]),
                elements=len(case["faces"]), point_m=case["nodes"][probe],
                c_uz_m=relative(frame["configuration"][6 * probe + 2], -1.8248e-5, 0.03))


def projected(samples):
    """Global x resultants per layer about the common attachment plane, in N."""
    points = []
    for sample in samples:
        x = np.asarray(sample["mid_axes"])[:2, 0]
        mx, my, mxy = sample["m_attachment_kn"]
        points.append((x @ sample["q_kn_m"], x @ [[mx, mxy], [mxy, my]] @ x))
    return 1000 * np.asarray(points)


def nodal_values(samples):
    matrix = [shape(*p["rs"])[0] for p in samples]
    return np.linalg.solve(matrix, projected(samples))


def force_pair(qm):
    return dict(qx_n_m=float(qm[0]), mxx_nm_m=float(qm[1]))


def section_traces(sections, nx):
    traces = {}
    for row in sections:
        key = (row["cut"], row["right_of_cut"], row["offset_m"])
        trace = traces.setdefault(key, [])
        # In either orientation, 0 and 2 are the ends of the vertical
        # edge and 5 is its midpoint (the same convention as the mesh).
        qm = nodal_values(row["samples"])[[0, 2, 5]]
        trace.append((row["y_m"], qm))
    report = []
    for (cut, right, offset), trace in sorted(traces.items()):
        integral = np.zeros(2)
        points = []
        for y, qm in trace:
            length = abs(y[1] - y[0])
            integral += length * (qm[0] + qm[1] + 4 * qm[2]) / 6
            points.extend(dict(y_m=coord, **force_pair(value)) for coord, value in zip(y, qm))
        report.append(dict(x_m=10 * cut / nx, trace_from="right" if right else "left",
                           offset_m=offset, mean=force_pair(integral / 5),
                           points=sorted(points, key=lambda p: p["y_m"])))
    return report


def plate(nx, ny, example):
    name = f"{nx}x{ny}"
    raw = json.loads(subprocess.check_output([example, name]))
    if raw["triangles"] != 4 * nx * ny:
        raise ValueError(f"{name}: unexpected layer mesh")
    a3 = dict(case_id=f"b02-a3-{name}", scope="A3 displacement only", nx=nx, ny=ny,
              elements=raw["triangles"], a3_dz_m=relative(raw["tip_dz_m"], -3.90625e-5, 0.005))
    layers = []
    for row in raw["layers"]:
        samples = row["root_samples"]
        q, m = nodal_values(samples)[0]
        # A1 is local vertex 0 of this one adjacent triangle in each layer.
        nearest = min(samples, key=lambda p: (p["rs"][0] * 10 / nx)**2 + (p["rs"][1] * 5 / ny)**2)
        layers.append(dict(offset_m=row["offset_m"], qx_n_m=relative(q, -500, 0.001),
                           mxx_nm_m=relative(m, 5000, 0.001),
                           nearest_gp_rs=nearest["rs"], nearest_gp=force_pair(projected([nearest])[0])))
    a1 = dict(grid=name, layers=layers, sections=section_traces(raw["sections"], nx),
              root_force_z_kn=raw["root_force_z_kn"], root_moment_y_kn_m=raw["root_moment_y_kn_m"])
    return a3, a1


def measure(release):
    suffix = ".exe" if os.name == "nt" else ""
    rows = {name: cylinder(name, release / ("soil-fem-3d" + suffix)) for name in CYLINDERS}
    diagnoses = []
    for nx, ny in PLATES:
        a3, a1 = plate(nx, ny, release / "examples" / ("t10_layered_cantilever" + suffix))
        rows[a3["case_id"]] = a3
        diagnoses.append(a1)
    rows["b02-a1"] = dict(case_id="b02-a1", scope="A1 layer QX and MXX", grids=diagnoses)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--release", type=Path, default=ROOT / "target/release")
    args = parser.parse_args()
    rows = measure(args.release.resolve())
    args.output.mkdir(parents=True, exist_ok=False)
    for name, row in rows.items():
        (args.output / (name + ".json")).write_text(
            json.dumps(row, ensure_ascii=False, allow_nan=False, indent=2) + "\n",
            encoding="utf-8", newline="\n")
