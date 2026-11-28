"""Compare the final *production* fields on the coarsest physical quadrature set."""

import argparse
import json
from pathlib import Path

import numpy as np
from t6_reference import shape
from t6_study import norms, symmetric


def fields(folder, targets=None):
    run = json.loads((folder / "run.json").read_bytes())
    g = run["frames"][-1]["active"]
    q = np.asarray(run["frames"][-1]["configuration"]).reshape(-1, 6)
    elems = []
    points = []
    for e in g["elements"]:
        x = np.array(e["positions_m"])
        dx = (x[1:3] - x[0]).T
        inverse = np.linalg.pinv(dx)
        rs = [s["rs"] for s in e["resultants"]]
        fit = np.linalg.inv(np.array([shape(*p)[0] for p in rs]))
        values = []
        for sec in e["resultants"]:
            axes = np.array(sec["mid_axes"])[:2]
            values.append(
                np.r_[
                    (axes.T @ symmetric(sec["n_kn_m"]) @ axes).ravel(),
                    (axes.T @ symmetric(sec["m_kn"]) @ axes).ravel(),
                    axes.T @ sec["q_kn_m"],
                ]
            )
        elems.append(
            (e["face"]["group"], x[0], inverse, q[e["face"]["nodes"]], fit @ values)
        )
        area = np.linalg.norm(np.cross(*dx.T)) / 2
        for j, p in enumerate(rs):
            points.append(
                (
                    e["face"]["group"],
                    x[0] + dx @ p,
                    area * (0.109951743655322 if j < 3 else 0.223381589678011),
                )
            )
    values = []
    for group, x, _ in points if targets is None else targets:
        for wing, o, inv, u, section in elems:
            if group != wing:
                continue
            r, s = inv @ (x - o)
            if min(r, s, 1 - r - s) >= -1e-12:
                h = shape(r, s)[0]
                values.append(np.r_[h @ u, h @ section])
                break
        else:
            raise ValueError("physical point outside declared wing")
    return points, np.array(values)


def estimate(a, b, c):
    d1 = float(a - b)
    d2 = float(b - c)
    if d1 * d2 <= 0 or abs(d2) >= abs(d1):
        return dict(order=None, extrapolated=None)
    p = float(np.log2(d1 / d2))
    return dict(order=p, extrapolated=float(c + (c - b) / (2**p - 1)))


def measure(root):
    rows = json.loads((root / "series.json").read_bytes())
    out = []
    for distorted in [False, True]:
        data = [r for r in rows if r["distorted"] == distorted]
        coords = None
        vals = []
        for n in [2, 4, 8]:
            pts, v = fields(root / f"n{n}-distort{int(distorted)}", coords)
            if coords is None:
                coords = pts
            vals.append(v)
        w = np.array([p[2] for p in coords])
        diff = []
        for first, second in zip(vals, vals[1:]):
            diff.append(
                {
                    k: norms(second[:, s] - first[:, s], w, first[:, s])
                    for k, s in [
                        ("u", slice(3)),
                        ("rotation", slice(3, 6)),
                        ("N", slice(6, 15)),
                        ("M", slice(15, 24)),
                        ("Q", slice(24, 27)),
                    ]
                }
            )
        # Refinement ratio exactly two for the regular affine family only.
        est = (
            {
                k: estimate(*(r["midpoint_global_displacement"][j] for r in data))
                for j, k in enumerate("xyz")
            }
            if not distorted
            else None
        )
        out.append(
            dict(
                distorted=distorted,
                common_points=len(coords),
                adjacent_field_differences=diff,
                joint_midpoint_xyz=est,
                compliance=estimate(*(r["compliance"] for r in data))
                if not distorted
                else None,
            )
        )
    (root / "convergence.json").write_text(
        json.dumps(out, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path)
    measure(p.parse_args().root)
