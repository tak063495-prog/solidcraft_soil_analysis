"""Run and compare degree, load-step, and face-quadrature series on 12 tetrahedra."""

import argparse
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np

REFERENCE = Path(__file__).resolve().parents[1] / "2026-11-25-approved-six"
ROOT = REFERENCE.parents[1]
sys.path.insert(0, str(REFERENCE))
from c02_geometry import lagrange
from c02_reference import reference


def integral(data):
    return np.einsum("tqc,q->tc", data["slip"], data["weight"])


def weighted_field(data):
    return data["slip"] * np.sqrt(data["weight"])[None, :, None]


def compare(first, second):
    for key in ("position", "parent", "weight"):
        np.testing.assert_array_equal(first[key], second[key])
    delta = weighted_field(first) - weighted_field(second)
    mismatch = first["status"] != second["status"]
    total = np.sum(delta**2)
    return dict(
        field_relative=np.linalg.norm(delta) / np.linalg.norm(weighted_field(second)),
        integral_relative=np.linalg.norm(integral(first) - integral(second))
        / np.linalg.norm(integral(second)),
        different_labels=int(mismatch.sum()),
        squared_difference_on_mismatched_labels=np.sum(delta[mismatch] ** 2) / total,
        squared_difference_by_stage=(
            np.sum(delta.reshape(7, 64, -1) ** 2, axis=(1, 2)) / total
        ).tolist(),
    )


def time_series(values):
    coarse, middle, fine = map(weighted_field, values)
    a, b = coarse - middle, middle - fine
    return dict(
        pairs=[compare(values[0], values[1]), compare(values[1], values[2])],
        observed_order=np.log2(np.linalg.norm(a) / np.linalg.norm(b)),
        cosine=np.vdot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)),
    )


def run(folder):
    mesh, ground, case = [
        json.loads((ROOT / "examples" / (name + ".json")).read_bytes())
        for name in (
            "joint-water-column-tet10-mesh",
            "ground-joint-water-mc",
            "stages-joint-water-cycle",
        )
    ]
    # Check quartic nodal interpolation and barycentric polynomial reproduction.
    alpha = np.array([a for a in product(range(5), repeat=4) if sum(a) == 4])
    nodes = alpha / 4
    phi, _ = lagrange(nodes, alpha, 4)
    np.testing.assert_allclose(phi, np.eye(len(alpha)), atol=1e-13)
    bary = np.array([[0.17, 0.22, 0.32, 0.29], [0.09, 0.58, 0.19, 0.14]])
    phi, deriv = lagrange(bary, alpha, 4)
    for powers in alpha:
        nodal = np.prod(nodes**powers, axis=1)
        expected = np.prod(bary**powers, axis=1)
        np.testing.assert_allclose(phi @ nodal, expected, rtol=1e-12, atol=1e-14)
        exact = expected[:, None] * powers / bary
        gr = np.einsum("qnd,n->qd", deriv, nodal)
        np.testing.assert_allclose(
            gr[:, 1:] - gr[:, :1], exact[:, 1:] - exact[:, :1], rtol=1e-11, atol=1e-14
        )
    grids = [(p, n, 6) for p in (2, 3) for n in (64, 128, 256)]
    grids += [(4, 128, 6), (4, 256, 6), (3, 256, 8), (3, 256, 10)]
    data, costs = {}, []
    folder.mkdir()
    for p, n, q in grids:
        path = folder / f"p{p}-n{n}-q{q}"
        path.mkdir()
        cost = reference(mesh, ground, case, p, n, q, path)
        assert cost["accepted"] == 7 * n <= 2065
        costs.append(cost)
        with np.load(path / "slip.npz", allow_pickle=False) as fields:
            data[p, n, q] = {k: fields[k] for k in fields}
        if p <= 3 and n <= 128:
            with np.load(
                REFERENCE / f"c02-p{p}-n{n}" / "slip.npz", allow_pickle=False
            ) as old:
                for key in old:
                    np.testing.assert_array_equal(old[key], data[p, n, q][key])
    face = [integral(data[3, 256, q]) for q in (6, 8, 10)]
    report = dict(
        costs=costs,
        time={p: time_series([data[p, n, 6] for n in (64, 128, 256)]) for p in (2, 3)},
        quartic_time=compare(data[4, 128, 6], data[4, 256, 6]),
        degree={
            n: [compare(data[p, n, 6], data[p + 1, n, 6]) for p in (2, 3)]
            for n in (128, 256)
        },
        face_integrals=[
            np.linalg.norm(a - b) / np.linalg.norm(b) for a, b in zip(face, face[1:])
        ],
    )
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_directory", type=Path)
    run(parser.parse_args().new_directory)
