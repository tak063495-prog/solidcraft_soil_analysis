"""Separate stress interpolation, polygonal boundary and equilibrium effects at the unchanged probes."""

import sys, json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

sys.path[:0] = [
    str(Path.cwd() / "interim/2026-11-24-t9d-benchmarks"),
    str(Path(__file__).resolve().parent),
]
from kirsch import geometry, boundary, traction, solution, A, R, K, G


def exterior_gauss8(xyz, faces, radius):
    abscissa, weight = np.polynomial.legendre.leggauss(8)
    t = (abscissa + 1) / 2
    w = weight / 2
    load = np.zeros_like(xyz)
    for _, nodes in faces.values():
        vertices = xyz[nodes]
        if not np.allclose(np.linalg.norm(vertices[:, :2], axis=1), radius):
            continue
        normal = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
        for i, u in enumerate(t):
            for j, v in enumerate(t):
                bary = np.array([1 - u, u * (1 - v), u * v])
                force = solution(bary @ vertices)[0] @ normal[:2]
                load[nodes, :2] += np.outer(bary, force) * u * w[i] * w[j]
    return load[:, :2].ravel()


class ProbeOnFacet(ValueError):
    """The archived strict-interior point locator did not select one cell."""


def run(nt, rs, *, interpolation_only=False):
    xyz, count, ring, _ = geometry(None, divisions=(nt, len(rs) - 1, 1), radii=rs)
    xyz = xyz[:count]
    cells = np.array(ring)
    verts = xyz[cells]
    dx = verts[:, 1:] - verts[:, :1]
    vol = np.abs(np.linalg.det(dx)) / 6
    inv = np.linalg.inv(dx.transpose(0, 2, 1))
    grad = np.concatenate((-inv.sum(axis=1, keepdims=True), inv), axis=1)
    lame = K - 2 * G / 3
    g = grad[:, :, :2]
    xy = np.array(
        [
            [r * np.cos(t), r * np.sin(t), 0.37]
            for r in (1.1, 1.5, 2.25)
            for t in np.deg2rad([7, 27, 52, 79])
        ]
    )
    ids = []
    incident = []
    for ix, x in enumerate(xy):
        local = np.einsum("eij,ej->ei", inv, x - verts[:, 0])
        bary = np.column_stack((1 - local.sum(axis=1), local))
        selected = np.flatnonzero((bary > 0).all(axis=1))
        if len(selected) != 1:
            raise ProbeOnFacet(
                f"probe {ix} has {len(selected)} strict-interior elements"
            )
        ids.append(selected.item())
        # Read both one-sided traces; this tolerance never enters the gate.
        adjacent = np.flatnonzero(bary.min(axis=1) >= -1e-12)
        incident.append(adjacent)
    ref = np.array([solution(x)[0] for x in xy])

    def stress(u, which=ids):
        du = np.einsum("eia,eib->eab", u[cells[which]], g[which])
        return lame * np.trace(du, axis1=1, axis2=2)[:, None, None] * np.eye(2) + G * (
            du + du.transpose(0, 2, 1)
        )

    from kirsch import S0

    interp = stress(np.array([solution(x)[1] for x in xyz])) + S0

    def rms(d):
        return np.linalg.norm(d) / (np.sqrt(len(d)) * 30000)

    if interpolation_only:
        return dict(theta=nt, radial=len(rs) - 1, exact_node_interp=rms(interp - ref))
    ke = (
        lame * np.einsum("eia,ejb->eiajb", g, g)
        + G * np.einsum("eib,eja->eiajb", g, g)
        + G * np.einsum("eic,ejc,ab->eiajb", grad, grad, np.eye(2))
    )
    ke = (ke * vol[:, None, None, None, None]).reshape(-1, 8, 8)
    dofs = (2 * cells[:, :, None] + np.arange(2)).reshape(-1, 8)
    mat = coo_matrix(
        (
            ke.ravel(),
            (
                np.broadcast_to(dofs[:, :, None], ke.shape).ravel(),
                np.broadcast_to(dofs[:, None, :], ke.shape).ravel(),
            ),
        ),
        shape=(2 * count, 2 * count),
    ).tocsr()
    free = xyz[:, :2].ravel() != 0
    assert sum(free) <= 2000 and len(cells) not in (768, 6144)
    f = boundary(ring)
    tr = lambda key, x: solution(x)[0]
    b = traction(xyz, f, R, tr)[:, :2].ravel()
    # Analytical traction on the polygonal hole isolates the error due to
    # representing a traction-free circular wall by straight facets.
    fc = traction(xyz, f, A, tr)[:, :2].ravel()
    b8 = exterior_gauss8(xyz, f, R)
    fields = {}
    for label, load in [("raw", b), ("polygon_exact", b + fc), ("outer_gauss8", b8)]:
        u = np.zeros(mat.shape[0])
        u[free] = spsolve(mat[free][:, free], load[free])
        fields[label] = stress(u.reshape(-1, 2))
        if label == "raw":
            trace_norms = [
                (
                    np.linalg.norm(stress(u.reshape(-1, 2), els) - target, axis=(1, 2))
                    / 30000
                ).tolist()
                for els, target in zip(incident, ref, strict=True)
            ]
            trace_jump = [
                float(
                    np.linalg.norm(np.diff(stress(u.reshape(-1, 2), els), axis=0))
                    / 30000
                )
                if len(els) == 2
                else None
                for els in incident
            ]
    row = dict(
        nt=nt,
        nr=len(rs) - 1,
        free=int(sum(free)),
        trace_rms_range=[
            np.sqrt(np.mean([fn(values) ** 2 for values in trace_norms]))
            for fn in (min, max)
        ],
        facet_probes=[
            dict(
                index=i,
                elements=els.tolist(),
                selected=ids[i],
                one_sided_errors=trace_norms[i],
                jump=trace_jump[i],
            )
            for i, els in enumerate(incident)
            if len(els) > 1
        ],
        elements=len(cells),
        outer_load_quadrature=np.linalg.norm(b8 - b) / np.linalg.norm(b8),
        outer_quadrature_stress_effect=rms(fields["outer_gauss8"] - fields["raw"]),
        actual=rms(fields["raw"] - ref),
        exact_hole=rms(fields["polygon_exact"] - ref),
        boundary_effect=rms(fields["raw"] - fields["polygon_exact"]),
        exact_node_interp=rms(interp - ref),
        solved_minus_interp=rms(fields["raw"] - interp),
        point_norms=(np.linalg.norm(fields["raw"] - ref, axis=(1, 2)) / 30000).tolist(),
    )
    return row


if __name__ == "__main__":
    rings = np.geomspace(A, R, 17)

    def band_grid(inner, outer):
        return np.r_[
            np.geomspace(A, 2.5, inner + 1), np.geomspace(2.5, R, outer + 1)[1:]
        ]

    rows = []
    for nt, k in [(24, 0), (24, 3), (20, 3), (20, 8), (20, 18)]:
        rs = (
            band_grid(18, 6)
            if k == 18
            else np.sort(np.r_[rings, np.sqrt(rings[:k] * rings[1 : k + 1])])
        )
        rows.append(run(nt, rs))
    old = json.loads(
        (Path.cwd() / "interim/2026-11-25-approved-six/kirsch.json").read_bytes()
    )
    for measured, archived in zip(rows, old, strict=True):
        np.testing.assert_allclose(
            measured["actual"],
            archived["points"]["rms"]["tensor_over_30MPa"],
            rtol=1e-8,
        )
    sweep = []
    for nt in (14, 16, 18, 20, 22, 24, 26, 28, 30, 32):
        nr = 500 // nt - 1
        inner = round(0.75 * nr)
        sweep.append(run(nt, band_grid(inner, nr - inner)))
    interpolation = [
        run(20 * f, band_grid(18 * f, 6 * f), interpolation_only=True)
        for f in (1, 2, 4)
    ]
    print(
        json.dumps(
            dict(
                original_grids=rows,
                fixed_dof_sweep=sweep,
                interpolation_only=interpolation,
            ),
            indent=2,
        )
    )
