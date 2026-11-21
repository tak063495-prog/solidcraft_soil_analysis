"""Input-only affine Mindlin shell reference, with explicit line multipliers.

Integrate the planar thickness polynomial analytically. The transverse field
uses the eight MITC6 tangential samples in natural coordinates. The solver
never reads a production stiffness matrix, state, or applied nodal vector.
"""

import numpy as np
from scipy import sparse
from scipy.linalg import null_space
from scipy.sparse.linalg import spsolve

EYE = np.eye(3)
GAUSS = np.polynomial.legendre.leggauss(3)


def skew(v):
    return np.cross(v, EYE).T


def shape(r, s):
    l = 1 - r - s
    return (
        np.array(
            [
                l * (2 * l - 1),
                r * (2 * r - 1),
                s * (2 * s - 1),
                4 * l * r,
                4 * r * s,
                4 * s * l,
            ]
        ),
        np.array(
            [
                [1 - 4 * l, 4 * r - 1, 0, 4 * (l - r), 4 * s, -4 * s],
                [1 - 4 * l, 0, 4 * s - 1, -4 * r, 4 * r, 4 * (l - s)],
            ]
        ),
    )


def interpolation(r, s):
    return np.array(
        [[1, r, s, 0, 0, 0, s * r, s * s], [0, 0, 0, 1, r, s, -r * r, -r * s]]
    )


class Element:
    def __init__(self, xyz, face, sec):
        self.ids = np.asarray(face["nodes"])
        self.dofs = (6 * self.ids[:, None] + np.arange(6)).ravel()
        self.e = face["offset_m"]
        dx = (xyz[self.ids[1:3]] - xyz[self.ids[0]]).T
        cross = np.cross(dx[:, 0], dx[:, 1])
        self.area = np.linalg.norm(cross) / 2
        self.d = cross / (2 * self.area)
        ex = dx[:, 0] / np.linalg.norm(dx[:, 0])
        self.basis = np.array([ex, np.cross(self.d, ex)])
        self.j = self.basis @ dx
        self.inv = np.linalg.inv(self.j).T
        self.t = sec["thickness_m"]
        E, nu = sec["young_kpa"], sec["poisson"]
        self.elastic = (
            E / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
        )
        self.g = E / (2 * (1 + nu))
        self.alpha = sec["drilling_factor"]

        def raw(r, s):
            h, dh = shape(r, s)
            out = np.zeros((2, 6, 6))
            out[:, :, :3] = dh[:, :, None] * self.d
            out[:, :, 3:] = h[None, :, None] * np.cross(self.d, dx.T)[:, None, :]
            return out.reshape(2, 36)

        mat, rhs = [], []
        for start, tangent in [((0, 0), (1, 0)), ((0, 0), (0, 1)), ((1, 0), (-1, 1))]:
            for a in (1 - 1 / np.sqrt(3), 1 + 1 / np.sqrt(3)):
                p = np.array(start) + np.array(tangent) * a / 2
                mat.append(tangent @ interpolation(*p))
                rhs.append(tangent @ raw(*p))
        p = (1 / 3, 1 / 3)
        self.coeff = np.linalg.solve(
            np.vstack([mat, interpolation(*p)]), np.vstack([rhs, raw(*p)])
        )
        self.samples = []
        self.k = np.zeros((36, 36))
        self.drill = np.zeros_like(self.k)
        for a, w in [
            (0.091576213509771, 0.109951743655322),
            (0.445948490915965, 0.223381589678011),
        ]:
            for r, s in [(a, a), (1 - 2 * a, a), (a, 1 - 2 * a)]:
                B, D, S, spin = self.strains(r, s)
                mid = B - self.e * D
                wt = w * self.area
                self.k += wt * (
                    self.t * (mid.T @ self.elastic @ mid + self.g * 5 / 6 * S.T @ S)
                    + self.t**3 / 12 * D.T @ self.elastic @ D
                )
                self.drill += wt * self.t * self.g * self.alpha * np.outer(spin, spin)
                self.samples.append((r, s, wt, B, D, S))

    def strains(self, r, s):
        h, dh = shape(r, s)
        dx, dy = self.inv @ dh
        ex, ey = self.basis

        def membrane(a, b):
            return np.stack(
                [dx[:, None] * a, dy[:, None] * b, dy[:, None] * a + dx[:, None] * b]
            )

        B, D = np.zeros((3, 6, 6)), np.zeros((3, 6, 6))
        B[:, :, :3] = membrane(ex, ey)
        D[:, :, 3:] = membrane(np.cross(self.d, ex), np.cross(self.d, ey))
        spin = np.zeros((6, 6))
        spin[:, :3] = (dy[:, None] * ex - dx[:, None] * ey) / 2
        spin[:, 3:] = h[:, None] * self.d
        return (
            B.reshape(3, 36),
            D.reshape(3, 36),
            self.inv @ interpolation(r, s) @ self.coeff,
            spin.ravel(),
        )


def system(case):
    shell, joint = case["shell"], case["joint"]
    xyz = np.array(shell["nodes"])
    size = 6 * len(xyz)
    section = {w["group"]: w.get("section", shell["section"]) for w in joint["wings"]}
    elems = [Element(xyz, f, section[f["group"]]) for f in shell["faces"]]

    def assemble(matrices):
        rr, cc, values = [], [], []
        for e, k in zip(elems, matrices):
            rr.extend(np.repeat(e.dofs, 36))
            cc.extend(np.tile(e.dofs, 36))
            values.extend(k.ravel())
        return sparse.coo_matrix((values, (rr, cc)), shape=(size, size)).tocsc()

    physical = assemble(e.k for e in elems)
    drill = assemble(e.drill for e in elems)
    points = np.array(joint["line_m"])
    t = points[-1] - points[0]
    t /= np.linalg.norm(t)
    s = (points - points[0]) @ t
    wings = [np.array(w["edge_nodes"]) for w in joint["wings"]]
    arms = [xyz[ids] - points for ids in wings]
    B = []
    law = joint["rotation"]["kind"]
    for i in range(len(s)):
        for side in range(1, len(wings)):
            axes = EYE if law == "rigid" else null_space(t[None, :]).T
            block = np.zeros((3 + len(axes), size))
            for a, sign in [(0, -1), (side, 1)]:
                k = wings[a][i] * 6
                block[:3, k : k + 3] = sign * EYE
                block[:3, k + 3 : k + 6] = sign * skew(arms[a][i])
                block[3:, k + 3 : k + 6] = sign * axes
            B.extend(block)
    C = sparse.csr_matrix(np.array(B))
    f = np.array(shell["steps"][1]["nodal_load"])
    spring = sparse.lil_matrix((size, size))
    for i in range(0, len(s) - 2, 2):
        L = s[i + 2] - s[i]
        if law == "spring":
            idx = np.concatenate(
                [(6 * ids[i : i + 3, None] + [3, 4, 5]).ravel() for ids in wings]
            )
            mass = L / 30 * np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]])
            block = joint["rotation"]["stiffness_kn_m_per_rad_m"] * np.kron(
                mass, np.outer(t, t)
            )
            spring[np.ix_(idx, idx)] += np.block([[block, -block], [-block, block]])
        for load in joint["loads"].get(2, joint["loads"].get("2", [])):
            side = next(
                a for a, w in enumerate(joint["wings"]) if w["group"] == load["wing"]
            )
            lo = max(s[i], load["interval_m"][0])
            hi = min(s[i + 2], load["interval_m"][1])
            if hi <= lo:
                continue

            # Antiderivatives of the three quadratic line basis functions.
            def primitive(x):
                return L * np.array(
                    [
                        x - 1.5 * x * x + 2 * x**3 / 3,
                        2 * x * x - 4 * x**3 / 3,
                        2 * x**3 / 3 - x * x / 2,
                    ]
                )

            weights = primitive((hi - s[i]) / L) - primitive((lo - s[i]) / L)
            for j, w in enumerate(weights):
                k = 6 * wings[side][i + j]
                F = w * np.array(load["force_kn_m"])
                f[k : k + 3] += F
                f[k + 3 : k + 6] += w * np.array(load["couple_kn_m_per_m"]) - np.cross(
                    arms[side][i + j], F
                )
    return elems, physical, drill, spring.tocsc(), C, f


def solve(case):
    elems, K, A, S, C, f = system(case)
    free = case["shell"]["steps"][1]["free_dofs"]
    k = (K + A + S).tocsc()
    sub = C[:, free]
    matrix = sparse.bmat([[k[free, :][:, free], sub.T], [sub, None]], format="csc")
    x = spsolve(matrix, np.r_[f[free], np.zeros(C.shape[0])])
    q = np.zeros(f.size)
    q[free] = x[: len(free)]
    transfer = -C.T @ x[len(free) :]
    return dict(
        q=q,
        transfer=transfer - S @ q,
        support=k @ q - f - transfer,
        load=f,
        physical=q @ (K @ q) / 2,
        artificial=q @ (A @ q) / 2,
        spring=q @ (S @ q) / 2,
        constraint=float(np.max(np.abs(C @ q))),
        elements=elems,
    )


def five(case):
    # Only membrane translations and tangent rotations on isolated faces.
    # Shared line variables are u_J, theta_J and (for releases) one relative
    # angle. Project the latter variables onto each physical tangent plane.
    elems, K, _, S, _, f = system(case)
    joint = case["joint"]
    xyz = np.array(case["shell"]["nodes"])
    free = set(case["shell"]["steps"][1]["free_dofs"])
    d = {node: e.d for e in elems for node in e.ids}
    seam = {node for w in joint["wings"] for node in w["edge_nodes"]}
    cols = []

    def new():
        cols.append(np.zeros(f.size))
        return cols[-1]

    for i in range(len(xyz)):
        if i in seam or 6 * i not in free:
            continue
        for v in EYE:
            new()[6 * i : 6 * i + 3] = v
        for v in null_space(d[i][None, :]).T:
            new()[6 * i + 3 : 6 * i + 6] = v
    t = np.array(joint["line_m"][-1]) - joint["line_m"][0]
    t /= np.linalg.norm(t)
    for i, j in enumerate(joint["line_m"]):
        nodes = [w["edge_nodes"][i] for w in joint["wings"]]
        for axis in EYE:
            col = new()
            for n in nodes:
                col[6 * n : 6 * n + 3] = axis
        axes = [(EYE[k], nodes) for k in range(3)]
        if joint["rotation"]["kind"] != "rigid":
            axes.append((t, nodes[1:]))
        for axis, members in axes:
            col = new()
            for n in members:
                col[6 * n : 6 * n + 3] = np.cross(axis, xyz[n] - j)
                col[6 * n + 3 : 6 * n + 6] = (
                    axis  # K observes only its tangent projection. The line couple uses the full spin.
                )
    L = np.array(cols).T
    k = L.T @ ((K + S) @ L)
    values, vectors = np.linalg.eigh(k)
    keep = values > max(k.shape) * np.finfo(float).eps * max(values)
    load = L.T @ f
    if np.linalg.norm(vectors[:, ~keep].T @ load) > 1e-10:
        raise ValueError("load on unobservable five-DOF mode")
    u = L @ (vectors[:, keep] @ ((vectors[:, keep].T @ load) / values[keep]))
    return dict(
        compliance=float(f @ u),
        nullity=int(sum(~keep)),
        physical=float(u @ (K @ u) / 2),
        spring=float(u @ (S @ u) / 2),
    )


def modes(case):
    """Eigenvalues on the exact constraint nullspace before/after support."""
    _, K, A, S, C, _ = system(case)
    k = (K + A + S).toarray()
    result = {}
    for name, ids in [
        ("free", np.arange(len(k))),
        ("supported", case["shell"]["steps"][1]["free_dofs"]),
    ]:
        basis = null_space(C[:, ids].toarray())
        reduced = basis.T @ k[np.ix_(ids, ids)] @ basis
        eig = np.linalg.eigvalsh(reduced)
        # Numerical rank diagnostic at the machine-precision matrix scale.
        tol = len(eig) * np.finfo(float).eps * max(abs(eig))
        result[name] = dict(
            nullity=int(sum(abs(eig) <= tol)),
            tolerance=float(tol),
            eigenvalues=eig.tolist(),
        )
    return result
