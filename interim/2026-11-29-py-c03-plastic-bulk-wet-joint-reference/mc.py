"""Independent small-strain, perfect non-associated Mohr-Coulomb reference.

Tension-positive symmetric tensors, kPa, tensor shear (not engineering gamma).
References and limitations: README.ja.md. No production solver imports.
"""
from dataclasses import dataclass
from itertools import combinations, permutations
import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Material:
    E: float = 100000.
    nu: float = .3
    c: float = 30.
    phi: float = 30.
    psi: float = 5.

    @property
    def G(self):
        return self.E / (2 * (1 + self.nu))

    @property
    def lam(self):
        return self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))

    @property
    def C(self):
        return 2 * self.G * np.eye(3) + self.lam * np.ones((3, 3))

    def elastic(self, e):
        return 2 * self.G * e + self.lam * np.trace(e) * np.eye(3)

    def compliance(self, s):
        return (s - self.nu / (1 + self.nu) * np.trace(s) * np.eye(3)) / (2 * self.G)

    def planes(self, angle):
        a = np.sin(np.deg2rad(angle))
        rows = []
        for i, j in permutations(range(3), 2):
            r = np.zeros(3)
            r[i], r[j] = 1 + a, -(1 - a)
            rows.append(r)
        return np.array(rows)

    @property
    def k(self):
        return 2 * self.c * np.cos(np.deg2rad(self.phi))


def yield_value(s, mat=Material()):
    return float(np.max(mat.planes(mat.phi) @ np.linalg.eigvalsh(s) - mat.k))


def implicit(old, de, mat=Material(), tangent=False):
    """Backward Euler: enumerate face, two-face edge and apex cones.

    Six principal-space inequalities, nonnegative multipliers and all inactive
    inequalities checked. At most three generators suffice in R^3. No tension
    cutoff, hardening, smoothing or production active-set implementation.
    """
    trial = old + mat.elastic(de)
    t, Q = np.linalg.eigh(trial)
    F, G, C = mat.planes(mat.phi), mat.planes(mat.psi), mat.C
    v = F @ t - mat.k
    active, mult = (), np.zeros(0)
    s, D = t.copy(), np.eye(3)
    if v.max() > 1e-9:
        found = False
        for count in (1, 2, 3):
            for ids in combinations(range(6), count):
                A, B = F[list(ids)], G[list(ids)]
                H = A @ C @ B.T
                if np.linalg.matrix_rank(H) < count:
                    continue
                lm = np.linalg.solve(H, v[list(ids)])
                candidate = t - C @ B.T @ lm
                if lm.min() >= -1e-13 and (F @ candidate - mat.k).max() < 1e-7:
                    if np.min(np.diff(candidate)) < -1e-7:
                        continue
                    s, active, mult = candidate, ids, lm
                    D = np.eye(3) - C @ B.T @ np.linalg.solve(H, A)
                    found = True
                    break
            if found:
                break
        if not found:
            raise RuntimeError("MC return: no admissible face/edge/apex")
    stress = (Q * s) @ Q.T
    dep = mat.compliance(trial - stress)
    branch = ('elastic', 'face', 'edge', 'apex')[len(active)]
    result = dict(stress=stress, dep=dep, branch=branch,
                  work=float(np.sum(stress * dep)), multipliers=mult)
    if tangent:
        def apply(direction):
            h = Q.T @ mat.elastic(direction) @ Q
            out = np.diag(D @ np.diag(h))
            for i, j in combinations(range(3), 2):
                if abs(t[i] - t[j]) < 1e-10:
                    # Repeated elastic eigenvalues have identity spectral map.
                    # Plastic branch crossings at repeated trial values have no
                    # unique classical derivative and must be tested one-sided.
                    if active:
                        raise ValueError("nondifferentiable repeated trial eigenvalue")
                    ratio = 1.
                else:
                    ratio = (s[i] - s[j]) / (t[i] - t[j])
                out[i, j] = out[j, i] = ratio * h[i, j]
            return Q @ out @ Q.T
        result['tangent'] = apply
    return result


def explicit(old, de, subdivisions, mat=Material()):
    """Forward Euler rate integration; never calls implicit or projects stress.

    Elastic yield entry is located by scalar bisection/Brent. Left-endpoint
    spectral normals are frozen for a linearized consistency correction solved
    by projected Gauss-Seidel. This is a stabilized explicit Euler scheme, not
    an end-stress return: no eigenvector iteration, no implicit call. Residual
    nonlinear yield drift is reported. The O(h) drift stabilization is included
    explicitly in the right hand side; see README for its limitations.
    """
    if not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("positive integer subdivisions required")
    s = old.copy()
    ep = np.zeros((3, 3))
    work, drift = 0., 0.
    F, G = mat.planes(mat.phi), mat.planes(mat.psi)
    for _ in range(subdivisions):
        h = de / subdivisions
        f0, f1 = yield_value(s, mat), yield_value(s + mat.elastic(h), mat)
        if f0 < -1e-8:
            if f1 <= 0:
                s += mat.elastic(h)
                continue
            a = brentq(lambda x: yield_value(s + mat.elastic(x*h), mat), 0., 1.)
            s += mat.elastic(a*h)
            h = (1-a)*h
        vals, Q = np.linalg.eigh(s)
        fs = F @ vals - mat.k
        A, B = F, G
        rate = np.diag(Q.T @ mat.elastic(h) @ Q)
        rhs, H = fs + A @ rate, A @ mat.C @ B.T
        lm = np.zeros(6)
        for _iteration in range(2000):
            previous = lm.copy()
            for j in range(6):
                lm[j] = max(0., lm[j] + (rhs[j] - H[j] @ lm) / H[j, j])
            if np.max(abs(lm-previous)) < 1e-16:
                break
        else:
            raise RuntimeError("explicit rate complementarity did not converge")
        if np.max(rhs-H@lm)>1e-7 or np.max(abs(lm*(rhs-H@lm)))>1e-9:
            raise RuntimeError("explicit frozen-normal complementarity residual")
        dp = (Q * (B.T @ lm)) @ Q.T
        s += mat.elastic(h - dp)
        ep += dp
        work += float(np.sum(s * dp))
        drift = max(drift, yield_value(s, mat))
    return dict(stress=s, dep=ep, work=work, max_yield_drift=drift)


def hydrostatic(z, head=.75, gamma_w=10.):
    return gamma_w * max(head-z, 0.)


def wet_update(old_effective, de, z, alpha=1., mat=Material()):
    result = implicit(old_effective, de, mat)
    p = hydrostatic(z)
    return dict(**result, effective=result['stress'],
                total=result['stress']-alpha*p*np.eye(3), pore=p)
