"""Independent beam diagnostics in the consistent units of each reference."""

from pathlib import Path
import argparse
import json

import numpy as np
from scipy.linalg import eigh
from scipy.integrate import solve_bvp
from continuation import Controls, Evaluation, serial, simple_point, trace


def euler(n):
    # Cubic Hermite displacement/rotation with the consistent initial-stress matrix.
    h = 1 / n
    k = np.zeros((2 * n + 2, 2 * n + 2))
    g = np.zeros_like(k)
    ke = (
        np.array(
            [
                [12, 6 * h, -12, 6 * h],
                [6 * h, 4 * h * h, -6 * h, 2 * h * h],
                [-12, -6 * h, 12, -6 * h],
                [6 * h, 2 * h * h, -6 * h, 4 * h * h],
            ]
        )
        / h**3
    )
    ge = np.array(
        [
            [36, 3 * h, -36, 3 * h],
            [3 * h, 4 * h * h, -3 * h, -h * h],
            [-36, -3 * h, 36, -3 * h],
            [3 * h, -h * h, -3 * h, 4 * h * h],
        ]
    ) / (30 * h)
    for i in range(n):
        ix = slice(2 * i, 2 * i + 4)
        k[ix, ix] += ke
        g[ix, ix] += ge
    k = k[2:, 2:]
    g = g[2:, 2:]
    loads, modes = eigh(k, g, subset_by_index=[0, 0])
    mode = modes[:, 0] / modes[-2, 0]
    transverse = np.r_[0, mode[::2]]
    exact = 1 - np.cos(np.pi * np.linspace(0, 1, n + 1) / 2)
    return dict(
        elements=n,
        Pcr=loads[0],
        exact=np.pi**2 / 4,
        relative_load_error=abs(loads[0] * 4 / np.pi**2 - 1),
        mode=transverse,
        mode_max_error=np.max(abs(transverse - exact)),
        eigen_residual=np.linalg.norm(k @ mode - loads[0] * g @ mode)
        / np.linalg.norm(k @ mode),
    )


# Abaqus/Standard 2025, snapbuckling_deep.inp. Coordinates, section and E
# determine the nondimensional BVP. Positive load points downward.
R = np.hypot(95.37169508, 30.07057995)
A = 2.289428295
EI = 1e6 * A**3 / 12
force_unit = EI / R**2
theta0 = np.arctan2(-30.07057995, 95.37169508)
span = np.pi - 2 * theta0
target = np.array([-2 * np.cos(theta0), 0.0])


def rhs(y, nx, ny):
    theta, m = y[2], y[3]
    c, s = np.cos(theta), np.sin(theta)
    stretch = 1 + (nx * c + ny * s) * force_unit / (1e6 * A)
    return np.array([stretch * c, stretch * s, 1 + m, stretch * (nx * s - ny * c)])


def integrate(z, n):
    # Integrate on the n equal arc intervals, splitting an odd central
    # interval at the actual concentrated load.
    nx, ny, m, p = z
    points = np.unique(np.r_[np.linspace(0, 1, n + 1), 0.5])
    y = np.array([0 * nx, 0 * nx, 0 * nx + theta0 + np.pi / 2, m])
    for left, right in zip(points[:-1], points[1:]):
        h = (right - left) * span
        fy = ny + (p if left >= 0.5 else 0)
        a = rhs(y, nx, fy)
        b = rhs(y + h * a / 2, nx, fy)
        c = rhs(y + h * b / 2, nx, fy)
        d = rhs(y + h * c, nx, fy)
        y = y + h * (a + 2 * b + 2 * c + d) / 6
        if right == 0.5:
            crown = y.copy()
    return np.array([y[0] - target[0], y[1], y[3]]), crown


def collocation(z):
    # Two half-arcs and the force jump solved simultaneously by collocation.
    # Unknown clamp force and moment are determined by end/continuity data.
    grid = np.linspace(0, 1, 33)
    angles = theta0 + np.array([grid, 1 + grid]) * span / 2
    seed = np.vstack(
        [
            block
            for ang in angles
            for block in (
                np.cos(ang) - np.cos(theta0),
                np.sin(ang) - np.sin(theta0),
                ang + np.pi / 2,
                np.zeros_like(grid),
            )
        ]
    )

    def ode(t, y, p):
        nx, ny = p
        return span / 2 * np.vstack((rhs(y[:4], nx, ny), rhs(y[4:], nx, ny + z[-1])))

    def bc(a, b, p):
        return np.r_[
            a[:3] - [0, 0, theta0 + np.pi / 2], a[4:] - b[:4], b[4:6] - target, b[7]
        ]

    sol = solve_bvp(ode, bc, grid, seed, p=z[:2], tol=1e-8, max_nodes=3000)
    return dict(
        converged=sol.success,
        load=z[-1] * force_unit,
        downward_crown=(1 - np.sin(theta0) - sol.y[5, 0]) * R if sol.success else None,
    )


def independent_curve():
    """Independent displacement-controlled collocation from the undeformed input."""
    grid = np.linspace(0, 1, 33)
    angles = theta0 + np.array([grid, 1 + grid]) * span / 2
    y = np.vstack(
        [
            v
            for ang in angles
            for v in (
                np.cos(ang) - np.cos(theta0),
                np.sin(ang) - np.sin(theta0),
                ang + np.pi / 2,
                np.zeros_like(grid),
            )
        ]
    )
    forces = np.zeros(3)
    curve = []
    for displacement in (5, 20, 40, 60, 80, 100, 110, 115, 120, 125, 130):

        def ode(t, y, p):
            nx, ny, load = p
            return span / 2 * np.vstack((rhs(y[:4], nx, ny), rhs(y[4:], nx, ny + load)))

        def bc(a, b, p):
            return np.r_[
                a[:3] - [0, 0, theta0 + np.pi / 2],
                a[4:] - b[:4],
                b[4:6] - target,
                b[7],
                a[5] - (1 - np.sin(theta0) - displacement / R),
            ]

        sol = solve_bvp(ode, bc, grid, y, p=forces, tol=1e-8, max_nodes=5000)
        curve.append(
            dict(
                downward_crown=displacement,
                converged=sol.success,
                load=sol.p[2] * force_unit if sol.success else None,
                max_rms_residual=float(max(sol.rms_residuals)),
            )
        )
        if not sol.success:
            break
        y = sol.sol(grid)
        forces = sol.p
    return curve


def arch(n):
    def evaluate(parent, z):
        batch = np.broadcast_to(z[:, None], (4, 5)).astype(complex).copy()
        batch[:, 1:] += 1j * 1e-24 * np.eye(4)
        values, crown = integrate(batch, n)
        jac = values[:, 1:].imag / 1e-24
        return Evaluation(
            values[:, 0].real,
            jac[:, :3],
            jac[:, 3],
            {},
            1.0,
            details=dict(crown=crown[:, 0].real),
        )

    # Solve the discrete zero-load starting equilibrium (RK quadrature error).
    from scipy.optimize import root

    seed = root(lambda x: evaluate({}, np.r_[x, 0.0]).residual, np.zeros(3)).x

    def stop(points):
        loads = [p["z"][-1] for p in points]
        return len(loads) > 6 and loads[-1] < max(loads) - 0.4

    result = trace(
        evaluate,
        {},
        np.r_[seed, 0.0],
        [0, 0, 0, 1],
        Controls(0.28, 165, 1e-10, 1e-10, 20),
        stop=stop,
    )
    zs = np.array([p["z"] for p in result["accepted"]])
    i = int(np.argmax(zs[:, -1]))
    fold = (
        simple_point(evaluate, {}, zs[i - 1], zs[i + 1], 1e-9)
        if 0 < i < len(zs) - 1
        else dict(kind="D")
    )
    curve = [
        dict(
            z=z,
            load=z[-1] * force_unit,
            downward_crown=(1 - np.sin(theta0) - integrate(z, n)[1][1]) * R,
        )
        for z in zs
    ]
    return dict(
        intervals=n,
        curve=curve,
        fold=fold,
        peak_sample=curve[i],
        termination=result["termination_reason"],
        estimated_peak_load=fold["z"][-1] * force_unit
        if fold.get("kind") == "simple fold"
        else None,
        independent_collocation=collocation(zs[min(5, len(zs) - 1)]),
        trials=result["trials"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False, parents=True)

    def save(name, obj):
        (args.output / name).write_text(
            json.dumps(obj, default=serial, allow_nan=False) + "\n", encoding="utf-8"
        )

    save("euler.json", [euler(n) for n in (4, 8, 16)])
    for n in (15, 30, 60):
        save(f"arch-{n}.json", arch(n))
    save("collocation.json", independent_curve())
