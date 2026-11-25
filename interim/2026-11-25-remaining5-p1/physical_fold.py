"""Opt-in physical crown-displacement fold check for a matched refinement pair."""

import argparse
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [
    str(HERE.parent / k)
    for k in (
        "2026-11-25-approved-six",
        "2026-11-24-t9d-benchmarks",
        "2026-11-22-impl-t7-physical-instability",
    )
]
from arch import beam_curve
from curved import AnalyticChordArch
from scipy.optimize import root
from continuation import direction, Evaluation, simple_point, serial
from references import force_unit, R

TOL = 1e-9
# Relative curvature spread between adjacent grids and independent bracket estimate.
SPREAD = 0.10


def local(model, bracket, z, axis, displacement_unit, load_unit):
    a, m, b = map(np.asarray, bracket)
    ea, em, eb = (model({}, point) for point in (a, m, b))
    ev = model({}, z)
    scale = max(np.linalg.norm(ea.k), np.linalg.norm(eb.k), 1.0)
    u, s, vh = np.linalg.svd(ev.k)
    w, v = u[:, -1], vh[-1]
    wq = w @ ev.q
    dx = displacement_unit * v[axis]
    c = w @ model.tangent_action_derivative(z, v)
    transverse = abs(wq) > TOL * np.linalg.norm(ev.q)
    coord = abs(v[axis]) > TOL
    curvature = -load_unit * c / (wq * dx**2) if transverse and coord else None
    x = np.array([p[axis] * displacement_unit for p in (a, m, b)])
    y = np.array([p[-1] * load_unit for p in (a, m, b)])
    monotone = np.prod(np.diff(x)) > 0
    secant = (
        2 * np.diff(np.diff(y) / np.diff(x)).item() / (x[-1] - x[0])
        if monotone
        else None
    )
    # Resample at two nearby physical displacements about the localized candidate.
    # Widely spaced route points can have substantial third/higher derivatives.
    raw_secant = secant
    force_ratios = []
    local_xy = []
    if (
        monotone
        and coord
        and min(x[0], x[-1]) < z[axis] * displacement_unit < max(x[0], x[-1])
    ):
        step = min(abs((z[axis] * displacement_unit) - x[[0, -1]])) / 8
        for sign in (-1, 0, 1):
            target = z[axis] + sign * step / displacement_unit

            def equations(y):
                return np.r_[model({}, y).residual, y[axis] - target]

            def jacobian(y):
                e = model({}, y)
                return np.vstack((np.column_stack((e.k, e.q)), np.eye(len(y))[axis]))

            seed = z + np.r_[v * ((target - z[axis]) / v[axis]), 0.0]
            equilibrium = root(equations, seed, jac=jacobian, tol=1e-12).x
            e = model({}, equilibrium)
            force_ratios.append(np.linalg.norm(e.residual) / e.denominator)
            if (
                not e.eligible
                or e.event != ev.event
                or abs(equilibrium[axis] - target) > TOL
            ):
                break
            local_xy.append(
                [equilibrium[axis] * displacement_unit, equilibrium[-1] * load_unit]
            )
        if len(local_xy) == 3:
            xs, ys = np.array(local_xy).T
            secant = 2 * np.diff(np.diff(ys) / np.diff(xs)).item() / (xs[-1] - xs[0])
        else:
            secant = None
    else:
        secant = None
    residual = max(np.linalg.norm(e.residual) / e.denominator for e in (ea, em, eb, ev))
    chord = b - a
    try:
        ta, tb = direction(ea, chord), direction(eb, chord)
        slopes = [ta[-1], tb[-1]]
        augmented = np.vstack((np.column_stack((ev.k, ev.q)), ta))
        bordered_full = np.linalg.matrix_rank(augmented) == len(z)
    except ArithmeticError:
        slopes = [0.0, 0.0]
        bordered_full = False
    same = all(e.eligible and e.event == ea.event for e in (ea, em, eb, ev))
    checks = dict(
        reversal=np.prod(slopes) < 0 and np.prod(np.diff(y)) < 0,
        bracket=0 <= (z - a) @ chord / (chord @ chord) <= 1,
        smooth_branch=same,
        original_force=max([residual, *force_ratios]) <= TOL,
        single_null=np.count_nonzero(s <= TOL * scale) == 1,
        transverse=transverse,
        bordered_full=bordered_full,
        physical_coordinate=coord and monotone,
        signed_curvature=secant is not None
        and curvature is not None
        and secant * curvature > 0
        and abs(curvature - secant) <= SPREAD * abs(secant),
    )
    return dict(
        checks=checks,
        physical_curvature=curvature,
        bracket_curvature=secant,
        raw_route_curvature=raw_secant,
        local_xy=local_xy,
        local_force=force_ratios,
        bracket_displacement=x,
        bracket_load=y,
        force_ratio=residual,
        signed_tangent_load=slopes,
        wq=wq,
        singular_tail=s[-2:],
        null_cut=TOL * scale,
    )


def classify(row, neighbour):
    a, b = row["physical_curvature"], neighbour["physical_curvature"]
    stable = (
        a is not None
        and b is not None
        and a * b > 0
        and abs(a - b) <= SPREAD * max(abs(a), abs(b))
    )
    return {
        "adjacent_stable": stable,
        "relative_spread": abs(a - b) / max(abs(a), abs(b))
        if a is not None and b is not None and max(abs(a), abs(b)) > 0
        else None,
        "candidate": bool(
            stable and all(row["checks"].values()) and all(neighbour["checks"].values())
        ),
    }


class Quartic:
    # Stable transverse, one null mode, load reverses, second derivative vanishes.
    def __call__(self, parent, z):
        x, lam = z
        return Evaluation(
            np.array([x**4 + lam]), np.array([[4 * x**3]]), np.ones(1), {}, 1.0
        )

    def tangent_action_derivative(self, z, v):
        return np.array([12 * z[0] ** 2 * v[0] ** 2])


class Pitchfork:
    def __call__(self, parent, z):
        x, lam = z
        return Evaluation(
            np.array([x * (x * x + lam)]),
            np.array([[3 * x * x + lam]]),
            np.array([x]),
            {},
            1.0,
        )

    def tangent_action_derivative(self, z, v):
        return np.array([6 * z[0] * v[0] ** 2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--physical-fold", action="store_true", required=True)
    parser.parse_args()
    rows = {}
    for n in (12, 24, 48, 96):
        model = AnalyticChordArch(n)
        route = beam_curve(
            n,
            extrapolate=True,
            model=model,
            fold_metrics=True,
            return_fold_bracket=True,
        )
        (axis,) = np.flatnonzero(model.free == 3 * model.apex + 1)
        row = local(
            model, route["fold_bracket"], route["fold"]["z"], axis, -R, force_unit
        )
        rows[n] = dict(
            row,
            legacy=route["fold"]["kind"],
            trials=len(route["trials"]),
            legacy_normalized=route["fold"]["curvature_over_scale"],
        )
    for n, other in ((12, 24), (24, 12), (48, 24), (96, 48)):
        rows[n].update(neighbour=other, **classify(rows[n], rows[other]))
    control = []
    for h in (0.125, 0.0625):
        bracket = np.array([[-h, -(h**4)], [0, 0], [h, -(h**4)]])
        model = Quartic()
        old = simple_point(model, {}, bracket[0], bracket[-1], TOL)
        # Exact equilibrium and zero mode are known for this manufactured control.
        control.append(
            dict(local(model, bracket, np.zeros(2), 0, 1.0, 1.0), legacy=old["kind"])
        )
    pitchfork = []
    for h in (0.125, 0.0625):
        bracket = np.array([[-h, -h * h], [0, 0], [h, -h * h]])
        model = Pitchfork()
        old = simple_point(model, {}, bracket[0], bracket[-1], TOL)
        pitchfork.append(
            dict(local(model, bracket, np.zeros(2), 0, 1.0, 1.0), legacy=old["kind"])
        )
    print(
        json.dumps(
            dict(
                beams=rows,
                quartic=control,
                quartic_pair=classify(*control),
                pitchfork=pitchfork,
                pitchfork_pair=classify(*pitchfork),
                relative_spread_limit=SPREAD,
            ),
            default=serial,
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
