"""Standalone linear pile references. No project solver imports or subprocesses.

Sources, signs, derivations and limits: README.ja.md. Engineering acceptance
is deliberately undefined; numerical verification is not P6 acceptance.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def parameters(length, stiffness, spring, load, tip=0.0):
    values = np.array([length, stiffness, spring, load, tip], dtype=float)
    if not np.isfinite(values).all() or min(length, stiffness, spring) <= 0 or tip < 0:
        raise ValueError("finite L, EA/EI, distributed k > 0 and tip k >= 0 required")


def axial_exact(s, length, ea, spring, load, tip):
    """EA*u''-k*u=0, N=-EA*u', N(0)=P, N(L)=Kt*u(L)."""
    parameters(length, ea, spring, load, tip)
    s = stations(s, length)
    a = np.sqrt(spring / ea)
    if a * length > 20:
        raise ValueError("aL > 20 outside this small-reference conditioning scope")
    r = tip / (ea * a)
    c = load / (ea * a * np.sinh(a * length) + tip * np.cosh(a * length))
    t = a * (length - s)
    u = c * (np.cosh(t) + r * np.sinh(t))
    n = ea * a * c * (np.sinh(t) + r * np.cosh(t))
    return np.array([u, n, -spring * u]).T


def stations(s, length):
    s = np.atleast_1d(np.asarray(s, dtype=float))
    if s.ndim != 1 or not np.isfinite(s).all() or np.any(s < 0) or np.any(s > length):
        raise ValueError("stations must lie in [0,L]")
    return s


def beam_exact(s, length, ei, spring, load):
    """Finite free-tip beam, M(0)=M(L)=V(L)=0, V(0)=-H.

    Matrix exponential of the constant-coefficient ODE in t=s/L;
    independent of the finite element assembly and interpolation.
    """
    parameters(length, ei, spring, load)
    s = stations(s, length)
    alpha = spring * length**4 / ei
    if alpha**0.25 > 5:
        raise ValueError("k*L^4/EI > 625 outside transfer-matrix conditioning scope")
    a = np.zeros((4, 4))
    a[:3, 1:] = np.eye(3)
    a[3, 0] = -alpha
    transfer = expm(a)
    known = np.array([0., 0., 0., load * length**3 / ei])
    initial = known.copy()
    initial[:2] = np.linalg.solve(transfer[2:, :2], -(transfer @ known)[2:])
    state = np.array([expm(a * (x / length)) @ initial for x in s])
    return np.column_stack((state[:, 0], state[:, 1] / length,
                            -ei * state[:, 2] / length**2,
                            -ei * state[:, 3] / length**3, -spring * state[:, 0]))


def beam_semi_infinite(s, ei, spring, load):
    """FHWA-HRT-06-106 ch.6 eq.(8), zero end moment; not finite L."""
    parameters(1., ei, spring, load)
    s = np.atleast_1d(np.asarray(s, dtype=float))
    if s.ndim != 1 or not np.isfinite(s).all() or np.any(s < 0):
        raise ValueError("nonnegative finite stations required")
    beta = (spring / (4 * ei))**0.25
    root = beta * (-1 + 1j)
    scale = load / (2 * ei * beta**3)
    derivatives = [scale * np.real(root**j * np.exp(root * s)) for j in range(4)]
    return np.column_stack((derivatives[0], derivatives[1], -ei * derivatives[2],
                            -ei * derivatives[3], -spring * derivatives[0]))


def hermite(t, h):
    """C1 cubic, physical slope DOFs; derivatives with respect to s."""
    shape = np.array([1-3*t*t+2*t**3, h*(t-2*t*t+t**3),
                      3*t*t-2*t**3, h*(-t*t+t**3)])
    first = np.array([(-6*t+6*t*t)/h, 1-4*t+3*t*t,
                      (6*t-6*t*t)/h, -2*t+3*t*t])
    second = np.array([(-6+12*t)/h**2, (-4+6*t)/h,
                       (6-12*t)/h**2, (-2+6*t)/h])
    third = np.array([12/h**3, 6/h**2, -12/h**3, 6/h**2])
    return shape, first, second, third


def finite_elements(kind, n, length, stiffness, spring, load, tip=0.0):
    parameters(length, stiffness, spring, load, tip)
    if kind not in ("axial", "beam") or isinstance(n, bool) or not isinstance(n, int) or not 2 <= n <= 512:
        raise ValueError("kind axial/beam and integer 2 <= n <= 512 required")
    if kind == "beam" and tip != 0:
        raise ValueError("beam tip is free; no implicit tip spring")
    h = length / n
    stride = 1 if kind == "axial" else 2
    rows, cols, values = [], [], []
    g, w = leggauss(4)  # exact polynomial integration of cubic foundation
    for e in range(n):
        ids = np.arange(stride*e, stride*(e+2))
        ke = np.zeros((2*stride, 2*stride))
        for t, weight in zip((g+1)/2, w*h/2):
            if kind == "axial":
                shape, derivative = np.array([1-t, t]), np.array([-1/h, 1/h])
            else:
                shape, _, derivative, _ = hermite(t, h)
            ke += weight * (stiffness*np.outer(derivative, derivative) + spring*np.outer(shape, shape))
        rows.extend(np.repeat(ids, len(ids)))
        cols.extend(np.tile(ids, len(ids)))
        values.extend(ke.ravel())
    size = stride*(n+1)
    matrix = coo_matrix((values, (rows, cols)), shape=(size, size)).tocsc()
    if kind == "axial":
        matrix[-1, -1] += tip
    force = np.zeros(size)
    force[0] = load
    solution = spsolve(matrix, force)
    residual = np.linalg.norm(matrix @ solution - force) / max(np.linalg.norm(force), 1.)
    if not np.isfinite(solution).all():
        raise ValueError("nonfinite numerical solution")
    return solution, float(residual)


def evaluate(kind, solution, s, length, stiffness, spring):
    stride = 1 if kind == "axial" else 2
    n = len(solution)//stride-1
    h = length/n
    result = []
    for x in stations(s, length):
        e = min(int(x/h), n-1)
        t = (x-e*h)/h
        local = solution[stride*e:stride*(e+2)]
        if kind == "axial":
            u = np.array([1-t, t]) @ local
            result.append([u, -stiffness * (local[1]-local[0])/h, -spring*u])
        else:
            shape, first, second, third = hermite(t, h)
            u = shape @ local
            result.append([u, first @ local, -stiffness*(second @ local),
                           -stiffness*(third @ local), -spring*u])
    return np.array(result)


def convergence(kind, length, stiffness, spring, load, tip=0., ns=(4, 8, 16, 32)):
    fields = ["u", "N", "shaft_reaction"] if kind == "axial" else ["w", "slope", "M", "V", "soil_reaction"]
    results = []
    previous = None
    for n in ns:
        solution, residual = finite_elements(kind, n, length, stiffness, spring, load, tip)
        g, w = leggauss(6)
        s = np.concatenate([(e+(g+1)/2)*length/n for e in range(n)])
        weights = np.tile(w*length/(2*n), n)
        exact = (axial_exact(s, length, stiffness, spring, load, tip) if kind == "axial"
                 else beam_exact(s, length, stiffness, spring, load))
        actual = evaluate(kind, solution, s, length, stiffness, spring)
        denominator = np.sqrt(np.sum(weights[:, None]*exact**2, axis=0))
        if np.any(denominator == 0):
            raise ValueError("zero reference norm; relative error undefined")
        errors = np.sqrt(np.sum(weights[:, None]*(actual-exact)**2, axis=0))/denominator
        orders = (None if previous is None else
                  dict(zip(fields, (np.log(previous[1]/errors)/np.log(n/previous[0])).tolist())))
        results.append(dict(n=n, h_m=length/n, relative_L2=dict(zip(fields, errors.tolist())),
                            observed_order=orders, algebraic_relative_residual=residual,
                            head_displacement_m=float(solution[0]), accepted=None))
        previous = n, errors
    return results


def report():
    # Explicit mathematical verification parameters. None is a fitted 3D coefficient.
    length, ea, ei = 6., 30_000_000., 2_500_000.
    ka, kt, kh = ea/length**2, ea/length, 4*ei/length**4
    return dict(schema="pile-independent-diagnostic/1", accepted=None,
                tolerance="未定義", p6_complete=False, original_3d_comparison="判定不能",
                parameter_role="manufactured verification only; no identified 3D spring mapping",
                verification_parameters=dict(length_m=length, ea_kn=ea, ei_kn_m2=ei,
                                             axial_k_kn_m2=ka, tip_k_kn_m=kt, lateral_k_kn_m2=kh,
                                             axial_head_load_kn=100., lateral_head_load_kn=10.),
                axial=convergence("axial", length, ea, ka, 100., kt),
                finite_beam=convergence("beam", length, ei, kh, 10.),
                semi_infinite_head_m=float(beam_semi_infinite([0.], ei, kh, 10.)[0, 0]),
                finite_beam_head_m=float(beam_exact([0.], length, ei, kh, 10.)[0, 0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = report()
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
