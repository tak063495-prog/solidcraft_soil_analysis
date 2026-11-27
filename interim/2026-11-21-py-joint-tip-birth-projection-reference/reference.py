"""Dry A-02 surface reference; imports no project solver. See README.ja.md.

Positive normal jump/traction is opening/tension. Units: m, kN, kPa.
The partial-domain projection is a diagnostic, never an authorization to BIRTH.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.linalg import lstsq, solve_banded

BASE = "f377e8eda41777c84c3bb849af556bf61683b4cb"
KN, KT, MU = 10000., 5000., .5
FORCE_ATOL = 1e-10  # Existing limited projection contract, NOT extended acceptance.
XY = np.array([[0, 0], [1, 0], [0, 1], [.5, 0], [.5, .5], [0, .5]])
COORDS = np.column_stack((np.vstack((XY, XY[1:])), np.zeros(11)))


def shape(xy):
    x, y = np.asarray(xy).T
    a = 1-x-y
    return np.array([a*(2*a-1), x*(2*x-1), y*(2*y-1),
                     4*a*x, 4*x*y, 4*y*a]).T


def quadrature(n=4, cut=.25):
    """s=x+y in [cut,1], v=x/s in [0,1], dA=s ds dv."""
    if type(n) is not int or not 1 <= n <= 32 or not 0 <= cut < 1:
        raise ValueError("bounded dry surface quadrature required")
    z, w = leggauss(n)
    s = cut+(1-cut)*(z+1)/2
    v = (z+1)/2
    ss, vv = np.meshgrid(s, v, indexing="ij")
    weights = np.outer(w*(1-cut)/2, w/2)*ss
    return np.column_stack((ss.ravel()*vv.ravel(), ss.ravel()*(1-vv.ravel()))), weights.ravel()


def assembly(xy, weights):
    """All 11 physical nodes; shared tip is ONE node with exactly zero row."""
    n = shape(xy)
    b = np.zeros((11, len(weights)))
    b[1:6] = -(n[:, 1:]*weights[:, None]).T
    b[6:] = (n[:, 1:]*weights[:, None]).T
    return b


def mul(a, b):
    out = {}
    for (i, j), u in a.items():
        for (k, l), v in b.items():
            out[i+k, j+l] = out.get((i+k, j+l), F(0))+u*v
    return out


# Polynomial coefficients independently expanded from barycentric Tri6 shapes.
POLYS = [
    {(0, 0): F(1), (1, 0): F(-3), (0, 1): F(-3),
     (2, 0): F(2), (1, 1): F(4), (0, 2): F(2)},
    {(2, 0): F(2), (1, 0): F(-1)},
    {(0, 2): F(2), (0, 1): F(-1)},
    {(1, 0): F(4), (2, 0): F(-4), (1, 1): F(-4)},
    {(1, 1): F(4)},
    {(0, 1): F(4), (0, 2): F(-4), (1, 1): F(-4)},
]
G = {(1, 0): F(3), (0, 1): F(3), (2, 0): F(-2),
     (1, 1): F(-4), (0, 2): F(-2)}  # 1-N_tip, in span of unshared shapes


def integral(poly, cut=F(1, 4)):
    # Integral x^i y^j on s>=cut: beta(i+1,j+1)*(1-cut^(i+j+2))/(i+j+2).
    return sum((c*F(math.factorial(i)*math.factorial(j), math.factorial(i+j+1))
                *(1-cut**(i+j+2))/F(i+j+2) for (i, j), c in poly.items()), F(0))


def exact_nodal(vector, cut=F(1, 4)):
    coeff = np.array([float(integral(mul(p, G), cut)) for p in POLYS[1:]])
    return np.vstack((np.zeros((1, 3)), -coeff[:, None]*vector, coeff[:, None]*vector))


def traction(xy, vector):
    s = np.sum(xy, axis=1)
    return (3*s-2*s*s)[:, None]*np.asarray(vector)


def project(xy, weights, nodal, ids=None):
    """Unique minimum area-weighted traction norm, equality constraints.

    Solve C z=f with C=B/sqrt(w), z=sqrt(w)*t. GELSD SVD is
    independent of the production projection implementation. Canonical GP IDs
    make repeated/permuted calls byte-identical on this runtime.
    """
    xy, weights, nodal = map(np.asarray, (xy, weights, nodal))
    ids = np.arange(len(weights)) if ids is None else np.asarray(ids)
    if (xy.shape != (len(weights), 2) or nodal.shape != (11, 3)
            or len(set(ids.tolist())) != len(weights) or len(ids) != len(weights)
            or np.any(weights <= 0) or not all(np.isfinite(x).all() for x in (xy, weights, nodal))):
        raise ValueError("invalid or nonfinite projection input")
    order = np.argsort(ids)
    b = assembly(xy[order], weights[order])
    root = np.sqrt(weights[order])
    z, _, rank, _ = lstsq(b/root, nodal, cond=None, lapack_driver="gelsd")
    t = z/root[:, None]
    residual = b@t-nodal
    feasible = bool(np.max(np.abs(residual)) <= FORCE_ATOL)
    admissible = bool(np.all(t[:, 2] <= 0) and
                      np.all(np.linalg.norm(t[:, :2], axis=1) <= MU*(-t[:, 2])))
    return dict(traction=t[np.argsort(order)], rank=int(rank), residual=residual,
                reproducible=feasible, admissible=admissible,
                outcome="representable_dry" if feasible and admissible else "rejected")


def balances(xy, weights, t, target):
    f = assembly(xy, weights)@t
    # Arbitrary nodal virtual field: includes affine and non-affine modes.
    du = np.arange(33).reshape(11, 3)*1e-5
    jump = assembly(xy, np.ones(len(weights))).T@du
    work_gp = float(np.sum(weights*np.sum(t*jump, axis=1)))
    return dict(nodal_max_kN=float(np.max(np.abs(f-target))),
                action_reaction_kN=np.sum(f, axis=0).tolist(),
                moment_kNm=np.sum(np.cross(COORDS, f), axis=0).tolist(),
                side_a_moment_defect_kNm=np.sum(np.cross(COORDS[:6], (f-target)[:6]), axis=0).tolist(),
                virtual_work_kNm=work_gp, target_virtual_work_kNm=float(np.sum(target*du)),
                virtual_work_defect_kNm=work_gp-float(np.sum(target*du)),
                assembly_virtual_work_defect_kNm=work_gp-float(np.sum(f*du)))


def clip(t):
    out = t.copy()
    out[:, 2] = np.minimum(out[:, 2], 0)
    r = np.linalg.norm(out[:, :2], axis=1)
    limit = MU*(-out[:, 2])
    out[:, :2] *= np.minimum(1, np.divide(limit, r, out=np.ones_like(r), where=r != 0))[:, None]
    return out


def projection_cases():
    xy, w = quadrature()
    cases = {}
    for name, vector in dict(compression=[0, 0, -10], inside=[2, 1, -10],
                             tension=[2, 1, 3], outside=[8, 6, -10]).items():
        target = exact_nodal(vector)
        p = project(xy, w, target)
        t = p.pop("traction")
        p.pop("residual")
        reverse = project(xy[::-1], w[::-1], target, np.arange(len(w))[::-1])["traction"][::-1]
        repeat = project(xy, w, target)["traction"]
        corrections = {}
        for label, corrected in [("clip", clip(t)), ("smooth_then_clip", clip(.5*t+.5*np.average(t, axis=0, weights=w)))]:
            jump0 = t/np.array([KT, KT, KN])
            corrections[label] = dict(**balances(xy, w, corrected, target),
                fixed_jump_work_change_kNm=float(np.sum(w*np.sum((corrected-t)*jump0, axis=1))),
                elastic_energy_change_kNm=float(.5*np.sum(w*np.sum((corrected**2-t**2)/[KT, KT, KN], axis=1)) ))
        cases[name] = dict(**p, **balances(xy, w, t, target),
            analytic_traction_max_kpa=float(np.max(np.abs(t-traction(xy, vector)))),
            repeat_exact=bool(np.array_equal(t, repeat)), order_exact=bool(np.array_equal(t, reverse)),
            corrections=corrections, accepted=None,
            scope="partial activity: diagnostic only; BIRTH remains rejected")
    impossible = exact_nodal([0, 0, -10])
    impossible[0, 2] = 1  # shared tip cannot carry a difference force
    p = project(xy, w, impossible)
    cases["unrepresentable_tip"] = dict(outcome=p["outcome"],
        residual_max_kN=float(np.max(np.abs(p["residual"]))), accepted=None)
    return cases


def birth_ledger():
    """Closed-form scalar spring/body-force bookkeeping, NOT a BIRTH solver."""
    # Nodes 0,1 lower bar; 2,3 upper bar. Joint 1/2. P=10kN.
    bulk = np.array([-10., 10., -10., 10.])
    contact = np.array([0., -10., 10., 0.])
    external = bulk+contact
    u = np.array([0., -.001, -.002, -.003])
    # Joint-only death keeps BOTH owners. Retained load balances the removed joint.
    retained = -contact
    # B reference starts with zero traction/slip/work, while owners are still present.
    joint_birth_residual = bulk-external-retained
    # Body birth additionally adds upper internal force and a 2kN body load.
    new_internal = np.array([0., 0., -10., 10.])
    new_body_load = np.array([0., 0., -1., -1.])
    delta = new_internal-new_body_load
    return dict(model="two elastic bars, scalar axial DOFs; P=10kN, upper weight=2kN",
        node_ids=[0, 1, 2, 3], virtual_displacement_m=u.tolist(),
        joint_only=dict(owner_active=[True, True], bulk_internal_kN=bulk.tolist(),
            old_contact_kN=contact.tolist(), external_kN=external.tolist(),
            retained_kN=retained.tolist(), death_residual_kN=joint_birth_residual.tolist(),
            rebirth_B_preparation_residual_kN=joint_birth_residual.tolist(),
            frozen_generation=1, new_generation=2, new_slip_m=0., new_work_kNm=0.,
            birth_jump_m=-.001, birth_reference_m=-.001),
        crossing_body_birth=dict(owner_active_before=[True, False], owner_active_after=[True, True],
            added_internal_kN=new_internal.tolist(), added_body_load_kN=new_body_load.tolist(),
            joint_only_retained_load_reuse_defect_kN=delta.tolist(),
            virtual_work_defect_kNm=float(delta@u), generation_policy="仕様未定義", accepted=None),
        partial_face_rebirth=dict(active_face_identity_policy="仕様未定義", generation_policy="仕様未定義", accepted=None),
        note="Preparation bookkeeping only. BIRTH equilibrium path, real energy transfer and material generations are not specified or accepted.")


def numerical_series(levels=(2, 4, 8, 16)):
    """Separate midpoint GP integration and P1 normal-gap approximation.

    Analytic smooth input t=(.2,.1,-1)*exp(s). Independent adaptive 1D
    integrals supply targets. P1 gap follows a separate weak penalty solve.
    No fitting to solver output, Richardson correction or acceptance threshold.
    """
    a = .25
    vector = np.array([.2, .1, -1.])
    # Analytic v-integrals of each unshared shape.
    ns = [lambda s: 2*s*s/3-s/2, lambda s: 2*s*s/3-s/2,
          lambda s: 2*s*(1-s), lambda s: 2*s*s/3, lambda s: 2*s*(1-s)]
    coeff = np.array([quad(lambda s: n(s)*math.exp(s)*s, a, 1, epsabs=1e-13, epsrel=1e-13)[0] for n in ns])
    target = np.vstack((np.zeros((1, 3)), -coeff[:, None]*vector, coeff[:, None]*vector))
    norm2 = quad(lambda s: s*math.exp(2*s), a, 1, epsabs=1e-13)[0]
    gp_rows, h_rows = [], []
    for n in levels:
        s = a+(np.arange(n)+.5)*(1-a)/n
        v = (np.arange(n)+.5)/n
        ss, vv = np.meshgrid(s, v, indexing="ij")
        xy = np.column_stack((ss.ravel()*vv.ravel(), ss.ravel()*(1-vv.ravel())))
        w = ss.ravel()*(1-a)/n**2
        f = assembly(xy, w)@(np.exp(ss.ravel())[:, None]*vector)
        gp_rows.append(dict(n=n, gp_count=n*n, surface_elements=1,
            force_l2_kN=float(np.linalg.norm(f-target)), force_max_kN=float(np.max(np.abs(f-target))),
            side_moment_error_kNm=float(np.linalg.norm(np.sum(np.cross(COORDS[:6], (f-target)[:6]), axis=0)))))
        grid = np.linspace(a, 1, n+1)
        # Independent P1 weak penalty equilibrium: integral KN*g_h*v = integral q*v.
        # The surface loading is q=-exp(s); no bulk or production solver is used.
        z, wz = leggauss(8)  # fixed GP rule for ALL h levels
        ss = (grid[:-1, None]+grid[1:, None])/2+np.diff(grid)[:, None]*z/2
        ww = np.diff(grid)[:, None]*wz/2
        basis = np.array([(1-z)/2, (1+z)/2])
        mass = np.zeros((n+1, n+1))
        rhs = np.zeros(n+1)
        for e in range(n):
            mass[e:e+2, e:e+2] += (basis*(ww[e]*ss[e]))@basis.T
            rhs[e:e+2] += basis@(ww[e]*ss[e]*np.exp(ss[e]))
        band = np.zeros((3, n+1))
        band[0, 1:] = np.diag(mass, 1)
        band[1] = np.diag(mass)
        band[2, :-1] = np.diag(mass, -1)
        nodal_pressure = solve_banded((1, 1), band, rhs)
        ih = nodal_pressure[:-1, None]*basis[0]+nodal_pressure[1:, None]*basis[1]
        err = np.sqrt(np.sum(ww*ss*(ih-np.exp(ss))**2)/norm2)
        c = np.array([np.sum(ww*ss*ih*fn(ss)) for fn in ns])
        fh = np.vstack((np.zeros((1, 3)), -c[:, None]*vector, c[:, None]*vector))
        h_rows.append(dict(n=n, strip_elements=n, gp_per_element=8,
            traction_relative_l2=float(err), gap_relative_l2=float(err),
            force_l2_kN=float(np.linalg.norm(fh-target)),
            weak_force_residual_kN=float(np.max(np.abs(mass@nodal_pressure-rhs))),
            work_balance_defect_kNm=float(abs(nodal_pressure@mass@nodal_pressure-nodal_pressure@rhs)/KN),
            minimum_pressure_kpa=float(np.min(nodal_pressure))))
    for rows, key in [(gp_rows, "force_l2_kN"), (h_rows, "traction_relative_l2")]:
        for i, row in enumerate(rows):
            row["observed_order"] = None if i == 0 else math.log(rows[i-1][key]/row[key])/math.log(row["n"]/rows[i-1]["n"])
    return dict(input="dry smooth penalty contact patch, exp(x+y); no bulk solver",
                gp_series=gp_rows, element_series=h_rows, accepted=None,
                tolerance="未定義", minimum_order="未定義")


def report(levels=(2, 4, 8, 16)):
    return dict(base_commit=BASE, accepted=None, A_02="未達", NC_01="未達", P6="未達",
                projection=projection_cases(), birth=birth_ledger(), numerical=numerical_series(levels))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--pc-series", action="store_true")
    args = p.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = report((16, 32, 64, 128) if args.pc_series else (2, 4, 8, 16))
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
