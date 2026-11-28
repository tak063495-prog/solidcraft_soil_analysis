"""Fixed-order observations of two decompositions and transported subspaces."""
import argparse
from decimal import Decimal as D, getcontext
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from b01eig_solve import load_inputs

# All dots accumulate in explicit DOF order. Decimal also resolves sines much
# smaller than sqrt(binary64 epsilon) after the eigensolver has returned.
PREC = 54


def dec(x):
    return D.from_float(x) if isinstance(x, float) else D(x)


def array(x):
    return np.vectorize(dec, otypes=[object])(x)


def dot(a, b):
    return sum((x*y for x, y in zip(a, b, strict=True)), D(0))


def product(a, b):
    return np.array([[dot(row, col) for col in b.T] for row in a], dtype=object)


def norm(a):
    flat = a.ravel(order="C")
    return dot(flat, flat).sqrt()


def orth(a):
    q = a.copy()
    for j in range(q.shape[1]):
        v = q[:, j]
        for i in range(j):
            v -= dot(q[:, i], v)*q[:, i]
        v /= norm(v)
    return q


def angles(a, b):
    qa, qb = orth(a), orth(b)
    diff = qb - product(qa, product(qa.T, qb))
    # Singular values of the projection residual give sin(theta). Thin SVD
    # avoids both 1-cos cancellation and squaring a very small singular value.
    residual = mp.matrix([[str(x) for x in row] for row in diff])
    return [mp.nstr(mp.asin(min(mp.mpf(1), s)), 44)
            for s in mp.svd_r(residual, compute_uv=False)]


def cross(a, b):
    return np.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
                     a[0]*b[1]-a[1]*b[0]], dtype=object)


def physical(c):
    x = array(c["points"])
    center = np.array([sum(col, D(0))/len(x) for col in x.T])
    axes = array(np.eye(3))
    modes = []
    for axis in axes:
        modes.append(np.tile(np.concatenate((axis, 0*axis)), len(x)))
    for axis in axes:
        modes.append(np.concatenate([np.concatenate((cross(axis,p-center),axis)) for p in x]))
    return np.column_stack(modes)


def component_norms(a):
    return {"translation_force_kn": [str(norm(v.reshape(-1,6)[:,:3])) for v in a.T],
            "rotation_moment_kn_m": [str(norm(v.reshape(-1,6)[:,3:])) for v in a.T]}


def drill(c, vector):
    points = array(c["points"])
    u, w = vector.reshape(-1,6)[:,:3], vector.reshape(-1,6)[:,3:]
    normal2, tangent2 = D(0), D(0)
    for ids in c["elements"]:
        a,b,d = (points[i] for i in ids[:3])
        normal = cross(b-a,d-a)
        normal /= norm(normal)
        for i in ids:
            comp = dot(w[i],normal)
            normal2 += comp**2
            tangent2 += norm(w[i]-comp*normal)**2
    return dict(translation_coordinate_l2=str(norm(u)), rotation_coordinate_l2=str(norm(w)),
                incidence_normal_rotation_squared=str(normal2),
                incidence_tangent_rotation_squared=str(tangent2),
                normal_fraction_of_incident_rotation_squared=str(normal2/(normal2+tangent2)))


def transport(v, q):
    return np.column_stack([product(col.reshape(-1,3),q.T).ravel() for col in v.T])


def describe(c,values,vectors,sym,rigid):
    v = vectors[:,:6]
    t = D.from_float(1e-10)*values[-1]
    a = [abs(x) for x in values]
    # Residual of the six-dimensional invariant subspace, independent of
    # arbitrary signs and rotations among the six near-zero eigenvectors.
    q = orth(v)
    hq = product(sym,q)
    invariant = hq-product(q,product(q.T,hq))
    first_positive = vectors[:,6:7]
    r7 = product(sym,first_positive)-values[6]*first_positive
    return dict(eigenvalues=[str(x) for x in values], threshold=str(t),
                eigensystem_residual_frobenius=str(norm(product(sym,vectors)-vectors*np.array(values))),
                orthogonality_defect_frobenius=str(norm(product(vectors.T,vectors)-array(np.eye(len(values))))),
                numerical_zero_count=sum(x<=t for x in values),
                count_abs_le_threshold=sum(x<=t for x in a),
                count_below_negative_threshold=sum(x < -t for x in values),
                gap_6_to_7=str(values[6]-values[5]),
                cluster_radius=str(max(a[:6])), gap_7_to_8=str(values[7]-values[6]),
                effective_condition_quotient=str(values[-1]/values[6]),
                full_algebraic_abs_ratio=str(max(a)/min(a)) if min(a) else None,
                rigid_principal_angles_rad=angles(q,rigid),
                invariant_residual_frobenius=str(norm(invariant)),
                mode7_residual_l2=str(norm(r7)),
                mode7=drill(c,first_positive))


def inspect(c, raw, base):
    k = array(c["stiffness"])
    sym = (k+k.T)/2
    k64 = np.array(c["stiffness"])
    sym64 = array((k64+k64.T)*0.5)
    free = c["free_dofs"]
    reduced = k64[np.ix_(free,free)]
    s = np.linalg.svd(reduced, compute_uv=False)
    r = physical(c)
    physical_res = product(k,r)
    report = dict(id=c["id"], shape=c["shape"], rotation=c["rotation"],
                  norm_K_minus_KT_frobenius=str(norm(k-k.T)),
                  norm_S64_minus_exact_symmetric_frobenius=str(norm(sym64-sym)),
                  constrained_Kff_sigma_max=float(s[0]), constrained_Kff_sigma_min=float(s[-1]),
                  constrained_Kff_condition=float(s[0]/s[-1]),
                  raw_Kr=component_norms(physical_res), symmetric_Kr=component_norms(product(sym,r)),
                  raw_Kr_relative=[str(norm(col)/(norm(k)*norm(mode)))
                                   for col,mode in zip(physical_res.T,r.T,strict=True)])
    basis = {}
    for key in ("binary64","mp"):
        v = array(raw[key+"_vectors"])
        values = [dec(x) for x in raw[key+"_values"]]
        report[key] = describe(c,values,v,sym,r)
        basis[key] = v[:,:7]
        if base:
            u = transport(base[key],array(c["rotation_matrix"]))
            report[key]["rotated_identity_cluster_angles_rad"] = angles(u[:,:6],v[:,:6])
            report[key]["rotated_identity_mode7_angle_rad"] = angles(u[:,6:7],v[:,6:7])
            report[key]["changes_from_identity"] = [str(x-y) for x,y in zip(values,base[key+"_values"],strict=True)]
        basis[key+"_values"] = values
    report["binary64_vs_mp_cluster_angles_rad"] = angles(basis["mp"][:,:6],basis["binary64"][:,:6])
    report["binary64_vs_mp_mode7_angle_rad"] = angles(basis["mp"][:,6:7],basis["binary64"][:,6:7])
    return report,basis


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    getcontext().prec = PREC
    mp.mp.dps = PREC
    cases = sorted(load_inputs(),key=lambda c:(c["shape"], c["rotation"]!="identity",c["rotation"]))
    args.output.mkdir(parents=True,exist_ok=False)
    for c in cases:
        raw = json.loads((args.raw/(c["id"]+".json")).read_text(encoding="utf-8"))
        if c["id"] != raw["id"]:
            raise ValueError("mislabeled eigenpairs")
        identity = c["rotation"]=="identity"
        report, basis = inspect(c, raw, None if identity else base)
        if identity:
            base = basis
        (args.output/(c["id"]+".json")).write_text(
            json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8",newline="\n")
        print(c["id"], flush=True)


if __name__=="__main__":
    main()
