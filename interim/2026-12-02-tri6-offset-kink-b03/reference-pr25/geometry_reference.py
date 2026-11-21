"""B-03 independent Euclidean geometry/rigid-arm diagnostics; no solver imports.

Coordinates in m; rotations in rad. This module never decides B-03 acceptance.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import lstsq, svdvals
from scipy.spatial.transform import Rotation

META = {
    "schema": "b03-geometry-v1",
    "length_unit": "m",
    "coordinate_basis": "global_cartesian_right_handed",
    "normal_sign": "oriented_face_winding",
    "offset_sign": "attachment=neutral+e*normal",
}
VECTORS = ("neutral_m", "normal", "attachment_m")
SCALARS = ("e_m", "t_m")
FIELDS = [*META, "shared_node_id", "face_id", *SCALARS,
          *(f"{v}_{k}" for v in VECTORS for k in "xyz")]


def skew(v):
    x, y, z = np.asarray(v, dtype=float)
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def arm(r):
    """[u(r), rotation] = B(r)[u(0), rotation]."""
    b = np.eye(6)
    b[:3, 3:] = -skew(r)
    return b


def mismatch(x1, n1, e1, x2, n2, e2):
    """Oriented point separation; not a finite-solid penetration classifier."""
    x1, x2, n1, n2 = map(lambda v: np.asarray(v, float), (x1, x2, n1, n2))
    p1, p2 = x1 + e1*n1, x2 + e2*n2
    d = p2 - p1
    return {"attachment1_m": p1.tolist(), "attachment2_m": p2.tolist(),
            "delta_m": d.tolist(), "norm_m": float(np.linalg.norm(d)),
            "delta_dot_n1_m": float(d@n1), "delta_dot_n2_m": float(d@n2),
            "tangent_to_face1_m": (d-(d@n1)*n1).tolist(),
            "signed_gap_to_face1_m": float(d@n1),
            "signed_gap_to_face2_m": float(-d@n2)}


def kink_closed(theta, e1, e2):
    # n1=(0,0,1), n2=(sin(theta),0,cos(theta)); common neutral origin.
    return np.array([e2*np.sin(theta), 0., e2*np.cos(theta)-e1])


def step_closed(h, s, e1, e2):
    # n1=n2=ez, x2-x1=(s,0,h).
    return np.array([s, 0., h+e2-e1])


def parallel_layer_gap(h, t1, t2):
    """Signed separation of normal thickness intervals, with tangent overlap assumed.

    Positive: gap, zero: touch, negative: overlap depth of normal intervals.
    This does not establish overlap of finite 3D plates.
    """
    return max(-t1/2, h-t2/2) - min(t1/2, h+t2/2)


def numerical_rank(a):
    # Floating-point diagnostic only (not an engineering acceptance tolerance).
    s = svdvals(a)
    threshold = 128*np.finfo(float).eps*max(1., float(s[0]) if len(s) else 0.)
    return int(np.count_nonzero(s > threshold)), threshold


def candidate_points(x, normals, offsets):
    x, normals, offsets = np.asarray(x, float), np.asarray(normals, float), np.asarray(offsets, float)
    p = x + offsets[:, None]*normals
    result = {"independent_rigid_arms": p.copy()}
    b = normals.sum(axis=0)
    if np.linalg.norm(b) > 128*np.finfo(float).eps:
        b /= np.linalg.norm(b)
        c = x.mean(axis=0) + offsets.mean()*b
        result["average_normal_common_point"] = np.tile(c, (2, 1))
    else:
        result["average_normal_common_point"] = None
    a = normals[1]-normals[0]
    if np.linalg.norm(a) > 128*np.finfo(float).eps:
        a /= np.linalg.norm(a)
        m = x.mean(axis=0)
        # Explicit choice of one angle-bisector plane: a.(r-m)=0.
        result["orthogonal_bisector_projection"] = p - ((p-m)@a)[:, None]*a
    else:
        result["orthogonal_bisector_projection"] = None
    return p, result


def diagnose_candidate(p, q):
    if q is None:
        return {"defined": False, "reason": "bisector/average normal is not unique",
                "accepted": None}
    bp = np.vstack([arm(r) for r in p])
    bq = np.vstack([arm(r) for r in q])
    # Eliminate a 6-DOF common rigid body. B(r)^-1 = B(-r).
    c = np.hstack([-arm(q[1])@arm(-q[0]), np.eye(6)])
    represented = c@bq
    physical = c@bp
    rank, threshold = numerical_rank(physical)
    qrank, _ = numerical_rank(bq)
    crank, _ = numerical_rank(c)
    # Six basis virtual motions and twelve basis loads: all bilinear work pairs.
    # f_common=Bq.T f. Reference physical work is f.T Bp dq.
    work = bq.T - bp.T
    # Independent check: resultant and moment from cross products.
    f = np.arange(1., 13.).reshape(2, 6)/13
    expected = np.r_[f[:, :3].sum(axis=0),
                     sum((np.cross(q[i], f[i, :3])+f[i, 3:] for i in range(2)))]
    conjugate = bq.T@f.ravel()
    cross_map = np.empty((6, 12))
    for j, basis_load in enumerate(np.eye(12)):
        loads = basis_load.reshape(2, 6)
        cross_map[:, j] = np.r_[loads[:, :3].sum(axis=0),
            sum((np.cross(q[i], loads[i, :3])+loads[i, 3:] for i in range(2)))]
    # Solve the constraints independently; residual for six physical rigid motions.
    recovered = lstsq(bq, bp, lapack_driver="gelsd")[0]
    return {"defined": True, "accepted": None,
            "geometry_error_vectors_m": (q-p).tolist(),
            "max_geometry_error_m": float(np.max(np.linalg.norm(q-p, axis=1))),
            "represented_rigid_map_rank": qrank, "constraint_nullity": 12-crank,
            "represented_six_mode_max_residual": float(np.max(np.abs(represented))),
            "original_six_mode_constraint_matrix": physical.tolist(),
            "original_six_mode_max_constraint_residual": float(np.max(np.abs(physical))),
            "original_rigid_subspace_dimension": 6-rank,
            "rank_roundoff_threshold": threshold,
            "original_field_fit_max_residual": float(np.max(np.abs(bq@recovered-bp))),
            "represented_work_cross_product_residual": float(max(np.max(np.abs(conjugate-expected)), np.max(np.abs(bq.T-cross_map)))),
            "original_work_defect_matrix": work.tolist(),
            "original_work_defect_max_coefficient": float(np.max(np.abs(work))),
            "common_force_moment_example": conjugate.tolist()}


def orientation_diagnostics(n, e, t, expected_n=None, expected_e=None):
    errors = []
    n = np.asarray(n, float)
    if n.shape != (3,) or not np.isfinite(n).all():
        return ["nonfinite_or_invalid_normal"]
    if not np.isfinite([e, t]).all() or t <= 0:
        errors.append("nonfinite_offset_or_nonpositive_thickness")
    elif abs(e) > t/2:
        errors.append("offset_outside_thickness")
    # Unit check uses a roundoff bound solely to validate a numeric unit-vector representation.
    if abs(float(n@n)-1.) > 128*np.finfo(float).eps:
        errors.append("normal_not_unit")
    if expected_n is None:
        errors.append("normal_orientation_unconfirmed")
    else:
        ref = np.asarray(expected_n, float)
        if n@ref < 0:
            errors.append("normal_orientation_reversed")
        if not np.array_equal(n, ref):
            errors.append("normal_differs_from_reference")
    if expected_e is not None and e != expected_e:
        errors.append("signed_offset_differs_from_reference")
        if e*expected_e < 0:
            errors.append("offset_sign_reversed")
    return errors


def validate_document(doc):
    if not isinstance(doc, dict) or set(doc) != {"metadata", "rows"}:
        raise ValueError("document keys must be metadata, rows")
    if doc["metadata"] != META:
        raise ValueError("unit/basis/normal/offset/schema metadata mismatch (no conversion)")
    if not isinstance(doc["rows"], list) or not doc["rows"]:
        raise ValueError("nonempty rows required")
    ids = set()
    faces = {}
    for row in doc["rows"]:
        if set(row) != {"shared_node_id", "face_id", *VECTORS, *SCALARS}:
            raise ValueError("missing/extra row fields")
        key = tuple(row[k] for k in ("shared_node_id", "face_id"))
        if not all(type(v) is str and v.strip() == v and v for v in key):
            raise ValueError("nonempty unambiguous string IDs required")
        if key in ids:
            raise ValueError("duplicate node/face ID")
        ids.add(key)
        for field in VECTORS:
            values = row[field]
            if (not isinstance(values, list) or len(values) != 3 or
                any(type(v) not in (int, float) for v in values) or
                not np.isfinite(values).all()):
                raise ValueError("invalid/nonfinite vector: " + field)
        for field in SCALARS:
            if type(row[field]) not in (int, float) or not np.isfinite(row[field]):
                raise ValueError("invalid/nonfinite scalar: " + field)
        errs = orientation_diagnostics(row["normal"], row["e_m"], row["t_m"], row["normal"])
        if errs:
            raise ValueError(",".join(errs))
        constants = (tuple(row["normal"]), row["e_m"], row["t_m"])
        if row["face_id"] in faces and faces[row["face_id"]] != constants:
            raise ValueError("face normal/e/t must be constant in this piecewise planar contract")
        faces[row["face_id"]] = constants
    return ids


def read_document(path):
    path = Path(path)
    if path.suffix.lower() == ".json":
        def no_duplicates(pairs):
            out = {}
            for k, v in pairs:
                if k in out:
                    raise ValueError("duplicate JSON key: " + k)
                out[k] = v
            return out
        doc = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=no_duplicates)
    elif path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames != FIELDS:
                raise ValueError("CSV headers/order mismatch")
            rows = []
            for row in reader:
                if set(row) != set(FIELDS) or any(row[k] != v for k, v in META.items()):
                    raise ValueError("CSV metadata/missing/extra field mismatch")
                rows.append({"shared_node_id": row["shared_node_id"], "face_id": row["face_id"],
                             **{k: float(row[k]) for k in SCALARS},
                             **{v: [float(row[f"{v}_{k}"]) for k in "xyz"] for v in VECTORS}})
        doc = {"metadata": META.copy(), "rows": rows}
    else:
        raise ValueError("only explicit CSV or JSON is supported")
    validate_document(doc)
    return doc


def write_document(doc, path):
    validate_document(doc)
    path = Path(path)
    if path.suffix == ".json":
        path.write_text(json.dumps(doc, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    elif path.suffix == ".csv":
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDS)
            writer.writeheader()
            for row in doc["rows"]:
                writer.writerow({**META, **{k: row[k] for k in ("shared_node_id", "face_id", *SCALARS)},
                                 **{f"{v}_{k}": row[v][i] for v in VECTORS for i, k in enumerate("xyz")}})
    else:
        raise ValueError("unsupported output suffix")


def compare_documents(reference, actual):
    expected_ids, actual_ids = validate_document(reference), validate_document(actual)
    if expected_ids != actual_ids:
        raise ValueError("missing/extra node/face ID")
    index = {(r["shared_node_id"], r["face_id"]): r for r in actual["rows"]}
    rows = []
    for ref in reference["rows"]:
        key = (ref["shared_node_id"], ref["face_id"])
        row = index[key]
        errors = orientation_diagnostics(row["normal"], row["e_m"], row["t_m"], ref["normal"], ref["e_m"])
        for k in (*VECTORS, *SCALARS):
            if not np.array_equal(row[k], ref[k]):
                errors.append("reference_mismatch:"+k)
        residual = np.asarray(row["attachment_m"])-(np.asarray(row["neutral_m"])+row["e_m"]*np.asarray(row["normal"]))
        # Undefined engineering tolerance -> no numerical discrepancy is accepted.
        if np.any(residual != 0):
            errors.append("attachment_identity_mismatch")
        rows.append({"shared_node_id": key[0], "face_id": key[1], "errors": errors,
                     "attachment_identity_residual_m": residual.tolist()})
    return {"accepted": None, "tolerance": "未定義", "contract_consistent": not any(r["errors"] for r in rows),
            "rows": rows, "production_acceptance": "B-03未達; synthetic roundtrip is not body output"}


def series():
    cases = []
    for t in (0.02, 0.1, 1.):
        for deg in (0., 15., 45., 90., 135., 179.):
            theta = np.deg2rad(deg)
            rotated_edges = Rotation.from_rotvec([0., theta, 0.]).apply(np.eye(3)[:2])
            numeric_normal = np.cross(*rotated_edges)
            numeric_normal /= np.linalg.norm(numeric_normal)
            n = np.array([[0., 0., 1.], numeric_normal])
            for a, b in ((-0.5, -0.5), (0., 0.), (0.25, 0.25), (0.5, 0.5), (-0.5, 0.5), (0.25, -0.5), (0.6, 0.6), (-0.6, -0.6)):
                e = np.array([a*t, b*t])
                p, candidates = candidate_points(np.zeros((2, 3)), n, e)
                calc = mismatch([0, 0, 0], n[0], e[0], [0, 0, 0], n[1], e[1])
                closed = kink_closed(theta, *e)
                # Stable equivalent norm formula, including almost flat equal offsets.
                if e[0]*e[1] >= 0:
                    length = np.sqrt((e[1]-e[0])**2 + 4*e[0]*e[1]*np.sin(theta/2)**2)
                else:
                    length = np.sqrt(e@e-2*e[0]*e[1]*np.cos(theta))
                cases.append({"kind": "kink", "t_m": t, "theta_deg": deg, "e_m": e.tolist(),
                              **calc, "closed_delta_m": closed.tolist(), "closed_norm_m": float(length),
                              "closed_vector_error_m": float(np.max(np.abs(np.asarray(calc["delta_m"])-closed))),
                              "closed_norm_error_m": abs(calc["norm_m"]-float(length)),
                              "range_rejected": bool(np.any(np.abs(e)>t/2)),
                              "candidates": {k: diagnose_candidate(p, q) for k, q in candidates.items()}, "accepted": None})
        for h in (-t, -t/2, 0., t/2, t):
            for s in (0., t/3):
                for a, b in ((-0.5, -0.5), (0., 0.5), (0.5, -0.5), (0.6, 0.6), (-0.6, -0.6)):
                    e = np.array([a*t, b*t])
                    x = np.array([[0., 0., 0.], [s, 0., h]])
                    n = np.array([[0., 0., 1.], [0., 0., 1.]])
                    p, candidates = candidate_points(x, n, e)
                    calc = mismatch(x[0], n[0], e[0], x[1], n[1], e[1])
                    closed = step_closed(h, s, *e)
                    cases.append({"kind": "step", "t_m": t, "h_m": h, "s_m": s, "e_m": e.tolist(),
                                  **calc, "closed_delta_m": closed.tolist(),
                                  "closed_vector_error_m": float(np.max(np.abs(np.asarray(calc["delta_m"])-closed))),
                                  "closed_norm_m": float(np.hypot(s, h+e[1]-e[0])),
                                  "closed_norm_error_m": abs(calc["norm_m"]-float(np.hypot(s, h+e[1]-e[0]))),
                                  "normal_layer_gap_m": parallel_layer_gap(h, t, t),
                                  "range_rejected": bool(np.any(np.abs(e)>t/2)),
                                  "candidates": {k: diagnose_candidate(p, q) for k, q in candidates.items()}, "accepted": None})
    return {"accepted": None, "tolerance": "未定義", "cases": cases,
            "case_count": len(cases), "scope": "Euclidean geometry and infinitesimal rigid-arm algebra only"}


def synthetic_document():
    normals = [[0., 0., 1.], [1., 0., 0.]]
    return {"metadata": META.copy(), "rows": [
        {"shared_node_id": "corner-1", "face_id": f"face-{i+1}", "neutral_m": [0., 0., 0.],
         "normal": n, "e_m": .05, "t_m": .1, "attachment_m": (.05*np.array(n)).tolist()}
        for i, n in enumerate(normals)]}


def dump(obj, path):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("series")
    s.add_argument("output")
    s = sub.add_parser("synthetic")
    s.add_argument("output")
    s = sub.add_parser("compare")
    s.add_argument("reference")
    s.add_argument("actual")
    s.add_argument("output")
    args = parser.parse_args()
    try:
        if args.command == "series":
            dump(series(), args.output)
            return 2
        if args.command == "synthetic":
            write_document(synthetic_document(), args.output)
            return 0
        report = compare_documents(read_document(args.reference), read_document(args.actual))
        dump(report, args.output)
        return 2 if report["contract_consistent"] else 1
    except (ValueError, TypeError, KeyError, OSError, OverflowError) as exc:
        dump({"accepted": None, "contract_consistent": False, "error": str(exc)}, args.output)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
