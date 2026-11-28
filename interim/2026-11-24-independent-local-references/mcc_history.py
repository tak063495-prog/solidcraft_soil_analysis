"""Constant-cell-pressure MCC rate equations with axial strain as independent variable."""

import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np


FIELDS = ("p_kpa", "q_kpa", "v", "axial_strain", "lateral_strain", "volumetric_strain")


def reference(spec, pc0, substeps):
    p0, m, k, lam = (spec[key] for key in ("p0_kpa", "M", "kappa", "lambda"))
    shear = 250.0
    v0 = (
        spec["v_lambda"] - lam * math.log(pc0 / spec["p1_kpa"]) + k * math.log(pc0 / p0)
    )
    # On the elastic branch q=3(p-p0); solve f=q²+M²p(p-pc0)=0.
    a, b, c = 9 + m * m, -18 * p0 - m * m * pc0, 9 * p0 * p0
    py = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
    qy = 3 * (py - p0)

    def elastic(q):
        ev = -math.log1p(-k / v0 * math.log1p(q / (3 * p0)))
        return ev / 3 + q / (3 * shear), ev

    ey, evy = elastic(qy)

    def rates(q, ev):
        p = p0 + q / 3
        v = v0 * math.exp(-ev)
        bulk = v * p / k
        pc = p + q * q / (m * m * p)
        # Associated flow dep_v=L*A, dep_s=L*B and logarithmic hardening.
        A, B = m * m * (2 * p - pc), 2 * q
        D = A / 3 + B
        H = m * m * p * pc * v * A / (lam - k)
        compliance = 1 / (9 * bulk) + 1 / (3 * shear)
        determinant = D * D + H * compliance
        if determinant <= 0:
            raise ValueError("singular drained axial tangent")
        dq = H / determinant
        return dq, dq / (3 * bulk) + D * A / determinant

    def advance(q, ev, h):
        k1, l1 = rates(q, ev)
        k2, l2 = rates(q + h * k1 / 2, ev + h * l1 / 2)
        k3, l3 = rates(q + h * k2 / 2, ev + h * l2 / 2)
        k4, l4 = rates(q + h * k3, ev + h * l3)
        return q + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6, ev + h * (
            l1 + 2 * l2 + 2 * l3 + l4
        ) / 6

    rows = []
    previous_axial = 0.0
    for i in range(1, 513):
        axial = spec["terminal_axial_strain"] * i / 512
        if axial <= ey:
            lo, hi = 0.0, qy
            for _ in range(55):
                mid = (lo + hi) / 2
                if elastic(mid)[0] < axial:
                    lo = mid
                else:
                    hi = mid
            q = (lo + hi) / 2
            _, ev = elastic(q)
        else:
            if previous_axial <= ey:
                q, ev = qy, evy
            h = (axial - max(ey, previous_axial)) / substeps
            for _ in range(substeps):
                q, ev = advance(q, ev, h)
        rows.append([p0 + q / 3, q, v0 * math.exp(-ev), axial, (ev - axial) / 2, ev])
        previous_axial = axial
    return np.array(rows)


def differences(candidate, ref):
    # All approved samples exclude the shared initial state and have nonzero reference values.
    if (
        np.any(ref == 0)
        or not np.isfinite(candidate).all()
        or not np.isfinite(ref).all()
    ):
        raise ValueError("non-finite history or zero pointwise denominator")
    errors = np.abs((candidate - ref) / ref)
    return {
        f: dict(relative=float(errors[:, j].max()), step=int(errors[:, j].argmax() + 1))
        for j, f in enumerate(FIELDS)
    }


def assess(spec, fem, archive, baseline):
    limit = spec["endpoint_relative_max_exclusive"]
    references, refinements, cases = {}, {}, []
    for pc in spec["pc_targets_kpa"]:
        levels = [reference(spec, pc, n) for n in (8, 16, 32)]
        refs = [differences(a, b) for a, b in zip(levels, levels[1:])]
        refinements[pc] = refs
        references[pc] = levels[-1]
    expected = {(pc, n, h) for pc in references for n in (12, 96) for h in (1e-4, 5e-5)}
    records = []
    old_cases = {
        (r["initial_pc_kpa"], r["elements"], r["material_substep"]): r for r in baseline
    }
    if set(old_cases) != expected or len(baseline) != len(expected):
        raise ValueError("archived FEM case inventory mismatch")
    for run in fem:
        pc = run["initial_pc_kpa"]
        expected.remove((pc, run["elements"], run["material_substep"]))
        history = run["history"]
        if (
            len(history) != 512
            or run["constant_shear_kpa"] != 250
            or any(
                row["stage"] != "shear" or row["fraction"] != i / 512
                for i, row in enumerate(history, 1)
            )
        ):
            raise ValueError("incomplete drained history")
        old = old_cases[pc, run["elements"], run["material_substep"]]["history"]
        if len(old) != len(history) or any(
            o["fraction"] != r["fraction"] for o, r in zip(old, history)
        ):
            raise ValueError("archived FEM sample coordinate mismatch")
        repeat_delta = {
            key: max(abs(o[key] - r[key]) for o, r in zip(old, history))
            for key in ("p_kpa", "q_kpa", "v")
        }
        vals = np.array([[row[f] for f in FIELDS] for row in history])
        diffs = differences(vals, references[pc])
        ref_ok = all(x["relative"] < limit / 100 for x in refinements[pc][-1].values())
        cases.append(
            dict(
                pc_kpa=pc,
                elements=run["elements"],
                material_substep=run["material_substep"],
                archived_fem_max_abs_difference=repeat_delta,
                maxima=diffs,
                passed=ref_ok and all(x["relative"] < limit for x in diffs.values()),
            )
        )
        records.append(
            dict(
                pc_kpa=pc,
                elements=run["elements"],
                material_substep=run["material_substep"],
                fem=vals.tolist(),
            )
        )
    if expected:
        raise ValueError("missing FEM cases")
    with gzip.GzipFile(str(archive), "wb", mtime=0) as out:
        out.write(
            json.dumps(
                dict(
                    fields=FIELDS,
                    fine_reference={k: v.tolist() for k, v in references.items()},
                    histories=records,
                ),
                allow_nan=False,
                separators=(",", ":"),
            ).encode()
        )
    return dict(
        limit_exclusive=limit,
        substeps_per_sample=[8, 16, 32],
        samples=512,
        refinement=refinements,
        cases=cases,
        passed=all(c["passed"] for c in cases),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("criteria", "fem", "baseline", "archive", "summary"):
        p.add_argument(arg, type=Path)
    args = p.parse_args()
    result = assess(
        json.loads(args.criteria.read_bytes())["mcc"],
        json.loads(args.fem.read_bytes()),
        args.archive,
        json.loads(args.baseline.read_bytes()),
    )
    args.summary.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
