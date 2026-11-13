"""Strict proposed interchange for future dry A-02 exports, not a solver adapter.

JSON = {metadata, rows}; CSV has identical scalar columns + .meta.json sidecar.
Two independent runs of the same physical input, one with reversed GP order,
are mandatory. A clean report exits 2 (undefined acceptance); rejection exits 1.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from reference import (BASE, KN, KT, FORCE_ATOL, COORDS, assembly, balances,
                       birth_ledger, quadrature, traction)

META = dict(schema="dry-joint-projection-diagnostic/1", units="m-kN-kPa",
            normal="a_to_b", jump="u_b-u_a", tension_positive=True,
            water=False, complete=True, accepted=None,
            protocol="提案・未承認", source_kind="synthetic_reference")
FIELDS = ["run", "frame", "kind", "id", "generation", "active", "x", "y", "z",
          "weight", "tx", "ty", "tn", "slip_x", "slip_y", "slip_z",
          "jump_x", "jump_y", "jump_z", "reference_x", "reference_y", "reference_z",
          "work", "energy", "fx", "fy", "fz", "status"]
VECTORS = [("tx", "ty", "tn"), ("slip_x", "slip_y", "slip_z"),
           ("jump_x", "jump_y", "jump_z"), ("reference_x", "reference_y", "reference_z"),
           ("fx", "fy", "fz")]
INPUT = dict(base=BASE, joint_count=1, water=False, cut=0., quadrature=4,
             vector=[2, 1, -10], kn=KN, kt=KT, mu=.5,
             events=["birth_A", "joint_death", "rebirth_B", "elastic_load"],
             owner_activity=[True, True], tip_node=0,
             separate_bar_ledger=dict(axial_force_kN=10., upper_weight_kN=2.,
                                      virtual_displacement_m=[0., -.001, -.002, -.003]))
INPUT_HASH = hashlib.sha256(json.dumps(INPUT, sort_keys=True).encode()).hexdigest()
REFUSALS = ["partial_face_rebirth", "crossing_body_birth", "tension_initial_force",
            "outside_friction_initial_force", "unrepresentable_tip_force", "water_combination"]


def row(**values):
    result = dict.fromkeys(FIELDS, 0.)
    result.update(run="forward", frame=0, kind="contact", id=0, generation=1,
                  active=1, status="diagnostic_only")
    result.update(values)
    return {k: v.item() if isinstance(v, np.generic) else v for k, v in result.items()}


def fixture():
    """Input-only analytic export, solely for testing the comparator."""
    xy, w = quadrature(cut=0.)
    initial = traction(xy, INPUT["vector"])
    # Representable by physical nodal displacements: zero jump at the shared tip.
    s = xy.sum(axis=1)
    refs = (3*s-2*s*s)[:, None]*[1e-12, -2e-12, -1e-12]
    rows = []
    for run in ["forward", "reversed"]:
        for frame in range(4):
            t = initial if frame < 2 else (np.zeros_like(initial) if frame == 2 else .2*initial)
            active = frame != 1
            generation = 1 if frame < 2 else 2
            energy = .5*np.sum(t*t/[KT, KT, KN], axis=1)
            jump = refs+t/[KT, KT, KN]
            for i in (range(len(w)) if run == "forward" else reversed(range(len(w)))):
                rows.append(row(run=run, frame=frame, id=i, generation=generation,
                    active=int(active), x=xy[i, 0], y=xy[i, 1], weight=w[i],
                    tx=t[i, 0], ty=t[i, 1], tn=t[i, 2],
                    jump_x=jump[i, 0], jump_y=jump[i, 1], jump_z=jump[i, 2],
                    reference_x=refs[i, 0], reference_y=refs[i, 1], reference_z=refs[i, 2],
                    energy=energy[i], work=energy[i] if frame == 3 else 0.))
            forces = assembly(xy, w)@t if active else np.zeros((11, 3))
            for i, force in enumerate(forces):
                rows.append(row(run=run, frame=frame, kind="node", id=i,
                    generation=generation, active=int(active), x=COORDS[i, 0],
                    y=COORDS[i, 1], z=COORDS[i, 2], fx=force[0], fy=force[1], fz=force[2]))
        for i, refusal in enumerate(REFUSALS):
            rows.append(row(run=run, frame=4, kind="refusal", id=i, status=refusal+":rejected"))
        # Separate full-node body/contact/external/retained-force reference.
        # Frames 5/6/7 are initial/death/B-preparation of the 4-node bar ledger.
        ledger = birth_ledger()["joint_only"]
        for frame in [5, 6, 7]:
            parts = {"ledger_internal": ledger["bulk_internal_kN"],
                     "ledger_contact": ledger["old_contact_kN"] if frame == 5 else [0.]*4,
                     "ledger_external": ledger["external_kN"],
                     "ledger_retained": [0.]*4 if frame == 5 else ledger["retained_kN"]}
            for kind, values in parts.items():
                for i, force in enumerate(values):
                    rows.append(row(run=run, frame=frame, kind=kind, id=i,
                        generation=1 if frame < 7 else 2, x=[0., 1., 1., 2.][i],
                        active=int(frame != 6) if kind == "ledger_contact" else 1, fx=force))
    return dict(metadata=dict(META, input_sha256=INPUT_HASH), rows=rows)


def unique_object(pairs):
    result = {}
    for k, v in pairs:
        if k in result:
            raise ValueError("duplicate JSON key: "+k)
        result[k] = v
    return result


def load(path):
    if path.suffix.lower() == ".csv":
        meta = json.loads(path.with_suffix(".meta.json").read_text(encoding="utf-8"), object_pairs_hook=unique_object)
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != FIELDS:
                raise ValueError("CSV columns must match the documented contract exactly")
            rows = list(reader)
        for r in rows:
            for k in set(FIELDS)-{"run", "kind", "status"}:
                if k in {"id", "generation", "active", "frame"}:
                    r[k] = int(r[k])
                else:
                    r[k] = float(r[k])
        return dict(metadata=meta, rows=rows)
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_object)


def save(data, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".csv":
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(data["rows"])
        path.with_suffix(".meta.json").write_text(json.dumps(data["metadata"], ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    else:
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def compare(data):
    expected = fixture()
    errors, differences = [], {}
    if (not isinstance(data, dict) or set(data) != {"metadata", "rows"}
            or not isinstance(data["metadata"], dict) or not isinstance(data["rows"], list)):
        raise ValueError("metadata and rows required")
    meta = data["metadata"]
    if set(meta) != set(expected["metadata"]):
        errors.append("metadata fields missing/unknown")
    for k, value in expected["metadata"].items():
        # Future native exports may identify themselves; synthetic exports never become native evidence.
        if k == "source_kind":
            if meta.get(k) not in {"native_export", "synthetic_reference"}:
                errors.append("source_kind invalid")
        elif type(meta.get(k)) is not type(value) or meta.get(k) != value:
            errors.append("metadata mismatch: "+k)
    index = {}
    for r in data["rows"]:
        if not isinstance(r, dict) or set(r) != set(FIELDS):
            raise ValueError("row fields missing/unknown")
        for k in set(FIELDS)-{"run", "kind", "status"}:
            if type(r[k]) not in (int, float) or not np.isfinite(r[k]):
                raise ValueError("nonfinite/non-numeric: "+k)
        for k in ["id", "frame", "generation", "active"]:
            if type(r[k]) is not int:
                raise ValueError("integer required: "+k)
        key = (r["run"], r["frame"], r["kind"], r["id"])
        if key in index:
            raise ValueError("duplicate physical sample: "+str(key))
        index[key] = r
    exp = {(r["run"], r["frame"], r["kind"], r["id"]): r for r in expected["rows"]}
    if set(index) != set(exp):
        errors.append("missing/extra frame, run, GP, node or refusal")
    for key in set(index)&set(exp):
        r, e = index[key], exp[key]
        for k in FIELDS:
            if k in {"run", "kind", "status", "id", "frame", "generation", "active"}:
                if r[k] != e[k]:
                    errors.append(f"lifecycle/scope mismatch {key}:{k}")
                continue
            delta = abs(r[k]-e[k])
            differences[k] = max(differences.get(k, 0.), delta)
            # Force gate is inherited only for the existing full-owner projection.
            # Other values use roundoff to detect mutation, NOT a physical acceptance tolerance.
            limit = FORCE_ATOL if k in {"fx", "fy", "fz"} else 64*np.finfo(float).eps*max(1., abs(e[k]))
            if delta > limit:
                errors.append(f"raw sample changed (clip/smooth/history/field) {key}:{k}")
        if key[0] == "reversed":
            other = index.get(("forward", *key[1:]))
            if other is not None and any(r[k] != other[k] for k in FIELDS if k != "run"):
                errors.append("nonreproducible projection/order dependence: "+str(key))
    audits = []
    xy, weights = quadrature(cut=0.)
    for run in ["forward", "reversed"]:
        for frame in range(4):
            gp_keys = [(run, frame, "contact", i) for i in range(len(weights))]
            node_keys = [(run, frame, "node", i) for i in range(11)]
            if not all(k in index for k in gp_keys+node_keys):
                continue
            gp = [index[k] for k in gp_keys]
            t = np.array([[r[k] for k in VECTORS[0]] for r in gp])
            if frame == 1:
                t *= 0  # frozen states do not assemble forces
            f = np.array([[index[k][a] for a in VECTORS[-1]] for k in node_keys])
            audit = balances(xy, weights, t, f)
            audit.update(run=run, frame=frame)
            audit["elastic_energy_kNm"] = float(sum(w*r["energy"] for w, r in zip(weights, gp)))
            audit["post_birth_work_kNm"] = float(sum(w*r["work"] for w, r in zip(weights, gp)))
            audit["new_generation_energy_work_defect_kNm"] = (
                audit["post_birth_work_kNm"]-audit["elastic_energy_kNm"] if frame >= 2 else None)
            audits.append(audit)
            if audit["nodal_max_kN"] > FORCE_ATOL:
                errors.append(f"GP/all-node assembly mismatch {run}:{frame}")
    ledger_audits = []
    for run in ["forward", "reversed"]:
        for frame in [5, 6, 7]:
            kinds = ["ledger_internal", "ledger_contact", "ledger_external", "ledger_retained"]
            if not all((run, frame, kind, i) in index for kind in kinds for i in range(4)):
                continue
            values = [np.array([[index[(run, frame, kind, i)][k] for k in VECTORS[-1]]
                                for i in range(4)]) for kind in kinds]
            residual = values[0]+values[1]-values[2]-values[3]
            ledger_audits.append(dict(run=run, frame=frame, all_node_residual_kN=residual.tolist(),
                virtual_work_defect_kNm=float(residual[:, 0]@INPUT["separate_bar_ledger"]["virtual_displacement_m"])))
            if np.max(np.abs(residual)) > FORCE_ATOL:
                errors.append(f"full-node body/contact/external/retained imbalance {run}:{frame}")
    return dict(accepted=None, A_02="未達", tolerance="拡張範囲は未定義",
                contract_valid=not errors, errors=sorted(set(errors)), max_absolute_differences=differences,
                force_and_virtual_work_audits=audits, full_node_ledger_audits=ledger_audits,
                source_kind=meta.get("source_kind"),
                specification="提案・未承認", exit_code=1 if errors else 2)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--write-synthetic", action="store_true")
    args = p.parse_args()
    if args.write_synthetic:
        save(fixture(), args.output)
        return 0
    try:
        if args.input is None:
            raise ValueError("--input required")
        result = compare(load(args.input))
    except (ValueError, TypeError, KeyError, OSError) as exc:
        result = dict(accepted=None, contract_valid=False, errors=[str(exc)], exit_code=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return result["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())
