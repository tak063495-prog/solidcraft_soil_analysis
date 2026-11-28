"""Strict CSV/JSON pile-field comparison; no implicit conversion/interpolation.

CSV requires an explicit JSON metadata sidecar. JSON uses {metadata, rows}.
Rows: id, quantity, x_m, y_m, z_m, value. See README.ja.md for the contract.
No engineering tolerance has been defined for these references: accepted=null.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


class ContractError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ContractError(message)


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def read_json(path):
    def bad(value):
        raise ContractError("nonfinite JSON constant: " + value)
    return json.loads(Path(path).read_text(encoding="utf-8-sig"),
                      object_pairs_hook=strict_object, parse_constant=bad)


def confirmed(value):
    require(value is not None, "unconfirmed/null physical metadata")
    if isinstance(value, str):
        require(value.strip() and value.strip().lower() not in
                {"unknown", "unconfirmed", "tbd", "未確認", "未定義"}, "unconfirmed metadata")
    elif isinstance(value, dict):
        require(bool(value), "empty metadata object")
        for item in value.values():
            confirmed(item)
    elif isinstance(value, list):
        require(bool(value), "empty metadata list")
        for item in value:
            confirmed(item)
    elif isinstance(value, (int, float)):
        require(math.isfinite(value), "nonfinite metadata")
    else:
        raise ContractError("unsupported metadata value")


def validate_metadata(meta):
    require(isinstance(meta, dict) and set(meta) ==
            {"schema", "case", "units", "basis", "signs", "frame", "quantities", "expected_ids", "provenance"},
            "metadata keys incomplete or unsupported")
    require(meta["schema"] == "pile-field-comparison/1", "unsupported schema")
    require(meta["units"] == "m-kN-kPa", "explicit m-kN-kPa units required; convert upstream")
    case = meta["case"]
    require(isinstance(case, dict) and set(case) ==
            {"id", "model", "geometry", "materials", "loads", "boundary", "interface", "water", "reference_state"},
            "incomplete physical conditions")
    require(case["model"] in {"3d-solid-elastic-bonded", "1d-axial-linear-springs", "1d-winkler-eb-free-tip"},
            "unsupported pile model")
    confirmed(case)
    confirmed(meta["frame"])
    basis = meta["basis"]
    require(isinstance(basis, dict) and set(basis) == {"origin_m", "axes"}, "explicit basis required")
    try:
        origin, axes = np.asarray(basis["origin_m"], float), np.asarray(basis["axes"], float)
    except (ValueError, TypeError) as exc:
        raise ContractError("nonnumeric coordinate basis") from exc
    require(origin.shape == (3,) and axes.shape == (3, 3) and np.isfinite(origin).all() and np.isfinite(axes).all(),
            "invalid coordinate basis")
    require(np.allclose(axes @ axes.T, np.eye(3), rtol=0., atol=1e-12) and
            abs(np.linalg.det(axes)-1.) < 1e-12, "orthonormal right-handed basis required")
    quantities, ids, signs = meta["quantities"], meta["expected_ids"], meta["signs"]
    require(isinstance(quantities, dict) and bool(quantities) and isinstance(ids, dict) and isinstance(signs, dict),
            "quantities, ID inventory and signs required")
    require(set(quantities) == set(ids) == set(signs), "quantity/sign/ID inventories differ")
    for name, spec in quantities.items():
        require(isinstance(spec, dict) and set(spec) == {"unit", "component", "location", "definition"},
                "quantity unit/component/location/definition required")
        require(spec["unit"] in {"m", "rad", "kN", "kN*m", "kN/m", "kPa"}, "unsupported quantity unit")
        confirmed(spec)
        confirmed(signs[name])
        require(isinstance(ids[name], list) and len(ids[name]) > 0 and
                all(isinstance(i, str) and i.strip() for i in ids[name]), "nonempty string IDs required")
        require(len(set(ids[name])) == len(ids[name]), "duplicate expected ID")
    provenance = meta["provenance"]
    require(isinstance(provenance, dict) and set(provenance) == {"source", "sha256"}, "source provenance required")
    require(isinstance(provenance["source"], str) and bool(provenance["source"].strip()) and
            isinstance(provenance["sha256"], str) and re.fullmatch(r"[0-9a-f]{64}", provenance["sha256"]),
            "source and SHA256 required")
    return meta


def index_rows(meta, rows):
    require(isinstance(rows, list) and bool(rows), "empty or invalid field rows")
    result = {}
    for row in rows:
        require(isinstance(row, dict) and set(row) == {"id", "quantity", "x_m", "y_m", "z_m", "value"},
                "row keys must be id,quantity,x_m,y_m,z_m,value")
        identifier, quantity = row["id"], row["quantity"]
        require(isinstance(identifier, str) and isinstance(quantity, str), "string ID/quantity required")
        require(quantity in meta["quantities"] and identifier in meta["expected_ids"][quantity],
                "unexpected quantity or ID")
        key = quantity, identifier
        require(key not in result, "duplicate field ID")
        try:
            require(not any(isinstance(row[k], bool) for k in ("x_m", "y_m", "z_m", "value")), "boolean field value")
            numbers = np.array([float(row[k]) for k in ("x_m", "y_m", "z_m", "value")])
        except (TypeError, ValueError) as exc:
            raise ContractError("nonnumeric field value") from exc
        require(np.isfinite(numbers).all(), "nonfinite field value")
        result[key] = numbers
    expected = {(q, i) for q, ids in meta["expected_ids"].items() for i in ids}
    require(set(result) == expected, "missing IDs; partial comparison refused")
    return result


def load(path, metadata_path=None):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        require(metadata_path is not None, "CSV requires --actual-metadata JSON sidecar")
        meta = read_json(metadata_path)
        with path.open(newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None and len(set(reader.fieldnames)) == len(reader.fieldnames),
                    "duplicate/absent CSV headers")
            rows = list(reader)
    else:
        require(metadata_path is None, "JSON embeds metadata; sidecar is not allowed")
        envelope = read_json(path)
        require(isinstance(envelope, dict) and set(envelope) == {"metadata", "rows"}, "invalid JSON envelope")
        meta, rows = envelope["metadata"], envelope["rows"]
    validate_metadata(meta)
    index_rows(meta, rows)
    return meta, rows


def compare(reference, actual):
    ref_meta, ref_rows = reference
    act_meta, act_rows = actual
    validate_metadata(ref_meta)
    validate_metadata(act_meta)
    for key in ref_meta:
        if key != "provenance":
            require(ref_meta[key] == act_meta[key], "metadata mismatch: " + key)
    refs, acts = index_rows(ref_meta, ref_rows), index_rows(act_meta, act_rows)
    results = {}
    for quantity in ref_meta["quantities"]:
        keys = [(quantity, i) for i in ref_meta["expected_ids"][quantity]]
        for key in keys:
            require(np.array_equal(refs[key][:3], acts[key][:3]), "coordinate/ID mismatch; no interpolation")
        r = np.array([refs[k][3] for k in keys])
        a = np.array([acts[k][3] for k in keys])
        with np.errstate(over="ignore", invalid="ignore"):
            delta = a-r
            norm = float(np.linalg.norm(r))
            delta_norm = float(np.linalg.norm(delta))
            rms = float(np.sqrt(np.mean(delta**2)))
        require(np.isfinite(delta).all() and all(math.isfinite(v) for v in (norm, delta_norm, rms)),
                "comparison arithmetic overflow; rescale explicitly upstream")
        relative = None if norm == 0 else delta_norm/norm
        require(relative is None or math.isfinite(relative), "relative norm overflow")
        results[quantity] = dict(count=len(keys), unit=ref_meta["quantities"][quantity]["unit"],
                                 max_absolute_error=float(np.max(np.abs(delta))),
                                 rms_error=rms,
                                 relative_l2=relative,
                                 accepted=None, tolerance="未定義")
    return dict(schema="pile-field-difference/1", case=ref_meta["case"], frame=ref_meta["frame"],
                reference=ref_meta["provenance"], actual=act_meta["provenance"], quantities=results,
                accepted=None, tolerance="未定義", p6_complete=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--actual", type=Path, required=True)
    parser.add_argument("--actual-metadata", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compare(load(args.reference), load(args.actual, args.actual_metadata))
        result["input_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [args.reference, args.actual] +
                                  ([] if args.actual_metadata is None else [args.actual_metadata])}
        with args.output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(2, f"comparison refused: {exc}\n")


if __name__ == "__main__":
    main()
