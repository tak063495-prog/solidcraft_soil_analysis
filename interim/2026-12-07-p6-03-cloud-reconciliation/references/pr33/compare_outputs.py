"""Strict diagnostic-only CSV/JSON exchange; no production schema change.

Exit 1 = invalid exchange, 2 = valid differences but acceptance undefined.
Never converts units/signs, guesses components, fills records or interpolates.
"""
import argparse
import csv
import json
import math
from pathlib import Path

SCHEMA = "p6-reference-exchange/1"
FIELDS = ("id","location","time_s","x_m","y_m","z_m","quantity","component","unit","value")
UNITS = {"stress":("Pa",{"xx","yy","zz","xy","yz","xz"}),
         "displacement":("m",{"x","y","z"}),
         "velocity":("m/s",{"x","y","z"}),
         "acceleration":("m/s^2",{"x","y","z"}),
         "period":("s",{"scalar"}), "frequency":("Hz",{"scalar"}),
         "energy":("J",{"scalar"})}
CONVENTIONS = {"stress_sign":"tension_positive", "shear":"tensor_components",
               "stress_basis":"dry_total_equals_effective", "axes":"right_handed_xyz",
               "vector_sign":"positive_along_named_global_axis"}


def unique(pairs):
    d = {}
    for k,v in pairs:
        if k in d:
            raise ValueError("duplicate JSON key: "+k)
        d[k] = v
    return d


def parse_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique,
                      parse_constant=lambda s: (_ for _ in ()).throw(ValueError("nonfinite JSON: "+s)))


def finite(x):
    return type(x) in (float,int) and math.isfinite(x)


def validate(d):
    if not isinstance(d,dict) or set(d) != {"schema","case","conventions","complete","records"}:
        raise ValueError("exact dataset keys required")
    if d["schema"] != SCHEMA or d["complete"] is not True or d["conventions"] != CONVENTIONS:
        raise ValueError("schema/completion/unit-sign-basis-axis conventions differ")
    case = d["case"]
    if not isinstance(case,dict) or set(case) != {"name","physical_conditions","time_meaning"}:
        raise ValueError("explicit case and physical conditions required")
    if not isinstance(case["name"],str) or not case["name"] or not isinstance(case["physical_conditions"],dict) or not case["physical_conditions"]:
        raise ValueError("empty case conditions")
    if case["time_meaning"] not in {"physical_seconds","static_state_label_not_physical_time"}:
        raise ValueError("time meaning required")
    def valid_tree(x):
        if isinstance(x,dict):
            return all(isinstance(k,str) and k and valid_tree(v) for k,v in x.items())
        if isinstance(x,list):
            return all(valid_tree(v) for v in x)
        return (isinstance(x,str) and bool(x)) or finite(x)
    if not valid_tree(case):
        raise ValueError("invalid/nonfinite case metadata")
    if not isinstance(d["records"],list) or not d["records"]:
        raise ValueError("records must be a nonempty list")
    indexed = {}
    for r in d["records"]:
        if not isinstance(r,dict) or set(r) != set(FIELDS):
            raise ValueError("exact record columns required")
        if not isinstance(r["id"],str) or not r["id"] or r["location"] not in {"node","gauss_point","physical_point","monitor"}:
            raise ValueError("explicit string ID/location required")
        if any(not finite(r[k]) for k in ("time_s","x_m","y_m","z_m","value")) or r["time_s"] < 0:
            raise ValueError("nonfinite/missing coordinate, time or value")
        rule = UNITS.get(r["quantity"])
        if rule is None or r["unit"] != rule[0] or r["component"] not in rule[1]:
            raise ValueError("quantity/component/SI unit mismatch")
        key = tuple(r[k] for k in ("id","location","time_s","quantity","component"))
        if key in indexed:
            raise ValueError("duplicate ID/time/component")
        indexed[key] = r
    return indexed


def load(path, manifest=None):
    path = Path(path)
    if path.suffix.lower() == ".json":
        if manifest is not None:
            raise ValueError("JSON carries its own metadata")
        d = parse_json(path)
    elif path.suffix.lower() == ".csv":
        if manifest is None:
            raise ValueError("CSV requires explicit --*-manifest; metadata never guessed")
        d = parse_json(manifest)
        if "records" in d:
            raise ValueError("CSV manifest must not contain records")
        with path.open(encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames != list(FIELDS):
                raise ValueError("exact ordered CSV header required")
            records = list(reader)
        for r in records:
            if set(r) != set(FIELDS) or any(v is None or v == "" for v in r.values()):
                raise ValueError("missing/extra CSV columns or values")
            for key in ("time_s","x_m","y_m","z_m","value"):
                r[key] = float(r[key])
        d["records"] = records
    else:
        raise ValueError("only .csv or .json accepted")
    validate(d)
    return d


def compare(reference, candidate):
    a,b = validate(reference),validate(candidate)
    for k in ("schema","case","conventions","complete"):
        # Numeric value types also match (True never equals 1 as metadata).
        if json.dumps(reference[k],sort_keys=True) != json.dumps(candidate[k],sort_keys=True):
            raise ValueError("metadata differs: "+k)
    if a.keys() != b.keys():
        raise ValueError("missing/extra ID, time, location or component; interpolation forbidden")
    groups = {}
    for key,r in a.items():
        s = b[key]
        if any(r[k] != s[k] for k in ("x_m","y_m","z_m","unit")):
            raise ValueError("coordinate/unit mismatch at "+str(key))
        group = r["quantity"]+"/"+r["component"]
        g = groups.setdefault(group,{"unit":r["unit"],"count":0,"sum_squared_difference":0.,"max_abs_difference":0.})
        difference = s["value"]-r["value"]
        if not math.isfinite(difference) or not math.isfinite(difference*difference):
            raise ValueError("difference overflow")
        g["count"] += 1
        g["sum_squared_difference"] += difference*difference
        g["max_abs_difference"] = max(g["max_abs_difference"],abs(difference))
    for g in groups.values():
        g["rms_difference"] = math.sqrt(g.pop("sum_squared_difference")/g["count"])
    return {"status":"diagnostic_only_tolerance_undefined","accepted":None,
            "tolerance":None,"tolerance_source":"未定義: no inherited gate for these independent elastic/dynamic cases",
            "records":len(a),"groups":groups,"interpolation":False}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference")
    p.add_argument("candidate")
    p.add_argument("output")
    p.add_argument("--reference-manifest")
    p.add_argument("--candidate-manifest")
    args = p.parse_args()
    try:
        out = compare(load(args.reference,args.reference_manifest),load(args.candidate,args.candidate_manifest))
        with Path(args.output).open("x",encoding="utf-8",newline="\n") as f:
            json.dump(out,f,ensure_ascii=False,indent=2,allow_nan=False)
            f.write("\n")
    except (ValueError,TypeError,OSError,KeyError) as exc:
        p.exit(1,str(exc)+"\n")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
