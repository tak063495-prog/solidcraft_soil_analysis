"""Strict CSV/JSON external-output diagnostics. Never invokes an application.

No production tolerances exist for these three transient benchmark adapters.
Success reading input is not numerical acceptance. Exit 2 means valid report
with undefined acceptance; exit 1 means invalid input; exit 0 is demo only.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from diagnose import errors
from reference import CASES, fields, metadata, storage


def finite(value):
    if isinstance(value, bool):
        raise ValueError("boolean is not a measurement")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError("non-finite measurement")
    return result


def read_input(path, meta_path=None):
    if path.suffix.lower() == ".csv":
        if meta_path is None:
            raise ValueError("CSV requires --metadata JSON (units/coordinates)")
        data = json.loads(meta_path.read_text(encoding="utf-8-sig"))
        with path.open(encoding="utf-8-sig", newline="") as f:
            data["samples"] = list(csv.DictReader(f))
        return data
    if path.suffix.lower() != ".json":
        raise ValueError("expected .csv or .json")
    if meta_path is not None:
        raise ValueError("JSON embeds metadata; do not also supply --metadata")
    return json.loads(path.read_text(encoding="utf-8-sig"))


def compare(data):
    if not isinstance(data, dict):
        raise ValueError("top-level JSON must be an object")
    name = data.get("case")
    if name not in CASES:
        raise ValueError("unknown benchmark case")
    for key, expected in metadata(name).items():
        if data.get(key) != expected or type(data.get(key)) is not type(expected):
            raise ValueError(f"metadata mismatch: {key}: expected {expected}")
    samples = data.get("samples")
    if not isinstance(samples, list) or len(samples) < 2:
        raise ValueError("at least two head/flux samples required")
    seen = set()
    actual_h, actual_q, exact_h, exact_q = [], [], [], []
    times = set()
    per_sample = []
    for row in samples:
        if not isinstance(row, dict):
            raise ValueError("sample must be an object")
        t, x, h, q = [finite(row[k]) for k in ("time", "x", "head", "flux")]
        if t < CASES[name].start or t > CASES[name].end:
            raise ValueError("sample outside documented smooth benchmark time window")
        if (t, x) in seen:
            raise ValueError("duplicate (time,x)")
        seen.add((t, x))
        times.add(t)
        hh, qq = [float(v) for v in fields(name, x, t)]
        actual_h.append(h)
        actual_q.append(q)
        exact_h.append(hh)
        exact_q.append(qq)
        per_sample.append({"time": t, "x": x, "head_error": h-hh, "flux_error": q-qq})
    if len(times) < 2:
        raise ValueError("time series must contain at least two distinct times")
    pde = {"head_m": errors(actual_h, exact_h), "flux_m_per_h": errors(actual_q, exact_q),
           "weighting": "equal supplied samples; no interpolation/extrapolation",
           "per_sample": per_sample, "tolerance": "未定義", "accepted": None}
    balance = compare_balance(name, data.get("balance", []))
    mechanics = compare_mechanics(name, data.get("mechanics", []))
    return {"schema": "transient-water-comparison-1", "case": name,
            "input_provenance": data.get("provenance", "未提供"),
            "pde": pde, "mass_balance": balance, "mechanics": mechanics,
            "acceptance": "未定義", "accepted": None,
            "C-04": "未達", "A-04": "未達", "P6": "未達"}


def compare_balance(name, rows):
    """Interval cumulative amounts per unit area, not trapezoids of sparse q."""
    if not isinstance(rows, list):
        raise ValueError("balance must be a list")
    if not rows:
        return {"status": "未評価: balance 欠落", "accepted": None}
    results = []
    previous = None
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("balance interval must be an object")
        t0, t1, fi, fo, s0, s1 = [finite(row[k]) for k in
            ("start", "end", "inflow_m", "outflow_m", "storage_start_m", "storage_end_m")]
        if t1 <= t0 or fi < 0 or fo < 0 or t0 < CASES[name].start or t1 > CASES[name].end:
            raise ValueError("invalid balance interval or cumulative flows")
        if previous is not None and (t0 != previous[0] or abs(s0-previous[1]) > 1e-12):
            raise ValueError("balance intervals must be contiguous with continuous storage")
        defect = s1-s0-fi+fo
        results.append({"start": t0, "end": t1, "defect_m": defect,
                        "relative": abs(defect)/max(abs(s1-s0), fi+fo, 1e-30)})
        previous = (t1, s1)
    return {"status": "計測のみ", "intervals": results, "tolerance": "未定義", "accepted": None}


def compare_mechanics(name, rows):
    if not isinstance(rows, list):
        raise ValueError("mechanics must be a list")
    if not rows:
        return {"status": "未評価: mechanics 欠落", "accepted": None}
    if name != "terzaghi":
        raise ValueError("mechanical reference exists only for linear 1-D Terzaghi")
    result = []
    seen = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("mechanics row must be an object")
        t = finite(row["time"])
        if not CASES[name].start <= t <= CASES[name].end or t in seen:
            raise ValueError("invalid/duplicate mechanics time")
        seen.add(t)
        if row.get("law") != "linear_1d_mv_0.001_per_kPa_gamma_9.81":
            raise ValueError("mechanical law/units not the documented fixed benchmark")
        x = np.array([finite(v) for v in row["x"]])
        hh = np.array([finite(v) for v in row["head"]])
        if len(x) < 3 or x.shape != hh.shape or x[0] != 0 or x[-1] != 1 or np.any(np.diff(x) <= 0):
            raise ValueError("mechanics requires ordered full-domain nodal head including 0 and 1")
        exact_nodal = fields(name, x, t)[0]
        # Same quadrature for the conditional split, plus an explicit sampling error.
        scale = .001*9.81
        conditional = scale*(1-np.trapezoid(hh, x))
        sampled_exact = scale*(1-np.trapezoid(exact_nodal, x))
        exact = scale*(1-quad(lambda z: float(fields(name, z, t)[0]), 0, 1,
                              epsabs=1e-13)[0])
        observed = finite(row["settlement_m"])
        result.append({"time": t, "mechanical_path_error_m": observed-conditional,
                       "pde_propagation_error_m": conditional-sampled_exact,
                       "sampling_quadrature_error_m": sampled_exact-exact,
                       "total_error_m": observed-exact})
    return {"status": "線形1Dのみ・湿潤JOINT経路は未評価", "rows": result,
            "tolerance": "未定義", "accepted": None}


def demo(name, directory):
    """Synthetic exact samples for adapter verification, never application output."""
    c = CASES[name]
    data = metadata(name)
    data["provenance"] = "synthetic analytical adapter demo; NOT production evidence"
    data["samples"] = []
    for t in np.linspace(c.start, c.end, 3):
        for x in np.linspace(0, 1, 5):
            h, q = fields(name, x, t)
            data["samples"].append({"time": float(t), "x": float(x),
                                    "head": float(h), "flux": float(q)})
    if name == "terzaghi":
        x = np.linspace(0, 1, 33)
        h = fields(name, x, c.end)[0]
        data["mechanics"] = [{"time": c.end, "law": "linear_1d_mv_0.001_per_kPa_gamma_9.81",
                              "x": x.tolist(), "head": h.tolist(),
                              "settlement_m": float(.001*9.81*(1-np.trapezoid(h, x)))}]
    directory.mkdir(parents=True, exist_ok=True)
    (directory/f"{name}.json").write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    meta = {k: v for k, v in data.items() if k != "samples"}
    (directory/f"{name}.metadata.json").write_text(json.dumps(meta, indent=2)+"\n", encoding="utf-8")
    with (directory/f"{name}.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["time", "x", "head", "flux"])
        writer.writeheader()
        writer.writerows(data["samples"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path, nargs="?")
    p.add_argument("--metadata", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--demo", type=Path)
    args = p.parse_args()
    if args.demo:
        for name in CASES:
            demo(name, args.demo)
        return 0
    if not args.input or not args.output:
        p.error("input and --output required (or --demo directory)")
    try:
        report = compare(read_input(args.input, args.metadata))
    except (ValueError, KeyError, TypeError, OSError) as e:
        p.exit(1, f"invalid input: {e}\n")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("診断を書き出しました。受入許容値は未定義（exit 2）。")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
