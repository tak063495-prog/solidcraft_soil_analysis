"""Measure the default and selected arch element on the same continuation path."""
import argparse
import json
from pathlib import Path
import resource
import sys
import time

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent / "2026-11-24-t9d-benchmarks"
sys.path[:0] = [str(BENCH), str(HERE.parent / "2026-11-22-impl-t7-physical-instability")]
from arch import beam_curve, summarize
from curved import CurvedArch, AnalyticChordArch
from references import force_unit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("elements", type=int, choices=(12, 96))
    parser.add_argument("method", choices=("baseline", "new"))
    args = parser.parse_args()
    baseline = args.method == "baseline"
    model = None if baseline else (CurvedArch if args.elements == 12 else AnalyticChordArch)(args.elements)
    band = json.loads((BENCH / "criteria.json").read_bytes())["arch"]["first_peak_load_reading"]
    start = time.monotonic()
    case = beam_curve(args.elements, extrapolate=True, model=model, fold_metrics=True)
    fold = case["fold"]
    row = dict(elements=args.elements, baseline=baseline, seconds=time.monotonic()-start,
               peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               **summarize(case, "fold_load", band),
               curvature_over_scale=fold.get("curvature_over_scale"),
               fd_curvature_over_scale=fold.get("fd_curvature_over_scale"),
               candidate_load=float(fold["z"][-1]*force_unit) if "z" in fold else None)
    print(json.dumps(row, allow_nan=False))


if __name__ == "__main__":
    main()
