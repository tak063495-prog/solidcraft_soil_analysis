"""Run small separated refinement series, or a larger Windows-local series."""
import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.integrate import quad

from finite_volume import backward_euler, faces, semidiscrete
from reference import CASES, fields, primary, storage


def errors(value, truth):
    d = np.asarray(value)-truth
    return {"rms": float(np.sqrt(np.mean(d*d))), "linf": float(np.max(abs(d)))}


def with_orders(rows, step):
    for previous, current in zip(rows, rows[1:]):
        for quantity in ("head", "flux"):
            e0, e1 = previous[quantity]["rms"], current[quantity]["rms"]
            current[quantity]["order"] = (float(np.log(e0/e1)/np.log(previous[step]/current[step]))
                                            if e0 > 1e-14 and e1 > 1e-14 else None)
    return rows


def analytic_audit(name):
    """Independent quadratures and finite-difference residuals of exact fields."""
    c = CASES[name]
    change = quad(lambda x: float(storage(name, primary(name, x, c.end))[0]
                                  -storage(name, primary(name, x, c.start))[0]),
                  0, 1, epsabs=1e-12)[0]
    net = quad(lambda t: float(fields(name, 0., t)[1]-fields(name, 1., t)[1]),
               c.start, c.end, epsabs=1e-12)[0]
    t = (c.start+c.end)/2
    x = np.array([.2, .5, .8])
    delta = 1e-5
    d_storage = (storage(name, primary(name, x, t+delta))[0]
                 -storage(name, primary(name, x, t-delta))[0])/(2*delta)
    d_flux = (fields(name, x+delta, t)[1]-fields(name, x-delta, t)[1])/(2*delta)
    result = {"continuum_delta_storage_m": change, "continuum_net_inflow_m": net,
              "continuum_mass_defect_m": change-net,
              "pde_finite_difference_max_residual_per_h": float(np.max(abs(d_storage+d_flux)))}
    if name != "moving_level":
        x = np.linspace(0, 1, 33)
        a, b = fields(name, x, c.start, 64), fields(name, x, c.start, 128)
        result["series_64_128_max_head_difference_m"] = float(np.max(abs(a[0]-b[0])))
        result["series_64_128_max_flux_difference_m_per_h"] = float(np.max(abs(a[1]-b[1])))
    result["checks_passed"] = (abs(change-net) < 1e-11
        and result["pde_finite_difference_max_residual_per_h"] < 2e-8
        and result.get("series_64_128_max_head_difference_m", 0.) < 1e-13
        and result.get("series_64_128_max_flux_difference_m_per_h", 0.) < 1e-13)
    return result


def study(name, large=False):
    c = CASES[name]
    grid = [20, 40, 80, 160] if large else [10, 20, 40, 80]
    n_time = 160 if large else 80
    counts = [40, 80, 160, 320] if large else [10, 20, 40, 80]
    exact_time, _ = semidiscrete(name, n_time, c.start, c.end)
    qtime = faces(name, exact_time, c.end)[0]
    temporal = []
    for count in counts:
        h, mass = backward_euler(name, n_time, c.start, c.end, count)
        temporal.append({**mass, "head": errors(h, exact_time),
                         "flux": errors(faces(name, h, c.end)[0], qtime)})
    spatial = []
    for n in grid:
        h, mass = semidiscrete(name, n, c.start, c.end)
        tighter, _ = semidiscrete(name, n, c.start, c.end, rtol=2e-12, atol=2e-14)
        x = (np.arange(n)+.5)/n
        truth = primary(name, x, c.end)
        _, qtruth = fields(name, np.arange(n+1)/n, c.end)
        spatial.append({"n": n, "dx": 1/n, **mass,
                        "head": errors(h, truth),
                        "flux": errors(faces(name, h, c.end)[0], qtruth),
                        "time_control_rms": errors(h, tighter)["rms"]})
    # Local reference implementation checks, never application acceptance gates.
    time_ok = all(.8 < r["head"]["order"] < 1.2 for r in with_orders(temporal, "dt")[1:])
    space_ok = all(1.7 < r["head"]["order"] < 2.3 for r in with_orders(spatial, "dx")[1:])
    mass_ok = all(abs(r["mass_defect"]) < 1e-10 for r in temporal+spatial)
    control_ok = all(r["time_control_rms"] < .01*r["head"]["rms"] for r in spatial)
    analytical = analytic_audit(name)
    return {"time": temporal, "space": spatial, "analytical_diagnostics": analytical,
            "reference_checks": {"time_order": time_ok, "space_order": space_ok,
                                 "mass": mass_ok, "time_control": control_ok,
                                 "analytical": analytical["checks_passed"]},
            "application_acceptance": "未定義・未評価", "C-04": "未達", "A-04": "未達"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--profile", choices=["smoke", "windows-local"], default="smoke")
    args = p.parse_args()
    begin = time.perf_counter()
    report = {"schema": "transient-water-study-1", "profile": args.profile,
              "python": platform.python_version(), "numpy": np.__version__,
              "scipy": scipy.__version__, "platform": platform.platform(),
              "cases": {name: study(name, args.profile == "windows-local") for name in CASES}}
    report["elapsed_seconds"] = time.perf_counter()-begin
    report["reference_checks_passed"] = all(all(c["reference_checks"].values())
                                                for c in report["cases"].values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+"\n",
                           encoding="utf-8")
    print(json.dumps({"output": str(args.output), "reference_checks_passed": report["reference_checks_passed"],
                      "elapsed_seconds": report["elapsed_seconds"]}, ensure_ascii=False))
    return 0 if report["reference_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
