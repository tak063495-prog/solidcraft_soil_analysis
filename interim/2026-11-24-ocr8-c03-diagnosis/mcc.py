"""Stress-parametrized MCC solution on the OCR 8 dilating branch."""
import csv
import gzip
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
SPEC = json.loads((HERE.parent / "2026-11-24-t9d-benchmarks/criteria.json").read_bytes())["mcc"]
SOURCE = HERE.parent / "2026-11-24-independent-local-references"
P0, M, K, L = (SPEC[k] for k in ("p0_kpa", "M", "kappa", "lambda"))
PC0, G = 40., 250.
V0 = SPEC["v_lambda"] - L * math.log(PC0 / SPEC["p1_kpa"]) + K * math.log(PC0 / P0)


def pc(p):
    return p + 9 * (p - P0)**2 / (M*M*p)


def volume(p, cap):
    return V0 - K * math.log(p/P0) - (L-K) * math.log(cap/PC0)


def rate(p):
    # d eps_v^p / dp from logarithmic hardening and the exact state surface.
    cap = pc(p)
    cap_prime = 1 + (9/M**2) * (1 - (P0/p)**2)
    return (L-K) * cap_prime / (volume(p, cap) * cap)


def integrate(fun, left, right):
    return quad(fun, left, right, epsabs=1e-14, epsrel=1e-13)[0]


def closed_path(sample):
    py = brentq(lambda p: pc(p)-PC0, P0, PC0)
    axial = sample * SPEC["terminal_axial_strain"] / 512
    q = lambda p: 3*(p-P0)
    ey = math.log(V0/volume(py, PC0))/3 + q(py)/(3*G)
    if axial <= ey:
        p = brentq(lambda p: math.log(V0/volume(p, PC0))/3 + q(p)/(3*G) - axial, P0, py)
        cap, ep = PC0, 0.
    else:
        def shear_rate(p):
            return rate(p) * 2*q(p)/(M*M*(2*p-pc(p)))

        def residual(p):
            return math.log(V0/volume(p,pc(p)))/3 + q(p)/(3*G) + integrate(shear_rate, py, p) - axial

        critical = P0 / (1-M/3)
        p = brentq(residual, (py+critical)/2, py, xtol=1e-14)
        cap, ep = pc(p), integrate(rate, py, p)
    v = volume(p, cap)
    ev = math.log(V0/v)
    return dict(p=p, q=q(p), pc=cap, v=v, axial=axial, lateral=(ev-axial)/2,
                ev=ev, ev_plastic=ep, ev_elastic=ev-ep, yield_axial=ey, yield_p=py)


def assess(path):
    data = json.loads(gzip.decompress((SOURCE / "mcc-history.json.gz").read_bytes()))
    ref = np.asarray(data["fine_reference"]["40"])
    keys = ("p", "q", "v", "axial", "lateral", "ev")
    exact = [closed_path(i) for i in range(1,52)]
    exact_array = np.array([[x[k] for k in keys] for x in exact])
    rows = list(csv.DictReader(path.open()))
    levels, control = [], []
    for row in rows:
        if int(row["sample"]) != 50:
            continue
        num = {k:float(v) for k,v in row.items()}
        num["ev_elastic"] = num["ev"] - num["ev_plastic"]
        for k in ("ev", "ev_plastic", "ev_elastic"):
            num["delta_"+k] = num[k] - exact[49][k]
        num["log_volume_residual"] = num["ev"] - math.log(V0/num["v"])
        veq = volume(num["p"],num["pc"])
        num["ev_minus_state_surface_at_fem_stress"] = math.log(veq/num["v"])
        num["state_surface_shift_due_to_fem_stress"] = math.log(exact[49]["v"]/veq)
        levels.append(num)
    for run in data["histories"]:
        if run["pc_kpa"] != PC0:
            continue
        fem = np.array(run["fem"])
        subset = np.array([[float(r[k]) for k in keys] for r in rows
                  if float(r["h"]) == run["material_substep"] and int(r["load_divisions"]) == 1])
        control.append(dict(elements=run["elements"], h=run["material_substep"],
                       max_point_replay_abs_delta=dict(zip(keys,np.max(abs(subset-fem[:51]),axis=0))),
                       max_volumetric_identity_error=float(np.max(abs(fem[:,5]-np.log(V0/fem[:,2])))),
                       max_axial_coordinate_error=float(np.max(abs(fem[:,3]-ref[:,3]))),
                       point_49_to_51=fem[48:51].tolist()))
    yield_window = []
    for row in rows:
        i = int(row["sample"])
        if i in (30, 31) and float(row["h"]) in (1e-4, 5e-5) and int(row["load_divisions"]) == 1:
            yield_window.append(dict(sample=i, h=float(row["h"]), reference=exact[i-1],
                                     native={k: float(row[k]) for k in (*keys, "pc", "ev_plastic", "plastic_steps")}))
    coarse, fine = levels[:2]
    return dict(fields=keys, v0=V0, first_51_rk4_minus_quadrature_max_abs=dict(
                zip(keys, np.max(abs(ref[:51]-exact_array),axis=0))),
                yield_window=yield_window, independent_samples_49_to_51=exact[48:51], archived_cases=control, levels=levels,
                first_order_extrapolation=2*fine["ev"]-coarse["ev"],
                extrapolation_error=2*fine["ev"]-coarse["ev"]-exact[49]["ev"])


if __name__ == "__main__":
    print(json.dumps(assess(Path(sys.argv[1])), indent=2, allow_nan=False))
