"""Attribute accumulated plastic flow at sample 50 using plastic shear as coordinate."""
import csv
import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

root = Path(sys.argv[1]) / "interim"
raw = Path(sys.argv[2])
sys.path[:0] = [str(root/"2026-11-24-ocr8-c03-diagnosis"),
                str(root/"2026-11-24-independent-local-references")]
from mcc import G, K, L, M, P0, V0, closed_path, pc, volume
from mcc_history import FIELDS, differences


def records(file):
    with (raw/file).open() as stream:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]


def normal_ratio(p):
    return M**2*(2*p-pc(p))/(6*(p-P0))


def rates(shear, state):
    p = state[0]
    slope = 1 + 9/M**2*(1-(P0/p)**2)
    ratio = normal_ratio(p)
    return [volume(p,pc(p))*pc(p)*ratio/((L-K)*slope), ratio]


exact = closed_path(50)
shear50 = exact["axial"]-exact["ev"]/3-exact["q"]/(3*G)
reference = solve_ivp(rates, [0,1.02*shear50], [exact["yield_p"],0.],
                      dense_output=True, method="DOP853",rtol=1e-13,atol=1e-15)
assert reference.success
ref_end = reference.sol(shear50)
assert abs(ref_end[0]-exact["p"]) < 1e-11 and abs(ref_end[1]-exact["ev_plastic"]) < 1e-13
samples, plastic = records("samples.csv"), records("plastic.csv")
ref_dir = root/"2026-11-24-independent-local-references"
archive = json.loads(gzip.decompress((ref_dir/"mcc-history.json.gz").read_bytes()))
stored = json.loads(gzip.decompress((root/"2026-11-24-ocr8-c03-split"/"mcc-histories.json.gz").read_bytes()))
ref = np.asarray(archive["fine_reference"]["40"])
assert tuple(archive["fields"]) == FIELDS
results = []
for h,div in sorted({(r["h"],int(r["divisions"])) for r in samples},reverse=True):
    rows = [r for r in samples if r["h"]==h and r["divisions"]==div]
    assert [r["sample"] for r in rows]==list(range(1,len(rows)+1))
    values = np.array([[r[k] for k in ("p","q","v","axial","lateral","ev")] for r in rows])
    end = rows[49]
    item = dict(h=h,load_divisions=div,samples=len(rows),
                maxima_fraction=differences(values,ref[:len(rows)]),
                delta_ev=end["ev"]-exact["ev"],delta_ep=end["ep"]-exact["ev_plastic"],
                log_volume_max_abs=max(abs(r["ev"]-math.log(V0/r["v"])) for r in rows),
                accepted_material_steps=int(rows[-1]["steps"]),control_steps=int(rows[-1]["control_steps"]))
    if div==1 and h>=5e-5:
        cross = {}
        for run in stored:
            if run["analytic_elastic"] and run["material_substep"]==h:
                hist = np.asarray([[r[k] for k in FIELDS] for r in run["history"]])
                cross[run["elements"]] = np.max(abs(values-hist),axis=0).tolist()
        assert set(cross)=={12,96} and np.max(list(cross.values()))<1e-9
        item["max_replay_minus_archive_abs_in_field_order"] = cross
        trace = [r for r in plastic if r["h"]==h]
        ds = np.array([r["sp"] for r in trace])
        coords = np.cumsum(ds)
        assert np.min(ds)>0 and coords[-1]<reference.t[-1]
        p_ref,ep_ref = reference.sol(coords)
        right = np.dot(normal_ratio(p_ref),ds)
        endpoint_flow = sum(M**2*(2*r["p"]-r["pc"])/(2*r["q"])*r["sp"] for r in trace)
        terms = dict(exact_path_endpoint_quadrature=right-ep_ref[-1],
                     actual_stress_direction=endpoint_flow-right,
                     terminal_plastic_shear=ep_ref[-1]-exact["ev_plastic"],
                     discrete_flow_closure=end["ep"]-endpoint_flow)
        assert abs(sum(terms.values())-item["delta_ep"])<2e-15
        item["plastic_error_terms"] = terms
        # Compare inverse hardening quadratures on the same accepted pressure increments.
        # Crossing substep can include an elastic fraction. Restrict this diagnostic
        # to substeps that start plastic to use v_before at the actual plastic parent.
        item["max_plastic_normalized_yield"] = max(abs((r["q"]/M)**2+r["p"]*(r["p"]-r["pc"]))/r["pc"]**2 for r in trace)
        tr=trace[1:]
        item["plastic_step_elastic_quadrature_defect"] = sum(
            r["total_ev"] - r["ep"] - K*math.log(r["p"]/r["p_before"])*(1/r["v_before"]+1/r["v_after"])/2 for r in tr)
        item["hardening_after_first_plastic_step"] = {
            "old_v_inverse_law_defect":sum(r["ep"]-(L-K)*r["log_pc_ratio"]/r["v_before"] for r in tr),
            "trapezoid_minus_old_v_same_pc_history":sum((L-K)*r["log_pc_ratio"]*0.5*(1/r["v_after"]-1/r["v_before"]) for r in tr)}
    results.append(item)
print(json.dumps({"field_order":FIELDS,"point50_reference":exact,"cases":results},indent=2,allow_nan=False))
