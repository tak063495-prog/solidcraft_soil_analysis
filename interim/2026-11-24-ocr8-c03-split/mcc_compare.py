"""Score all four OCR 8 cases and decompose the sample-50 strain error."""
import argparse
import gzip
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "2026-11-24-independent-local-references"
sys.path[:0] = [str(SOURCE), str(HERE.parent / "2026-11-24-ocr8-c03-diagnosis")]
from mcc_history import FIELDS, differences
from mcc import closed_path, SPEC, K, M, G, PC0


def intersection(rows, h):
    # The specimen solves for one lateral increment per axial interval. Rebuild
    # the elastic path on that actual affine interval from the preceding sample.
    start, end = rows[29:31]
    ev = end["volumetric_strain"] - start["volumetric_strain"]
    dz = end["axial_strain"] - start["axial_strain"]
    dr = end["lateral_strain"] - start["lateral_strain"]
    def f(t):
        p = start["p_kpa"] * math.exp(-start["v"] * math.expm1(-ev*t)/K)
        q = start["q_kpa"] + 2*G*(dz-dr)*t
        return ((q/M)**2 + p*(p-PC0))/PC0**2
    t = brentq(f, 0., 1., xtol=1e-15)
    subdivisions = math.ceil((2*abs(dr)+abs(dz))/h)
    k = math.floor(t*subdivisions)
    lo, hi = k/subdivisions, (k+1)/subdivisions
    assert f(lo) < 0 < f(hi)
    for _ in range(52):
        mid = (lo+hi)/2
        if f(mid)<=0:
            lo=mid
        else:
            hi=mid
    return dict(fraction=t, substep=k+1, subdivisions=subdivisions,
                bisection_fraction=lo, bisection_width=hi-lo,
                bisection_residual=f(lo), bisection_minus_brent=lo-t,
                axial=start["axial_strain"]+dz*t,
                normalized_yield=f(t), initial=f(0.), uncorrected_end=f(1.))


def compare(runs):
    archive = json.loads(gzip.decompress((SOURCE/"mcc-history.json.gz").read_bytes()))
    ref = np.array(archive["fine_reference"]["40"])
    controls = json.loads((SOURCE/"mcc.json").read_bytes())
    limit = SPEC["endpoint_relative_max_exclusive"]
    ref_ok = all(x["relative"] < limit/100 for x in controls["refinement"]["40"][-1].values())
    assert ref_ok and tuple(archive["fields"]) == FIELDS
    previous = {(s["elements"],s["material_substep"]): s["fem"]
                for s in archive["histories"] if s["pc_kpa"] == PC0}
    exact = closed_path(50)
    cases = []
    expected = {(e,h,b) for e,h in previous for b in (False,True)}
    for run in runs:
        e,h,b = (run[k] for k in ("elements","material_substep","analytic_elastic"))
        expected.remove((e,h,b))
        history = run["history"]
        assert run["ocr"] == 8 and run["constant_shear_kpa"] == G and len(history) == 512
        assert all(row["stage"]=="shear" and row["fraction"]==i/512 for i,row in enumerate(history,1))
        values = np.array([[s[k] for k in FIELDS] for s in history])
        if not b:
            np.testing.assert_array_equal(values, previous[e,h])
        maxima = differences(values,ref)
        sample = history[49]
        ev, ep = sample["volumetric_strain"], sample["plastic_volumetric_strain"]
        details = dict(total=ev, elastic=ev-ep, plastic=ep, delta_total=ev-exact["ev"],
                       delta_elastic=ev-ep-exact["ev_elastic"], delta_plastic=ep-exact["ev_plastic"])
        case = dict(elements=e, h=h, analytic=b, maxima=maxima,
                    passed=ref_ok and all(x["relative"]<limit for x in maxima.values()),
                    point50=details, cell_steps=run["cell_substeps"], control_steps=run["control_substeps"],
                    seconds=run["seconds"])
        if b:
            case["intersection"] = intersection(history,h)
            case["elastic_sample30_abs_error"] = abs(history[29]["volumetric_strain"]-closed_path(30)["ev"])
        cases.append(case)
    assert not expected and len(cases)==8
    return dict(limit_exclusive=limit, reference_point50=exact, cases=cases)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("histories",type=Path)
    p.add_argument("output",type=Path)
    a = p.parse_args()
    a.output.write_text(json.dumps(compare(json.loads(a.histories.read_bytes())),indent=2,
                                   allow_nan=False)+"\n",encoding="utf8")
