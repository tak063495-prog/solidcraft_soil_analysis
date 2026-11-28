"""Score 512 sample histories against the pre-existing independent OCR 8 quadrature."""
import argparse
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-24-ocr8-c03-split"))
from mcc_compare import FIELDS, G, SOURCE, SPEC, differences


def score(runs):
    archive = json.loads(gzip.decompress((SOURCE/"mcc-history.json.gz").read_bytes()))
    assert tuple(archive["fields"]) == FIELDS
    ref = np.array(archive["fine_reference"]["40"])
    limit = SPEC["endpoint_relative_max_exclusive"]
    baseline = json.loads((SOURCE/"mcc.json").read_bytes())
    assert all(d["relative"] < limit/100 for d in baseline["refinement"]["40"][-1].values())
    results = []
    keys = set()
    for run in runs:
        e, width = run["elements"], run["material_substep"]
        mode = ("plastic_refinement" if "plastic_material_width" in run else
                "analytic_elastic" if run["analytic_elastic"] else "default")
        key = mode,e,width
        assert key not in keys
        keys.add(key)
        history = run["history"]
        assert run["ocr"] == 8 and run["constant_shear_kpa"] == G and len(history) == 512
        assert all(row["stage"] == "shear" and row["fraction"] == i/512
                   for i,row in enumerate(history,1))
        maxima = differences(np.array([[row[k] for k in FIELDS] for row in history]),ref)
        results.append(dict(mode=mode,elements=e,material_substep=width,
                            **{k:run[k] for k in ("control_substeps","cell_substeps","seconds",
                               "max_free_relative","max_lateral_error_kpa","endpoint_passed")},
                            maxima=maxima,passed=bool(run["endpoint_passed"] and
                                                     all(x["relative"] < limit for x in maxima.values()))))
    assert keys == {(m,e,h) for m in ("default","analytic_elastic","plastic_refinement")
                    for e in (12,96) for h in (1e-4,5e-5)}
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histories", nargs="+", type=Path)
    inputs = []
    for path in parser.parse_args().histories:
        raw = path.read_bytes()
        inputs.extend(json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw))
    print(json.dumps(score(inputs),indent=2,allow_nan=False))
