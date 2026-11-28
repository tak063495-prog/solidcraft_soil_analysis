"""Read archived full-field and near-hole observations, preserving their original scales."""
import json
import math
import sys
from pathlib import Path

interim=Path(sys.argv[1])/"interim"
def read(path):
    return json.loads((interim/path).read_bytes())

slip={}
for name,ref in read("2026-11-24-approved-four/c01-fields.json")["refinements"].items():
    s=ref["fields"]["plastic_slip"]
    parts={stage: row["coarse_norm"]*row["fine_norm"]*(row["cosine"] or 0.)
           for stage,row in s["by_stage"].items()}
    total=sum(parts.values())
    common=s["common"]
    assert math.isclose(total,common["coarse_norm"]*common["fine_norm"]*common["cosine"],rel_tol=1e-12)
    slip[name]={"inner_product_by_stage":parts,
                "reverse_reload_return_fraction_of_signed_total":sum(parts[k] for k in ("reverse","reload","return"))/total,
                "coarse_to_medium_max_location":s["adjacent"]["256"]["maximum_location"],
                "medium_to_fine_max_location":s["adjacent"]["512"]["maximum_location"]}
hole=[]
for case in read("2026-11-24-t9d-benchmarks/kirsch-results.json"):
    end=case["stages"][-1]
    assert end["release"]==1
    points=end["near_field"]
    total=sum(p["weight"] for p in points)
    hole.append({"tet4":case["retained_tets"],"near_gp":len(points),
                 "weighted_rms_over_30_MPa":math.sqrt(sum(p["weight"]*p["relative_to_30MPa"]**2 for p in points)/total) if total else None,
                 "max_over_30_MPa":max((p["relative_to_30MPa"] for p in points),default=None),
                 "wall_displacement_m":end["wall_axis_displacement_m"],
                 "wall_reference_m":end["wall_axis_reference_m"]})
print(json.dumps({"c01":slip,"kirsch":hole},indent=2,allow_nan=False))
