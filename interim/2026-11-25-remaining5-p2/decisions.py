"""Apply the user-approved T7 classification and registered Kirsch grid."""
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def decide(folder):
    result=json.loads((folder/"t7.json").read_bytes())
    limit=result["relative_spread_limit"]
    assert limit==0.10
    assert not any(result[k]["candidate"] for k in ("quartic_pair","pitchfork_pair"))
    original=json.loads((HERE.parent/"2026-11-24-t9d-benchmarks/criteria.json").read_bytes())
    low,high=original["arch"]["first_peak_load_reading"]
    beams={}
    for n,row in result["beams"].items():
        neighbour=result["beams"][str(row["neighbour"])]
        ca,cb=row["physical_curvature"],neighbour["physical_curvature"]
        simple=(all(row["checks"].values()) and all(neighbour["checks"].values())
                and ca*cb>0 and abs(ca-cb)<=limit*max(abs(ca),abs(cb)))
        assert simple==row["candidate"]
        peak=row["local_xy"][1][1]
        beams[n]=dict(official="simple fold" if simple else "unconfirmed",
                      peak_load=peak,within_original_load_reading=low<=peak<=high,
                      fold_and_load_pass=simple and low<=peak<=high,
                      relative_curvature_spread=row["relative_spread"],
                      previous_class=row["legacy"],previous_normalized=row["legacy_normalized"])
    k=json.loads((folder/"kirsch.json").read_bytes())
    raw=json.loads((folder/"kirsch-raw.json").read_bytes())
    assert (k["theta"],k["radial"],k["solves"],k["elements"],k["points"]["count"],len(raw))==(22,21,4,2772,12,12)
    np.testing.assert_allclose(1/np.array(k["radii_m"]),np.linspace(1,.1,22),atol=2e-15,rtol=0)
    for row,(radius,theta) in zip(raw,((r,np.deg2rad(t)) for r in (1.1,1.5,2.25) for t in (7,27,52,79)),strict=True):
        np.testing.assert_allclose(row["x_m"],[radius*np.cos(theta),radius*np.sin(theta),.37],atol=1e-14,rtol=0)
        assert min(row["barycentric"])>0 and len(row["raw_stress_tension_positive_kpa"])==6
    rms=k["points"]["rms"]["tensor_over_30MPa"]
    return dict(authority="user decision for T7 and Kirsch",T7=dict(legacy_threshold=1e-9,
        curvature_relative_width=limit,load_reading=[low,high],chord_beams=beams),
        Kirsch=dict(registered_grid="inverse-radius-22x21",original_raw_point_count=len(raw),
                    original_rms_limit=.03,tensor_rms_over_30MPa=rms,passed=rms<=.03))


if __name__=="__main__":
    print(json.dumps(decide(HERE),indent=2,allow_nan=False))
