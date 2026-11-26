import json
from pathlib import Path
import sys

import numpy as np

from face_bounds import P8, p6


def run(out):
    prior=json.loads((P8/"c02/comparison.json").read_text())
    b,_=p6.face.triangle_rule(1)
    child=np.where(b.max(axis=1)>.5,b.argmax(axis=1),3)
    summary={}
    for sign in ("plus","minus"):
        z=np.load(out/f"{sign}.npz")
        q=np.load(P8/"c02"/sign/"slip.npz")
        residual=np.load(P8/"c02"/sign/"residual.npz")
        classes=z["class_code"][:,:,child]
        observed=q["status"].reshape(classes.shape)
        assert not np.any((classes<3)&(classes!=observed))
        gap=np.linalg.norm(z["slip_center_integral"]-residual["child"][:,:,:2],axis=2)
        assert np.all(gap<=z["slip_radius_known"]+z["slip_radius_uncertain"])
        assert np.all(np.linalg.norm(residual["signed"],axis=2)<=z["slip_rule_error_upper"])
        area=z["uncertain_area"].sum(axis=1)
        row=dict(uncertain_area_mean=float(area.mean()),uncertain_area_max=float(area.max()),
                 fully_unresolved_anchors=int(sum(area>1-1e-12)),
                 area_per_stage_mean=area.reshape(7,64).mean(1).tolist(),
                 class_counts=np.bincount(z["class_code"].ravel(),minlength=4).tolist(),
                 checked_resolved_active_labels=int((classes<3).sum()))
        for name in ("slip_rule_error_upper","slip_radius_known","slip_radius_uncertain",
                     "traction_radius_known","traction_radius_uncertain","traction_leaf_oscillation_upper"):
            v=z[name]
            row[name]=dict(max_leaf=float(v.max()),max_sum=float(v.sum(1).max()),
                           sum_time_l2=float(np.linalg.norm(v.sum(1))))
        bands=prior[sign]["bands"]["all"]
        row["fixed_path_integral_bound_over_p8_reference"]=row["slip_rule_error_upper"]["sum_time_l2"]/bands["own_active_signed_integral"]["denominator"]
        row["p8_field_difference"]=bands["slip"]["relative"]
        row["p8_integral_difference"]=bands["own_active_signed_integral"]["relative"]
        summary[sign]=row
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
