"""Compare enclosures at identical anchors and inventory the per-step certificates."""
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P8=HERE.parent/"2026-11-26-remaining-p8/c02"
P9=HERE.parent/"2026-11-26-remaining-p9/c02"


def run(out):
    summary={}
    source=json.loads((P8/"comparison.json").read_text())
    for sign in ("plus","minus"):
        old=np.load(P9/(sign+".npz"))
        bounds=np.load(out/(sign+".npz"))
        child=np.load(P8/sign/"residual.npz")
        idx=np.arange(448)
        center=bounds["center_slip"]
        assert np.all(np.linalg.norm(center-child["child"][:,:,:2],axis=2)<=bounds["slip_radius"])
        assert np.all(np.linalg.norm(child["signed"],axis=2)<=bounds["slip_error"])
        h=[json.loads(l) for l in (out/f"{sign}-intervals.jsonl").read_text().splitlines()]
        new=bounds["traction_radius"].sum(1)
        before=(old["traction_radius_known"]+old["traction_radius_uncertain"]).sum(1)
        rows={}
        for phase in range(7):
            sub=idx[64*phase:64*(phase+1)]
            steps=[r for r in h if r["phase"]==phase]
            row=dict(steps=len(steps),inverse_positive=sum(r["inverse_candidate"] is not None for r in steps),
                     history_matrix_available=sum(r["history_matrix"] is not None for r in steps),
                     count_crossing=sum(r["count_cross_area"]>0 for r in steps),
                     inconclusive_tangent=sum(r["coercivity_lower"] is None for r in steps),
                     radius_p9_kn_rms=float(np.linalg.norm(before[sub])/8),
                     radius_p10_kn_rms=float(np.linalg.norm(new[sub])/8),
                     reduced_anchors=int((new[sub]<before[sub]).sum()))
            rows[phase]=row
        ref=source[sign]["bands"]["all"]
        failed=next((r for r in h if r["history_matrix"] is None),None)
        summary[sign]=dict(phases=rows,intervals=len(h),active_tractions_checked=len(h)*2660,
            signed_slips_checked=448*2660,leaf_child_integrals_checked=448*380,
            lower_positive=sum(r["inverse_candidate"] is not None for r in h),
            conditional_history_matrices=sum(r["history_matrix"] is not None for r in h),
            first_unresolved=dict(phase=failed["phase"],t=failed["t"]) if failed else None,
            p8_field_relative=ref["slip"]["relative"],p8_integral_relative=ref["own_active_signed_integral"]["relative"],
            fixed_path_upper=float(np.linalg.norm(bounds["slip_error"].sum(1))/ref["own_active_signed_integral"]["denominator"]),
            old_fixed_path_upper=float(np.linalg.norm(old["slip_rule_error_upper"].sum(1))/ref["own_active_signed_integral"]["denominator"]),
            continuous_reequilibrated_upper=None,passed=False)
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
