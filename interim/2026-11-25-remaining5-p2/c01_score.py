"""Use the original complete-grid, signed-direction and categorical gates."""
from pathlib import Path
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"2026-11-24-c01-event4-retry"))
from c01_retry_fields import CATEGORIES,STAGES
from c01_retry_metrics import History,align,compare_pair,observed
from collect import check_dissipation_anchors


def score(raw,out):
    out.mkdir()
    histories, routes=[],{}
    for n in (256,512,1024):
        path=raw/f"n{n}"
        case=json.loads((path/"stages.json").read_bytes())
        control=case["joint_system"]["global_step_control"]["path_error"]
        assert control["reequilibrate_endpoints"] and control["slip_scale_m"]==1e-4
        result=json.loads((path/"result.json").read_bytes())
        rows=result["stages"]
        assert [r["stage"] for r in rows]==list(STAGES)
        assert result["complete"]==all(r["last_accepted"]==1 for r in rows)
        assert max(r["attempts"] for r in rows)<=2065 and result["seed_queries"]==result["arc_reviews"]==0
        h=History(path/"fields","native",n,rows)
        anchors=check_dissipation_anchors(h,path/"dissipation.jsonl")
        histories.append(h)
        routes[n]=dict(result,anchors=anchors,slip_budget_m=control["slip_scale_m"]/n**2)
        (out/f"manifest-{n}.json").write_text(json.dumps(h.meta,indent=2)+"\n")
        (out/f"case-{n}.json").write_bytes((path/"stages.json").read_bytes())
    fields={}
    for field in histories[0].meta["fields"]:
        category=field in CATEGORIES
        adjacent={}
        for i,n in enumerate((256,512)):
            a,b=align(histories[i:i+2],field,n)
            adjacent[n]=dict(times=len(a),**compare_pair(a,b,n,category))
        a,b,c=align(histories,field,256)
        metric=(dict(coarse_mismatches=int(np.count_nonzero(a!=b)),fine_mismatches=int(np.count_nonzero(b!=c)))
                if category else observed(a,b,c))
        fields[field]=dict(adjacent=adjacent,common=metric)
        if field=="plastic_slip":
            np.savez_compressed(out/"slip-common.npz",n256=a,n512=b,n1024=c)
            fields[field]["by_stage"]={name:observed(*(arr[k*256:(k+1)*256] for arr in (a,b,c))) for k,name in enumerate(STAGES)}
    complete=all(row["complete"] for row in routes.values()) and len(a)==1792
    passed=complete and all(v["common"]["trend"] in ("decreasing_aligned","identical") for k,v in fields.items() if k not in CATEGORIES)
    passed=passed and all(v["common"]["coarse_mismatches"]==v["common"]["fine_mismatches"]==0 for k,v in fields.items() if k in CATEGORIES)
    (out/"summary.json").write_text(json.dumps(dict(routes=routes,fields=fields,common_times=len(a),
        full_domain=complete,finite_three_level_gate=passed),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    score(Path(sys.argv[1]),Path(sys.argv[2]))
