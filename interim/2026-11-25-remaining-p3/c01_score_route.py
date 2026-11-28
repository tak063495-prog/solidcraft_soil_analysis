"""Apply the P2 signed, categorical and all-field gates at the common anchors."""
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P2=HERE.parent/"2026-11-25-remaining5-p2"
sys.path.insert(0,str(P2))
from c01_score import History, STAGES, CATEGORIES, check_dissipation_anchors, align, compare_pair, observed


def score(raw, out, fixed):
    out.mkdir(parents=True)
    series, routes = [], {}
    for n in (256,512,1024):
        path = raw/f"n{n}"
        case = json.loads((path/"stages.json").read_bytes())
        original = json.loads((P2/f"c01/case-{n}.json").read_bytes())
        targets = {"reverse":-0.075,"reload":0.06}
        for stage in original["stages"]:
            if stage["name"] in targets:
                stage["target_factor"] = targets[stage["name"]]
        if fixed:
            original["joint_system"]["global_step_control"]["path_error"]["slip_scale_m"] = n**2*(1e-4/1024**2)
        assert case == original
        control = case["joint_system"]["global_step_control"]["path_error"]
        result = json.loads((path/"result.json").read_bytes())
        rows = result["stages"]
        assert [r["stage"] for r in rows] == list(STAGES)[:len(rows)]
        assert result["seed_queries"] == result["arc_reviews"] == 0
        assert max(r["attempts"] for r in rows) <= 2065
        h = History(path/"fields","native",n,rows)
        anchors = check_dissipation_anchors(h,path/"dissipation.jsonl")
        series.append(h)
        routes[n] = dict(result, anchors=anchors, slip_budget_m=control["slip_scale_m"]/n**2)
        for name,value in ((f"manifest-{n}",h.meta),(f"case-{n}",case)):
            (out/f"{name}.json").write_text(json.dumps(value,indent=2)+"\n")
    # The same compare_pair and observed calls used by P2 for the common
    # prefix. The two categorical streams use exact integer comparisons.
    fields={}
    for name in series[0].meta["fields"]:
        cat = name in CATEGORIES
        adjacent={}
        for i,n in enumerate((256,512)):
            a,b=align(series[i:i+2],name,n)
            adjacent[n]=dict(times=len(a),**compare_pair(a,b,n,cat))
        a,b,c=align(series,name,256)
        if cat:
            common=dict(coarse_mismatches=int(np.count_nonzero(a!=b)),fine_mismatches=int(np.count_nonzero(b!=c)))
        else:
            common=observed(a,b,c)
        fields[name]=dict(adjacent=adjacent,common=common)
        if name=="plastic_slip":
            np.savez_compressed(out/"slip-common.npz",n256=a,n512=b,n1024=c)
            fields[name]["by_stage"]={stage:observed(*(x[i*256:(i+1)*256] for x in (a,b,c)))
                                      for i,stage in enumerate(STAGES)}
    complete=len(a)==1792 and all(r["complete"] for r in routes.values())
    passed=complete and all(v["common"]["trend"] in ("identical","decreasing_aligned")
                            for k,v in fields.items() if k not in CATEGORIES)
    passed=passed and all(v["common"]["coarse_mismatches"]==v["common"]["fine_mismatches"]==0
                          for k,v in fields.items() if k in CATEGORIES)
    (out/"summary.json").write_text(json.dumps(dict(routes=routes,fields=fields,common_times=len(a),
        full_domain=complete,finite_three_level_gate=passed),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    score(Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3]=="fixed")
