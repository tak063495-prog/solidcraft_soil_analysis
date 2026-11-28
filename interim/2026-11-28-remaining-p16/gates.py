"""Reapply the p11 gate to all archived fields at the same 1,792 loads."""
import json
from pathlib import Path
import sys

import numpy as np

from attribution import LEVELS,P14,P15,STAGES,observed,read,restore
from alternate_score import CATEGORIES


def measure(key,vals):
    if key in CATEGORIES:
        ab,bc = [int((a!=b).sum()) for a,b in zip(vals[:-1],vals[1:],strict=True)]
        return dict(mismatch=[ab,bc],passes=(ab==bc==0))
    r = observed(*vals)
    return dict(**r,passes=r["trend"] in ("identical","decreasing_aligned"))


def run(out):
    results = {}
    for name,p in (("outer",P15),("middle",P14)):
        d = p/"c01/alternate"/name
        ref = read(d/"score.json")
        fields = {}
        for key in ref["fields"]:
            vals = [restore(d,n,[key])[key] for n in LEVELS]
            assert all(len(a)==1792 for a in vals)
            stages = {}
            for j,s in enumerate(STAGES):
                stages[s] = measure(key,[a[j*256:(j+1)*256] for a in vals])
            fields[key] = dict(full=measure(key,vals),stages=stages)
        for k,sc in fields.items():
            for a,b in [(sc["full"],ref["fields"][k]["full"]),*[
                    (sc["stages"][s],ref["fields"][k]["stages"][s]) for s in STAGES]]:
                for metric,val in a.items():
                    if isinstance(val,float):
                        np.testing.assert_allclose(val,b[metric],rtol=1e-12,atol=1e-30)
                    else:
                        assert val == b[metric]
        count = sum(r["full"]["passes"] for r in fields.values())
        results[name] = dict(common_coordinates=1792,passed_fields=count,fields=fields,
                             passed=count==len(fields)==22)
    (out/"gates.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
