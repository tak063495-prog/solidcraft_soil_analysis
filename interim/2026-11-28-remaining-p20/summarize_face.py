"""Attribute the measured 833 image to the original target group and its complement."""
import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
P19=HERE.with_name("2026-11-28-remaining-p19")
sys.path.insert(0,str(P19))
from leaf_score import score


def read_table(path):
    with gzip.open(path,"rt") as f:
        return {(int(r["block"]),int(r["leaf"])):r for r in csv.DictReader(f)}


def run(out):
    results={}
    for sign in ("minus","plus"):
        before=P19/"c02/baseline"
        old=read_table(before/f"{sign}-833-by-leaf.csv.gz")
        rec=json.loads((out/f"{sign}.json").read_text())
        last=rec["steps"][-1]
        summary=dict(consecutive=rec["consecutive"],total=rec["total"],first_failed=last["step"],
                     full_history_relative_upper=rec["full_history_relative_upper"],threshold=rec["threshold"])
        assert not last["included"]
        if last["step"]==833:
            s=score(out,sign)
            new=read_table(out/f"{sign}-833-by-leaf.csv.gz")
            assert old.keys()==new.keys()
            groups={}
            for label,target in (("target",True),("complement",False)):
                keys=[k for k,r in old.items() if ((r["refined"]=="True")==(sign=="plus"))==target]
                groups[label]=dict(leaves=len(keys),**{phase:{col:sum(float(table[k][col]) for k in keys)
                    for col in ("nominal_nm","projected_history_share_nm","total_nm")}
                    for phase,table in (("before",old),("after",new))})
            summary.update(image_nm=s["image_nm"],excess_nm=s["excess_nm"],fixed_nm=s["fixed_nm"],
                groups=groups,nominal_and_history=s["all"])
            np.testing.assert_allclose(s["fixed_nm"]+sum(g["after"]["total_nm"] for g in groups.values()),s["image_nm"],rtol=3e-13)
        results[sign]=summary
    (out/"summary.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
