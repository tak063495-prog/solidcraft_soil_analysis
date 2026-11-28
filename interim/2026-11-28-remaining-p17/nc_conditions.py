"""Tabulate each failed line-search candidate against the original MCC predicates."""
import csv
import json
from pathlib import Path
import sys


def run(out):
    summary=json.loads((out/"summary.json").read_text())
    with (out/"candidates.csv").open("w",newline="") as stream:
        writer=csv.writer(stream,lineterminator="\n")
        writer.writerow(["path","parent","point","evaluation","iteration","k","alpha","finite","log_pc",
            "log_admissible","multiplier","multiplier_nonnegative","pressure","pressure_positive",
            "original_norm","candidate_norm","strict_decrease","all_admissible","accepted"])
        for path,result in summary["paths"].items():
            rows=list(map(json.loads,(out/f"{path}-material.jsonl").read_text().splitlines()))
            start=max(i for i,r in enumerate(rows) if "iteration" in r)
            initial=rows[start]["norm"]
            assert initial==result["residual"]
            candidates=[r for r in rows[start+1:] if "alpha" in r]
            assert len(candidates)==result["candidate_count"]
            for k,r in enumerate(candidates):
                assert r["alpha"]==2.**-k
                finite=r["finite"]
                hardening=abs(r["log_pc"])<100
                multiplier=r["multiplier"]>=0
                compression=r["pressure"]>0
                admissible=finite and hardening and multiplier and compression
                decrease=r["candidate_norm"] is not None and r["candidate_norm"]<initial
                assert not(admissible and decrease)
                writer.writerow([path,result["parent"],351,result["material_calls"],rows[start]["iteration"],
                    k,r["alpha"],finite,r["log_pc"],hardening,r["multiplier"],multiplier,r["pressure"],compression,
                    initial,r["candidate_norm"],decrease,admissible,admissible and decrease])


if __name__=="__main__":
    run(Path(sys.argv[1]))
