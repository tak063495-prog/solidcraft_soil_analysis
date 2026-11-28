"""Report archived positive bounds, split matrices, and residual-inclusive balls."""
import json
from pathlib import Path
import sys

def rows(p):
    return [json.loads(line) for line in p.read_text().splitlines()]

def run(folder):
    p10=Path(__file__).resolve().parent.parent/"2026-11-26-remaining-p10/c02"
    before=json.loads((p10/"summary.json").read_text())
    summary={}
    for sign in ("plus","minus"):
        old=rows(p10/(sign+"-intervals.jsonl"))
        values={}
        for radius in ("0","1e-08"):
            seq=[rows(path/f"{sign}-{radius}.jsonl") for path in (folder,folder/"weighted")]
            by_phase={}
            chain=True
            for i,(a,b,c) in enumerate(zip(old,*seq,strict=True)):
                assert a["phase"]==b["phase"]==c["phase"] and a["t"]==b["t"]==c["t"]
                assert b["residual_upper_kn"]==c["residual_upper_kn"]
                lower=max(b["lower"] or 0,c["lower"] or 0)
                include=bool(lower and float(radius)*lower>b["residual_upper_kn"])
                chain=chain and include
                r=by_phase.setdefault(a["phase"],dict(steps=0,p10=0,unweighted=0,weighted=0,positive=0,includes=0,consecutive=0))
                for k,v in (("steps",1),("p10",bool(a["inverse_candidate"])),("unweighted",bool(b["lower"])),
                            ("weighted",bool(c["lower"])),("positive",bool(lower)),("includes",include),("consecutive",chain)):
                    r[k]+=v
            values[radius]=dict(phases=by_phase,total={k:sum(r[k] for r in by_phase.values()) for k in by_phase[0]},
                                first={k:v for k,v in seq[1][0].items() if k in ("phase","t","lower","residual_upper_kn","schur_margin")})
        summary[sign]=dict(radii_m=values,fixed_saved_path_slip_relative_upper=before[sign]["fixed_path_upper"],
            reequilibrated_continuous_history_relative_upper=None,threshold=.001,
            passed=False)
    (folder/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")

if __name__=="__main__":
    run(Path(sys.argv[1]))
