"""Match counterfactuals by attempt and full parent, and compare equal-budget paths."""
from collections import Counter
import gzip
import json
import shutil
import tarfile
from pathlib import Path
import sys
from freeze import HERE, P7, digest, save


def read_trace(path):
    with (gzip.open(path,"rt") if path.suffix==".gz" else path.open()) as f:
        for line in f:
            if line.startswith("{"):
                yield json.loads(line)


def compress(src,dest):
    with src.open("rb") as f,gzip.open(dest,"wb",compresslevel=6) as z:
        shutil.copyfileobj(f,z)


def run(work,out):
    trained=json.loads((HERE/"c01/predictor.json").read_text())
    previous=[r for r in read_trace(P7/"c01/original1024/trace.jsonl.gz") if r.get("stage")=="reload"]
    baseline=[r for r in read_trace(work/"baseline.trace") if r.get("stage")=="reload"]
    assert len(baseline)==len(previous)==1787
    for r,s in zip(previous,baseline,strict=True):
        assert (r["from"],r["to"],r["accepted"],r["reason"])==(s["from"],s["to"],s["accepted"],s["reason"])
    replay=[r for r in read_trace(work/"paired.trace") if "nominal" in r]
    assert len(replay)==len(trained["training_rejections"])==163
    by_reason={}
    counts=Counter()
    for r,t in zip(replay,trained["training_rejections"],strict=True):
        i=t["attempt"]
        name=f"{i:04}.json"
        assert r["file"]==name and r["unchanged"] and r["from"]==t["parent"] and r["nominal"]==t["end"]
        src=work/"parents"/name
        parent=json.loads(src.read_bytes())
        assert parent[:3]==[i,r["from"],r["nominal"]]
        cap=trained["caps"][min(255,int(256*r["from"]))]
        assert r["to"]==min(r["nominal"],r["from"]+cap) or (r["to"]==r["nominal"] and r["nominal"]-r["from"]<=cap)
        reason="+".join(t["causes"])
        group=by_reason.setdefault(reason,Counter())
        for c in (counts,group):
            c["probes"]+=1
            c["narrower"]+=int(r["to"]<r["nominal"])
            c["accepted"]+=int(r["accepted"])
            c["rejected"]+=int(not r["accepted"])
        r.update(parent_sha256=digest(src),archived_causes=t["causes"])
    with gzip.open(out/"paired-observations.jsonl.gz","wt") as f:
        for row in replay:
            print(json.dumps(row,allow_nan=False),file=f)
    with tarfile.open(out/"rejected-parents.tar.gz","w:gz",compresslevel=6) as tar:
        for p in sorted((work/"parents").glob("*.json")):
            tar.add(p,arcname=p.name,recursive=False)
    compress(work/"prefix.json",out/"reverse-end-checkpoint.json.gz")
    observed={}
    for name in ("prefix","baseline","predicted","paired"):
        shutil.copyfile(work/(name+".out"),out.parent/"checks"/(name+"-resource.log"))
        compress(work/(name+".trace"),out/(name+"-trace.jsonl.gz"))
        info=json.loads((out/(name+".json")).read_text())
        observed[name]={s["stage"]:dict(accepted=s["accepted"],attempts=s["attempts"],rejected=s["attempts"]-s["accepted"],
                progress=s["last_accepted"]) for s in info["stages"]}
        assert all(s["attempts"]<=1787 for s in info["stages"])
    # Recontact quarter-solve charges belong to the shared physical prefix.
    prefix_rows=list(read_trace(work/"prefix.trace"))
    charges={}
    for idx,(stage,row) in enumerate(observed["prefix"].items()):
        extra=sum(r.get("stage")==idx and "charge" in r for r in prefix_rows)
        charges[stage]=dict(outer=row["attempts"],extra=extra,total=row["attempts"]+extra)
        assert charges[stage]["total"]<=1787
    controlled=[r for r in read_trace(work/"predicted.trace") if r.get("stage")=="reload"]
    assert all(r["parent_unchanged"] for r in controlled if not r["accepted"])
    completed=min(baseline[-1]["from"],observed["predicted"]["reload"]["progress"])
    last_grid=int(completed*1024)/1024
    common={}
    for label,rows in (("baseline",baseline),("predicted",controlled)):
        last=next(i for i,r in enumerate(rows) if r["accepted"] and r["to"]==last_grid)
        common[label]=dict(coordinate=last_grid,attempts=last+1,rejected=sum(not r["accepted"] for r in rows[:last+1]))
    save(out/"reload-comparison.json",dict(training_set=True,parent_observations=counts,by_archived_cause=by_reason,
        frozen_predictor=digest(HERE/"c01/predictor.json"),prefix_sha256=digest(work/"prefix.json"),
        prefix_charged=charges,paths=observed,matched_progress=common,baseline_matches_p7=1787,
        same_parent_trials=len(replay),route_stage_cap=1787,mesh_cap=2065,unchanged_default=2000))


if __name__=="__main__":
    run(Path(sys.argv[1]),Path(sys.argv[2]))
