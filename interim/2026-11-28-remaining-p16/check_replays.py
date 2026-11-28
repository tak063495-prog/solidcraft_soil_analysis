"""Verify each observation-only replay's charged attempts and retained failure records."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")
P15 = HERE.with_name("2026-11-28-remaining-p15")
STAGES = ("geostatic","shear","open","recontact","reverse","reload","return")


def audit(work,out):
    runs = {}
    for route,p in (("outer",P15),("middle",P14)):
        for n in (256,512,1024):
            key = f"{route}-{n}"
            raw = work/(key+".trace")
            src = p/"c01/alternate"/route/str(n)
            with raw.open("rb") as a,gzip.open(src/"trace.jsonl.gz","rb") as b:
                measured,previous = [hashlib.file_digest(f,"sha256").hexdigest() for f in (a,b)]
            if measured!=previous:
                with raw.open("rb") as a,gzip.open(out/(key+"-replay.jsonl.gz"),"wb") as b:
                    shutil.copyfileobj(a,b)
            rows = {s:dict(outer=0,rejected=0,extra=0) for s in STAGES}
            with raw.open() as f:
                for line in f:
                    if not line.startswith("{"):
                        continue
                    r = json.loads(line)
                    s = r.get("stage")
                    if s in rows and "accepted" in r:
                        rows[s]["outer"]+=1
                        if not r["accepted"]:
                            assert r["parent_unchanged"]
                            rows[s]["rejected"]+=1
                    if "charge" in r:
                        rows[STAGES[s]]["extra"]+=1
            charges = json.loads((src/"charges.json").read_text())
            for s,row in rows.items():
                row["total"]=row["extra"]+row["outer"]
                assert row["total"] <= 1787 <= (2065*12+11)//12
                assert all(row[k]==charges[s][k] for k in ("outer","extra","total"))
            runs[key] = dict(sha256=measured,archived_sha256=previous,trace_identical=measured==previous,
                archive=str((src/"trace.jsonl.gz").relative_to(HERE.parents[1])),stages=rows)
    (out/"replay-cost.json").write_text(json.dumps(runs,indent=2)+"\n")


if __name__=="__main__":
    audit(*map(Path,sys.argv[1:]))
