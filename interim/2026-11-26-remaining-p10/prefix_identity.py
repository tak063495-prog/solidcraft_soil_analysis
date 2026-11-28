"""Reuse the archive's strict recursive comparison and name both checkpoints."""
import json
from pathlib import Path
import sys

p9=Path("interim/2026-11-26-remaining-p9")
sys.path.insert(0,str(p9.resolve()))
from prefix_compare import check

if __name__=="__main__":
    new,trace,out=map(Path,sys.argv[1:])
    check(p9/"c01/reload/reverse-end.json.gz",new,trace,out)
    row=json.loads(out.read_text())
    old_hash=row.pop("p8_checkpoint_sha256")
    new_hash=row.pop("p9_checkpoint_sha256")
    row.update(p9_checkpoint_sha256=old_hash,p10_checkpoint_sha256=new_hash)
    out.write_text(json.dumps(row,indent=2)+"\n")
