import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"2026-11-26-remaining-p9"))
from prefix_compare import check
new,trace,out=map(Path,sys.argv[1:])
check(HERE.parent/"2026-11-26-remaining-p10/c01/reload/reverse-end.json.gz",new,trace,out)
row=json.loads(out.read_text())
row["p10_checkpoint_sha256"]=row.pop("p8_checkpoint_sha256")
row["p11_checkpoint_sha256"]=row.pop("p9_checkpoint_sha256")
out.write_text(json.dumps(row,indent=2)+"\n")
