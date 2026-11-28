"""Apply P5's fixed norms at the same 448 times and 72 passive points."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P5 = HERE.parent/"2026-11-26-remaining-p5/c02"
NATIVE = HERE.parent/"2026-11-25-remaining-p3/c02/p4"


def run(out):
    dirs=dict(native=NATIVE,d4=P5/"d4",d8=P5/"d8",local1=out/"local1",focused=out/"focused",local2=out/"local2")
    fields={k:np.load(v/"principal.npz") for k,v in dirs.items()}
    w=fields["native"]["weight"]
    pos=fields["native"]["position"][:,:2]
    edge=np.minimum(pos,1-pos).min(axis=1)<.1
    pairs={}
    for left,right in (("d4","d8"),("native","d8"),("d8","local1"),("local1","local2"),
                       ("local1","focused"),("focused","local2"),("native","focused"),("d8","local2"),("native","local1"),("native","local2")):
        row={}
        for k in ("parent","weight","position"):
            np.testing.assert_array_equal(fields[left][k],fields[right][k])
        for k in ("slip","signed_area_integral"):
            weight=np.sqrt(w)[None,:,None] if k=="slip" else 1.
            a,b=(fields[t][k]*weight for t in (left,right))
            diff=(fields[left][k]-fields[right][k])*weight
            assert len(a)==len(b)==448
            rel=np.linalg.norm(diff)/np.linalg.norm(b)
            sq=diff**2
            row[k]=dict(numerator=np.linalg.norm(diff),denominator=np.linalg.norm(b),relative=rel,
                within_width=bool(rel<=.001),stage_squared_share=(sq.reshape(7,64,-1).sum(axis=(1,2))/sq.sum()).tolist())
            if k=="slip":
                row[k]["edge_squared_share"]=sq[:,edge].sum()/sq.sum()
        pairs[left+"-"+right]=row
    prev=json.loads((P5/"comparisons.json").read_text())
    for key in ("d4-d8","native-d8"):
        for k in pairs[key]:
            assert pairs[key][k]["relative"]==prev["pairs"][key][k]["relative"]
    files=[v/"principal.npz" for v in dirs.values()]+[NATIVE/"increments.npz"]
    hashes={str(x.relative_to(HERE.parents[1])) if x.is_relative_to(HERE.parents[1]) else str(x):
            hashlib.sha256(x.read_bytes()).hexdigest() for x in files}
    (out/"summary.json").write_text(json.dumps(dict(width=.001,pairs=pairs,sha256=hashes),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve())
