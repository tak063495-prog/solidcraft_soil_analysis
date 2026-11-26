"""Verify a regenerated prefix without modifying either checkpoint."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

from collections import Counter


def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def check(old_path,new_path,trace_path,out):
    raw=gzip.decompress(old_path.read_bytes()),new_path.read_bytes()
    a,b=map(json.loads,raw)
    counts=Counter()
    def compare(x,y,path):
        if type(x) is dict:
            assert x.keys()==y.keys(),path
            for k,v in x.items(): compare(v,y[k],path+(k,))
        elif type(x) is list:
            assert len(x)==len(y),path
            for v,w in zip(x,y,strict=True): compare(v,w,path+("[]",))
        elif path[-1]=="elapsed_seconds" or path==("provenance","solver_binary_sha256"):
            if x!=y: counts[".".join(path)]+=1
        else:
            assert type(x)==type(y) and repr(x)==repr(y),path
    compare(a,b,())
    assert a["last_committed"]["stage_progress"]==1
    def canonical(value):
        return hash_bytes(json.dumps(value,sort_keys=True,separators=(",",":")).encode())
    totals=Counter()
    for line in trace_path.read_text().splitlines():
        if not line.startswith("{"): continue
        row=json.loads(line)
        if "charge" in row:
            totals[("geostatic","shear","open","recontact","reverse")[row["stage"]]]+=1
        if "accepted" in row: totals[row["stage"]]+=1
    assert max(totals.values())<=1787
    record=dict(p8_checkpoint_sha256=hash_bytes(raw[0]),p9_checkpoint_sha256=hash_bytes(raw[1]),
                solver_sha256=[v["provenance"]["solver_binary_sha256"] for v in (a,b)],
                common_canonical_sha256={k:canonical(a[k]) for k in ("last_committed","continuation")},
                different_leaves=dict(counts),prefix_charges=dict(totals))
    out.write_text(json.dumps(record,indent=2)+"\n")


if __name__=="__main__":
    check(*map(Path,sys.argv[1:]))
