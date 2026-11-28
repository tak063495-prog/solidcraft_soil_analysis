"""Require independent nominal integration and propagated histories to be byte-identical."""
import hashlib
import json
from pathlib import Path
import sys


def check(original,replay,dest):
    names={p.name for p in original.glob("*.npz")}|{"plus.json","minus.json"}
    assert names=={p.name for p in replay.iterdir()}
    hashes={}
    for name in sorted(names):
        raw=(replay/name).read_bytes()
        assert raw==(original/name).read_bytes(),name
        hashes[name]=hashlib.sha256(raw).hexdigest()
    dest.write_text(json.dumps(hashes,indent=2)+"\n")


if __name__=="__main__":
    check(*map(Path,sys.argv[1:]))
