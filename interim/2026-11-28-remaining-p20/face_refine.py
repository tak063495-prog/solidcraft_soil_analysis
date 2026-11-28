"""Recertify the p19 quadrature with one further four-way split of the ranked cohort."""
import csv
import gzip
import importlib.util
import inspect
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P19 = HERE.with_name("2026-11-28-remaining-p19")
spec = importlib.util.spec_from_file_location("p19face", P19/"face.py")
p19 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p19)
previous = p19.previous


def changed(source, old, new):
    assert source.count(old)==1, old
    return source.replace(old,new)


def run(sign, out):
    with gzip.open(P19/f"c02/baseline/{sign}-833-by-leaf.csv.gz","rt") as stream:
        rows=list(csv.DictReader(stream))
    assert all(float(a["total_nm"])>=float(b["total_nm"]) for a,b in zip(rows,rows[1:]))
    keys=lambda seq:[(int(r["block"]),int(r["leaf"])) for r in seq]
    ranked=keys(r for r in rows if (r["refined"]=="True")==(sign=="plus"))
    assert len(ranked)=={"minus":16047,"plus":4230}[sign]
    chosen=set(ranked)
    # The 16-cell nominal integral is the common starting point in p18.
    # Retain the minus top 328 at 64 cells, then take the ranked complementary
    # leaves to 64. The previously refined positive group goes directly to 256.
    retained=keys(r for r in rows if r["refined"]=="True") if sign=="minus" else []
    source=inspect.getsource(previous.integrate)
    anchor='    with np.load(out) as data:'
    source=changed(source,anchor,'''    if int(step)==833:
        index = {key:i for i,key in enumerate(keys)}
        selected = np.array([index[k] for k in retained+ranked if k in index],dtype=int)
        cohort = len(selected)
        if sign=="plus":
            finer=(mats[:,None] @ finer[None]).reshape(-1,15,15)
'''+anchor)
    scope=dict(vars(previous),negative_group=328,ranked=ranked,retained=retained)
    exec(compile(source,str(HERE/"face_refine.py"),"exec"),scope)
    previous.p15.integrate=scope["integrate"]
    # Keep every support and guard of the predecessor. Only the restriction
    # matrices of the selected history leaves change, in the same paired integral.
    source=inspect.getsource(p19.measured)
    anchor='    replace("    for block, mesh'
    pos=source.index(anchor)
    source=source[:pos]+'''    if int(dest.stem.split("-")[1])==833:
        replace("    norm_a =", "    finer = (mats[:,None] @ restrict[None]).reshape(-1,15,15)\\n    norm_a =")
        replace("        for k in active:", "        active = sorted(active,key=lambda k: ranks.get((block_id,k),len(ranks)))\\n        for k in active:")
        replace("cells = restrict @ weights[k]", "cells = (finer if (block_id,k) in ranks else restrict) @ weights[k]")
'''+source[pos:]
    source=changed(source,'scope = dict(vars(previous.p16))','scope = dict(vars(previous.p16), ranks=ranks)')
    scope=dict(vars(p19),ranks={key:i for i,key in enumerate(ranked)})
    exec(compile(source,str(HERE/"face_refine.py"),"exec"),scope)
    previous.p16.paired=scope["measured"]
    previous.p16.run(sign,out,"projected")


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
