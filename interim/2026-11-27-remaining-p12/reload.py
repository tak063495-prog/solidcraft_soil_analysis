"""Freeze the entrance-only policy or compare three independently accepted histories."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "2026-11-26-remaining-p11/c01/reload"


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def freeze(out):
    plan=read(BASE/"plan.json")
    controls=[r for r in plan["controls"] if r["transition"]==["slip","open"]]
    assert controls and all(r["target"]==r["span"][0] for r in controls)
    save(out/"plan.json",dict(bands=plan["bands"],controls=controls,
         p11_plan_sha256=hashlib.sha256((BASE/"plan.json").read_bytes()).hexdigest()))


def verify_prefix(out, work):
    sys.path.insert(0,str(HERE.parent/"2026-11-26-remaining-p9"))
    from prefix_compare import check
    path=out/"prefix-identity.json"
    old=HERE.parent/"2026-11-26-remaining-p10/c01/reload/reverse-end.json.gz"
    check(old,work/"prefix.json",work/"prefix.trace",path)
    audit=read(path)
    prior=read(BASE/"prefix-identity.json")
    assert audit["common_canonical_sha256"]==prior["common_canonical_sha256"]
    audit["p10_checkpoint_sha256"]=audit.pop("p8_checkpoint_sha256")
    audit["p11_checkpoint_sha256"]=prior["p11_checkpoint_sha256"]
    audit["p12_checkpoint_sha256"]=audit.pop("p9_checkpoint_sha256")
    save(path,audit)


def score(out):
    plan=read(out/"plan.json")
    paths={}
    for label,file in (("ordinary",BASE/"ordinary.jsonl.gz"),("p11",BASE/"cross.jsonl.gz"),("p12",out/"entrance.jsonl.gz")):
        with gzip.open(file,"rt") as stream:
            rows=[json.loads(s) for s in stream if s.startswith("{")]
        assert not any(r.get("stage")==5 and "charge" in r for r in rows)
        path=paths[label]=[r for r in rows if r.get("stage")=="reload"]
        assert len(path)<=1787
        for r in path:
            if r["accepted"]:
                a=r["audit"]
                assert a["event"] is None and a["max_slip"]<=a["slip_budget"] and a["max_work"]<=a["work_budget"] and r["force"]<=1e-10
            else:
                assert r["parent_unchanged"]
    actions=[]
    increment=1/1024
    for r in paths["p12"]:
        parent=r["from"]
        grid=(math.floor(parent*1024)+1)/1024
        if grid<=parent+16*sys.float_info.epsilon:
            grid+=1/1024
        nominal=min(1.,parent+increment,grid)
        targets=[] if any(l<=parent<h for l,h in plan["bands"]) else [
            a["target"] for a in plan["controls"] if parent<a["span"][0] and a["span"][1]<=nominal]
        assert r["to"]==min(nominal,*targets) if targets else r["to"]==nominal
        if r["to"]<nominal:
            actions.append([parent,nominal,r["to"],r["accepted"]])
        audit=r["audit"]
        increment=audit["proposed_increment"] if audit and audit["proposed_increment"] is not None else (r["to"]-parent)/2

    def counts(path):
        a=[r for r in path if r["accepted"]]
        return dict(accepted=len(a),rejected=len(path)-len(a),attempts=len(path),progress=a[-1]["to"])
    shared=set.intersection(*[{r["to"] for r in path if r["accepted"] and (1024*r["to"]).is_integer()}
                              for path in paths.values()])
    def at(t):
        return {label:counts(path[:next(i+1 for i,r in enumerate(path) if r["accepted"] and r["to"]==t)])
                for label,path in paths.items()}
    save(out/"comparison.json",dict(actions=actions,common=at(max(shared)),p11_common=at(.900390625) if .900390625 in shared else None,
                                   budget={label:counts(path) for label,path in paths.items()}))


if __name__=="__main__":
    if sys.argv[1]=="prefix":
        verify_prefix(Path(sys.argv[2]),Path(sys.argv[3]))
    else:
        {"freeze":freeze,"score":score}[sys.argv[1]](Path(sys.argv[2]))
