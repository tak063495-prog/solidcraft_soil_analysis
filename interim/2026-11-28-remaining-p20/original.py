"""Compare p19's fine control route to the completed original 256/512 routes."""
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
P19=HERE.with_name("2026-11-28-remaining-p19")
sys.path.insert(0,str(HERE.with_name("2026-11-28-remaining-p16")))
from gates import measure
from attribution import P7, STAGES
from contact import Fields, ROOT, save


def field_bytes(n, route, name):
    key=name+".bin"
    folder=P7/f"c01/{route}{n}"
    with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as a, zipfile.ZipFile(folder/"fields-tail.zip") as b:
        raw=a.read(key)+b.read(key)
    meta=json.loads((folder/"manifest.json").read_text())["fields"][name]
    assert hashlib.sha256(raw).hexdigest()==meta["sha256"]
    return raw


def solve(work,out):
    case=P7/"c01/original1024"
    cmd=[work/"probe/target/release/examples/c01_methods",ROOT/"examples/joint-water-column-tet10-mesh.json",
         ROOT/"examples/ground-joint-water-mc.json",case/"input.json",out/"result.json","--methods",case/"options.json"]
    clean={k:v for k,v in os.environ.items() if not re.match(r"P[0-9]+_",k)}
    env=dict(clean,P7_SWITCH_HALF="1",P10_LABELS="1",P14_GAP="1",P15_RELOAD="1",P16_SLIP="1",
        P17_WORK="1",P18_OBSERVE="1",P18_CONTROL="1",P19_CONTROL="1",P20_SLIP_OBSERVE="1",
        P11_CROSS=str(HERE.with_name("2026-11-27-remaining-p12")/"c01/reload/plan.json"))
    with (work/"original.trace").open("xb") as trace, Fields(work/"fields",256,"native") as sink, \
         subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=trace,env=env) as proc:
        for line in proc.stdout:
            sink.native(*json.loads(line))
        if proc.wait():
            raise RuntimeError(proc.returncode)
    assert sink.count==1792


def score(work,out):
    manifest=json.loads((work/"fields/manifest.json").read_text())
    fields={}
    prefix_digests={}
    with zipfile.ZipFile(out/"fields-xor-p7-alternate1024.zip","w",zipfile.ZIP_DEFLATED) as archive:
        for key,m in manifest["fields"].items():
            raw=(work/f"fields/{key}.bin").read_bytes()
            assert hashlib.sha256(raw).hexdigest()==m["sha256"]
            prefix_len=5*256*int(np.prod(m["shape"]))*np.dtype(m["dtype"]).itemsize
            prefix=raw[:prefix_len]
            archived=field_bytes(1024,"original",key)[:prefix_len]
            current,old=(np.frombuffer(b,m["dtype"]) for b in (prefix,archived))
            prefix_digests[key]=dict(sha256=hashlib.sha256(prefix).hexdigest(),
                archived_sha256=hashlib.sha256(archived).hexdigest(),different=int((current!=old).sum()),
                max_absolute_difference=float(np.max(np.abs(current.astype(float)-old))))
            mask=np.frombuffer(raw,"u1")^np.frombuffer(field_bytes(1024,"alternate",key),"u1")
            if mask.any():
                archive.writestr(key+".bin",mask.tobytes())
            vals=[np.frombuffer(b,m["dtype"]).reshape(1792,*m["shape"]) for b in
                  [field_bytes(256,"original",key),field_bytes(512,"original",key),raw]]
            fields[key]=dict(full=measure(key,vals),stages={stage:measure(key,[v[j*256:(j+1)*256] for v in vals])
                for j,stage in enumerate(STAGES)})
    save(out/"first-five-sha256.json",prefix_digests)
    save(out/"manifest.json",manifest)
    count=sum(r["full"]["passes"] for r in fields.values())
    save(out/"gates.json",dict(common_coordinates=1792,levels=[256,512,1024],passed_fields=count,
                               passed=count==22,fields=fields))


def slip(work,out):
    with (work/"original.trace").open() as stream:
        rows=[json.loads(s) for s in stream if s.startswith("{")]
    existing=json.loads((P19/"c01/reload/control/result.json").read_text())
    actual=json.loads((out/"result.json").read_text())
    assert actual==existing
    charges={}
    for i,s in enumerate(actual["stages"]):
        stage=s["stage"]
        notes=[r for r in rows if r.get("stage")==stage and "accepted" in r]
        assert len(notes)==s["attempts"] and sum(r["accepted"] for r in notes)==s["accepted"]
        assert all(r["parent_unchanged"] for r in notes if not r["accepted"])
        charges[stage]=len(notes)+sum("charge" in r and r["stage"]==i for r in rows)
    assert max(charges.values())<=1787<=2065
    with gzip.open(P19/"c01/reload/control/trace.jsonl.gz","rt") as f:
        ref=[json.loads(x) for x in f if x.startswith("{")]
    trials=lambda data:[r for r in data if "accepted" in r and r.get("stage") in ("reload","return")]
    assert trials(rows)==trials(ref)
    with gzip.open(out/"observations.jsonl.gz","wt") as f:
        for r in rows:
            if "p20_predict" in r or "p20_actual" in r:
                print(json.dumps(r,separators=(",",":")),file=f)
    matched=[]
    # Forecasts and material observations are emitted during each trial.
    # The archived global verdict list is emitted after the solver returns.
    live=[r for r in rows if "p20_predict" in r or "p20_actual" in r]
    verdicts=[(i,r) for i,r in enumerate((r for r in trials(rows) if r["stage"]=="reload"),1)
              if .5625<=r["from"]<.625]
    for predict,observed,(i,r) in zip(live[::2],live[1::2],verdicts,strict=True):
        assert predict["p20_predict"]==i and "p20_actual" in observed
        assert predict["parent"]==observed["p20_actual"]==r["from"] and predict["to"]==observed["to"]==r["to"]
        threshold=predict["budget_m"]
        assert observed["budget_m"]==threshold==r["audit"]["slip_budget"]
        p,a=[np.abs(v["error_m"]).max(axis=1) for v in (predict,observed)]
        assert max(a)==r["audit"]["max_slip"]
        matched.append(dict(attempt=predict["p20_predict"],parent=r["from"],to=r["to"],
            accepted=r["accepted"],predicted=np.flatnonzero(p>threshold).tolist(),
            actual=np.flatnonzero(a>threshold).tolist(),predicted_argmax=int(p.argmax()),
            actual_argmax=int(a.argmax()),predicted_ratio=(p/threshold).tolist(),actual_ratio=(a/threshold).tolist()))
    assert len(matched)==101 and sum(not r["accepted"] for r in matched)==15
    save(out/"comparison.json",dict(charges=charges,matched_p19_reload_return_trials=True,
        guard=1/65536,points=32,slip_budget_m=threshold,trials=matched))


if __name__=="__main__":
    globals()[sys.argv[1]](Path(sys.argv[2]).resolve(),Path(sys.argv[3]).resolve())
