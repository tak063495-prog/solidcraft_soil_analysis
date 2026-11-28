"""Compare only shared accepted loads, retaining exactly the p10 22-field rule."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
P7=HERE.parent/"2026-11-26-remaining-p7"
P10=HERE.parent/"2026-11-26-remaining-p10/c01/alternate"
sys.path.insert(0,str(P7))
from contact import CATEGORIES,History,STAGES,observed

def save(path,x):
    path.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")

def p10_bytes(n,meta,field,common,xor,tail):
    f=field+".bin"
    size=3*256*int(np.prod(meta["shape"]))*np.dtype(meta["dtype"]).itemsize
    with common.open(f) as stream:
        prefix=stream.read(size)
    if f in xor.namelist():
        prefix=np.bitwise_xor(np.frombuffer(prefix,"u1"),np.frombuffer(xor.read(f),"u1")).tobytes()
    data=prefix+tail.read(f)
    assert hashlib.sha256(data).hexdigest()==meta["sha256"]
    return data

def score(raw,out):
    histories=[]
    old_slip=[]
    for n in (256,512,1024):
        folder=raw/f"alternate{n}"
        sub=out/str(n)
        sub.mkdir(exist_ok=True)
        for f in ("charges.json","result.json","trace.jsonl.gz","input.json"):
            shutil.copyfile(folder/f,sub/f)
        shutil.copyfile(folder/"fields/manifest.json",sub/"manifest.json")
        h=History(folder/"fields","native",256,json.loads((sub/"result.json").read_text())["stages"])
        histories.append(h)
        baseline=json.loads((P10/str(n)/"manifest.json").read_text())["fields"]
        with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as common, zipfile.ZipFile(P10/str(n)/"prefix-xor.zip") as xor, \
                zipfile.ZipFile(P10/str(n)/"fields-after-open.zip") as tail, zipfile.ZipFile(sub/"fields-xor-p10.zip","w",zipfile.ZIP_DEFLATED) as delta:
            for name,m in baseline.items():
                old=p10_bytes(n,m,name,common,xor,tail)
                if name=="plastic_slip":
                    old_slip.append(np.frombuffer(old,m["dtype"]).reshape(1792,32,3))
                new=(folder/"fields"/(name+".bin")).read_bytes()
                assert len(new)<=len(old)
                old=np.frombuffer(old,"u1")[:len(new)]
                x=np.bitwise_xor(old,np.frombuffer(new,"u1"))
                if x.any():
                    delta.writestr(name+".bin",x.tobytes())
                assert hashlib.sha256(np.bitwise_xor(old,x)).hexdigest()==h.meta["fields"][name]["sha256"]
    times=min(h.count for h in histories)
    def metric(name,values):
        if name in CATEGORIES:
            a,b,c=values
            counts=[int((a!=b).sum()),int((b!=c).sum())]
            return dict(mismatch=counts,passes=counts==[0,0])
        r=observed(*values)
        return dict(**r,passes=r["trend"] in ("identical","decreasing_aligned"))
    fields={}
    for name in histories[0].meta["fields"]:
        vals=[h.grid(name,256)[:times] for h in histories]
        fields[name]=dict(full=metric(name,vals),stages={s:metric(name,[a[j*256:(j+1)*256] for a in vals])
            for j,s in enumerate(STAGES) if j*256<times})
    assert len(fields)==22
    score=dict(common_coordinates=times,passed_fields=sum(r["full"]["passes"] for r in fields.values()),
        by_stage={s:sum(r["stages"][s]["passes"] for r in fields.values()) for j,s in enumerate(STAGES) if j*256<times},
        fields=fields)
    score["passed"]=times==1792 and score["passed_fields"]==22
    save(out/"score.json",score)
    # Signed point contributions at every common load and sums through each later stage.
    v=np.stack([h.grid("plastic_slip",256)[:times].reshape(times,32,3) for h in histories])
    def dots(x):
        a,b=-np.diff(x,axis=0)
        return (a*b).sum(-1)
    current=dots(v)
    before=dots(np.stack(old_slip)[:,:times])
    np.savez_compressed(out/"signed-32.npz",progress=np.tile(np.arange(1,257)/256,7)[:times],p10=before,p11=current)
    rows={}
    for j,s in enumerate(STAGES):
        sl=slice(j*256,min((j+1)*256,times))
        if sl.start==sl.stop or sl.start>=times: continue
        rows[s]={}
        for label,x in (("p10",before),("p11",current)):
            y=x[sl]
            rows[s][label]=dict(coordinates=len(y),sum_m2=float(y.sum()),by_point_m2=y.sum(0).tolist(),
                endpoint=y[-1].tolist(),endpoint_total=float(y[-1].sum()),endpoint_negative=int((y[-1]<0).sum()))
    save(out/"contributions.json",rows)

if __name__=="__main__":
    score(Path(sys.argv[1]),Path(sys.argv[2]))
