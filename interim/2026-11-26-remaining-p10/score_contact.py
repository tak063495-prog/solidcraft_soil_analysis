"""Reconstruct exact accepted endpoints, four-way work, and the original 22-field gate."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
P7 = HERE.parent/"2026-11-26-remaining-p7"
P8 = HERE.parent/"2026-11-26-remaining-p8"
P9 = HERE.parent/"2026-11-26-remaining-p9"
sys.path.insert(0,str(P7))
import contact
from contact import observed, STAGES, CATEGORIES, History


def read(path):
    with gzip.open(path,"rt") as f:
        return [json.loads(l) for l in f if l.startswith("{")]


def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")


def reload(folder):
    series={}
    for name in ("baseline","predicted","windows","narrow"):
        z=read((folder if name=="narrow" else P9/"c01/reload")/(name+"-trace.jsonl.gz"))
        series[name]=[r for r in z if r.get("stage")=="reload"]
    new=series["narrow"]
    table=json.loads((folder/"plan.json").read_text())
    widths=json.loads((P8/"c01/predictor.json").read_text())["caps"]
    bands=table["windows"]
    inc=1/1024
    inside=short=0
    for r in new:
        p=r["from"]
        step=(math.floor(p*1024)+1)/1024
        if step<=p+16*np.finfo(float).eps:
            step+=1/1024
        normal=min(1.,step,p+inc)
        inband=any(a<=p<b for a,b in bands)
        cap=widths[int(p*256)]
        target=min(normal,p+cap) if inband and normal-p>cap else normal
        assert target==r["to"]
        inside+=inband
        short+=target<normal
        au=r["audit"]
        inc=au["proposed_increment"] if au and au["proposed_increment"] is not None else (target-p)/2
    def report(path):
        assert all(r["parent_unchanged"] for r in path if not r["accepted"])
        for r in path:
            if r["accepted"]:
                a=r["audit"]
                assert a["event"] is None and a["max_slip"]<=a["slip_budget"] and a["max_work"]<=a["work_budget"]
                assert r["force"]<=1e-10
        return dict(attempts=len(path),accepted=sum(r["accepted"] for r in path),rejected=sum(not r["accepted"] for r in path),
                    end=max(r["to"] for r in path if r["accepted"]))
    shared=set.intersection(*[{r["to"] for r in rs if r["accepted"] and (r["to"]*1024).is_integer()} for rs in series.values()])
    p9t=json.loads((P9/"c01/reload/comparison.json").read_text())["matched"]["t"]
    def matched(t):
        counts={k:report(v[:next(i for i,r in enumerate(v) if r["accepted"] and r["to"]==t)+1]) for k,v in series.items()}
        for row in counts.values():
            row["extra_accepted"]=row["accepted"]-counts["baseline"]["accepted"]
            row["rejections_saved"]=counts["baseline"]["rejected"]-row["rejected"]
            row["net_extra_trials"]=row["extra_accepted"]-row["rejections_saved"]
        return dict(t=t,counts=counts)
    report4={k:report(v) for k,v in series.items()}
    assert max(r["attempts"] for r in report4.values())<=1787
    save(folder/"comparison.json",dict(full_budget=report4,common=matched(max(shared)),p9_coordinate=matched(p9t) if p9t in shared else None,
        new_window_parents=inside,new_shortened=short,point_events=len(table["events"]),distinct_windows=len(bands),
        controlled_span=sum(b-a for a,b in bands),previous_span=len(json.loads((P9/"windows.json").read_text())["cells"])/256))


def alternate(raw,folder):
    old=json.loads((P7/"c01/alternate-score.json").read_text())
    plan=json.loads((folder/"plan.json").read_text())
    histories=[]
    contacts=[]
    prefix_diffs={}
    for n in (256,512,1024):
        path=raw/f"alternate{n}"
        dest=folder/str(n)
        dest.mkdir(exist_ok=True)
        result=json.loads((path/"result.json").read_text())
        h=History(path/"fields","native",256,result["stages"])
        histories.append(h)
        for file in ("trace.jsonl.gz","charges.json","result.json","input.json"):
            shutil.copyfile(path/file,dest/file)
        shutil.copyfile(path/"fields/manifest.json",dest/"manifest.json")
        rows=read(dest/"trace.jsonl.gz")
        cand={(r["contact_span"][0],r["contact_span"][2]):r for r in rows if "contact_span" in r}
        accepted={}
        for r in rows:
            if r.get("stage")=="recontact" and r["accepted"]:
                c=cand[r["from"],r["to"]]
                accepted[r["to"]]=c["points"][-1]
        assert set(plan["knots"])<=accepted.keys() and 1. in accepted
        contacts.append(accepted)
        # Store each changed field stream losslessly; first three stages are
        # already present in P7. Verify identity before referencing those bytes.
        prefix_diffs[n]={}
        with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as common, zipfile.ZipFile(dest/"prefix-xor.zip","w",zipfile.ZIP_DEFLATED) as delta, zipfile.ZipFile(dest/"fields-after-open.zip","w",zipfile.ZIP_DEFLATED,compresslevel=6) as tail:
            for name,meta in h.meta["fields"].items():
                f=name+".bin"
                count=3*256*math.prod(meta["shape"])*np.dtype(meta["dtype"]).itemsize
                with (path/"fields"/f).open("rb") as src,common.open(f) as prev:
                    new,prior=src.read(count),prev.read(count)
                    if new!=prior:
                        xor=np.bitwise_xor(np.frombuffer(new,"u1"),np.frombuffer(prior,"u1"))
                        delta.writestr(f,xor.tobytes())
                        x,y=[np.frombuffer(v,meta["dtype"]) for v in (new,prior)]
                        prefix_diffs[n][name]=dict(max_absolute=float(np.max(abs(x-y))),
                                                  changed_components=int((x!=y).sum()))

                    tail.writestr(f,src.read())
    save(folder/"prefix-differences.json",prefix_diffs)
    records=[]
    for t in [*plan["knots"],1.]:
        detail=[]
        for j,point in enumerate((3,5,13)):
            data=[c[t][j] for c in contacts]
            v=np.array([p[1] for p in data])
            a,b=v[0]-v[1],v[1]-v[2]
            detail.append(dict(point=point,status=[p[0] for p in data],slip=v.tolist(),
                coarse_medium=a.tolist(),medium_fine=b.tolist(),components=(a*b).tolist(),dot=float(a@b)))
        records.append(dict(t=t,points=detail))
    save(folder/"matched.json",records)
    endpoint=np.stack([h.grid("plastic_slip",256)[4*256-1].reshape(32,3) for h in histories])
    old_endpoint=np.load(P8/"c01/alternate/signed.npz")["endpoint"]
    all_points=[]
    for i in range(32):
        item=dict(point=i)
        for name,e in (("p10",endpoint),("previous",old_endpoint)):
            a,b=-np.diff(e[:,i],axis=0)
            item[name]=dict(slip=e[:,i].tolist(),dot=float(a@b),components=(a*b).tolist())
        all_points.append(item)
    for row in records[-1]["points"]:
        np.testing.assert_array_equal(row["slip"],endpoint[:,row["point"]])
    save(folder/"endpoint-32.json",all_points)
    fields={}
    times=min(h.count for h in histories)
    def score(name,v):
        if name in CATEGORIES:
            a,b,c=v
            counts=[int((a!=b).sum()),int((b!=c).sum())]
            return dict(mismatch=counts,passes=counts==[0,0])
        m=observed(*v)
        return dict(**m,passes=m["trend"] in ("identical","decreasing_aligned"))
    for name in histories[0].meta["fields"]:
        v=[h.grid(name,256)[:times] for h in histories]
        stages={}
        for j,stage in enumerate(STAGES):
            subset=[a[256*j:256*(j+1)] for a in v]
            if len(subset[0]):
                stages[stage]=dict(times=len(subset[0]),**score(name,subset))
        fields[name]=dict(full=score(name,v),by_stage=stages,prior_passed=old["fields"][name]["passes_on_common_coordinates"])
    assert len(fields)==22
    save(folder/"score.json",dict(common_coordinates=times,counts=[h.count for h in histories],fields=fields,
        by_stage={s:sum(f["by_stage"].get(s,{}).get("passes",False) for f in fields.values()) for s in STAGES},
        passed_fields=sum(f["full"]["passes"] for f in fields.values()),prior_passed=sum(f["prior_passed"] for f in fields.values()),
        full_domain=times==1792,passed=times==1792 and all(f["full"]["passes"] for f in fields.values())))


if __name__=="__main__":
    if sys.argv[1]=="reload":
        reload(Path(sys.argv[2]))
    else:
        alternate(Path(sys.argv[2]),Path(sys.argv[3]))
