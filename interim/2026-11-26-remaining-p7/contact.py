"""Run each path from its case file and rescore only measured common accepted coordinates."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
P4 = HERE.parent/"2026-11-25-remaining-p4"
P3 = HERE.parent/"2026-11-25-remaining-p3"
P2 = HERE.parent/"2026-11-25-remaining5-p2"
ROOT = HERE.parents[1]
sys.path.insert(0,str(P2))
from c01_native import Fields
from c01_score import CATEGORIES, STAGES, observed
from c01_retry_metrics import History


def save(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def run(binary,folder,n,route):
    folder.mkdir(parents=True)
    if route=="original":
        case = json.loads((P4/f"c01/original{n}-input.json").read_text())
    else:
        case = json.loads((P3/f"c01/fixed/case-{n}.json").read_text())
    assert case["joint_system"]["global_step_control"]["path_error"]["slip_scale_m"]/n**2 == 1e-4/1024**2
    save(folder/"input.json",case)
    save(folder/"options.json",dict(full_field_grid=256))
    cmd = [binary,ROOT/"examples/joint-water-column-tet10-mesh.json",ROOT/"examples/ground-joint-water-mc.json",
           folder/"input.json",folder/"result.json","--methods",folder/"options.json"]
    with (folder/"trace.jsonl").open("xb") as log, Fields(folder/"fields",256,"native") as fields, subprocess.Popen(
            cmd,stdout=subprocess.PIPE,stderr=log,env=dict(os.environ,P7_SWITCH_HALF="1")) as proc:
        for line in proc.stdout:
            fields.native(*json.loads(line))
        if proc.wait():
            raise RuntimeError(proc.returncode)
    result = json.loads((folder/"result.json").read_text())
    rows = list(map(json.loads,(folder/"trace.jsonl").read_text().splitlines()))
    charges = {}
    for k,stage in enumerate(result["stages"]):
        name = stage["stage"]
        actual = [x for x in rows if x.get("stage")==name and "accepted" in x]
        extra = [x for x in rows if x.get("stage")==k and "charge" in x]
        assert len(actual)==stage["attempts"] and sum(x["accepted"] for x in actual)==stage["accepted"]
        assert all(x["parent_unchanged"] for x in actual if not x["accepted"])
        assert len(actual)+len(extra)<=1787<=2065
        charges[name] = dict(outer=len(actual),extra=len(extra),total=len(actual)+len(extra),
            extra_rebalances=sum(x["charge"]=="quarter_rebalance" for x in extra))
    save(folder/"charges.json",charges)
    with (folder/"trace.jsonl").open("rb") as src,gzip.open(folder/"trace.jsonl.gz","wb") as dest:
        shutil.copyfileobj(src,dest)
    (folder/"trace.jsonl").unlink()


def diagnose(folder):
    from collections import Counter
    with gzip.open(folder/"trace.jsonl.gz","rt") as stream:
        rows = list(map(json.loads,stream))
    half = [r for r in rows if "half" in r]
    extra_rejected = {(r["from"],r["to"]) for r in half if r["outside"]}
    reasons = {}
    for stage in STAGES:
        counts = Counter()
        for r in rows:
            if r.get("stage")!=stage or r["accepted"]:
                continue
            a = r["audit"] or {}
            why = []
            if a.get("event") is not None:
                why.append("friction_event")
            if a and a["max_slip"]>a["slip_budget"]:
                why.append("slip")
            if a and a["max_work"]>a["work_budget"]:
                why.append("dissipation")
            if stage=="recontact" and (r["from"],r["to"]) in extra_rejected:
                why.append("switch_half")
            counts["+".join(why) or r["reason"]] += 1
        reasons[stage] = counts
    return dict(candidate_halves=len(half),outside_halves=sum(r["outside"] for r in half),
                additional_rejected_attempts=len(extra_rejected),rejections_by_combination=reasons)


def archive_fields(folder, dest, common, manifest):
    """Share the identical first four stages losslessly; verify the second route byte for byte."""
    create = not common.exists()
    with zipfile.ZipFile(common,"x" if create else "r",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as prefix, \
         zipfile.ZipFile(dest/"fields-tail.zip","w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as tail:
        for name,f in manifest["fields"].items():
            filename = name+".bin"
            width = int(np.prod(f["shape"]))*np.dtype(f["dtype"]).itemsize
            with (folder/"fields"/filename).open("rb") as stream:
                data = stream.read(4*256*width)
                assert len(data)==4*256*width
                if create:
                    prefix.writestr(filename,data)
                else:
                    assert prefix.read(filename)==data
                tail.writestr(filename,stream.read())


def score(raw,out,routes=("original","alternate")):
    out.mkdir(exist_ok=True,parents=True)
    prior_orig = json.loads((P4/"c01/summary.json").read_text())["original_fixed_fields"]
    prior_alt = {k:v["common"] for k,v in json.loads((P3/"c01/fixed/summary.json").read_text())["fields"].items()}
    for route in routes:
        prior = {"original":prior_orig,"alternate":prior_alt}[route]
        histories = []
        for n in (256,512,1024):
            folder = raw/f"{route}{n}"
            dest = out/f"{route}{n}"
            dest.mkdir(exist_ok=True)
            result = json.loads((folder/"result.json").read_text())
            h = History(folder/"fields","native",256,result["stages"])
            histories.append(h)
            for file in ("input.json","options.json","result.json","charges.json","trace.jsonl.gz"):
                shutil.copyfile(folder/file,dest/file)
            shutil.copyfile(folder/"fields/manifest.json",dest/"manifest.json")
            save(dest/"cutbacks.json",diagnose(folder))
            archive_fields(folder,dest,out/f"first-four-n{n}.zip",h.meta)
        times = min(h.count for h in histories)
        pairs = {}
        for name in histories[0].meta["fields"]:
            a,b,c = (h.grid(name,256)[:times] for h in histories)
            if name in CATEGORIES:
                counts = [int(np.count_nonzero(a!=b)),int(np.count_nonzero(b!=c))]
                metric = dict(mismatch=counts)
                ok = counts==[0,0]
            else:
                metric = observed(a,b,c)
                ok = metric["trend"] in ("decreasing_aligned","identical")
            by_stage = {}
            for i,stage in enumerate(STAGES):
                sl = slice(i*256,(i+1)*256)
                if not len(a[sl]):
                    continue
                if name in CATEGORIES:
                    local = dict(mismatch=[int((a[sl]!=b[sl]).sum()),int((b[sl]!=c[sl]).sum())])
                else:
                    local = observed(a[sl],b[sl],c[sl])
                by_stage[stage] = dict(coordinates=len(a[sl]),**local)
            pairs[name]=dict(measured=metric,by_stage=by_stage,passes_on_common_coordinates=ok,
                             archived_full_domain=prior[name])
        assert len(pairs)==22
        save(out/f"{route}-score.json",dict(stages=STAGES,common_coordinates=times,
             counts=[h.count for h in histories],full_domain=times==1792,
             finite_three_level_gate=times==1792 and all(v["passes_on_common_coordinates"] for v in pairs.values()),
             fields=pairs))


if __name__=="__main__":
    if sys.argv[1]=="score":
        score(Path(sys.argv[2]),Path(sys.argv[3]))
    else:
        run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),int(sys.argv[3]),sys.argv[4])
