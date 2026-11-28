"""Repeat all field comparisons using the same extraction environment at the three levels."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile

import numpy as np

from original import P7, STAGES, field_bytes, measure, save


def run(work,out):
    base=json.loads((out/"gates.json").read_text())
    info={}
    for n in (256,512):
        raw=work/f"recheck/original{n}"
        archived=P7/f"c01/original{n}"
        dest=out/f"same-environment/{n}"
        dest.mkdir(parents=True,exist_ok=True)
        manifest=json.loads((raw/"fields/manifest.json").read_text())
        assert manifest["count"]==1792
        for name in ("result.json","charges.json"):
            assert json.loads((raw/name).read_text())==json.loads((archived/name).read_text())
            shutil.copyfile(raw/name,dest/name)
        shutil.copyfile(raw/"fields/manifest.json",dest/"manifest.json")
        shutil.copyfile(raw/"trace.jsonl.gz",dest/"trace.jsonl.gz")
        differences={}
        with zipfile.ZipFile(dest/"fields-xor-p7-original.zip","w",zipfile.ZIP_DEFLATED) as z:
            for k,m in manifest["fields"].items():
                data=(raw/f"fields/{k}.bin").read_bytes()
                assert hashlib.sha256(data).hexdigest()==m["sha256"]
                prior=field_bytes(n,"original",k)
                x=np.frombuffer(data,"u1")^np.frombuffer(prior,"u1")
                if x.any():
                    z.writestr(k+".bin",x.tobytes())
                    a,b=(np.frombuffer(s,m["dtype"]) for s in (data,prior))
                    differences[k]=dict(count=int(np.sum(a!=b)),max_abs=float(np.abs(a-b).max()))
        info[str(n)]=differences
    fields={}
    for k,m in json.loads((out/"manifest.json").read_text())["fields"].items():
        vals=[np.fromfile(d/f"{k}.bin",m["dtype"]).reshape(1792,*m["shape"])
              for d in (work/"recheck/original256/fields",work/"recheck/original512/fields",work/"fields")]
        fields[k]=dict(full=measure(k,vals),stages={s:measure(k,[v[j*256:(j+1)*256] for v in vals])
            for j,s in enumerate(STAGES)})
    count=sum(row["full"]["passes"] for row in fields.values())
    save(out/"same-environment/gates.json",dict(common_coordinates=1792,levels=[256,512,1024],
        passed_fields=count,passed=count==22,fields=fields,extraction_differences=info,
        changed_full_verdicts=[k for k in fields if fields[k]["full"]["passes"]!=base["fields"][k]["full"]["passes"]],
        changed_stage_verdicts=[[k,s] for k in fields for s in STAGES
            if fields[k]["stages"][s]["passes"]!=base["fields"][k]["stages"][s]["passes"]]))


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
