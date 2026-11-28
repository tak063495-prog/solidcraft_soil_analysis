"""Score independent references and the shared discrepancy of the two archived solvers."""
from contextlib import ExitStack
from pathlib import Path
import argparse
import hashlib
import io
import json
import sys
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/"2026-11-23-impl-root-causes-c"),str(HERE.parent/"2026-12-17-c02-spatial-cause-diagnosis")]
from c02 import common_indices
from c02diag_extract import SplitZip, OLD
from c02diag_mapping import ContactPlan
from c02diag_metrics import difference
import c02diag_driver as diag


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs",type=Path)
    args=parser.parse_args()
    legacy=diag.legacy()
    index=json.loads((OLD/"archive-sha256.json").read_bytes())
    whole=hashlib.sha256()
    parts=[]
    for part in index["parts"]:
        path=OLD/part["file"]
        local=hashlib.sha256()
        assert path.stat().st_size==part["bytes"]
        with path.open("rb") as f:
            while block:=f.read(2**20):
                whole.update(block)
                local.update(block)
        assert local.hexdigest()==part["sha256"]
        parts.append(path)
    assert whole.hexdigest()==index["zip_sha256"]
    result={"archive_sha256":whole.hexdigest(),"costs":{},"time":{},"degree":{},"shared_discrepancy":{}}
    with ExitStack() as stack:
        fields={}
        states={}
        for p in (2,3):
            for n in (64,128):
                path=args.runs/f"c02-p{p}-n{n}"
                data=stack.enter_context(np.load(path/"slip.npz",allow_pickle=False))
                fields[p,n]=data["slip"]
                states[p,n]=data["status"]
                if (p,n)==(2,64):
                    area,pos,parents=data["weight"],data["position"],data["parent"]
                else:
                    for key,exact in [("weight",area),("position",pos),("parent",parents)]:
                        np.testing.assert_array_equal(data[key],exact)
                result["costs"][f"p{p}-n{n}"]=json.loads((path/"cost.json").read_bytes())
        def integral(hist,w):
            return np.einsum("tqc,q->tc",hist,w)
        def measures(a,b):
            return dict(area= difference(a,b,area[None,:,None]),integral=difference(integral(a,area),integral(b,area)))
        def contact_disagreement(first,second):
            different=states[first]!=states[second]
            return dict(disagreeing_time_points=int(different.sum()),total_time_points=different.size,
                        largest_disagreeing_area_m2=float(np.max(different@area)))
        for p in (2,3):
            result["time"][str(p)]=measures(fields[p,64],fields[p,128])
            result["time"][str(p)]["contact_classification"]=contact_disagreement((p,64),(p,128))
        for n in (64,128):
            result["degree"][str(n)]=measures(fields[2,n],fields[3,n])
            result["degree"][str(n)]["contact_classification"]=contact_disagreement((2,n),(3,n))
        ref=fields[3,128]
        z=stack.enter_context(zipfile.ZipFile(stack.enter_context(SplitZip(parts))))
        def read(name):
            raw=z.read("paths/"+name)
            entry=index["entries"]["paths/"+name]
            assert len(raw)==entry["bytes"] and hashlib.sha256(raw).hexdigest()==entry["sha256"]
            return raw
        ground,case=[json.loads(read(k+".json")) for k in ("ground","stages")]
        for h in (0,1):
            mesh=json.loads(read(f"mesh-h{h}.json"))
            system=legacy.old.Water(mesh,ground,case).contact
            transfer=ContactPlan(mesh,system,h,pos,parents,"cubic")
            sides={}
            for side in ("native","independent"):
                folder=f"h{h}-g128/{side}/"
                manifest=json.loads(read(folder+"manifest.json"))
                values,rows=[],[]
                for chunk in manifest["chunks"]:
                    b=read(folder+chunk["file"])
                    assert hashlib.sha256(b).hexdigest()==chunk["sha256"] and len(rows)==chunk["start"]
                    with np.load(io.BytesIO(b),allow_pickle=False) as data:
                        raw=data["plastic_slip"]
                        values.append(legacy.old.codec.decode(raw,manifest["encoding"]) if manifest.get("encoding") else raw)
                    rows.extend(chunk["rows"])
                keys,indices=common_indices(rows,64)
                assert len(keys)==448
                raw=np.concatenate(values)[indices]
                sides[side]=transfer.history(raw)
                record=measures(sides[side],ref)
                record["direct_integral"]=difference(integral(raw,system.weight),integral(ref,area))
                result["shared_discrepancy"][f"h{h}-{side}"]=record
            result["shared_discrepancy"][f"h{h}-mean"]=measures((sides["native"]+sides["independent"])/2,ref)
            result["shared_discrepancy"][f"h{h}-solver-difference"]=measures(sides["native"],sides["independent"])
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=="__main__":
    main()
