"""Signed adjacent differences, with and without the recontact endpoint."""
from pathlib import Path
import hashlib
import json
import zipfile

import numpy as np
from freeze import HERE, P7, save


def metric(x,y):
    dot=float(np.sum(x*y))
    a,b=map(float,(np.linalg.norm(x),np.linalg.norm(y)))
    return dict(dot=dot,coarse_norm=a,fine_norm=b,cosine=dot/(a*b) if a*b else None)


def run(out):
    out.mkdir(exist_ok=True,parents=True)
    slips=[]
    for n in (256,512,1024):
        root=P7/f"c01/alternate{n}"
        field=json.loads((root/"manifest.json").read_text())["fields"]["plastic_slip"]
        with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as pre, zipfile.ZipFile(root/"fields-tail.zip") as tail:
            raw=pre.read("plastic_slip.bin")+tail.read("plastic_slip.bin")
        assert hashlib.sha256(raw).hexdigest()==field["sha256"]
        slips.append(np.frombuffer(raw,field["dtype"]).reshape(7,256,*field["shape"]))
    s=np.stack(slips)
    endpoint=s[:,3,-1]
    post=s[:,4:]
    increments=post-endpoint[:,None,None]
    data=dict(levels=[256,512,1024],endpoint=endpoint,increments=increments)
    report={}
    for name,values in (("raw",post),("endpoint_removed",increments)):
        x,y=values[0]-values[1],values[1]-values[2]
        data[name+"_signed_adjacent"]=np.stack((x,y))
        stages={}
        for j,stage in enumerate(("reverse","reload","return")):
            by_point=[]
            for i in range(x.shape[2]):
                u,v=x[j,:,i],y[j,:,i]
                by_point.append(dict(point=i,**metric(u,v),
                    signed_sum_coarse=u.sum(axis=0).tolist(),signed_sum_fine=v.sum(axis=0).tolist(),
                    dot_xyz=np.sum(u*v,axis=0).tolist()))
            stages[stage]=dict(**metric(x[j],y[j]),points=sorted(by_point,key=lambda r:r["dot"]))
        report[name]=dict(**metric(x,y),stages=stages)
    # Expand (endpoint error + increment error) for each adjacent pair.
    b0,b1=endpoint[0]-endpoint[1],endpoint[1]-endpoint[2]
    inc0,inc1=data["endpoint_removed_signed_adjacent"]
    report["dot_decomposition_m2"]={}
    for j,stage in enumerate(("reverse","reload","return")):
        pieces=dict(inherited=float(256*np.sum(b0*b1)),
                    endpoint_coarse_times_increment_fine=float(np.sum(b0*inc1[j])),
                    increment_coarse_times_endpoint_fine=float(np.sum(inc0[j]*b1)),
                    increment_only=float(np.sum(inc0[j]*inc1[j])))
        # Evaluate the operation-order remainder directly from the signed arrays.
        x,y=data["raw_signed_adjacent"][:,j]
        reconstructed0,reconstructed1=b0+inc0[j],b1+inc1[j]
        pieces["roundoff_from_endpoint_subtractions"]=float(np.sum(
            (x-reconstructed0)*y + reconstructed0*(y-reconstructed1)))
        np.testing.assert_allclose(sum(pieces.values()),report["raw"]["stages"][stage]["dot"],rtol=2e-13)
        report["dot_decomposition_m2"][stage]=pieces
    all_raw=metric(s[0]-s[1],s[1]-s[2])
    prior=json.loads((P7/"c01/alternate-score.json").read_text())["fields"]["plastic_slip"]["measured"]
    np.testing.assert_allclose(all_raw["cosine"],prior["cosine"],atol=1e-14)
    report["whole_original"]=all_raw
    np.savez_compressed(out/"signed.npz",**data)
    save(out/"stages-points.json",report)


if __name__=="__main__":
    import sys
    run(Path(sys.argv[1]) if len(sys.argv)>1 else HERE/"c01/alternate")
