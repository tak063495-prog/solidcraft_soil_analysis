"""Compare on common passive points whose triangles conform to the physical bands."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from face_pair import geometry, p6
from face_score import masks, integral, score
from freeze import HERE, save


def band_geometry():
    old=json.loads((p6.P5.parent/"2026-11-26-remaining-p6/c02/focused/tiles.json").read_text())["faces"]
    corners=[np.array(f["corners"])[:,:2] for f in old]
    basis=[np.column_stack((c[1]-c[0],c[2]-c[0])) for c in corners]
    tiles=[[] for c in corners]
    rings=[0,.025,.05,.1,.2,.3,.4,.5]
    def xy(d,s,side):
        x=d+(1-2*d)*s
        return ((x,d),(1-d,x),(1-x,1-d),(d,1-x))[side]
    for side in range(4):
        for lo,hi in zip(rings,rings[1:]):
            for left,right in ((0.,.5),(.5,1.)):
                a,b,c,d=xy(lo,left,side),xy(lo,right,side),xy(hi,right,side),xy(hi,left,side)
                for tri in ((a,b,c),(c,d,a)):
                    pts=np.array(tri)
                    if np.abs(np.linalg.det(pts[1:]-pts[:1])) < 1e-16:
                        continue
                    matches=[]
                    for i,base in enumerate(basis):
                        bary=np.linalg.solve(base,(pts-corners[i][0]).T).T
                        if bary.min()>=-2e-15 and bary.sum(axis=1).max()<=1+2e-15:
                            matches.append((i,bary))
                    (i,bary),=matches
                    tiles[i].append(bary)
    g=geometry([dict(parent=f["parent"],leaves=ts) for f,ts in zip(old,tiles,strict=True)])
    regions=masks(g.positions)
    for k,area in zip(("west","east","south","north","interior"),(.09,.09,.09,.09,.64)):
        np.testing.assert_allclose(g.area[regions[k]].sum(),area,atol=3e-15,rtol=0)
    return g


def run(work,out):
    g=band_geometry()
    model=json.loads((work/"native-raw/operators.json").read_text())
    cols=np.any(g.scalar_c!=0,axis=0)
    assert len(model["trace"][0])==cols.sum()
    n_common=len(model["observation"])
    model["observation"]+=g.scalar_c[:,cols].tolist()
    model["observation_normal"]+=g.system.normal.tolist()
    opfile,raw=work/"band-operators.json",work/"band-steps.jsonl"
    opfile.write_text(json.dumps(model))
    subprocess.run([work/"native-probe/target/release/soil-fem-3d","study-joint-p-cycle",opfile,raw],check=True)
    samples=[]
    with raw.open() as stream,(work/"native-raw/steps.jsonl").open() as baseline:
        for row,prior in zip(map(json.loads,stream),map(json.loads,baseline),strict=True):
            assert row["x"]==prior["x"] and row["t"]==prior["t"] and row["traction"]==prior["traction"]
            if "integral" in row:
                assert row["slip"][:n_common]==prior["slip"] and row["integral"]==prior["integral"]
                samples.append(row["slip"][n_common:])
    native=np.array(samples)
    del samples
    assert len(native)==448
    op=np.ascontiguousarray(g.scalar_c[:,cols])
    regions=masks(g.positions)
    report=dict(passive_points=len(g.area),native_clock_and_active_traction_rows_matched=4156,
                native_principal_anchors_matched=448,
                regions={k:dict(points=int(m.sum()),area=g.area[m].sum()) for k,m in regions.items()})
    histories={}
    for name in ("plus","minus"):
        archive=np.load(out/name/"kinematics.npz")
        old=np.zeros((len(g.area),3))
        traction,flow=old.copy(),old.copy()
        samples=[]
        for (phase,t,width),u in zip(archive["clock"],archive["u"],strict=True):
            target=op @ u.reshape(-1,3)[cols]
            counts=np.maximum(1,np.ceil(np.linalg.norm(target-old,axis=1)/width)).astype(int)
            q,ds,_=p6.face.ref.project(g.system,old,traction,target[:,None,:],counts)
            old,traction=target,q[:,0]
            flow+=ds[:,0]
            if t*64==int(t*64):
                samples.append(flow.copy())
        ref=np.array(samples)
        del samples
        rows={}
        for k,mask in regions.items():
            a,b=integral(native,g.area,mask),integral(ref,g.area,mask)
            histories["native_"+k]=a
            histories[name+"_"+k]=b
            rows[k]=dict(slip=score(native[:,mask],ref[:,mask],np.sqrt(g.area[mask])[None,:,None]),
                         signed_integral=score(a,b))
        report[name]=rows
    np.savez_compressed(out/"physical-band-integrals.npz",**histories,position=g.positions,area=g.area,parent=g.parents)
    save(out/"physical-bands.json",report)


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
