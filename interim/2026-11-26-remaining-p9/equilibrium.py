"""Weak contact-force perturbation at one frozen trajectory, plus re-equilibrated displacement differences."""
import json
from pathlib import Path
import sys

import numpy as np

from face_bounds import geometry, P8, ALPHA, HERE
from face_pair import p6


def rule(folder):
    blocks,law=geometry(folder)
    faces=json.loads((folder/"tiles.json").read_text())["faces"]
    bary,_=p6.face.triangle_rule(1)
    rows=[]
    for face,(ids,_,_,_,_,sgn) in zip(faces,blocks,strict=True):
        xy=(bary@np.array(face["leaves"])).reshape(-1,2)
        l=np.column_stack((1-xy.sum(axis=1),xy))
        shape=np.stack([np.prod([(4*l[:,j]-k)/(k+1) for j in range(3) for k in range(a[j])],axis=0)
                        for a in ALPHA],axis=1)
        rows.append((ids,shape,sgn))
    weights=np.load(folder/"slip.npz")["weight"]
    n=np.repeat([r[-1] for r in rows],[len(r[1]) for r in rows])
    from types import SimpleNamespace
    system=SimpleNamespace(normal=np.column_stack((n*0,n*0,n)),
                           kn=np.full(len(n),law["normal_stiffness_kpa_m"]),
                           kt=np.full(len(n),law["tangential_stiffness_kpa_m"]),
                           mu=np.full(len(n),law["friction_coefficient"]))
    return rows,weights,system


def forces(folder,archive):
    rows,weights,syslaw=rule(folder)
    old=np.zeros((len(weights),3))
    traction=old.copy()
    result=[]
    for (_,t,width),u in zip(archive["clock"],archive["u"],strict=True):
        blocks=u.reshape(-1,3)
        target=np.concatenate([s@(blocks[ids[1]]-blocks[ids[0]]) for ids,s,_ in rows])
        counts=np.maximum(1,np.ceil(np.linalg.norm(target-old,axis=1)/width)).astype(int)
        q,_,_=p6.face.ref.project(syslaw,old,traction,target[:,None,:],counts)
        old,traction=target,q[:,0]
        if t*64!=int(t*64): continue
        force=np.zeros_like(blocks)
        start=0
        for ids,s,_ in rows:
            end=start+len(s)
            f=s.T@(weights[start:end,None]*traction[start:end])
            np.add.at(force,ids[0],-f)
            np.add.at(force,ids[1],f)
            start=end
        result.append(force)
    return np.array(result)


def run(out):
    keys=("plus","minus")
    folders=[P8/"c02"/k for k in keys]
    k=[np.load(folder/"kinematics.npz") for folder in folders]
    np.testing.assert_array_equal(k[0]["clock"],k[1]["clock"])
    du=k[0]["u"]-k[1]["u"]
    cross=[]
    for archive in k:
        cross.append(forces(folders[0],archive)-forces(folders[1],archive))
    out.mkdir(exist_ok=True,parents=True)
    np.savez_compressed(out/"equilibrium.npz",clock=k[0]["clock"],u_plus_minus=du,
                        fixed_path_plus_minus_force=cross)
    (out/"equilibrium.json").write_text(json.dumps(dict(
        clocks=len(du),max_nodal_m=float(np.linalg.norm(du.reshape(len(du),-1,3),axis=2).max()),
        history_relative_to_minus=float(np.linalg.norm(du)/np.linalg.norm(k[1]["u"])),
        force_l2_max_kn=[float(np.linalg.norm(f.reshape(len(f),-1),axis=1).max()) for f in cross],
        uniform_inverse_or_history_gain=None),indent=2)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
