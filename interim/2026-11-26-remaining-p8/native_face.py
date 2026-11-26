"""Recover the unsaved native active-point slip using output-only instrumentation."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from freeze import HERE, save

P3=HERE.parent/"2026-11-25-remaining-p3"
sys.path.insert(0,str(P3))
import c02_run as original


def run(work,out):
    probe=work/"native-probe"
    probe.mkdir()
    root=HERE.parents[1]
    for name in ("Cargo.lock","Cargo.toml","rust-toolchain.toml"):
        shutil.copyfile(root/name,probe/name)
    for name in ("src","examples"):
        shutil.copytree(root/name,probe/name)
    source=probe/"src/joint_p_cycle.rs"
    marker='                row["integral"] = json!(integral.as_slice());'
    s=source.read_text()
    assert s.count(marker)==1
    source.write_text(s.replace(marker,marker+'\n                row["active_slip"] = json!(old.iter().map(|s| s.plastic_slip_m).collect::<Vec<_>>());'))
    subprocess.run(["cargo","build","--release","--locked","--bin","soil-fem-3d"],cwd=probe,check=True)
    raw=work/"native-raw"
    dest=work/"native-repeat"
    original.run(raw,dest,probe/"target/release/soil-fem-3d",4,512,2)
    prior=P3/"c02/p4"
    repeat={}
    for name in ("increments","principal"):
        with np.load(prior/(name+".npz")) as a, np.load(dest/(name+".npz")) as b:
            repeat[name]={k:float(np.max(np.abs(a[k]-b[k]))) for k in a.files}
    inputs=[json.loads((root/"examples"/(s+".json")).read_text()) for s in
            ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
    g=original.geometry(inputs,4,10)
    samples=[]
    with (raw/"steps.jsonl").open() as stream:
        for row in map(json.loads,stream):
            if "active_slip" in row:
                samples.append(row["active_slip"])
    assert len(samples)==448
    np.savez_compressed(out/"native-active.npz",slip=samples,position=g.positions,weight=g.area,parent=g.parents)
    save(out/"native-replay.json",repeat)


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
