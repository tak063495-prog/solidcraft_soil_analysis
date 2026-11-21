"""Prepare a separate, observation-only source copy. Never edit production."""
import argparse, hashlib, json, shutil, subprocess
from pathlib import Path
BASE='09ff926d82e574561393be2ac877cf9f5627d750'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def prepare(out):
    out.mkdir(parents=True,exist_ok=False)
    paths=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    fingerprints={}
    for n in paths:
        if n.startswith(('src/','examples/')) or n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):
            b=subprocess.check_output(['git','show',f'{BASE}:{n}'],cwd=ROOT)
            p=out/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
            fingerprints[n]=hashlib.sha256(b).hexdigest()
    p=out/'src/stages.rs';s=p.read_text()
    old='"internal": eval.internal, "target": target,'
    new='''"joint_state": eval.joints.as_ref().map(|j| &j.state),
                    "joint_report": eval.joints.as_ref().map(|j| &j.report),
                    "joint_internal": eval.joints.as_ref().map(|j| &j.internal_kn),
                    "dof_rows": dofs.rows,
                    "internal": eval.internal, "target": target,'''
    assert s.count(old)==1;s=s.replace(old,new)
    old='"alpha": alpha, "displacement_m": candidate, "states": trial.states,'
    new='''"alpha": alpha, "displacement_m": candidate, "states": trial.states,
                                "joint_state": trial.joints.as_ref().map(|j| &j.state),
                                "joint_report": trial.joints.as_ref().map(|j| &j.report),'''
    assert s.count(old)==1;s=s.replace(old,new)
    s+='\n'+(HERE/'diagnostic_hook.rs').read_text();p.write_text(s)
    for name in ('nc_stop_probe.rs','nc_midpoint_map_probe.rs'):
        shutil.copyfile(HERE/name,out/'examples'/name)
    (out/'original-source-sha256.json').write_text(json.dumps(fingerprints,indent=2))
    return fingerprints
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);args=a.parse_args();prepare(args.output.resolve())
