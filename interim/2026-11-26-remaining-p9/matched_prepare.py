"""Add diagnostic equilibrium solves inside previously accepted halves, with P7 charges."""
import json
from pathlib import Path
import sys

import numpy as np

from alternate import LEVELS


def plan(folder):
    z=np.load(folder/"accepted-halves.npz")
    grids=[z[f"recontact_{n}_progress"] for n in LEVELS]
    union=np.unique(np.concatenate(grids))
    windows=[]
    for i in range(3):
        ends=[]
        for n,t in zip(LEVELS,grids):
            status=z[f"recontact_{n}_classes"][:,i]
            j,=np.flatnonzero(status[1:]!=status[:-1])
            ends.extend(t[j:j+2])
        knots=union[(union>=min(ends))&(union<=max(ends))]
        windows.append(sorted(set(np.r_[knots,(knots[1:]+knots[:-1])/2].tolist())))
    times=np.unique(np.concatenate(windows))
    groups={}
    for n,t in zip(LEVELS,grids):
        outer=t[::2]
        groups[n]=[[float(l),float(r),times[(times>l)&(times<=r)].tolist()]
                   for l,r in zip(outer,outer[1:]) if np.any((times>l)&(times<=r))]
    (folder/"matched-plan.json").write_text(json.dumps(dict(points=[3,5,13],grids=windows,
        levels=groups,stop=float(times[-1])),indent=2)+"\n")


def patch(src,folder):
    data=json.loads((folder/"matched-plan.json").read_text())
    path=src/"stages.rs"
    s=path.read_text()
    key='(std::env::var_os("P9_HALVES").is_some() && stage.name()=="recontact" && next==1.)'
    assert s.count(key)==1
    path.write_text(s.replace(key,key+' || (std::env::var_os("P9_MATCHED").is_some() && stage.name()=="recontact" && next>='+repr(data["stop"])+')'))
    path=src/"stages_joint_path_error.rs"
    s=path.read_text()
    key='        let trial = if control.resolve_slip_path'
    assert s.count(key)==1
    extra=(Path(__file__).parent/"matched.rs.inc").read_text()
    path.write_text(s.replace(key,extra+"\n"+key))


if __name__=="__main__":
    folder=Path(sys.argv[1])
    plan(folder)
    patch(Path(sys.argv[2]),folder)
