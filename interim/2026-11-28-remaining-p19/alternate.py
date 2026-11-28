"""Resolve the last crossing using the unchanged p16 equilibrium operators."""
import csv
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P18 = HERE.with_name("2026-11-28-remaining-p18")
sys.path.insert(0, str(P18))
from reload_split import FIELDS, POINTS, prior


def run(out):
    _, operators = prior.operators()
    rows = []
    history = {}
    for route in ("outer", "middle"):
        with np.load(P18/f"c01/alternate/{route}-reload.npz") as data:
            time, n = data["progress"], data["normals"]
            assert len(time) == 641 and time[0] == 0 and time[-1] == 1
            dt = data["traction_difference"]
            tau = dt - (dt*n).sum(-1)[..., None]*n
            for field in FIELDS:
                field_diff = data[field+"_difference"]
                projected = (field_diff @ operators[field].reshape(96,-1).T).reshape(tau.shape)
                dot = lambda a,b: np.einsum("...d,...d->...", a,b)
                # First separate the fine/coarse product difference at every time.
                t_level = dot(projected.mean(0), tau[1]-tau[0])
                f_level = dot(projected[1]-projected[0], tau.mean(0))
                shares = dot(projected, tau)
                growth = shares[1]-shares[0]
                atol = 1e-25 if field == "reaction" else 1e-32
                np.testing.assert_allclose(growth,data[field+"_growth"][...,1],rtol=3e-11,atol=atol)
                np.testing.assert_allclose(t_level+f_level,growth,rtol=3e-11,atol=atol)
                # Then separate the change during the crossing interval at each level.
                traction_step = dot((projected[:,1:]+projected[:,:-1])/2, np.diff(tau,axis=1))
                field_step = dot(np.diff(projected,axis=1), (tau[:,1:]+tau[:,:-1])/2)
                tr, fi = np.diff(traction_step,axis=0)[0], np.diff(field_step,axis=0)[0]
                np.testing.assert_allclose(tr+fi,np.diff(growth,axis=0),rtol=1e-9,atol=atol)
                for name,a in dict(time=time,traction_level=t_level,field_level=f_level,
                                   growth=growth,traction_step=traction_step,field_step=field_step).items():
                    history[f"{route}_{field}_{name}"] = a
                for p in POINTS:
                    j = np.flatnonzero(growth[:,p]<=0)[-1]
                    assert j+1<len(time) and (growth[j+1:,p]>0).all()
                    before,after = growth[j:j+2,p]
                    # Evaluate both orders of replacing the factors. This explicitly
                    # checks whether one factor alone crosses zero from the left end.
                    def crossed(field_index, traction_index):
                        a = dot(projected[:,field_index,p],tau[:,traction_index,p])
                        return a[1]-a[0]
                    rows.append(dict(route=route,field=field,point=p,start=time[j],end=time[j+1],
                        before=before,after=after,traction_level_before=t_level[j,p],
                        traction_level_after=t_level[j+1,p],field_level_before=f_level[j,p],
                        field_level_after=f_level[j+1,p],traction_change=tr[j,p],field_change=fi[j,p],
                        traction_only_end=crossed(j,j+1),field_only_end=crossed(j+1,j),
                        at_reload_end=growth[-1,p]))
    with (out/"last-crossings.csv").open("w") as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(out/"factor-terms.npz",**history)


if __name__ == "__main__":
    run(Path(sys.argv[1]))
