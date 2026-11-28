"""Factor the p19 endpoint product into field, traction length and unit direction."""
import csv
import itertools
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P19 = HERE.with_name("2026-11-28-remaining-p19")
sys.path.insert(0, str(P19))
from alternate import FIELDS, P18, prior


def run(out):
    _, operators = prior.operators()
    with (P19/"c01/alternate/last-crossings.csv").open() as stream:
        targets = list(csv.DictReader(stream))
    rows = []
    for route in ("outer", "middle"):
        with np.load(P18/f"c01/alternate/{route}-reload.npz") as data:
            t, n, dt = (data[k] for k in ("progress", "normals", "traction_difference"))
            assert len(t) == 641
            tau = dt - (dt*n).sum(-1)[..., None]*n
            for field in FIELDS:
                projected = (data[field+"_difference"] @ operators[field].reshape(96,-1).T).reshape(tau.shape)
                for ref in targets:
                    p = int(ref["point"])
                    if (ref["route"], ref["field"]) != (route,field) or not (p == 21 or p == 25 and field == "reaction"):
                        continue
                    j, = np.flatnonzero(t==float(ref["start"]))
                    assert t[j+1] == float(ref["end"])
                    f = projected[:,j:j+2,p]
                    tr = tau[:,j:j+2,p]
                    m = np.linalg.norm(tr,axis=-1)
                    assert (m>0).all()
                    u = tr/m[:,:,None]
                    def gap(fi, mi, ui):
                        pair = m[:,mi]*np.sum(f[:,fi]*u[:,ui],axis=-1)
                        return float(pair[1]-pair[0])
                    swaps = {"".join(map(str,s)):gap(*s) for s in itertools.product(range(2),repeat=3)}
                    for k, col in [("000","before"),("111","after"),("100","field_only_end"),("011","traction_only_end")]:
                        np.testing.assert_allclose(swaps[k],float(ref[col]),rtol=3e-11,atol=1e-39)
                    # Exact midpoint product identity, applied first to F·T, then T=m*u.
                    fm, um, mm = f.mean(1), u.mean(1), m.mean(1)
                    dm, du = m[:,1]-m[:,0], u[:,1]-u[:,0]
                    terms = dict(field=np.sum((f[:,1]-f[:,0])*tr.mean(1),axis=-1),
                        magnitude=dm*np.sum(fm*um,axis=-1), direction=mm*np.sum(fm*du,axis=-1))
                    contributions = {k:float(v[1]-v[0]) for k,v in terms.items()}
                    np.testing.assert_allclose(sum(contributions.values()),swaps["111"]-swaps["000"],rtol=2e-11)
                    rows.append(dict(route=route,field=field,point=p,interval=t[j:j+2].tolist(),
                        length_by_level_pair=m.tolist(),direction_by_level_pair=u.tolist(),projected_field=f.tolist(),
                        contributions=contributions,swaps_field_magnitude_direction=swaps,
                        necessary_at_right=[k for k,mask in [("field","011"),("magnitude","101"),("direction","110")] if swaps[mask]<=0]))
    assert len(rows)==8
    (out/"factors.json").write_text(json.dumps(rows,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
