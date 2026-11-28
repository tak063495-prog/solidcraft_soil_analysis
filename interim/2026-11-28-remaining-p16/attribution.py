"""Attribute measured squared norm changes to contact tractions on the elastic bodies."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")
P15 = HERE.with_name("2026-11-28-remaining-p15")
sys.path.insert(0,str(HERE.with_name("2026-11-26-remaining-p11")))
from alternate_score import P7,P10,p10_bytes,STAGES,observed
from c01_native import context

LEVELS = (256,512,1024)
NAMES = ("u","jump","reaction","plastic_slip","accumulated","traction")


def read(path):
    return json.loads(path.read_text())


def restore(route, n, names=NAMES):
    base = P10/str(n)
    meta = read(base/"manifest.json")["fields"]
    d = route/str(n)
    manifest = read(d/"manifest.json")
    assert manifest["count"] == 1792
    with zipfile.ZipFile(P7/f"c01/first-four-n{n}.zip") as shared, zipfile.ZipFile(base/"prefix-xor.zip") as xor, \
         zipfile.ZipFile(base/"fields-after-open.zip") as tail, zipfile.ZipFile(d/"fields-xor-p10.zip") as patch:
        vals = {}
        for key in names:
            m = manifest["fields"][key]
            raw = np.frombuffer(p10_bytes(n,meta[key],key,shared,xor,tail),"u1")
            if key+".bin" in patch.namelist():
                raw = raw ^ np.frombuffer(patch.read(key+".bin"),"u1")
            assert hashlib.sha256(raw).hexdigest() == m["sha256"]
            vals[key] = raw.view(m["dtype"]).reshape(-1,*m["shape"])
    return vals


def operators():
    study,mesh,ground,case,_ = context(256)
    op = study.prior.operators(mesh,ground,case,"cut")
    C,K = op.system.c, op.k
    # Prescribed components coincide across all three levels at a matched load.
    P = np.stack([op.reduce(e) for e in np.eye(len(K))],axis=1)
    B = C * op.system.weight[:,None,None]
    kred = P @ K @ P.T
    assert np.linalg.eigvalsh(kred).min() > 0
    U = (-P.T @ np.linalg.solve(kred,P @ B.reshape(-1,len(K)).T)).T
    L = {"u":U, "jump": U @ C.reshape(-1,len(K)).T,
         "reaction": U @ K.T + B.reshape(-1,len(K))}
    return op,{s:v.reshape(32,3,-1) for s,v in L.items()}


def measure(levels, L):
    diffs = {k:-np.diff(np.stack([x[k] for x in levels]),axis=0) for k in NAMES}
    result = {}
    for k in ("plastic_slip","accumulated","jump"):
        a,b = diffs[k].reshape(2,-1,32,levels[0][k].size//(len(levels[0][k])*32))
        result[k+"_dot"] = (a*b).sum(-1)
        result[k+"_receiver_growth"] = ((b-a)*(b+a)).sum(-1)
    for k,op in L.items():
        d = diffs[k].reshape(2,len(levels[0][k]),-1)
        traction = diffs["traction"]
        # Project measured differences onto each point's equilibrium response.
        # This retains cross terms, including interference between contact points.
        projected = (d @ op.reshape(96,-1).T).reshape(2,len(d[0]),32,3)
        shares = (projected*traction).sum(-1)
        leftover = d - traction.reshape(2,len(d[0]),96) @ op.reshape(96,-1)
        res = (d*leftover).sum(-1)
        sq = (d*d).sum(-1)
        np.testing.assert_allclose(shares.sum(-1)+res,sq,rtol=1e-10,atol=1e-40)
        result[k+"_point_growth"] = shares[1]-shares[0]
        result[k+"_residual_growth"] = res[1]-res[0]
        result[k+"_squared_norms"] = sq.T
    return result


def summary(terms,sl):
    return {k:np.asarray(v[sl].sum(axis=0)).tolist() for k,v in terms.items()}


def run(out, extra):
    op,L = operators()
    saved = {}
    for name,p in (("outer",P15),("middle",P14)):
        route = p/"c01/alternate"/name
        levels = [restore(route,n) for n in LEVELS]
        terms = measure(levels,L)
        score = read(route/"score.json")
        for k in NAMES:
            actual = observed(*[v[k] for v in levels])
            for stat,val in actual.items():
                ref = score["fields"][k]["full"][stat]
                if isinstance(val,float):
                    np.testing.assert_allclose(val,ref,rtol=1e-12,atol=1e-30)
                else:
                    assert val == ref
        assert (score["common_coordinates"],score["passed_fields"],len(score["fields"])) == (1792,17,22)
        np.savez_compressed(out/(name+"-1792.npz"),progress=np.tile(np.arange(1,257)/256,7),**terms)
        saved[name] = dict(stages={s:summary(terms,slice(j*256,(j+1)*256)) for j,s in enumerate(STAGES)},
                           common=1792,passes=17,fields=22)
        if extra:
            dense = []
            for n,grid in zip(LEVELS,levels,strict=True):
                with np.load(out/f"{name}-{n}.npz") as f:
                    vals = {k:f[k] for k in f.files}
                t = vals.pop("t")
                kt = np.array([op.k@u.ravel() for u in vals["u"]])
                # Apply the identical archived native field extraction order per load.
                vals["reaction"] = np.array([(internal + op.seed +
                    np.einsum("qag,qa,q->g",op.system.c,q,op.system.weight) - load.ravel()).reshape(-1,3)
                    for internal,q,load in zip(kt,vals["traction"],vals.pop("load"),strict=True)])
                ix = np.searchsorted(t,np.arange(1,257)/256)
                np.testing.assert_array_equal(t[ix],np.arange(1,257)/256)
                for k,v in vals.items():
                    np.testing.assert_array_equal(v[ix],grid[k][3*256:4*256])
                dense.append(vals)
            m = measure(dense,L)
            with gzip.open(P15/"c01/alternate/common-middle.jsonl.gz","rt") as f:
                records = list(map(json.loads,f))
            key = "p15_outer" if name=="outer" else "p14_middle"
            np.testing.assert_array_equal(t,[r["t"] for r in records])
            np.testing.assert_array_equal(m["plastic_slip_dot"],[r["points_m2"][key] for r in records])
            np.savez_compressed(out/(name+"-382.npz"),progress=t,**m)
            saved[name]["recontact_end"] = summary(m,slice(-1,None))
            saved[name]["shared_recontact_ends"] = len(t)
    (out/"attribution.json").write_text(json.dumps(saved,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]),"--dense" in sys.argv)
