"""Compare signed histories at the same 448 coordinates without rescaling inputs."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
NATIVE = HERE.parent / "2026-11-25-remaining-p3/c02/p4"
STAGES = ("geostatic", "shear", "open", "recontact", "reverse", "reload", "return")


def run(out):
    folders = {"native":NATIVE, "d4":out/"d4", "d8":out/"d8"}
    data = {label:np.load(path/"principal.npz") for label,path in folders.items()}
    fixed = np.load(HERE.parent / "2026-11-25-remaining5-p2/c02/p3-d8-n256/principal.npz")
    native = data["native"]
    w = native["weight"]
    pos = native["position"][:,:2]
    edge = np.minimum(pos,1-pos).min(axis=1) < .1
    for d in [*data.values(),fixed]:
        for key in ("parent","weight","position"):
            np.testing.assert_array_equal(d[key],native[key])
        assert d["slip"].shape == (448,72,3) and d["signed_area_integral"].shape == (448,3)
    report = dict(relative_width=.001, edge_width_m=.1, edge_area_fraction=w[edge].sum()/w.sum(),
        input_sha256={k:hashlib.sha256((v/"principal.npz").read_bytes()).hexdigest() for k,v in folders.items()},
        native_grid_sha256=hashlib.sha256((NATIVE/"increments.npz").read_bytes()).hexdigest(),pairs={})
    for left,right in (("d4","d8"),("native","d4"),("native","d8")):
        row = {}
        for field in ("slip","signed_area_integral"):
            weight = np.sqrt(w)[None,:,None] if field=="slip" else 1.
            diff=(data[left][field]-data[right][field])*weight
            sq=diff**2
            numerator=np.linalg.norm(diff)
            denom=np.linalg.norm(data[right][field]*weight)
            shares=sq.reshape(7,64,-1).sum(axis=(1,2))/sq.sum()
            row[field]=dict(norm=numerator,denominator_right_norm=denom,relative=numerator/denom,
                within_width=bool(numerator/denom<=.001),
                relative_to_fixed_p3_d8_n256_norm=numerator/np.linalg.norm(fixed[field]*weight),
                stage_squared_shares=dict(zip(STAGES,shares)),reverse_reload_share=shares[4:6].sum())
            if field=="slip":
                row[field].update(edge_share=sq[:,edge].sum()/sq.sum(),
                    edge_and_reverse_reload_share=sq[256:384,edge].sum()/sq.sum())
        report["pairs"][f"{left}-{right}"]=row
    (out/"comparisons.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]))
