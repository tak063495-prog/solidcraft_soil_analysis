"""Compare complete signed histories at the identical 448 prescribed coordinates."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P2 = HERE.parent / "2026-11-25-remaining5-p2"
P1 = HERE.parent / "2026-11-25-remaining5-p1"
sys.path.insert(0,str(P2))
from c02_compare import compare


def metrics(a,b):
    row = compare(a,b)
    return dict(slip_relative=row["slip"],signed_area_relative=row["signed_area_integral"])


def digest(folder):
    return hashlib.sha256((folder/"principal.npz").read_bytes()).hexdigest()


def main(folder):
    reference = P2/"c02/p3-d8-n256"
    coarse = P2/"c02/p3-d4-n256"
    delta = metrics(coarse,reference)
    archived = json.loads((P2/"c02/summary.json").read_bytes())["comparisons"]["face_4_8"]
    np.testing.assert_allclose(list(delta.values()), [archived["slip"],archived["signed_area_integral"]],rtol=1e-13)
    body = [folder / f"p{p}" for p in (3,4,5)]
    with np.load(reference/"principal.npz") as data:
        normalizers = dict(slip_area_weighted=float(np.linalg.norm(data["slip"]*np.sqrt(data["weight"])[None,:,None])),
                           signed_area_integral=float(np.linalg.norm(data["signed_area_integral"])))
    summary = dict(reference_history_norms=normalizers,decision_by="user",decision_time="2026-11-25T20:03:00+09:00",
                   adopted_reference_relative_width=0.001,comparison_scope="face subdivision d4 to d8, P3 N256",
                   reference_route=str(reference.relative_to(HERE.parents[1])),reference_sha256=digest(reference),
                   d4_sha256=digest(coarse),face_refinement=delta,
                   both_face_differences_within_width=all(x<=0.001 for x in delta.values()),
                   native_vs_archived_same_P3_Q10_N256=metrics(body[0],P1/"c02/p3-q10-switch1"),
                   body_to_fixed_reference={p.name:metrics(p,reference) for p in body},
                   diagnostic_P4_vs_archived_P4_d8_N256=metrics(body[1],P2/"c02/p4-d8-n256"),
                   adjacent_body={a.name+"_"+b.name:metrics(a,b) for a,b in zip(body,body[1:])},
                   body_certification="not awarded")
    assert summary["both_face_differences_within_width"]
    (folder/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main(Path(sys.argv[1]) if len(sys.argv)>1 else HERE/"c02")
