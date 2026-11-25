"""Compare both signed principal quantities on the original 448 load coordinates."""
import json
from pathlib import Path
import sys
import numpy as np
from c02_face import P1, PROPOSED_RELATIVE_WIDTH


def compare(first, second):
    with np.load(first/"principal.npz") as a, np.load(second/"principal.npz") as b:
        for field in ("parent","weight","position"):
            np.testing.assert_allclose(a[field], b[field], atol=2e-15,rtol=0)
        w = np.sqrt(b["weight"])[None,:,None]
        result = {}
        for name, factor in (("slip", w), ("signed_area_integral",1.)):
            x,y = a[name]*factor,b[name]*factor
            assert len(x)==len(y)==448
            result[name] = float(np.linalg.norm(x-y)/np.linalg.norm(y))
        result["within_proposed_width"] = max(result.values()) <= PROPOSED_RELATIVE_WIDTH
        return result


def main(raw):
    inputs = {p.name:p for p in raw.iterdir() if p.is_dir()}
    inputs.update({p.name:p for p in (P1/"c02").iterdir() if p.is_dir()})
    pairs = [
        ("face_2_4","p3-d2-n256","p3-d4-n256"),
        ("face_4_8","p3-d4-n256","p3-d8-n256"),
        ("time","p3-d8-n256","p3-d8-n512"),
        ("degree","p3-d8-n256","p4-d8-n256"),
        ("duffy_Q6","p3-q6-switch1","p3-d8-n256"),
        ("duffy_Q8","p3-q8-switch1","p3-d8-n256"),
        ("duffy_Q10","p3-q10-switch1","p3-d8-n256"),
    ]
    print(json.dumps(dict(proposed_relative_width=PROPOSED_RELATIVE_WIDTH,
        adoption="pending_user_decision",comparisons={k:dict(first=x,second=y,**compare(inputs[x],inputs[y]))
        for k,x,y in pairs}),indent=2,allow_nan=False))


if __name__=="__main__":
    main(Path(sys.argv[1]))
