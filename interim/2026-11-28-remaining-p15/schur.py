"""Compensate the two independent Schur solves and contractions separately."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P14 = HERE.with_name("2026-11-27-remaining-p14")
P13 = HERE.with_name("2026-11-27-remaining-p13")


def dot(a, b):
    return sum((F(x)*F(y) for x, y in zip(a, b, strict=True)), F(0))


def refined_solve(block, rhs):
    hi = np.linalg.solve(block, rhs)
    residual = [float(F(y)-dot(row, hi)) for row, y in zip(block, rhs, strict=True)]
    return hi, np.linalg.solve(block, residual)


def inspect(J, rr, lu, exact):
    A, row = J[:7,:7], J[7,:7]
    v, dv = refined_solve(A, -J[:7,7])
    z, dz = refined_solve(A, rr[:7])
    split = J[7,7] + row[:6]@v[:6] + row[6]*v[6]
    rhs = -rr[7] + row@z
    rows = {}
    for solve, inner, name in ((False,False,"existing"),(True,False,"solve_only"),
                               (False,True,"dot_only"),(True,True,"both")):
        if inner:
            d = F(J[7,7]) + dot(row,v)
            n = -F(rr[7]) + dot(row,z)
            if solve:
                d += dot(row,dv)
                n += dot(row,dz)
            d,n = float(d),float(n)
        else:
            d,n = split,rhs
            if solve:
                d += row[:6]@dv[:6] + row[6]*dv[6]
                n += row@dz
        x = float(n/d)
        err = abs(x-lu)
        rows[name] = dict(rhs=float(n),denominator=float(d),direction=x,error=err,
                          relative=err/abs(lu),relative_to_exact=abs(x-exact)/abs(exact),
                          passes=bool(err<=1e-14+5e-9*abs(lu)))
    return dict(variants=rows,v=v.tolist(),v_low=dv.tolist(),z=z.tolist(),z_low=dz.tolist(),
                lu=lu,exact=exact,lu_relative_to_exact=abs(lu-exact)/abs(exact),
                allowance=1e-14+5e-9*abs(lu))


def run(out):
    systems=json.loads((P13/"nc/direction-check.json").read_text())
    history=json.loads((P14/"nc/local-lu.json").read_text())["paths"]
    rows={}
    for name,(r,) in systems.items():
        p=history[name]
        assert r["point"]==p["point"]==351 and r["width"]==p["width_m"]
        row=inspect(np.array(r["matrix"]),np.array(r["rhs_vector"]),
                    p["one_local_refinement"]["dx"][-1],p["exact_direction"])
        assert row["variants"]["existing"]["direction"]==p["existing_binary64_schur"]["direction"]
        row.update(point=351,width_m=r["width"])
        rows[name]=row
    out.write_text(json.dumps(dict(paths=rows,atol=1e-14,rtol=5e-9,parent_updates=0,
                             new_global_trials=0),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
