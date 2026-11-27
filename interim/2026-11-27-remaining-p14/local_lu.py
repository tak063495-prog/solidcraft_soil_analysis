"""Exact residual and partial-pivot audit on the two archived binary64 systems."""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from plans import P13, save


def native(binary, matrix, b):
    raw = json.dumps(dict(matrix=matrix, b=[float(v) for v in b]))
    return json.loads(subprocess.check_output([binary], input=raw, text=True))


def elimination(A,b):
    a = [row[:] for row in A]
    y = b[:]
    pivots, terms = [], []
    perm = list(range(len(b)))
    for i in range(len(b)):
        j = max(range(i,len(b)),key=lambda k: abs(a[k][i]))
        a[i],a[j] = a[j],a[i]
        y[i],y[j] = y[j],y[i]
        perm[i],perm[j] = perm[j],perm[i]
        pivots.append(a[i][i])
        for k in range(i+1,len(b)):
            w = a[k][i]/a[i][i]
            if k == len(b)-1:
                terms.append(w*a[i][-1])
            for col in range(i,len(b)):
                a[k][col] -= w*a[i][col]
            y[k] -= w*y[i]
    x = [F(0)]*len(b)
    for i in reversed(range(len(x))):
        x[i] = (y[i]-sum(a[i][j]*x[j] for j in range(i+1,len(x))))/a[i][i]
    return x, pivots, perm, terms


def audit(A,b,values,exact):
    x = list(map(F,values))
    residual = [v-sum(c*z for c,z in zip(row,x,strict=True)) for row,v in zip(A,b,strict=True)]
    denom = [abs(v)+sum(abs(c*z) for c,z in zip(row,x,strict=True)) for row,v in zip(A,b,strict=True)]
    eta = max(abs(v)/d for v,d in zip(residual,denom,strict=True) if d)
    norm = max(map(abs,residual))/(max(sum(map(abs,row)) for row in A)*max(map(abs,x))+max(map(abs,b)))
    delta = abs(exact[-1]-x[-1])
    allowance = F(1e-14)+F(5e-9)*abs(x[-1])
    return dict(dx=values,residual=[float(v) for v in residual],normwise_backward_error=float(norm),
                componentwise_backward_error=float(eta),direction_absolute_error=float(delta),
                direction_relative_error=float(delta/abs(x[-1])),allowance=float(allowance),
                passes=delta<=allowance),residual


def run(binary,out):
    rows = json.loads((P13/"nc/direction-check.json").read_text())
    result = {}
    for name, (r,) in rows.items():
        M = r["matrix"]
        A = [list(map(F,row)) for row in M]
        b = [-F(v) for v in r["rhs_vector"]]
        exact,pivots,perm,terms = elimination(A,b)
        obs = native(binary,M,b)
        assert obs["dx"][-1] == r["native_dx"] and obs["permutation"] == perm
        assert float(exact[-1]) == r["exact_direction"]
        original,res = audit(A,b,obs["dx"],exact)
        correction = native(binary,M,res)["dx"]
        refined,_ = audit(A,b,[x+y for x,y in zip(obs["dx"],correction,strict=True)],exact)
        # Reproduce p8.point's split stress/hardening evaluation order as well
        # as p13's independent exact Schur reference.
        J,rr = np.array(M),np.array(r["rhs_vector"])
        v = -np.linalg.solve(J[:7,:7],J[:7,7])
        stress,hard = J[7,:6]@v[:6],J[7,6]*v[6]
        den = J[7,7]+stress+hard
        rhs = -rr[7]+J[7,:7]@np.linalg.solve(J[:7,:7],rr[:7])
        direction = float(rhs/den)
        split = dict(stress=float(stress),hardening=float(hard),denominator=float(den),rhs=float(rhs),direction=direction)
        for stage,d in (("native",original),("one_local_refinement",refined)):
            x = d["dx"][-1]
            split[stage] = dict(relative_error=abs(direction-x)/abs(x),
                                passes=abs(direction-x)<=1e-14+5e-9*abs(x))
        result[name] = dict(existing_binary64_schur=split,point=r["point"],width_m=r["width"],permutation=perm,
            native_pivots=obs["pivots"],exact_pivots=[float(x) for x in pivots],
            final_column_terms=[float(x) for x in terms],
            final_pivot_cancellation=float(sum(map(abs,terms))/abs(pivots[-1])),
            final_pivot_relative_error=float(abs(F(obs["pivots"][-1])/pivots[-1]-1)),
            exact_direction=float(exact[-1]),native=original,one_local_refinement=refined)
    save(out/"local-lu.json",dict(rtol=5e-9,atol=1e-14,global_attempts_added=0,paths=result))


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]))
