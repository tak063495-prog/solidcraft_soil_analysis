"""Cross-check the live 8-vector before returning it to the native line search."""
from fractions import Fraction as F
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.with_name("2026-11-28-remaining-p15")))
import schur

j = json.load(sys.stdin)
M, rhs, dx = j["matrix"], j["rhs"], j["native"]
res = [float(-F(y)-schur.dot(row,dx)) for row,y in zip(M,rhs,strict=True)]
sol = json.loads(subprocess.check_output([os.environ["P16_LU"]],input=json.dumps(dict(matrix=M,b=res)),text=True))
updated = [x+y for x,y in zip(dx,sol["dx"],strict=True)]
r = schur.inspect(np.array(M),np.array(rhs),updated[-1],updated[-1])["variants"]["both"]
report = dict(**j,refined=updated,schur=r["direction"],denominator=r["denominator"],
              allowance=1e-14+5e-9*abs(updated[-1]),passes=r["passes"])
with open(os.environ["P16_RECORD"],"a") as out:
    print(json.dumps(report,allow_nan=False),file=out)
assert r["passes"]
print(json.dumps(updated,allow_nan=False))
