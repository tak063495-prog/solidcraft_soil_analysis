"""Postmortem of the matrix, RHS and displacement correction as saved at rejection."""

from fractions import Fraction
import json
from pathlib import Path
import sys
import numpy as np

run = json.loads(Path(sys.argv[1]).read_bytes())
reason = run["failure_reason"]
res = float(reason.split("true linear residual ")[1].split()[0])
match = [
    f
    for f in Path(sys.argv[2]).glob("reject-*.json")
    if json.loads(f.read_bytes())["relative"] == res
]
assert len(match) == 1
record = json.loads(match[0].read_bytes())
a = np.array(record["a_column_major"]).reshape((348, 348), order="F")
b, x = map(np.array, (record["b"], record["x"]))
s = np.linalg.svd(a, compute_uv=False)
eig = np.linalg.eigvals(a)
e = np.linalg.eigvalsh((a + a.T) / 2)
# Each saved double is a dyadic rational. Compare the residual of exactly
# these saved operands. No correction or load step is computed here.
xf = [Fraction(float(v)) for v in x]
r = np.array(
    [
        float(
            sum(
                (Fraction(float(v)) * x0 for v, x0 in zip(row, xf, strict=True) if v),
                -Fraction(float(b0)),
            )
        )
        for row, b0 in zip(a, b, strict=True)
    ]
)
normb = np.linalg.norm(b)
print(
    json.dumps(
        dict(
            reason=reason,
            accepted=run["last_committed"]["stage_progress"],
            attempts=len(run["attempts"]),
            source=match[0].name,
            rank=int(np.linalg.matrix_rank(a)),
            smin=s[-1],
            smax=s[0],
            cond=s[0] / s[-1],
            rayleigh_min=e[0],
            negative_energy_modes=int(np.sum(e < 0)),
            eigen_min_real=eig.real.min(),
            eigen_max_imag=np.max(abs(eig.imag)),
            rhs_norm=normb,
            correction_norm=np.linalg.norm(x),
            native_absolute=record["absolute"],
            native_relative=res,
            exact_operands_absolute=np.linalg.norm(r),
            exact_operands_relative=np.linalg.norm(r) / normb,
            cancellation=np.linalg.norm(abs(a) @ abs(x)) / normb,
            backward=np.linalg.norm(r) / (s[0] * np.linalg.norm(x) + normb),
        )
    ),
    flush=True,
)
