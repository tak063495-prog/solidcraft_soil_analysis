import json
from pathlib import Path
import sys
import numpy as np

root = Path.cwd()
src = root / "interim"
sys.path[:0] = [
    str(src / x)
    for x in (
        "2026-11-25-approved-six",
        "2026-11-24-t9d-benchmarks",
        "2026-11-22-impl-t7-physical-instability",
    )
]
from arch import beam_curve
from curved import AnalyticChordArch
from references import force_unit, R

for n in (12, 24, 48, 96):
    arch = AnalyticChordArch(n)
    result = beam_curve(n, extrapolate=True, model=arch, fold_metrics=True)
    f = result["fold"]
    z = f["z"]
    ev = arch({}, z)
    left, singular, right = np.linalg.svd(ev.k)
    v = right[-1]
    crown = np.flatnonzero(arch.free == 3 * arch.apex + 1).item()
    curvature = float(left[:, -1] @ arch.tangent_action_derivative(z, v))
    wq = float(left[:, -1] @ ev.q)
    print(
        json.dumps(
            dict(
                elements=n,
                classification=f["kind"],
                load=z[-1] * force_unit,
                normalized=f["curvature_over_scale"],
                fd_normalized=f["fd_curvature_over_scale"],
                numerator=abs(curvature),
                scale=abs(f["curvature"]) / f["curvature_over_scale"],
                wq=abs(wq),
                physical_second_derivative=-force_unit
                * curvature
                / (wq * (v[crown] * R) ** 2),
                sv_last=singular[-2:].tolist(),
                rank=int(f["bordered_rank"]),
                dim=len(z),
                trials=len(result["trials"]),
                force=f["force_ratio"],
            )
        ),
        flush=True,
    )
