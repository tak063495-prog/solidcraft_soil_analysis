import json
import numpy as np
import sys

d = np.load(sys.argv[1])
for key in ("native_k", "independent_k"):
    k = d[key]
    eig = np.linalg.eigvals(k)
    s = np.linalg.svd(k, compute_uv=False)
    sym = (k + k.T) / 2
    w = np.linalg.eigvalsh(sym)
    print(
        json.dumps(
            dict(
                matrix=key,
                skew=np.linalg.norm(k - k.T) / np.linalg.norm(k),
                eigen_min_real=float(np.min(eig.real)),
                eigen_max_imag=float(np.max(abs(eig.imag))),
                negative_energy_count=int(sum(w < 0)),
                rayleigh_min=w[0],
                rayleigh_max=w[-1],
                singular_min=s[-1],
                singular_max=s[0],
                kappa=s[0] / s[-1],
                eps_kappa=(s[0] / s[-1]) * np.finfo(float).eps,
            )
        )
    )
