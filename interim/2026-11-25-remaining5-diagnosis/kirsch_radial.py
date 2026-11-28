"""Test radius-only grading at the existing free-DOF cap, independently of probe locations."""

import json
import numpy as np
from kirsch_split import run, ProbeOnFacet
from kirsch import A, R

results = []
for exponent in (0, 1, 2):
    for nt in range(12, 33, 2):
        nr = 500 // nt - 1
        # Uniform steps in log(r), a/r or (a/r)^2. The analytic Kirsch
        # stress corrections contain inverse second and fourth powers of r.
        if exponent == 0:
            radii = np.geomspace(A, R, nr + 1)
        else:
            radii = A * np.linspace(1.0, (A / R) ** exponent, nr + 1) ** (-1 / exponent)
            radii[[0, -1]] = [A, R]
        try:
            result = run(nt, radii)
        except ProbeOnFacet as failure:
            result = dict(nt=nt, nr=nr, locator_failure=str(failure))
        results.append(dict(exponent=exponent, **result))
print(json.dumps(results, indent=2))
