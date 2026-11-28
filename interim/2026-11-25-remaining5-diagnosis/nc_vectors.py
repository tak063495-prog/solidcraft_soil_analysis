"""Generate full-rank binary64 comparison vectors from one saved linear system."""

import gzip
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import qr, solve, solve_triangular, svd

source = Path(sys.argv[1])
with gzip.open(source, "rt") as stream:
    record = json.load(stream)
rhs = np.array(record["b"])
matrix = np.array(record["a_column_major"]).reshape((len(rhs), len(rhs)), order="F")
q, r = qr(matrix, mode="economic")
u, s, vh = svd(matrix, full_matrices=False)
answers = {
    "saved": record["x"],
    "LAPACK_LU": solve(matrix, rhs, assume_a="general").tolist(),
    "Householder_QR": solve_triangular(r, q.T @ rhs).tolist(),
    "full_SVD": (vh.T @ ((u.T @ rhs) / s)).tolist(),
}
print(json.dumps(answers, allow_nan=False))
