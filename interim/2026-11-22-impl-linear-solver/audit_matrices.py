"""Inspect the original reduced matrices, including the explicitly rejected solution."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import lu_factor, svdvals

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from diagnose_nc_dense_lu import decimal_residual


def extremes(values):
    return [float(values.min()), float(values.max())]


def measure(raw, witness):
    item = json.loads(raw)
    a = np.asarray(item['matrix'])
    b = np.asarray(item['rhs'])
    x = np.asarray(item['solution'])
    singular = svdvals(a)
    lu, _ = lu_factor(a)
    high = decimal_residual(a, b, x)
    b_norm = float(np.linalg.norm(b))
    x_norm = float(np.linalg.norm(x))
    if witness is not None:
        y = np.asarray(witness['solution'])
        if y.shape != b.shape or not np.isfinite(y).all():
            raise ValueError('invalid correction witness')
        corrected = decimal_residual(a, b, y)
    record = dict(
        sha256=hashlib.sha256(raw).hexdigest(), dimension=len(x),
        numerical_rank=int(sum(singular > len(x) * np.finfo(float).eps * singular[0])),
        sigma_max=float(singular[0]), sigma_penultimate=float(singular[-2]),
        sigma_min=float(singular[-1]), kappa2=float(singular[0] / singular[-1]),
        row_max_abs_range=extremes(abs(a).max(axis=1)),
        column_max_abs_range=extremes(abs(a).max(axis=0)),
        lapack_abs_U_pivots=extremes(abs(lu.diagonal())),
        lapack_U_growth=float(abs(np.triu(lu)).max() / abs(a).max()),
        rhs_norm2=b_norm, solution_norm2=x_norm,
        native_relative_residual=item['native_true_relative_residual'],
        decimal80_relative_residual=high,
        decimal80_backward_error_2=high * max(b_norm, 1e-30) / (singular[0] * x_norm + b_norm),
    )
    if witness is not None:
        record.update(corrected_decimal80_relative_residual=corrected,
                      corrected_native_relative_residual=witness['corrected_native_relative_residual'],
                      native_abs_U_pivots=witness['native_abs_U_pivots'])
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('systems', type=Path, help='observer JSONL or a one-line fixture')
    parser.add_argument('--solutions', type=Path, help='saved corrected solutions keyed by system SHA-256')
    args = parser.parse_args()
    witnesses = None
    if args.solutions:
        with args.solutions.open('rb') as stream:
            witnesses = {row['system_sha256']: row for row in map(json.loads, stream)}
    with args.systems.open('rb') as stream:
        records = [measure(line, None if witnesses is None else witnesses[hashlib.sha256(line).hexdigest()])
                   for line in stream]
    print(json.dumps(records, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
