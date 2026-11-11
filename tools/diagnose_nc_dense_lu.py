"""Audit captured rejected dense systems and offline LU iterative refinement.

All observations refer to the saved, unmodified A and b. Decimal arithmetic
measures cancellation independently of BLAS and Windows long-double support.
No result here accepts a solver step, edits a checkpoint, or certifies NC-01.
"""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

TOLERANCE = 1e-10


def decimal_residual(matrix, rhs, solution):
    with localcontext() as context:
        context.prec = 80
        d = Decimal.from_float
        b = [d(float(v)) for v in rhs]
        x = [d(float(v)) for v in solution]
        residual = [sum((d(float(a))*v for a, v in zip(row, x)), Decimal(0))-bv
                    for row, bv in zip(matrix, b)]
        numerator = sum((v*v for v in residual), Decimal(0)).sqrt()
        denominator = max(sum((v*v for v in b), Decimal(0)).sqrt(), Decimal('1e-30'))
        return float(numerator/denominator)


def run(systems, binary, output, refinements=8):
    if not 0 <= refinements <= 32:
        raise ValueError('refinements must be in 0..32')
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for index, line in enumerate(systems.read_text(encoding='utf-8').splitlines()):
        item = json.loads(line)
        a = np.asarray(item['matrix'], dtype=float)
        b = np.asarray(item['rhs'], dtype=float)
        if (a.ndim != 2 or a.shape != (len(b), len(b)) or len(b) == 0
                or not np.isfinite(a).all() or not np.isfinite(b).all()):
            raise ValueError('invalid captured system')
        request = dict(matrix=item['matrix'], rhs=item['rhs'], refinement_steps=refinements)
        child = subprocess.run([str(binary.resolve())], input=json.dumps(request, allow_nan=False),
                               text=True, capture_output=True, check=True)
        rows = json.loads(child.stdout)
        for row in rows:
            row['decimal80_relative_residual'] = decimal_residual(item['matrix'], item['rhs'], row['solution'])
            row['passes_unchanged_native_gate'] = row['true_relative_residual'] <= TOLERANCE
        singular = np.linalg.svd(a, compute_uv=False)
        record = dict(system=index, dimension=len(b), rejected=True, numerical_acceptance=False,
                      condition_number_2=float(singular[0]/singular[-1]),
                      singular_minimum=float(singular[-1]), singular_maximum=float(singular[0]),
                      native_saved_relative_residual=item['native_true_relative_residual'],
                      saved_solution_decimal80_relative_residual=decimal_residual(item['matrix'], item['rhs'], item['solution']),
                      replay_solution_exact=rows[0]['solution'] == item['solution'],
                      observations=rows,
                      minimum_native_relative_residual=min(r['true_relative_residual'] for r in rows),
                      any_refinement_passes_original_gate=any(r['passes_unchanged_native_gate'] for r in rows[1:]))
        records.append(record)
    if not records:
        raise ValueError('no captured rejected dense systems')
    report = dict(diagnostic_only=True, numerical_acceptance=False, tolerance=TOLERANCE,
                  original_systems_sha256=hashlib.sha256(systems.read_bytes()).hexdigest(),
                  diagnostic_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  decimal_precision=80, matrix_or_rhs_changed=False, refinement_applied_to_production=False,
                  systems=records)
    (output/'diagnosis.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--systems', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--refinements', type=int, default=8)
    args = parser.parse_args()
    report = run(args.systems, args.binary, args.output, args.refinements)
    print(json.dumps(dict(systems=len(report['systems']),
                         refinement_passes=sum(r['any_refinement_passes_original_gate'] for r in report['systems']),
                         numerical_acceptance=False), indent=2))
