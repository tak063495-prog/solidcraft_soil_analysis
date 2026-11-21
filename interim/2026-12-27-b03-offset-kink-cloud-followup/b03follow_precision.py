"""Decimal solves of the SAME rounded matrices; no new element or physics."""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import time
import numpy as np
from b03follow_math import dot, norm, write


def decimal_saddle(k, r, f, digits):
    n, nr = len(f), r.shape[1]
    with localcontext() as ctx:
        ctx.prec = digits
        number = lambda x: Decimal.from_float(float(x))
        a = [[number(v) for v in row]+[number(v) for v in rr] for row, rr in zip(k, r, strict=True)]
        a.extend([[number(v) for v in col]+[Decimal(0)]*nr for col in r.T])
        b = [number(v) for v in f]+[Decimal(0)]*nr
        size = n+nr
        for j in range(size):
            pivot = max(range(j, size), key=lambda i: abs(a[i][j]))
            if a[pivot][j] == 0:
                raise ValueError('singular decimal saddle matrix')
            a[j], a[pivot] = a[pivot], a[j]
            b[j], b[pivot] = b[pivot], b[j]
            row = a[j]
            for i in range(j+1, size):
                target = a[i]
                factor = target[j]/row[j]
                target[j] = Decimal(0)
                if factor:
                    for c in range(j+1, size):
                        target[c] -= factor*row[c]
                    b[i] -= factor*b[j]
        x = [Decimal(0)]*size
        for i in range(size-1, -1, -1):
            x[i] = (b[i]-sum((a[i][j]*x[j] for j in range(i+1, size)), Decimal(0)))/a[i][i]
        return x


def run(rotation, output):
    output.mkdir(parents=True, exist_ok=False)
    summary = json.loads((rotation/'summary.json').read_text())
    rows, cache = [], {}
    for side in ('native', 'independent'):
        worst = max(summary['cases'], key=lambda c: c['rotation_covariance'][side]['rotation_max_abs_rad'])
        target = json.loads((rotation/'assemblies'/(worst['id']+'.json')).read_text())
        base_id = f"{len(target['elements'])}-e{target['e']}-identity"
        base = json.loads((rotation/'assemblies'/(base_id+'.json')).read_text())
        solutions = []
        for record in (base, target):
            k = np.array(record[side+'_stiffness']); r = np.array(record[side]['rigid_basis']); f = np.array(record['force'])
            key = (side, record['id'])
            started = time.monotonic()
            solved = []
            for digits in (50, 80):
                x = decimal_saddle(k, r, f, digits)
                solved.append(x)
                write(output/(side+'-'+record['id']+'-'+str(digits)+'.json'), dict(
                    input='exact decimal conversion of saved float64 K, R, F', digits=digits,
                    displacement=[str(v) for v in x[:len(f)]], lagrange_multipliers=[str(v) for v in x[len(f):]],
                    elapsed_seconds=time.monotonic()-started))
            d50, d80 = (np.array([float(v) for v in x[:len(f)]]) for x in solved)
            cache[key] = dict(d=d80, same_float64_roundtrip=bool(np.array_equal(d50, d80)),
                digit50_80_difference_max_abs=max(float(abs(a-b)) for a, b in zip(*solved, strict=True)),
                ordinary_translation_difference_max_abs_m=float(np.max(np.abs((d80-np.array(record[side]['displacement'])).reshape(-1, 6)[:, :3]))),
                ordinary_rotation_difference_max_abs_rad=float(np.max(np.abs((d80-np.array(record[side]['displacement'])).reshape(-1, 6)[:, 3:]))),
                elapsed_seconds=time.monotonic()-started)
            solutions.append(d80)
            print(side, record['id'], digits, cache[key]['same_float64_roundtrip'], flush=True)
        rotation_matrix = np.array(target['rotation'])
        delta = (solutions[1].reshape(-1, 3)@rotation_matrix-solutions[0].reshape(-1, 3)).reshape(-1, 6)
        row = dict(side=side, worst_case=target['id'], base_case=base['id'],
            ordinary_covariance=worst['rotation_covariance'][side],
            decimal_translation_covariance_max_abs_m=float(np.max(np.abs(delta[:, :3]))),
            decimal_rotation_covariance_max_abs_rad=float(np.max(np.abs(delta[:, 3:]))),
            decimal_work_difference_kn_m=dot(target['force'], solutions[1])-dot(base['force'], solutions[0]),
            rounded_matrix_covariance_relative=worst['rotation_covariance'][side]['relative_stiffness_difference'],
            condition_number=worst[side]['reduced_condition_number'],
            arithmetic_scale_condition_times_eps=worst[side]['reduced_condition_number']*np.finfo(float).eps,
            decimal_solves=4,
            limitation='Decimal changes only the solution arithmetic of the saved rounded operators, not kernel arithmetic or the physical model.')
        rows.append(row)
    write(output/'summary.json', dict(diagnostic_only=True, passed=False, b03_complete=False,
        new_physical_equilibrium_cases=0, decimal_saddle_solves=8,
        all_50_80_digit_roundtrips_identical=all(v['same_float64_roundtrip'] for v in cache.values()),
        cases=rows, per_matrix={side+'-'+name: {k: v for k, v in result.items() if k != 'd'} for (side, name), result in cache.items()}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rotation', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.rotation.resolve(), a.output.resolve())
