"""Fixed-order analysis of the observed rejected rounded linear system only."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
import scipy.linalg as la


def norm(x):return math.sqrt(math.fsum(float(v)*float(v) for v in np.asarray(x).ravel(order='C')))


def residual(a,x,b):
    # Physical row/column order, independent of BLAS reduction order.
    return np.array([math.fsum([*(float(a[i,j])*float(x[j]) for j in range(len(x))),-float(b[i])]) for i in range(len(b))])


def analyze(path):
    value=json.loads(path.read_text(encoding='utf-8'));n=value['rows']
    a=np.asarray(value['column_major_kqq']).reshape((n,n),order='F');b=np.asarray(value['rhs']);x=np.asarray(value['interior'])
    denominator=max(norm(b),1.)
    eigen=la.eigvalsh(a,subset_by_index=[0,n-1])
    fixed=norm(residual(a,x,b))/denominator
    # Observation of the SAME stored matrix/RHS, never a replacement solution.
    factor=la.cho_factor(a,lower=True);direct=la.cho_solve(factor,b)
    rows=[];refined=x.copy()
    for i in range(4):
        rr=residual(a,refined,b)
        rows.append(dict(iteration=i,fixed_order_relative=norm(rr)/denominator,
                         maximum_solution_change=float(np.max(np.abs(refined-x)))))
        if i<3:refined=refined+la.cho_solve(factor,-rr)
    return dict(schema='pcfail-beam-linear-observation/1',passed=False,b01_complete=False,diagnostic_only=True,
        rejected_system_only=True,n=n,factor=value['factor'],original_gate=value['tolerance'],
        native_relative=value['relative'],fixed_order_original_relative=fixed,
        matrix_symmetry_max=float(np.max(np.abs(a-a.T))),minimum_eigenvalue=float(eigen[0]),
        maximum_eigenvalue=float(eigen[-1]),mixed_dof_algebraic_condition=float(eigen[-1]/eigen[0]),
        independent_double_cholesky_fixed_order_relative=norm(residual(a,direct,b))/denominator,
        rounded_matrix_refinement_observations=rows,production_policy_changed=False,
        interpretation='Same rejected rounded matrix/RHS; mixed-DOF algebra only. Refinement is not applied to the solver or accepted as a completed roof path.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    with a.output.open('x',encoding='utf-8') as f:json.dump(analyze(a.input),f,indent=2);f.write('\n')
