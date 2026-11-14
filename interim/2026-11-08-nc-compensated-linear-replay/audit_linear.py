"""Independent Decimal(80) elimination/residuals on the exact binary64 inputs."""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
D = Decimal.from_float
TOL = Decimal('1e-10')


def norm(x):
    return sum((v*v for v in x), Decimal(0)).sqrt()


def solve(a, b):
    n = len(b)
    a = [list(row)+[b[i]] for i, row in enumerate(a)]
    for j in range(n):
        p = max(range(j,n), key=lambda i: abs(a[i][j]))
        a[p], a[j] = a[j], a[p]
        if not a[j][j]: raise ValueError('singular Decimal system')
        for i in range(j+1,n):
            f = a[i][j]/a[j][j]
            a[i][j] = Decimal(0)
            for k in range(j+1,n+1): a[i][k] -= f*a[j][k]
    x = [Decimal(0)]*n
    for i in range(n-1,-1,-1):
        x[i] = (a[i][n]-sum((a[i][j]*x[j] for j in range(i+1,n)),Decimal(0)))/a[i][i]
    return x


def residual(a, b, x):
    return norm([sum((v*w for v,w in zip(row,x)),Decimal(0))-rhs for row,rhs in zip(a,b)])/max(norm(b),Decimal('1e-30'))


def audit(requests, responses):
    assert len(requests) == len(responses)
    rows = []
    with localcontext() as c:
        c.prec = 80
        for req, res in zip(requests, responses):
            assert req['label'] == res['label']
            n = len(req['reduced_residual'])
            a = [[D(float(req['matrix_columns'][j][i])) for j in range(n)] for i in range(n)]
            b = [D(float(x)) for x in req['reduced_residual']]
            exact = solve(a,b)
            reference_residual = residual(a,b,exact)
            assert reference_residual < Decimal('1e-65')
            row = dict(label=req['label'], reference_relative_residual=str(reference_residual),
                       baseline_original_relative_residual=res['baseline_relative_residual'],
                       candidate_accepted=res['accepted'])
            row['baseline_decimal_relative_residual'] = str(residual(a,b,[D(float(x)) for x in res['baseline_x']]))
            if res['accepted']:
                x = [D(float(x)) for x in res['candidate_x']]
                r = residual(a,b,x)
                row.update(candidate_decimal_relative_residual=str(r),
                           candidate_reference_relative_error=str(norm([u-v for u,v in zip(x,exact)])/norm(exact)),
                           original_relative_residual=res['original_relative_residual'],
                           compensated_relative_residual=res['compensated_relative_residual'],
                           corrections=res['corrections'], independent_linear_gate=bool(r<=TOL and res['original_relative_residual']<=float(TOL)))
            else: row['independent_linear_gate'] = False
            rows.append(row)
    gate = all(r['independent_linear_gate'] for r in rows)
    return dict(precision_decimal_digits=80, exact_binary64_input_conversion=True,
                reference_uses_rust_solution=False, matrix_count=len(rows), rows=rows,
                all_saved_matrices_pass_independent_linear_gate=gate,
                baseline_original_gate_failures=sum(r['baseline_original_relative_residual']>1e-10 for r in rows),
                passed=False, history_converged=False,
                interpretation='linear-system evidence only; no nonlinear or removal acceptance')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--requests',type=Path,required=True);p.add_argument('--responses',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=audit(json.loads(a.requests.read_text()),json.loads(a.responses.read_text()))
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
