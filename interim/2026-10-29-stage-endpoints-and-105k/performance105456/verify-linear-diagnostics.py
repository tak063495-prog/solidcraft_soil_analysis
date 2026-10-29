"""Check the solver's recorded original-matrix residuals; no matrix re-solve."""
from pathlib import Path
import argparse,hashlib,json,math
p=argparse.ArgumentParser();p.add_argument('progress',type=Path);p.add_argument('record',type=Path);p.add_argument('new_report',type=Path);a=p.parse_args()
events=[json.loads(line) for line in a.progress.read_text('utf-8').splitlines() if line.strip()]
record=json.loads(a.record.read_text('utf-8'))
assert record['passed'] and [e['sequence'] for e in events]==list(range(1,len(events)+1))
assert len({e['job_id'] for e in events})==1
assert hashlib.sha256(a.progress.read_bytes()).hexdigest()==record['progress_sha256']
reports=[e['progress']['linear_solve'] for e in events if e['progress']['kind']=='linear_solved']
assert len(reports)==record['linear_solves']==80
assert all(r['method']=='csr_gmres_ilu0' and math.isfinite(r['true_relative_residual']) and r['true_relative_residual']<=record['requested_relative_tolerance'] for r in reports)
value=dict(passed=True,linear_solves=len(reports),maximum_reported_true_relative_residual=max(r['true_relative_residual'] for r in reports),
    maximum_dimension=max(r['dimension'] for r in reports),maximum_iterations=max(r['iterations'] for r in reports),
    reported_solver_diagnostics_only=True,independent_matrix_residual_reconstruction=False)
assert value['maximum_reported_true_relative_residual']==record['maximum_original_matrix_true_relative_residual']
with a.new_report.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2)
print('PASS recorded true-residual budget',value['maximum_reported_true_relative_residual'])
