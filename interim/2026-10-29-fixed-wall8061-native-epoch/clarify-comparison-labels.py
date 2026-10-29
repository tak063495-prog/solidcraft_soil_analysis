"""Retain original evidence and correct copied report labels against actual inputs."""
from pathlib import Path
import hashlib,json,datetime
p=Path(__file__).resolve().parent
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
r=read(p/'comparison.json')
counts=[];bindings={}
for n in [5144,8061]:
    case=p/f'mesh{n}'
    meshpath=case/'workflow/analysis/ground/mesh.json'
    mesh=read(meshpath);counts.append(len(mesh['tetrahedra']))
    assert len(mesh['tetrahedra'])==n and read(case/'study-terminal.json')['passed']
    preflight=read(p/('left-preflight.json' if n==5144 else 'right-preflight.json'))
    assert preflight['passed']
    for q in [meshpath,case/'study-terminal.json',case/'field-verification.json']:
        bindings[str(q)]=sha(q)
assert [r['left_elements'],r['right_elements']]==[2865,5144]
assert read(p/'terminal.json')['passed']
source=(p/'compare.py').read_text('utf-8')
assert "==[5144,8061]" in source and "left=p/'mesh5144'" in source and "right=p/'mesh8061'" in source
corrected=dict(r,left_elements=counts[0],right_elements=counts[1])
corrected['label_correction'] = dict(original_sha256=sha(p/'comparison.json'),
    original_declared_counts=[2865,5144],actual_input_counts=counts,
    explanation='Copied summary labels were stale. Input paths, asserted mesh sizes and all numerical comparisons used 5144/8061.',
    numerical_rows_unchanged=True,independent_nonlinear_solution_comparison_complete=False)
out=p/'comparison-corrected-labels.json'
with out.open('x',encoding='utf-8') as f:json.dump(corrected,f,indent=2,allow_nan=False)
bindings[str(p/'comparison.json')]=sha(p/'comparison.json')
bindings[str(p/'compare.py')]=sha(p/'compare.py')
with (p/'comparison-label-correction.json').open('x',encoding='utf-8') as f:
    json.dump(dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        actual_counts=counts,original_report_unchanged=True,original_failure_retained=True,
        corrected_report_sha256=sha(out),bindings=bindings,p6_complete=False),f,indent=2)
print('PASS report labels corrected in a separate file; numerical results and original evidence retained')
