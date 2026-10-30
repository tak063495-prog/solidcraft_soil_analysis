"""Bind new mechanics results to exact original preparations and independently audited fields."""
from pathlib import Path
import hashlib,json,sys
case,output=map(Path,sys.argv[1:]);runtime=case/'runtime'
sys.path.insert(0,str(runtime/'tools'));import mixed_mcc_io
out=Path(__file__).resolve().parent;p=out.parent

import atexit,tempfile
_cache=tempfile.TemporaryDirectory(prefix='source-blocks-',dir=out)
atexit.register(_cache.cleanup)
sys.path.insert(0,str(p/'python-blocked-integration-v1/runtime/tools'))
from streamed_source_v2 import index_source
_runs=[]
def read(q):
    if q.name=='stage-results.json':
        run=index_source(q,Path(_cache.name)/('result-%d'%len(_runs)))
        _runs.append(run)
        return run
    return json.loads(q.read_text(encoding='utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def verify_sources():
    for run in _runs:run.verify_source()

study=read(case/'study-terminal.json');assert study['passed'] and study['original_three_inputs_byte_exact'] and study['independent_seed_load_and_cutback_audit_passed']
prep=read(case/'preparation.json');assert all(sha(Path(n))==h for n,h in prep['sha256'].items())
original=Path(prep['original_preparation']);old=read(original/'preparation.json')
assert all(sha(Path(n))==h for n,h in old['sha256'].items())
analysis=case/'workflow/analysis';source=original/'workflow/analysis'
for n in ['ground/mesh.json','ground/ground.json','stages.json']:assert sha(analysis/n)==sha(source/n)
assert all(sha(runtime/n)==h for n,h in read(case/'audit-source-manifest.json').items())
assert read(case/'field-verification.json')['seed_load_and_recorded_cutback_gates_passed']
mesh=read(analysis/'ground/mesh.json');result=read(analysis/'stage-results.json');assert result['complete'] and len(result['frames'])==5
mixed_mcc_io.validate(mesh,result)
paths=[Path(__file__),case/'preparation.json',case/'study-terminal.json',case/'field-verification.json',case/'audit-source-manifest.json',original/'preparation.json',original/'wall-spatial-topology-proof.json']
paths.extend(analysis/n for n in ['ground/mesh.json','ground/ground.json','stages.json','stage-results.json'])
paths.extend(runtime/'tools'/n for n in ['tet10_overlap.py','fem_reference.py','mixed_mcc_io.py'])
binding={str(q.resolve()):sha(q) for q in paths}
verify_sources()
with output.open('x') as f:json.dump(dict(passed=True,original_prepared_inputs_byte_exact=True,own_frozen_independent_fields_audited=True,elements=len(mesh['tetrahedra']),sha256=binding,p6_complete=False),f,indent=2)
print('PASS retained preparation and independent result audit',len(mesh['tetrahedra']),flush=True)
