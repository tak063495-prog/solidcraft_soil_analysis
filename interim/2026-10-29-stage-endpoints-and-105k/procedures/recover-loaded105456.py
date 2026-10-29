"""Recompute unloading from the retained fully plastic loaded105456 state.

The newest slot is already at final unloading and would test no remaining
solve. Select the preceding sealed slot (load=1) without changing its payload.
"""
from pathlib import Path
import datetime
import gc
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

p = Path(__file__).resolve().parent
main = p.parent
root = main.parents[1]
original = main/'n26'
source = original/'solver-job'
out = p/'n26'
read = lambda q: json.loads(q.read_text('utf-8-sig'))


def sha(q):
    h = hashlib.sha256()
    with q.open('rb') as f:
        for block in iter(lambda:f.read(4<<20),b''):
            h.update(block)
    return h.hexdigest()


def write(path,value):
    with path.open('x',encoding='utf-8') as f:
        json.dump(value,f,indent=2,allow_nan=False)


def binding(folder):
    return {q.relative_to(folder).as_posix():sha(q) for q in folder.rglob('*') if q.is_file()}


def stable(v):
    if isinstance(v,dict):
        return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):
        return [stable(x) for x in v]
    return v


write(p/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
assert read(original/'terminal.json')['passed'] and read(original/'analytical-verification.json')['passed']
assert read(original/'process-terminal.json')['exit_code']==0
solver = source/'solver/soil-fem-3d.exe'
expected = 'd55ffe42c3b61e0891fec513489c9f5881830733ac9e003a01a0cf69f495905d'
assert sha(solver)==expected
before = binding(source)
pointer = read(source/'recovery.json')
assert pointer['generation']==4 and pointer['checkpoint']=='recovery-checkpoint-0.json'
slot = source/'recovery-checkpoint-1.json'
saved = read(slot)
c = saved['run']['continuation']
assert c['stage_index']==0 and c['progress']==1. and not saved['run']['complete']
assert saved['run']['last_committed']['stage_progress']==1.
assert [f['name'] for f in saved['run']['frames']]==['initial_equilibrated']
chosen = dict(pointer,generation=3,checkpoint=slot.name,checkpoint_sha256=sha(slot),run_sha256=saved['run_sha256'])
write(p/'origin.json',dict(passed=True,generation=3,checkpoint=slot.name,checkpoint_sha256=sha(slot),
    selected_fully_loaded_plastic_state=True,remaining_unloading_must_be_recomputed=True,latest_final_state_not_used=True,
    continuation_stage_index=c['stage_index'],continuation_progress=c['progress'],solver_sha256=expected))
del saved,c
gc.collect()
immutable = p/'immutable-source'


def ignore(folder,names):
    excluded = ['result.json','status.json','progress.jsonl','run.claim','cancel.request','cancel.json','recovery-checkpoint-0.json']
    return [name for name in names if name in excluded]


shutil.copytree(source,immutable,ignore=ignore)
(immutable/'recovery.json').write_text(json.dumps(chosen,indent=2)+'\n','utf-8')
immutable_hashes = binding(immutable)
out.mkdir()
for name in ['mesh.json','ground.json','stages.json','preparation.json']:
    shutil.copy2(original/name,out/name)
shutil.copy2(main/'verify.py',p/'verify.py')
write(p/'frozen-procedures.json',dict(verifier_sha256=sha(p/'verify.py'),driver_sha256=sha(Path(__file__)),original_job_hashes=before,immutable_origin_hashes=immutable_hashes))
env = dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
records = []


def call(label,args):
    if (main/'pause-requested').exists():
        raise SystemExit(2)
    command = list(map(str,args))
    started = time.monotonic()
    with (p/(label+'.log')).open('x',encoding='utf-8') as log:
        child = subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(p/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code = child.wait()
    row = dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-started)
    records.append(row)
    write(p/(label+'-terminal.json'),row)
    print(json.dumps(row),flush=True)
    assert code==0,row


job = out/'solver-job'
call('prepare',[solver,'prepare-stages-recovery',immutable,job])
call('solve',[job/'solver/soil-fem-3d.exe','run-stages-job',job])
call('analytical',[sys.executable,'-X','utf8',p/'verify.py','26'])
actual = read(job/'result.json')
expected_result = read(source/'result.json')
assert set(actual)==set(expected_result)
keys = {name:stable(actual[name])==stable(expected_result[name]) for name in actual}
write(p/'comparison.json',dict(passed=all(keys.values()),keys_exact_except_elapsed_seconds=keys,original_result_sha256=sha(source/'result.json'),recovered_result_sha256=sha(job/'result.json')))
assert all(keys.values()),[name for name,ok in keys.items() if not ok]
assert binding(source)==before and binding(immutable)==immutable_hashes
write(p/'terminal.json',dict(passed=True,tetrahedra=105456,solver_sha256=sha(solver),loaded_generation3_used=True,
    unloading_recomputed=True,all_report_fields_exact_except_elapsed_seconds=True,all_nodes_and_points_analytical_audit_passed=True,
    original_parent_and_immutable_origin_unchanged=True,records=records,
    scope='Uniform hydrostatic nonlinear MCC patch; no CAD/excavation/fluid/structures.',p6_complete=False))
print('PASS105456 loaded plastic state recovery, recomputed unload, every report field and theory',flush=True)
