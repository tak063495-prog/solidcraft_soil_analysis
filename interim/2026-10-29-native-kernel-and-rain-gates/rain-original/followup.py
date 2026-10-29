"""Actual fresh-process recovery controls, then predeclared FV/RT0 gates."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time

p = Path(__file__).resolve().parent
root = p.parents[1]
read = lambda q: json.loads(q.read_text('utf-8-sig'))
sha = lambda q: hashlib.sha256(q.read_bytes()).hexdigest()
env = dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')


def write(path,value):
    with path.open('x',encoding='utf-8') as f:
        json.dump(value,f,indent=2,allow_nan=False)


def wait(path):
    while not path.exists():
        if (p/'pause-requested').exists():
            raise SystemExit(2)
        for q in p.rglob('*-terminal.json'):
            r = read(q)
            if 'exit_code' in r:
                assert r['exit_code']==r.get('expected_exit_code',0),(str(q),r)
        time.sleep(5)
    return read(path)


def call(label,args,folder,expected=0):
    if (p/'pause-requested').exists():
        raise SystemExit(2)
    started = time.monotonic()
    command = list(map(str,args))
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child = subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(folder/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code = child.wait()
    row = dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-started)
    write(folder/(label+'-terminal.json'),row)
    print(json.dumps(row),flush=True)
    assert code==expected,row


write(p/'followup-driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
procedures = {name:sha(p/name) for name in ['compare.py','followup.py']}
write(p/'followup-procedures.json',procedures)
assert wait(p/'rt0-4-dt100/terminal.json')['passed']
sys.path.insert(0,str(p/'runtime/tools'))
from verify_richards_cases import normalized
import numpy as np
# Check interpolation independently on analytic affine fields and exact centres.
node = next(n for n in ast.parse((p/'compare.py').read_text('utf-8')).body if isinstance(n,ast.FunctionDef) and n.name=='interpolate')
namespace = dict(np=np)
exec(compile(ast.Module(body=[node],type_ignores=[]),str(p/'compare.py'),'exec'),namespace)
interpolate = namespace['interpolate']
z,y,x = np.meshgrid(*([(np.arange(4)+.5)/4]*3),indexing='ij')
field = 2*z-3*y+5*x+.7
np.testing.assert_allclose(interpolate(field,(np.arange(4)+.5)/4),field,rtol=0,atol=1e-14)
t = (np.arange(8)+.5)/8
zc,yc,xc = np.meshgrid(*([np.clip(t,.125,.875)]*3),indexing='ij')
np.testing.assert_allclose(interpolate(field,t),2*zc-3*yc+5*xc+.7,rtol=0,atol=1e-14)
np.testing.assert_allclose(interpolate(np.full((4,4,4),3.2),t),3.2,rtol=0,atol=1e-14)
controls = p/'checkpoint-controls'
controls.mkdir()
solver = p/'solver/soil-fem-3d.exe'
folder = p/'rt0-4-dt100'
original = {name:sha(folder/name) for name in ['mesh.json','case.json','result.json']}
solver_sha = read(p/'scope.json')['solver_sha256']
assert sha(solver)==solver_sha
baseline = read(folder/'result.json')
counts = [1] + [max(i for i,f in enumerate(baseline['frames']) if f['name']==period)+1 for period in ['rain','pause']]
rows = []
for number,count in enumerate(counts):
    assert 0<count<len(baseline['frames'])-1
    cp = controls/f'checkpoint-{number}.json'
    dest = controls/f'resumed-{number}.json'
    call(f'checkpoint-{number}',[solver,'checkpoint-richards',folder/'mesh.json',folder/'case.json',count,cp],controls)
    saved = read(cp)
    assert saved['run']['runtime'] and not saved['run']['complete']
    digest = sha(cp)
    call(f'resume-{number}',[solver,'resume-richards',folder/'mesh.json',folder/'case.json',cp,dest],controls)
    assert normalized(read(dest))==normalized(baseline)
    assert sha(cp)==digest and sha(solver)==solver_sha
    rows.append(dict(count=count,time_s=saved['run']['last_committed']['time_s'],phase=saved['run']['frames'][-1]['name'],all_states_frames_attempts_exact_except_elapsed_seconds=True,checkpoint_unchanged=True))
bad = read(controls/'checkpoint-0.json')
bad['run']['last_committed']['cell_storage_offset_m3'][0] += 1e-6
write(controls/'tampered.json',bad)
out = controls/'rejected.json'
call('tamper-rejected',[solver,'resume-richards',folder/'mesh.json',folder/'case.json',controls/'tampered.json',out],controls,expected=1)
assert not out.exists()
assert original=={name:sha(folder/name) for name in original}
write(controls/'terminal.json',dict(passed=True,solver_sha256=solver_sha,restarts=rows,tampered_payload_rejected_without_output=True,original_inputs_and_continuous_result_unchanged=True,interpolation_affine_constant_exact_centres_checks=True))
assert wait(p/'raw-terminal.json')['passed']
assert procedures=={name:sha(p/name) for name in procedures}
comparison = subprocess.Popen([sys.executable,'-X','utf8',str(p/'compare.py')],cwd=root,env=env,
    stdout=(p/'comparison.log').open('x',encoding='utf-8'),stderr=subprocess.STDOUT,
    creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
write(p/'comparison-child.json',dict(pid=comparison.pid,parent_pid=os.getpid(),executable=sys.executable,arguments=[str(p/'compare.py')],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
code = comparison.wait()
assert code in (0,1)
report = read(p/'comparison.json')
assert (code==0)==report['passed']
write(p/'study-terminal.json',dict(passed=report['passed'],raw_runs_and_discrete_audits_passed=True,
    independent_reference_tests_passed=3,local_rain_checkpoint_controls_passed=True,physical_refinement_gates_passed=report['passed'],
    comparison_sha256=sha(p/'comparison.json'),comparison_exit_code=code,failed_gates_preserved=True,production_unmodified=True,p6_complete=False))
print('Local-rain study finished; refinement accepted:',report['passed'],flush=True)
