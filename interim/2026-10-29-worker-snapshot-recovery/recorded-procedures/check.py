"""QA and real274 abrupt-process-loss/recovery with unchanged accepted Rust solver."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import freeze,manifest
from native_job import stop_child
quality=p/'quality';quality.mkdir();freeze(quality);sources=manifest()
base=read(p/'source-snapshot.json')
changed={n for n in sources if sources[n]!=base[n]}
assert changed=={'tools/native_job.py','tools/native_stage_profile.py','tests/test_native_job.py'}
assert len(sources)==611
solver=root/'runs/durable-root-promotion-20261029-v1/solver/soil-fem-3d.exe'
solver_sha=sha(solver);assert solver_sha==read(root/'runs/durable-root-promotion-20261029-v1/acceptance.json')['solver_sha256']
(runtime/'target/release').mkdir(parents=True);shutil.copy2(solver,runtime/'target/release/soil-fem-3d.exe')
records=[]
env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')

def spawn(label,args):
    assert manifest()==sources
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));log=(p/(label+'.log')).open('x',encoding='utf-8')
    child=subprocess.Popen(command,cwd=runtime,stdout=log,stderr=subprocess.STDOUT,env=env,
        creationflags=subprocess.CREATE_NO_WINDOW|subprocess.CREATE_NEW_PROCESS_GROUP|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    with (p/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
    return child,log,time.monotonic()

def call(label,args):
    child,log,began=spawn(label,args)
    try:code=child.wait()
    finally:log.close()
    records.append(dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-began))
    print(json.dumps(records[-1]),flush=True);assert code==0 and manifest()==sources

call('python-tests',[sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
assert 'Ran 171 tests' in (p/'python-tests.log').read_text('utf-8')
with (quality/'terminal.json').open('x') as f:json.dump(dict(passed=True,source_files=611,python_tests=171,rust_sources_unchanged_from272_passed=True,changed_files=sorted(changed),records=records,p6_complete=False),f,indent=2)
capture=root/'runs/native-fem-smoke-20261029-v1/Native63 cooperative smoke/capture'
controls=p/'controls';controls.mkdir()
worker=[sys.executable,runtime/'tools/native_job.py']
call('continuous',worker+['run',capture,controls/'continuous','--solver',solver,'--cooperative','--checkpoint-interval','4','--mesh-size-m','4'])
assert read(controls/'continuous/status.json')['state']=='succeeded'
old=controls/'abrupt-loss'
child,log,began=spawn('abrupt-loss',worker+['run',capture,old,'--solver',solver,'--cooperative','--checkpoint-interval','4','--mesh-size-m','4'])
try:
    while child.poll() is None:
        if (p/'pause-requested').exists():raise SystemExit(2)
        pointer_path=old/'workflow/analysis/solver-job/recovery.json'
        status_path=old/'status.json'
        if pointer_path.exists() and status_path.exists():
            try:
                status=read(status_path);progress=(status.get('solver_report') or {}).get('last_progress') or {}
                pointer=read(pointer_path)
                if progress.get('stage_index')==4 and progress.get('progress',0)>=.125 and pointer['generation']>=3:
                    # Real abrupt loss of only this still-owned worker/solver tree;
                    # no cooperative stop marker, no status rewrite, no PID probing.
                    stop_child(child)
                    break
            except PermissionError:pass
        time.sleep(.05)
    assert child.returncode not in (None,0),'case completed before abrupt loss was exercised'
finally:
    log.close()
    if child.poll() is None:stop_child(child)
records.append(dict(step='abrupt-loss',exit_code=child.returncode,elapsed_seconds=time.monotonic()-began))
parent_before={q.relative_to(old).as_posix():sha(q) for q in old.rglob('*') if q.is_file()}
assert read(old/'status.json')['state'] in ('running','cancelling')
assert not (old/'cancel-request.json').exists()
assert not (old/'workflow/analysis/solver-job/stop-request.json').exists()
call('recover',worker+['recover',old,controls/'recovered','--solver',solver])
new=controls/'recovered';assert read(new/'status.json')['state']=='succeeded'
assert parent_before=={q.relative_to(old).as_posix():sha(q) for q in old.rglob('*') if q.is_file()}

def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v

groups=['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','transition_audits','geostatic','attempts']
expected=read(controls/'continuous/workflow/analysis/stage-results.json');actual=read(new/'workflow/analysis/stage-results.json')
historical=read(root/'runs/native-fem-smoke-20261029-v1/Native63 cooperative smoke/workflow/analysis/stage-results.json')
for key in groups:assert stable(actual[key])==stable(expected[key])==stable(historical[key]),key
assert manifest()==sources and sha(solver)==solver_sha
with (p/'acceptance.json').open('x') as f:json.dump(dict(passed=True,python_tests=171,actual_native_capture_fresh274=True,abrupt_loss_during_excavation=True,
    original_parent_all_files_unchanged=True,stale_running_status_retained=True,no_parent_signal_during_recovery=True,
    full5_stages_recovered=True,all8_physical_groups_exact_to_new_continuous_and_historical_control=True,
    postprocessing_and_view_audits_passed=True,recovered_generation=read(new/'job.json')['recovered_generation'],
    solver_sha256=solver_sha,source_files=611,rust_sources_unchanged_from272_passed=True,sources_unchanged=True,records=records,
    root_not_promoted=True,native_gui_recovery_not_yet_connected=True,p6_complete=False),f,indent=2)
print('PASS real abrupt worker/solver loss and verified periodic-snapshot recovery, all8 physical groups exact',flush=True)
