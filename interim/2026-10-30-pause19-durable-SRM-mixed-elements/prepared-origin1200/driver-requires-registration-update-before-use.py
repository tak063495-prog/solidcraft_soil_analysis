"""Run only a fully qualified, exactly registered immutable new origin job."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];study=p/'slope1200-epoch-SRM-origin-v1';out=study/'run-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    policy=read(root/'runs/goal-progress-20261027-v1/active-pause-policy.json')
    if (p/'pause-requested').exists() or policy['pause_trigger_observed'] or policy['latest_observed_weekly_usage_percent']>80:
        raise SystemExit('pause requested; no new solver phase')
    return policy
active();quality=p/'epoch-SRM-job-v1/quality-v2';gate=quality/'terminal.json'
if not gate.exists():
    print('Full380Rust219Python qualification is not terminal; no solver started and no run directory created',file=sys.stderr)
    raise SystemExit(78)
qualified=read(gate)
assert qualified['passed'] and (qualified['Rust_tests'],qualified['Python_tests'])==(380,219)
prep=read(study/'preparation.json');job=Path(prep['job']);solver=job/'solver/soil-fem-3d.exe'
assert prep['passed'] and read(job/'status.json')['state']=='prepared' and not (job/'run.claim').exists()
assert sha(job/'job.json')==prep['manifest_sha256'] and sha(solver)==prep['solver_sha256']
assert sha(quality/'source-snapshot.json')==prep['candidate_source_snapshot_sha256']
frozen=read(quality/'source-snapshot.json');runtime=quality.parent/'runtime'
assert all(sha(runtime/n)==h for n,h in frozen.items())
baseline=read(quality.parent/'Root657-source-snapshot.json');assert all(sha(root/n)==h for n,h in baseline.items())
policy=active();registration=read(study/'pause-registration-v2/receipt.json')
assert registration['passed'] and sha(Path(policy['script']))==policy['script_sha256']==registration['new_guard_sha256']
expected=dict(path=str(job),job_id=prep['job_id'],manifest_sha256=prep['manifest_sha256'],solver_sha256=prep['solver_sha256'])
assert policy['exact_registered_epoch_origin_stage_job']==expected
assert all(sha(Path(n))==h for n,h in prep['source_binding'].items())
memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
    'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
    encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024<16*(1<<30) or shutil.disk_usage(study).free<48*(1<<30):
    print('Original resource gate not met; no solver started',file=sys.stderr);raise SystemExit(75)
out.mkdir();driver_sha=sha(Path(__file__))
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out/'preflight.json',dict(passed=True,qualification_receipt_sha256=sha(gate),registration=registration,
    immutable_binding=expected,memory=memory,free_disk_bytes=shutil.disk_usage(study).free,
    prepared_job_has_no_recovery_snapshot_until_first_accepted_save=True,driver_sha256=driver_sha,p6_complete=False))
active();assert sha(Path(__file__))==driver_sha
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
command=[str(solver),'run-stages-job',str(job)];start=time.monotonic()
with (out/'solver.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(command,cwd=job,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out/'solver-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=str(solver),arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(json.dumps(dict(started=True,pid=child.pid,job_id=prep['job_id'])),flush=True);code=child.wait()
status=read(job/'status.json')
assert sha(solver)==prep['solver_sha256'] and sha(job/'job.json')==prep['manifest_sha256']
dump(out/'process-terminal.json',dict(exit_code=code,elapsed_seconds=time.monotonic()-start,status=status,immutable_solver_and_job_unchanged=True,p6_complete=False))
if code==2:
    assert status['state']=='checkpointed';print('Origin paused; inspect saved accepted recovery state',flush=True);raise SystemExit(2)
if code:raise SystemExit(code)
assert status['state']=='succeeded' and status['artifact']=='result.json'
assert sha(job/'result.json')==status['artifact_sha256']
result=read(job/'result.json');assert result['complete'] and result['provenance']['solver_binary_sha256']==prep['solver_sha256']
assert all(sha(runtime/n)==h for n,h in frozen.items()) and all(sha(root/n)==h for n,h in baseline.items())
dump(out/'terminal.json',dict(passed=True,execution_and_frozen_complete_origin_only=True,
    independent_all_GP_force_yield_and_original_field_audit_still_required=True,
    original_2D5percent_gate_and_fresh_replay_reference_not_replaced=True,origin_result_sha256=sha(job/'result.json'),
    before_epoch_search_require_separate_exact_epoch_guard_registration=True,p6_complete=False))
print('New1200 origin execution complete; independent mechanics audit and epoch job registration still required',flush=True)
