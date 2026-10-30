"""Run one registered immutable slope job; independent acceptance is separate."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent
root=p.parents[1]
study=p/'plane-strain-slope-comparison-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
assert len(sys.argv)==2 and sys.argv[1] in ['xy1-ny1-global64','xy1-ny1-local-half']
case_name=sys.argv[1]
folder=study/case_name
job=folder/'solver-job'
out=folder/'solve-v1';out.mkdir()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
registration=read(study/'increment-sensitivity-preparation-v1/pause-registration.json')
policy=read(root/'runs/goal-progress-20261027-v1/active-pause-policy.json')
assert registration['passed'] and registration['new_guard_sha256']==policy['script_sha256']
prior=(study/'increment-sensitivity-preparation-v1/previous-guard.py').read_text('utf-8')
guard=Path(policy['script']).read_text('utf-8')
added=''.join("    jobs.append(root/'"+name+"')"+chr(10) for name in registration['paths'])
assert guard.count(added)==1 and guard.replace(added,'',1)==prior
assert read(study/'increment-sensitivity-preparation-v1/terminal.json')['passed']
assert sha(Path(policy['script']))==policy['script_sha256']
assert str(job.relative_to(root)).replace('\\','/') in registration['paths']
assert read(study/'job-preparation-v2/terminal.json')['passed']
assert read(study/'reference-2d-v2/terminal.json')['passed']
assert read(study/'native-input-validation-v1/terminal.json')['passed']
manifest=read(job/'job.json');status=read(job/'status.json')
assert manifest['binding']['checkpoint_interval']==4
assert status['state']=='prepared' and status['process_id'] is None and status['progress_events']==0
assert not (job/'progress.jsonl').exists()
solver=job/'solver/soil-fem-3d.exe'
expected='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
assert sha(solver)==expected==manifest['binding']['identity']['solver_binary_sha256']
binding={str(job/name):h for name,h in manifest['binding']['files_sha256'].items()}
prepared=read(folder/'prepared.json')
binding.update({str(folder/name):h for name,h in prepared['input_sha256'].items()})
binding.update({str(q):sha(q) for q in [job/'job.json',study/'protocol.json',Path(__file__)]})
assert all(sha(Path(q))==h for q,h in binding.items())
dump(out/'source-binding.json',binding)
# Release one actual solve slot before either new increment study starts.
gate=study/'xy1-ny1/recovery-equality-v1/resume-terminal.json'
while not gate.exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    time.sleep(5)
assert read(gate)['passed'],'original accepted-state replay must finish first'
while True:
    if (p/'pause-requested').exists():raise SystemExit(2)
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    required=(10 if case_name=='xy1-ny1' else 16)*(1<<30)
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=required and shutil.disk_usage(folder).free>=32*(1<<30):break
    time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,free_physical_bytes=memory['FreePhysicalMemory']*1024,
    free_virtual_bytes=memory['FreeVirtualMemory']*1024,minimum_each_free_bytes=required,
    free_disk_bytes=shutil.disk_usage(folder).free,simultaneous_exclusive_machine_or_performance_not_claimed=True))
command=[str(solver),'run-stages-job',str(job)]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='1')
start=time.monotonic();requested=False
with (out/'solve.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,
        creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out/'solver-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
        arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    while child.poll() is None:
        if (p/'pause-requested').exists() and not requested:
            stop=[str(solver),'stop-stages-job',str(job)]
            with (out/'cooperative-stop.log').open('x',encoding='utf-8') as stoplog:
                code=subprocess.call(stop,cwd=folder,stdout=stoplog,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
            dump(out/'stop-request-terminal.json',dict(exit_code=code,matching_job_cooperative_stop_requested=True))
            assert code==0
            requested=True
        time.sleep(.5)
    code=child.returncode
status=read(job/'status.json')
assert all(sha(Path(q))==h for q,h in binding.items())
report=dict(passed=code==0,exit_code=code,elapsed_seconds=time.monotonic()-start,status=status,
    case=case_name,solver_sha256=expected,exact_inputs_and_reference_factor_unchanged=True,
    cooperative_pause_requested=requested,checkpointed=code==2,
    numerical_failure_is_not_a_physical_Fs_upper_bound=True,
    independent_force_VTU_all661_nodes_all200_prism_stresses_and_criteria_comparison_pending=True,p6_complete=False)
if code==0:
    assert status['state']=='succeeded' and status['artifact']=='result.json'
    result=job/'result.json';assert sha(result)==status['artifact_sha256']
    output=folder/'stage-results.json';assert not output.exists()
    shutil.copyfile(result,output);assert sha(output)==sha(result)
    report.update(result=str(output),result_sha256=sha(output),result_bytes=output.stat().st_size)
elif code==2:
    assert status['state']=='checkpointed'
dump(out/'terminal.json',report)
print(json.dumps(dict(case=case_name,passed=report['passed'],exit_code=code,status=status['state'],
    elapsed_seconds=report['elapsed_seconds'],independent_acceptance_pending=True)),flush=True)
if code:raise SystemExit(code)
