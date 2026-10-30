"""Register and seal the finer inputs; no nonlinear solve is started."""
from pathlib import Path
import ast,datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
out=p/'fixed-wall-finer-study-v1';out.mkdir()
analysis=out/'mesh34350/workflow/analysis';analysis.mkdir(parents=True)
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
inputs=p/'fixed-wall-finer-input-v1'
capacity=p/'fixed-wall-finer-capacity-preflight-v1'
assert read(inputs/'terminal.json')['passed']
report=read(capacity/'terminal.json');assert report['passed']
assert report['mixed_primary_system_dimension']==149006
assert report['original_prepared_case_unchanged'] and report['only_maximum_free_dofs_configuration_changed']
assert sha(capacity/'stages-budget200k.json')==report['preview_case_sha256']
assert sha(inputs/'workflow/analysis/stages.json')==report['original_prepared_case_sha256']
assert sha(inputs/'workflow/analysis/ground/mesh.json')==report['original_mesh_sha256']
solver=root/'runs/mcc-mixed-assembly-development-20261029-v1/solver/soil-fem-3d.exe'
assert sha(solver)=='f2ab7b5a1ea8684d2f05f6edc1e35e51258b9f9890ac3bc9b0d849a370e00f27'
policy_path=root/'runs/goal-progress-20261027-v1/active-pause-policy.json'
policy=read(policy_path);guard=Path(policy['script']);before=guard.read_bytes()
assert hashlib.sha256(before).hexdigest()==policy['script_sha256']=='c03a656d6e74c2fcaffca291ea8717b1317c1e36186560e5b4477d9c95092f67'
assert policy['goal_active'] and not policy['explicit_user_pause'] and not policy['pause_trigger_observed']
assert policy['threshold_percent']==80 and policy['threshold_comparison']=='strictly_greater_than'
assert policy['remaining_percent_trigger']==20 and policy['remaining_percent_comparison']=='strictly_less_than'
assert read(p/'Root635-default-owned-workflow-v1/pause-registration.json')['passed']
job=analysis/'solver-job';relative=job.relative_to(root).as_posix()
text=before.decode('utf-8')
anchor="    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/Root635-default-owned-workflow-v1/run/owned-case/job/workflow/analysis/solver-job')\n"
line="    jobs.append(root/'"+relative+"')\n"
assert text.count(anchor)==1 and relative not in text and not job.exists()
updated=text.replace(anchor,anchor+line,1);assert updated.replace(line,'',1)==text;ast.parse(updated)
(out/'previous-guard.py').write_bytes(before)
(out/'previous-policy.json').write_bytes(policy_path.read_bytes())
pending=guard.with_suffix('.fixed34350-pending.py')
with pending.open('x',encoding='utf-8',newline='\n') as f:f.write(updated)
pending.replace(guard)
policy.update(script_sha256=sha(guard),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
pending=policy_path.with_suffix('.fixed34350-pending.json')
with pending.open('x',encoding='utf-8') as f:json.dump(policy,f,indent=2,allow_nan=False)
pending.replace(policy_path)
dump(out/'pause-registration.json',dict(passed=True,job=relative,old_guard_sha256=hashlib.sha256(before).hexdigest(),
    new_guard_sha256=sha(guard),only_one_exact_job_path_added=True,
    all_previous_scope_preservers_and_user_threshold_unchanged=True,registered_before_job_creation=True))
(analysis/'ground').mkdir()
for name in ['mesh.json','ground.json']:
    src=inputs/'workflow/analysis/ground'/name;dst=analysis/'ground'/name
    shutil.copyfile(src,dst);assert sha(src)==sha(dst)
shutil.copyfile(capacity/'stages-budget200k.json',analysis/'stages.json')
bindings={str(q):sha(q) for q in [solver,inputs/'terminal.json',inputs/'wall-spatial-topology-proof.json',
    capacity/'terminal.json',capacity/'stages-budget200k.json',Path(__file__)]}
bindings.update({str(analysis/n):sha(analysis/n) for n in ['ground/mesh.json','ground/ground.json','stages.json']})
dump(out/'source-binding.json',bindings)
while True:
    active()
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=20*(1<<30) and shutil.disk_usage(out).free>=64*(1<<30):break
    time.sleep(5)
dump(out/'preparation-resource-preflight.json',dict(passed=True,minimum_each_free_bytes=20*(1<<30),
    free_physical_bytes=memory['FreePhysicalMemory']*1024,free_virtual_bytes=memory['FreeVirtualMemory']*1024,
    free_disk_bytes=shutil.disk_usage(out).free,nonlinear_solve_resource_preflight_is_separate=True))
active()
command=[str(solver),'prepare-durable-stages-job',str(analysis/'ground/mesh.json'),
    str(analysis/'ground/ground.json'),str(analysis/'stages.json'),'4',str(job)]
start=time.monotonic()
with (out/'prepare-job.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,
        creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out/'prepare-job-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
        arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code=child.wait()
dump(out/'prepare-job-terminal.json',dict(exit_code=code,elapsed_seconds=time.monotonic()-start))
assert code==0
manifest=read(job/'job.json');status=read(job/'status.json')
assert manifest['binding']['checkpoint_interval']==4
assert all(sha(job/n)==h for n,h in manifest['binding']['files_sha256'].items())
assert status['state']=='prepared' and status['process_id'] is None and status['progress_events']==0
assert not (job/'progress.jsonl').exists()
assert all(sha(Path(q))==h for q,h in bindings.items())
dump(out/'preparation.json',dict(passed=True,prepared=True,elements=34350,integration_points=927450,
    nodes=46375,mixed_primary_dimension_static_inventory=149006,configured_primary_budget=200000,
    same_F2_solver_as_completed14909_and8061=True,checkpoint_interval=4,
    all440_wall_nodes_and832_oriented_plates_exact=True,
    all_material_load_boundary_and_numerical_controls_retained_except_size_budget=True,
    all_sealed_native_input_hashes_verified=True,exact_job_registered_before_creation=True,
    no_nonlinear_solve_started=True,no_history_interpolation=True,
    future_solve_requires_completed_large_phases_and_fresh_memory_gate=True,
    full_solve_full_independent_audit_and_same_original_spatial_gates_still_required=True,p6_complete=False))
print('PASS sealed34350 job with exact fixed wall and only larger size budget; nonlinear solve not started',flush=True)
