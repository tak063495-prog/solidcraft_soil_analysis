"""Execute the sealed larger study after serialized phases and verified IO."""
from pathlib import Path,PurePosixPath
import datetime,hashlib,json,os,shutil,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1]
study=p/'fixed-wall-finer-study-v1';meshroot=study/'mesh34350'
out=study/'run-v2';out.mkdir()
analysis=meshroot/'workflow/analysis';job=analysis/'solver-job'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
 if (p/'pause-requested').exists():raise SystemExit(2)
def gate(q,require_pass=True):
 while not q.exists():active();time.sleep(5)
 result=read(q)
 if require_pass:assert result['passed'],str(q)
 return result
active()
prep=read(study/'preparation.json');assert prep['passed'] and prep['no_nonlinear_solve_started']
binding=read(study/'source-binding.json');assert all(sha(Path(q))==h for q,h in binding.items())
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
 script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out/'protocol.json',dict(prepared=True,elements=34350,integration_points=927450,stages=5,
 mixed_primary_static_inventory=149006,configured_primary_budget=200000,
 only_size_budget_changed_from_original_inputs=True,all_original_physics_and_numerical_accuracy_controls_retained=True,
 all440_wall_nodes_and832_oriented_plates_fixed=True,checkpoint_interval=4,
 original_F2_solver_not_recompiled=True,full_independent_package_and_original_seed_load_history_audits_required=True,
 prior_large_phases_serialized=True,storage_adapter_must_pass_full8061_report_exactness=True,
 spatial_comparison_and_Native_display_and_general100k_and_allP6_separate=True,p6_complete=False))
# Ordering barriers may have failed field-accuracy criteria. Their failures must
# remain visible and cannot block this independent spatial investigation.
for q in [p/'global-volume-comparison-v1/manager-v1/terminal.json',
 p/'plane-strain-slope-comparison-v1/xy2-ny2/field-comparison-v1/terminal.json',
 p/'plane-strain-slope-comparison-v1/xy1-ny1-local-half/field-comparison-v1/terminal.json']:
 gate(q,False)
gate(p/'geometry-audit-real14909-v1/terminal.json')
gate(p/'geometry-audit-real14909-v1/worker-passed.json')
gate(p/'storage-backed-fields8061-v1/adapter-terminal.json')
adapter=p/'storage-backed-field-audit-v1'
adapter_prep=read(adapter/'preparation.json')
assert sha(adapter/'audit-fields.py')==adapter_prep['storage_driver_sha256']
validation=p/'storage-backed-fields8061-v1'
validated=read(validation/'field-verification.json')
assert validated['sha256'][str(validation/'audit-fields.py')]==adapter_prep['storage_driver_sha256']
assert sha(validation/'audit-fields.py')==adapter_prep['storage_driver_sha256']
policy=read(root/'runs/goal-progress-20261027-v1/active-pause-policy.json')
assert sha(Path(policy['script']))==policy['script_sha256']
assert "jobs.append(root/'"+job.relative_to(root).as_posix()+"')" in Path(policy['script']).read_text('utf-8')
status=read(job/'status.json');manifest=read(job/'job.json')
assert status['state']=='prepared' and status['process_id'] is None and status['progress_events']==0
assert manifest['binding']['checkpoint_interval']==4
assert all(sha(job/n)==h for n,h in manifest['binding']['files_sha256'].items())
solver=job/'solver/soil-fem-3d.exe'
assert sha(solver)=='f2ab7b5a1ea8684d2f05f6edc1e35e51258b9f9890ac3bc9b0d849a370e00f27'
while True:
 active()
 memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
  'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
  encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
 if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=40*(1<<30) and shutil.disk_usage(out).free>=96*(1<<30):break
 time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,minimum_each_free_bytes=40*(1<<30),minimum_free_disk_bytes=96*(1<<30),
 free_physical_bytes=memory['FreePhysicalMemory']*1024,free_virtual_bytes=memory['FreeVirtualMemory']*1024,
 free_disk_bytes=shutil.disk_usage(out).free,earlier_large_phase_ordering_barriers_completed=True))
runtime=meshroot/'runtime';runtime.mkdir()
quality=p/'Root-source635-promotion-v1/quality-v1'
sources=read(quality/'source-snapshot.json')
assert len(sources)==635
with zipfile.ZipFile(quality/'source-snapshot.zip') as z:
 assert z.testzip() is None and len(z.namelist())==635 and set(z.namelist())==set(sources)
 for name,h in sources.items():
  q=PurePosixPath(name);assert not q.is_absolute() and '..' not in q.parts
  raw=z.read(name);assert hashlib.sha256(raw).hexdigest()==h
  target=runtime/name;target.parent.mkdir(parents=True,exist_ok=True)
  with target.open('xb') as f:f.write(raw)
dump(out/'runtime635-source.json',sources)
shutil.copyfile(adapter/'audit-fields.py',meshroot/'audit-fields.py')
assert sha(meshroot/'audit-fields.py')==adapter_prep['storage_driver_sha256']
dump(meshroot/'preparation.json',dict(prepared=True,elements=34350,nodes=46375,points=927450,
 all440_wall_nodes_and832_plates_exact=True,same_native_input_sha256=read(p/'fixed-wall-finer-input-v1/preparation.json')['same_native_input_sha256'],
 configured_primary_budget=200000,all_other_controls_retained=True,no_history_interpolation=True,p6_complete=False))
env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',
 OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',RAYON_NUM_THREADS='1')
for key in ('PYTHONPATH','PYTHONHOME'):env.pop(key,None)
records=[]
def call(label,command,cooperative=False):
 active();assert all(sha(Path(q))==h for q,h in binding.items())
 assert all(sha(runtime/n)==h for n,h in sources.items())
 start=time.monotonic();requested=False
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(list(map(str,command)),cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
   creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(command[0]),arguments=list(map(str,command[1:])),
   created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  while child.poll() is None:
   if cooperative and (p/'pause-requested').exists() and not requested:
    with (out/'cooperative-stop.log').open('x',encoding='utf-8') as stoplog:
     code=subprocess.call([str(solver),'stop-stages-job',str(job)],cwd=runtime,stdout=stoplog,stderr=subprocess.STDOUT,
      creationflags=subprocess.CREATE_NO_WINDOW)
    assert code==0;requested=True
   time.sleep(.5)
  code=child.returncode
 row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start,cooperative_stop_requested=requested)
 records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
 assert all(sha(Path(q))==h for q,h in binding.items())
 assert all(sha(runtime/n)==h for n,h in sources.items())
 if cooperative and code==2:
  assert read(job/'status.json')['state']=='checkpointed'
  dump(out/'terminal.json',dict(passed=False,checkpointed=True,unaccepted_trial_not_preserved_as_result=True,records=records,p6_complete=False))
  raise SystemExit(2)
 assert code==0,label
call('solve',[solver,'run-stages-job',job],True)
status=read(job/'status.json');assert status['state']=='succeeded' and status['artifact']=='result.json'
assert sha(job/'result.json')==status['artifact_sha256']
shutil.copyfile(job/'result.json',analysis/'stage-results.json')
assert sha(analysis/'stage-results.json')==status['artifact_sha256']
call('export',[sys.executable,'-X','utf8',runtime/'tools/export_stages.py',analysis/'ground/mesh.json',analysis/'stage-results.json',analysis/'results'])
call('verify',[sys.executable,'-X','utf8',runtime/'tools/verify_stage_package.py',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stage-results.json',analysis/'results','--solver',solver])
dump(meshroot/'terminal.json',dict(passed=True,scope='completed solve and independent full mechanical/rawVTU package; seed/load diagnostic pending',records=records,p6_complete=False))
call('independent-fields',[sys.executable,'-I','-X','utf8',meshroot/'audit-fields.py'])
fields=read(meshroot/'field-verification.json')
assert fields['seed_load_and_recorded_cutback_gates_passed'] and fields['all_seed_states_independently_evaluated']
assert fields['seed_state_count']==927450 and fields['rejected_trials_restored']>=3
assert len(fields['frames'])==5 and fields['floor_prescribed_pressure_kpa']==0.
assert sha(meshroot/'audit-fields.py')==adapter_prep['storage_driver_sha256']
assert all(sha(Path(q))==h for q,h in fields['sha256'].items())
final=dict(passed=True,elements=34350,integration_points=927450,records=records,
 same_frozen_F2_solver_as_original14909=True,all_original_numerical_accuracy_controls_retained=True,
 all440_wall_nodes_and832_oriented_plates_fixed=True,configured_primary_budget_only_changed_to200000=True,
 all5_stages_and_full_independent_mechanical_raw_VTU_seed_history_load_and_surface_checks_passed=True,
 storage_backed_diagnostic_validated_against_all_original8061_values=True,
 periodic_checkpoint_interval=4,spatial_comparison_pending=True,
 native_large_GUI_general100k_and_allP6_separate=True,p6_complete=False)
dump(meshroot/'study-terminal.json',final);dump(out/'terminal.json',final);dump(study/'terminal.json',final)
print('PASS34350 full study; unchanged original spatial comparison criteria still pending',flush=True)
