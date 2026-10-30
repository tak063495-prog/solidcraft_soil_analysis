"""Separate loading-increment and material-integration errors at fixed physics."""
from pathlib import Path
import ast,datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent
root=p.parents[1]
study=p/'plane-strain-slope-comparison-v1'
base=study/'xy1-ny1'
out=study/'increment-sensitivity-preparation-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
if (p/'pause-requested').exists():raise SystemExit(2)
assert read(base/'solve-v1/terminal.json')['passed']
assert read(base/'field-comparison-v1/full-independent-audit-terminal.json')['passed']
assert read(base/'field-comparison-v1/assessment.json')['passed'] is False
assert read(base/'field-comparison-v1/assessment.json')['metrics']['displacement_sample_relative_L2_per_depth']['0.5']>.05
cases={'xy1-ny1-global64':dict(initial_increment=1/64,minimum_increment=1/8192),
       'xy1-ny1-local-half':dict(material_maximum_substep=5e-5)}
policy_path=root/'runs/goal-progress-20261027-v1/active-pause-policy.json'
policy=read(policy_path);guard=Path(policy['script']);before=guard.read_bytes()
assert hashlib.sha256(before).hexdigest()==policy['script_sha256']
assert policy['goal_active'] and not policy['explicit_user_pause']
assert policy['threshold_comparison']=='strictly_greater_than' and policy['remaining_percent_comparison']=='strictly_less_than'
text=before.decode('utf-8')
anchor="    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy1-ny1/recovery-probe-v1')\n"
paths=[f'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/{name}/solver-job' for name in cases]
lines=''.join(f"    jobs.append(root/'{name}')\n" for name in paths)
assert text.count(anchor)==1 and all(path not in text for path in paths)
updated=text.replace(anchor,anchor+lines,1);ast.parse(updated)
assert updated.replace(lines,'',1)==text
(out/'previous-guard.py').write_bytes(before);(out/'previous-policy.json').write_bytes(policy_path.read_bytes())
guard.write_text(updated,'utf-8',newline='\n')
policy.update(script_sha256=sha(guard),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
temporary=policy_path.with_suffix('.increment-sensitivity-pending.json')
with temporary.open('x',encoding='utf-8') as f:json.dump(policy,f,indent=2,allow_nan=False)
temporary.replace(policy_path)
dump(out/'pause-registration.json',dict(passed=True,paths=paths,only2_exact_new_job_paths_added=True,
    old_guard_sha256=hashlib.sha256(before).hexdigest(),new_guard_sha256=sha(guard),
    all_old_scope_and_preservers_and_user_threshold_unchanged=True))
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
binding={str(base/name):sha(base/name) for name in ['mesh.json','ground.json','stages.json','reference-map.json','prepared.json']}
binding[str(study/'protocol.json')]=sha(study/'protocol.json')
solver=base/'solver-job/solver/soil-fem-3d.exe'
assert sha(solver)=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
original=read(base/'stages.json')
rows=[]
for name,changes in cases.items():
    if (p/'pause-requested').exists():raise SystemExit(2)
    folder=study/name;folder.mkdir()
    for source in ['mesh.json','ground.json','reference-map.json']:
        shutil.copyfile(base/source,folder/source);assert sha(base/source)==sha(folder/source)
    case=json.loads(json.dumps(original))
    case['settings'].update(changes)
    actual_changes={k for k in case['settings'] if case['settings'][k]!=original['settings'][k]}
    assert actual_changes==set(changes)
    restored=json.loads(json.dumps(case));restored['settings']=original['settings']
    assert restored==original
    assert case['settings']['relative_force_tolerance']==1e-8 and case['settings']['relative_linear_tolerance']==1e-10
    assert case['settings']['tangent_perturbation']==original['settings']['tangent_perturbation']
    dump(folder/'stages.json',case)
    row=read(base/'prepared.json')
    row.update(case=name,mesh=str(folder/'mesh.json'),parent_case='xy1-ny1',only_changed_numerical_settings=changes,
        input_sha256={q:sha(folder/q) for q in ['mesh.json','ground.json','stages.json','reference-map.json']},
        loading_pseudo_increment_not_physical_time_or_consolidation=True)
    dump(folder/'prepared.json',row)
    job=folder/'solver-job'
    assert str(job.relative_to(root)).replace('\\','/') in paths and not job.exists()
    command=[str(solver),'prepare-durable-stages-job',*[str(folder/q) for q in ['mesh.json','ground.json','stages.json']],'4',str(job)]
    start=time.monotonic()
    with (out/(name+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(name+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    terminal=dict(case=name,exit_code=code,elapsed_seconds=time.monotonic()-start,changed_numerical_settings=changes,
        material_load_boundary_and_FSS_unchanged=True)
    dump(out/(name+'-terminal.json'),terminal);assert code==0
    manifest=read(job/'job.json')
    assert manifest['binding']['checkpoint_interval']==4
    assert all(sha(job/q)==h for q,h in manifest['binding']['files_sha256'].items())
    rows.append(terminal)
assert all(sha(Path(q))==h for q,h in binding.items())
dump(out/'terminal.json',dict(passed=True,rows=rows,case_names=list(cases),
    all_original_physical_model_and_predeclared_2d_comparison_limits_unchanged=True,
    mesh_and_reference_map_byte_exact=True,individual_numerical_controls_separated=True,
    fresh_full_gravity_replay_required_no_foreign_checkpoint_injection=True,
    all_exact_future_job_paths_registered_before_creation=True,
    no_new_solve_started=True,driver_sha256=sha(Path(__file__)),source_binding=binding,p6_complete=False))
print(json.dumps(dict(prepared=True,case_names=list(cases),same_physics_and_thresholds=True,solves_not_started=True)),flush=True)
