"""Select a byte-exact retained loaded snapshot and actually recompute unloading at10k scale."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
read=lambda q:json.loads(q.read_text(encoding='utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
with (p/'loaded-recovery-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
while not (p/'n12/terminal.json').exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    time.sleep(5)
assert read(p/'n12/terminal.json')['passed']
parent=p/'n12/solver-job';assert read(p/'n12/analytical-verification.json')['passed']
latest=read(parent/'recovery.json');owner=read(parent/'job.json');selected=[]
for q in sorted(parent.glob('recovery-checkpoint-*.json')):
    r=read(q);c=r['run'].get('continuation') or {}
    if r['run']['last_committed']['stage_index']==1 and c.get('stage_index')==0 and c.get('progress')==1.:selected.append((q,r))
assert len(selected)==1
slot,saved=selected[0];generation=latest['generation']-(0 if latest['checkpoint']==slot.name else 1)
assert generation>0 and slot.name==f'recovery-checkpoint-{generation%2}.json'
out=p/'loaded-recovery';out.mkdir();source=out/'accepted-snapshot-source';source.mkdir()
allowed={'inputs/mesh.json','inputs/ground.json','inputs/case.json','solver/soil-fem-3d.exe'};assert set(owner['binding']['files_sha256'])==allowed
for n,h in owner['binding']['files_sha256'].items():
    target=source/n;target.parent.mkdir(parents=True,exist_ok=True);assert sha(parent/n)==h;shutil.copy2(parent/n,target);assert sha(target)==h
shutil.copy2(parent/'job.json',source/'job.json');shutil.copy2(slot,source/slot.name);assert sha(source/slot.name)==sha(slot)
pointer=dict(schema='soil-fem-stages-recovery/1',job_id=owner['job_id'],generation=generation,checkpoint=slot.name,checkpoint_sha256=sha(slot),run_sha256=saved['run_sha256'])
with (source/'recovery.json').open('x') as f:json.dump(pointer,f,indent=2)
case=out/'n12';case.mkdir()
for n in ['mesh.json','ground.json','stages.json','preparation.json']:shutil.copy2(p/'n12'/n,case/n)
shutil.copy2(p/'verify.py',out/'verify.py');solver=p/'runtime/target/release/soil-fem-3d.exe';assert sha(solver)==owner['binding']['files_sha256']['solver/soil-fem-3d.exe']
records=[]
def call(label,args):
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));began=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,env={**os.environ,'PYTHONUTF8':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OMP_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (out/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-began);records.append(row);print(json.dumps(row),flush=True);assert code==0
job=case/'solver-job';call('prepare-loaded-recovery',[solver,'prepare-stages-recovery',source,job])
assert sha(job/'inputs/checkpoint.json')==sha(slot)
with (out/'prepared.json').open('x') as f:json.dump(dict(prepared=True,retained_original_slot_byte_exact=True,no_numeric_checkpoint_content_modified=True,original_parent_and_original_pointer_untouched=True,
    selected_generation=generation,selected_slot=slot.name,source_checkpoint_sha256=sha(slot),snapshot_stage_index=1,snapshot_continuation_stage=0,snapshot_continuation_progress=1.,
    full_unloading_recomputation_pending=True,records=records,p6_complete=False),f,indent=2)
call('solve-loaded-recovery',[job/'solver/soil-fem-3d.exe','run-stages-job',job]);call('analytical-loaded-recovery',[sys.executable,out/'verify.py','12'])
assert read(case/'analytical-verification.json')['passed']
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
continuous=read(parent/'result.json');resumed=read(job/'result.json');assert stable(continuous)==stable(resumed),[k for k in continuous if stable(continuous[k])!=stable(resumed.get(k))]
assert sha(slot)==pointer['checkpoint_sha256'] and read(parent/'recovery.json')==latest
with (out/'terminal.json').open('x') as f:json.dump(dict(passed=True,elements=10368,actual_unloading_recomputed_from_loaded_plastic_snapshot=True,parallel_material_solver=True,solver_sha256=sha(solver),source_snapshot_byte_exact=True,
    both_analytic_full_field_audits_passed=True,all_result_fields_except_elapsed_seconds_exact=True,original_parent_not_modified=True,selected_generation=generation,records=records,
    full_excavation_performance_acceptance=False,p6_complete=False),f,indent=2)
print('PASS10k: retained loaded plastic snapshot, recomputed unloading, full theory and exact whole result',flush=True)
