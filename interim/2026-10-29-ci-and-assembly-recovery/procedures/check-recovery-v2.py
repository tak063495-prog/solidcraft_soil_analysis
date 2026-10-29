"""Recompute remaining stress release using only the accepted candidate's own snapshot."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'recovery';out.mkdir(exist_ok=True);assert not any(out.iterdir())
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
with (p/'recovery-v2-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
while not (p/'terminal.json').exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    for q in (p/'quality-v1').glob('*-terminal.json'):
        if read(q).get('exit_code',0)!=0:raise RuntimeError('Candidate quality failed before recovery: '+q.name)
    time.sleep(5)
assert read(p/'terminal.json')['passed'] and read(p/'quality-v1/terminal.json')['rust_tests']==278
source=p/'workflow/analysis/solver-job';binary=p/'solver/soil-fem-3d.exe';assert sha(binary)==read(p/'terminal.json')['solver_sha256']
before={q.relative_to(source).as_posix():sha(q) for q in source.rglob('*') if q.is_file()}
pointer=read(source/'recovery.json');selected=[]
for q in source.glob('recovery-checkpoint-*.json'):
    data=read(q);run=data['run'];c=run.get('continuation')
    if c and c.get('stage_index')==3 and 0<c.get('progress',0)<1:
        selected.append((q,data))
assert len(selected)==1,[(q.name) for q,_ in selected]
slot,data=selected[0];run=data['run'];c=run['continuation']
assert len(run['frames'])==4 and c['progress']==.875 and sum(run['last_committed']['active_elements'])==250
immutable=out/'immutable-source'
def ignore(folder,names):return [n for n in names if n in ['result.json','status.json','progress.jsonl','run.claim','cancel.request','cancel.json']]
shutil.copytree(source,immutable,ignore=ignore)
generation=int(slot.stem.rsplit('-',1)[1]);assert generation in [0,1]
# The latest pointer's generation is one newer than the other alternating slot.
generation=pointer['generation']-1;assert generation%2==int(slot.stem.rsplit('-',1)[1])
chosen={**pointer,'generation':generation,'checkpoint':slot.name,'checkpoint_sha256':sha(slot),'run_sha256':data['run_sha256']}
(immutable/'recovery.json').write_text(json.dumps(chosen,indent=2)+'\n','utf-8')
immutable_binding={q.relative_to(immutable).as_posix():sha(q) for q in immutable.rglob('*') if q.is_file()}
records=[];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
def call(label,args):
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));started=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (out/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-started);records.append(row)
    with (out/(label+'-terminal.json')).open('x') as f:json.dump(row,f,indent=2)
    print(json.dumps(row),flush=True);assert code==0
analysis=out/'workflow/analysis';(analysis/'ground').mkdir(parents=True)
for src,dst in [('inputs/mesh.json','ground/mesh.json'),('inputs/ground.json','ground/ground.json'),('inputs/case.json','stages.json')]:shutil.copy2(source/src,analysis/dst)
job=analysis/'solver-job';call('prepare',[binary,'prepare-stages-recovery',immutable,job]);call('solve',[job/'solver/soil-fem-3d.exe','run-stages-job',job])
shutil.copy2(job/'result.json',analysis/'stage-results.json')
call('export',[sys.executable,p/'runtime/tools/export_stages.py',analysis/'ground/mesh.json',analysis/'stage-results.json',analysis/'results'])
call('verify',[sys.executable,p/'runtime/tools/verify_stage_package.py',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stage-results.json',analysis/'results','--solver',binary])
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
expected=read(source/'result.json');actual=read(job/'result.json')
for key in ['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','transition_audits','geostatic','attempts']:assert stable(actual[key])==stable(expected[key]),key
assert before=={q.relative_to(source).as_posix():sha(q) for q in source.rglob('*') if q.is_file()}
assert immutable_binding=={q.relative_to(immutable).as_posix():sha(q) for q in immutable.rglob('*') if q.is_file()}
with (out/'terminal.json').open('x') as f:json.dump(dict(passed=True,private_scatter_candidate=True,source_and_immutable_copy_every_file_unchanged=True,selected_generation=generation,excavation_progress=c['progress'],remaining_excavation_stress_release_recomputed=True,all8_physical_groups_exact_to_continuous=True,view_and_mechanical_audits_passed=True,solver_sha256=sha(binary),records=records,production_not_modified=True,p6_complete=False),f,indent=2)
print('PASS actual remaining excavation after cached scatter snapshot; all8 exact and both parents unchanged',flush=True)
