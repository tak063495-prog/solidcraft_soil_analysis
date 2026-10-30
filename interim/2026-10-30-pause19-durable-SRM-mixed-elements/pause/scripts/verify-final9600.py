"""Deep-verify the final cooperative checkpoint by preparing, never running, a fork."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,time
import ijson
p=Path(__file__).resolve().parent;root=p.parents[1];pause=root/'runs/goal-progress-20261027-v1/pause-remaining20-20261030-183721'
out=pause/'final9600-resume-verification-v1';out.mkdir();job=p/'plane-strain-slope-comparison-v1/xy2-ny2/solver-job'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
assert (p/'pause-requested').exists()
status=read(job/'status.json');manifest=read(job/'job.json');solver=job/'solver/soil-fem-3d.exe'
assert status['state']=='checkpointed' and status['artifact']=='checkpoint.json'
assert manifest['job_id']=='9c8d90dd603569e430762ee2654f4d38b536f6072bdb5e4d0ce491df62ab006e'
assert sha(solver)==manifest['binding']['identity']['solver_binary_sha256']=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
assert sha(job/'checkpoint.json')==status['artifact_sha256']=='b5660a0cb8bcd20093cf289a99f8bf32568b26223aa16e06cc3ed49191955488'
before={n:sha(job/n) for n in [*manifest['binding']['files_sha256'],'job.json','status.json','checkpoint.json']}
for name,h in manifest['binding']['files_sha256'].items():assert before[name]==h
with (job/'checkpoint.json').open('rb') as f:progress=float(next(ijson.items(f,'run.geostatic.progress')))
assert progress==0.607421875
dump(out/'source-binding.json',dict(source_job=str(job),files_sha256=before,accepted_geostatic_progress=progress))
resume=out/'prepared-resume-job';command=[str(solver),'prepare-stages-resume',str(job),str(resume)];start=time.monotonic()
with (out/'prepare.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out/'prepare-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=str(solver),arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code=child.wait()
assert code==0,'full typed checkpoint/history verification failed; inspect prepare.log'
copied=read(resume/'job.json');saved=read(resume/'status.json')
assert copied['binding']['parent_job_id']==manifest['job_id'] and copied['binding']['identity']==manifest['binding']['identity']
assert copied['binding']['files_sha256']['inputs/checkpoint.json']==before['checkpoint.json']==sha(resume/'inputs/checkpoint.json')
assert sha(resume/'solver/soil-fem-3d.exe')==sha(solver)
assert saved['state']=='prepared' and not (resume/'run.claim').exists()
assert all(sha(job/n)==h for n,h in before.items())
result=dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_job=str(job),source_job_id=manifest['job_id'],
    accepted_geostatic_progress=progress,checkpoint_sha256=before['checkpoint.json'],
    solver_sha256=sha(solver),prepared_resume_job=str(resume),resume_job_id=copied['job_id'],resume_manifest_sha256=sha(resume/'job.json'),
    all_complete_checkpoint_bytes_exact_and_Rust_context_material_history_verified=True,
    no_original_file_changed=True,no_Newton_or_resume_run_started=True,unaccepted_trial_not_claimed_saved=True,
    independent_final_equilibrium_not_claimed=True,elapsed_seconds=time.monotonic()-start,p6_complete=False)
dump(out/'terminal.json',result);print(json.dumps(result),flush=True)
