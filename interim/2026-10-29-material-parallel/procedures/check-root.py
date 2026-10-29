"""Actual production build, whole-field theory and exact old controls after promotion."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert read(p/'promotion.json')['promoted']
with (p/'root-build-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
out=p/'root-control';out.mkdir();records=[]
env=dict(os.environ,CARGO_BUILD_JOBS='2',CARGO_TARGET_DIR=str(root/'target'),PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
def call(label,args,folder=out):
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));began=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (folder/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    records.append(dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-began));print(json.dumps(records[-1]),flush=True);assert code==0
cargo='C:/Users/link_/.cargo/bin/cargo.exe';call('format',[cargo,'fmt','--all','--check']);call('build',[cargo,'build','--locked','--release'])
solver=out/'solver/soil-fem-3d.exe';solver.parent.mkdir();shutil.copy2(root/'target/release/soil-fem-3d.exe',solver);solver_sha=sha(solver)
case=out/'n6';case.mkdir()
for name in ['mesh.json','ground.json','stages.json','preparation.json']:shutil.copy2(p/'n6'/name,case/name);assert sha(case/name)==sha(p/'n6'/name)
shutil.copy2(p/'verify.py',out/'verify.py')
job=case/'solver-job';call('prepare',[solver,'prepare-durable-stages-job',case/'mesh.json',case/'ground.json',case/'stages.json','4',job]);call('solve',[job/'solver/soil-fem-3d.exe','run-stages-job',job]);call('theory',[sys.executable,'-X','utf8',out/'verify.py','6'])
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
expected=read(p/'n6/solver-job/result.json');actual=read(job/'result.json')
for key in ['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','geostatic','attempts']:assert stable(actual[key])==stable(expected[key]),key
tested=read(p/'full-quality/source-snapshot.json');assert all(sha(root/n)==h for n,h in tested.items() if not n.endswith('.md'))
assert sha(solver)==sha(root/'target/release/soil-fem-3d.exe')==solver_sha
with (p/'actual-root-acceptance.json').open('x') as f:json.dump(dict(passed=True,actual_promoted_root_build=True,rust_tests=274,python_tests=171,
    all_production_non_markdown_matches_tested_candidate=True,actual_root1296_theory_and_all7_groups_exact=True,
    previously_validated10368_theory_loaded_recovery_and_native274_all8_source_exact=True,solver_sha256=solver_sha,
    private_accepted_solver_sha256=read(p/'acceptance.json')['solver_sha256'],maximum_material_workers=4,records=records,
    old_jobs_require_their_original_frozen_solver=True,mixed_backend_performance_acceptance_separate=True,p6_complete=False),f,indent=2)
print('PASS actual Root build and unchanged1296 fields/theory; accepted274/171 source exact',flush=True)
