"""Require QA and actual excavation recovery before a larger unchanged mixed case."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'large2865';out.mkdir(exist_ok=True);assert not any(out.iterdir())
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    h=hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
with (p/'large-v2-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
print('Waiting for full278/171, real274 all8 and retained excavation stress-release recovery before large2865',flush=True)
while not (p/'recovery/terminal.json').exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    for folder in [p/'quality-v1',p/'recovery']:
        for q in folder.glob('*-terminal.json'):
            if read(q).get('exit_code',0)!=0:raise RuntimeError('Candidate gate failed before large-case execution: '+q.name)
    time.sleep(5)
assert read(p/'recovery/terminal.json')['passed'] and read(p/'terminal.json')['passed'] and read(p/'quality-v1/terminal.json')['rust_tests']==278
candidate_sources=read(p/'quality-v1/source-snapshot.json');assert all(sha(p/'runtime'/n)==h for n,h in candidate_sources.items())
binary=p/'solver/soil-fem-3d.exe';assert sha(binary)==read(p/'terminal.json')['solver_sha256']
baseline=root/'runs/mcc-fixed-wall-durable-pair-20261029-v1/mesh2865';assert read(baseline/'study-terminal.json')['passed']
assert read(root/'runs/mcc-fixed-wall-durable-pair-20261029-v1/comparison.json')['wall_spatial_topology_and_sections_and_loads_preserved']
analysis=out/'workflow/analysis';(analysis/'ground').mkdir(parents=True)
paths=['ground/mesh.json','ground/ground.json','stages.json']
for n in paths:shutil.copy2(baseline/'workflow/analysis'/n,analysis/n);assert sha(analysis/n)==sha(baseline/'workflow/analysis'/n)
sources=read(baseline/'audit-source-manifest.json');assert len(sources)==600
for n,h in sources.items():
    assert sha(baseline/'runtime'/n)==h
    q=out/'runtime'/n;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(baseline/'runtime'/n,q);assert sha(q)==h
shutil.copy2(baseline/'audit-fields.py',out/'audit-fields.py');shutil.copy2(baseline/'audit-source-manifest.json',out/'audit-source-manifest.json')
binding={n:sha(analysis/n) for n in paths};records=[]
with (out/'preparation.json').open('x') as f:json.dump(dict(prepared=True,elements=2865,points=77355,plates=832,structural_nodes=440,original3_inputs_byte_exact=True,input_sha256=binding,solver_sha256=sha(binary),candidate_source_manifest_sha256=sha(p/'quality-v1/source-snapshot.json'),independent_auditor_source_manifest_sha256=sha(out/'audit-source-manifest.json'),p6_complete=False),f,indent=2)

env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
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
job=analysis/'solver-job';call('prepare',[binary,'prepare-durable-stages-job',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stages.json','4',job]);call('solve',[job/'solver/soil-fem-3d.exe','run-stages-job',job]);shutil.copy2(job/'result.json',analysis/'stage-results.json')
call('export',[sys.executable,p/'runtime/tools/export_stages.py',analysis/'ground/mesh.json',analysis/'stage-results.json',analysis/'results'])
call('verify',[sys.executable,p/'runtime/tools/verify_stage_package.py',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stage-results.json',analysis/'results','--solver',binary])
with (out/'terminal.json').open('x') as f:json.dump(dict(passed=True,scope='solver and mechanical package only; independent seed/load and all8 comparison pending',elements=2865,records=list(records),p6_complete=False),f,indent=2)
# Independent field auditor is the exact frozen600-source reader from the old
# reference study, separate from the new278/171-tested numerical candidate.
call('independent-fields',[sys.executable,out/'audit-fields.py']);assert read(out/'field-verification.json')['seed_load_and_recorded_cutback_gates_passed']
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
expected=read(baseline/'workflow/analysis/stage-results.json');actual=read(analysis/'stage-results.json')
for key in ['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','transition_audits','geostatic','attempts']:assert stable(actual[key])==stable(expected[key]),key
assert all(sha(analysis/n)==h for n,h in binding.items()) and all(sha(p/'runtime'/n)==h for n,h in candidate_sources.items()) and all(sha(out/'runtime'/n)==h for n,h in sources.items())
with (out/'study-terminal.json').open('x') as f:json.dump(dict(passed=True,private_scatter_candidate=True,elements=2865,integration_points=77355,plates=832,structural_nodes=440,original3_inputs_byte_exact=True,all5_stages_and_all8_physical_groups_exact=True,view_mechanical_and_independent_field_audits_passed=True,constitutive_equations_and_tolerances_unchanged=True,solver_sha256=sha(binary),records=records,input_sha256=binding,auditor_source_files=600,candidate_source_files=611,production_not_promoted=True,concurrent_other_jobs=True,isolated_throughput_benchmark=False,spatial_convergence_accepted=False,p6_complete=False),f,indent=2)
print('PASS larger2865 cached mixed assembly; all8 histories exact and independent frozen audit passed',flush=True)
