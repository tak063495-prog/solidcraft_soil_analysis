"""Full QA first, then unchanged-input field/performance controls after a measured gate."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime'
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert read(p/'preflight/terminal.json')['passed']
with (p/'validation-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import freeze,manifest
quality=p/'full-quality';quality.mkdir();freeze(quality);sources=manifest()
assert sources==read(p/'preflight/source-snapshot.json') and len(sources)==611
env=dict(os.environ,CARGO_BUILD_JOBS='2',CARGO_TARGET_DIR=str(runtime/'target'),RUST_TEST_THREADS='2',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
records=[]
def call(label,args,folder=quality):
    assert manifest()==sources
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));began=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (folder/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-began,sources_unchanged=manifest()==sources);records.append(row)
    with (folder/(label+'-terminal.json')).open('x') as f:json.dump(row,f,indent=2)
    print(json.dumps(row),flush=True);assert code==0 and row['sources_unchanged']
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
call('format',[cargo,'fmt','--all','--check'])
call('rust-tests',[cargo,'test','--locked','--release'])
log=(quality/'rust-tests.log').read_text('utf-8')
counts=[tuple(map(int,m)) for m in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',log)]
assert sum(n[0] for n in counts)==274 and all(n[1]==n[2]==0 for n in counts)
call('build',[cargo,'build','--locked','--release'])
solver=runtime/'target/release/soil-fem-3d.exe';solver_sha=sha(solver)
call('python-tests',[sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
assert 'Ran 171 tests' in (quality/'python-tests.log').read_text('utf-8')
with (quality/'terminal.json').open('x') as f:json.dump(dict(passed=True,rust_tests=274,python_tests=171,source_files=611,only_private_stages_rs_changed=True,
    production_not_promoted=True,solver_sha256=solver_sha,records=records,fields_and_performance_acceptance_pending=True,p6_complete=False),f,indent=2)
print('Full274 Rust/171 Python passed; waiting for the measured material-kernel gate',flush=True)
while not (p/'measurement-gate.json').exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    if (p/'measurement-declined.json').exists():raise RuntimeError('Material optimization gate declined; preserve QA without starting performance controls')
    time.sleep(5)
gate=read(p/'measurement-gate.json');assert gate['passed'] and gate['profile_full_field_and_equivalence_passed']
assert sha(solver)==solver_sha and manifest()==sources
controls=p/'controls';controls.mkdir()
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
groups=['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','geostatic','attempts']
field_reports=[]
for n in [6,12]:
    baseline=root/('runs/mcc-hydrostatic-scaling-20261029-v1/n6' if n==6 else 'runs/mcc-hydrostatic-durable-restart-20261029-v1/n12')
    out=p/('n'+str(n));out.mkdir()
    for name in ['mesh.json','ground.json','stages.json','preparation.json']:shutil.copy2(baseline/name,out/name);assert sha(baseline/name)==sha(out/name)
    call('prepare',[solver,'prepare-durable-stages-job',out/'mesh.json',out/'ground.json',out/'stages.json','4',out/'solver-job'],out)
    call('solve',[out/'solver-job/solver/soil-fem-3d.exe','run-stages-job',out/'solver-job'],out)
    verifier=root/'runs/mcc-hydrostatic-durable-restart-20261029-v1/verify.py'
    # The unchanged verifier locates n6/n12 beside itself.
    if not (p/'verify.py').exists():shutil.copy2(verifier,p/'verify.py')
    assert sha(p/'verify.py')==sha(verifier)
    call('theory',[sys.executable,'-X','utf8',p/'verify.py',str(n)],out)
    expected=read(baseline/'solver-job/result.json');actual=read(out/'solver-job/result.json')
    for key in groups:assert stable(expected[key])==stable(actual[key]),(n,key)
    field_reports.append(dict(subdivisions=n,elements=6*n**3,all7_physical_groups_exact=True,full_field_theory_and_rollback_passed=True,original3_inputs_byte_exact=True))
    with (out/'terminal.json').open('x') as f:json.dump(dict(passed=True,verification=field_reports[-1],solver_sha256=solver_sha,concurrent_other_jobs=True,isolated_throughput_benchmark=False,p6_complete=False),f,indent=2)
capture=root/'runs/native-fem-smoke-20261029-v1/Native63 cooperative smoke/capture'
call('native274',[sys.executable,runtime/'tools/native_job.py','run',capture,controls/'native274','--solver',solver,'--cooperative','--checkpoint-interval','4','--mesh-size-m','4'],controls)
expected=read(root/'runs/native-recovery-worker-development-20261029-v1/controls/continuous/workflow/analysis/stage-results.json')
actual=read(controls/'native274/workflow/analysis/stage-results.json')
for key in [*groups,'transition_audits']:assert stable(expected[key])==stable(actual[key]),key
assert read(controls/'native274/status.json')['state']=='succeeded' and read(controls/'native274/status.json')['solver_report']['state']=='succeeded'
assert manifest()==sources and sha(solver)==solver_sha
with (p/'acceptance.json').open('x') as f:json.dump(dict(passed=True,rust_tests=274,python_tests=171,only_private_stages_rs_changed=True,root_not_promoted=True,
    hydrostatic1296_and10368_full_fields_theory_and_all7_groups_exact=True,actual_native274_all8_groups_and_view_audits_exact=True,
    material_only_parallelism=True,serial_global_assembly_and_tolerances_unchanged=True,field_reports=field_reports,solver_sha256=solver_sha,
    records=records,general_mixed_backend_performance_not_claimed=True,isolated_throughput_benchmark=False,p6_complete=False),f,indent=2)
print('PASS private material-parallel full QA, unchanged1296/10368 fields and real native274 all8 groups; production not promoted',flush=True)
