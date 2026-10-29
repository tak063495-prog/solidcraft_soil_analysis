"""Require unchanged numerical fields after bounded pattern caching."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime';quality=p/'quality-v1';quality.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
with (p/'driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
env=dict(os.environ,CARGO_BUILD_JOBS='2',CARGO_TARGET_DIR=str(runtime/'target'),RUST_TEST_THREADS='2',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
records=[];sources=None
def call(label,args,folder=quality):
    if (p/'pause-requested').exists():raise SystemExit(2)
    if sources:assert manifest()==sources
    command=list(map(str,args));began=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (folder/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-began);records.append(row)
    with (folder/(label+'-terminal.json')).open('x') as f:json.dump(row,f,indent=2)
    print(json.dumps(row),flush=True);assert code==0
    if sources:assert manifest()==sources
cargo='C:/Users/link_/.cargo/bin/cargo.exe';call('format',[cargo,'fmt','--all'])
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest,freeze
freeze(quality);sources=manifest();original=read(p/'source-snapshot.json')
assert len(sources)==611 and {n for n,h in sources.items() if original[n]!=h}=={'src/linalg.rs','src/stages_mixed.rs'}
call('clippy',[cargo,'clippy','--locked','--all-targets','--','-D','warnings'])
call('focused-tests',[cargo,'test','--locked','--release','compiled_'])
focused=(quality/'focused-tests.log').read_text('utf-8');assert '4 passed; 0 failed; 0 ignored;' in focused
call('rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,m)) for m in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;', (quality/'rust-tests.log').read_text('utf-8'))]
assert sum(n[0] for n in counts)==278 and all(n[1]==n[2]==0 for n in counts)
call('build',[cargo,'build','--locked','--release'])
solver=p/'solver/soil-fem-3d.exe';solver.parent.mkdir();shutil.copy2(runtime/'target/release/soil-fem-3d.exe',solver);solver_sha=sha(solver)
call('python-tests',[sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
assert 'Ran 171 tests' in (quality/'python-tests.log').read_text('utf-8')
with (quality/'terminal.json').open('x') as f:json.dump(dict(passed=True,rust_tests=278,python_tests=171,source_files=611,solver_sha256=solver_sha,production_not_modified=True,all8_field_acceptance_pending=True,records=records,p6_complete=False),f,indent=2)
baseline=root/'runs/mcc-parallel-material-development-20261029-v1/controls/native274/workflow/analysis'
analysis=p/'workflow/analysis';(analysis/'ground').mkdir(parents=True)
for n in ['ground/mesh.json','ground/ground.json','stages.json']:shutil.copy2(baseline/n,analysis/n);assert sha(baseline/n)==sha(analysis/n)
job=analysis/'solver-job'
call('prepare',[solver,'prepare-durable-stages-job',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stages.json','4',job])
call('solve',[job/'solver/soil-fem-3d.exe','run-stages-job',job])
shutil.copy2(job/'result.json',analysis/'stage-results.json')
call('export',[sys.executable,runtime/'tools/export_stages.py',analysis/'ground/mesh.json',analysis/'stage-results.json',analysis/'results'])
call('verify',[sys.executable,runtime/'tools/verify_stage_package.py',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stage-results.json',analysis/'results','--solver',solver])
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
expected=read(baseline/'stage-results.json');actual=read(analysis/'stage-results.json')
for key in ['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','transition_audits','geostatic','attempts']:assert stable(expected[key])==stable(actual[key]),key
assert manifest()==sources and sha(solver)==solver_sha
with (p/'terminal.json').open('x') as f:json.dump(dict(passed=True,private_candidate_only=True,rust_tests=278,python_tests=171,all5_stages_and_all8_groups_exact=True,view_and_mechanical_audits_passed=True,scalar_floating_arithmetic_and_addition_order_unchanged=True,cache_limit_bytes=64<<20,dense_and_unprojected_and_overbudget_fallback=True,records=records,solver_sha256=solver_sha,production_not_modified=True,large_mixed_field_performance_and_recovery_pending=True,p6_complete=False),f,indent=2)
print('PASS private scatter candidate full278/171, actual274 all8 exact and view/mechanics audits',flush=True)
