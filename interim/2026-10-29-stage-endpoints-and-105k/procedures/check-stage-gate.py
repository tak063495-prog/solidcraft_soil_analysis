"""QA the new endpoint API/CLI without altering the production workflow."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime';quality=p/'quality-v2';quality.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
with (p/'v2-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
env=dict(os.environ,CARGO_BUILD_JOBS='2',CARGO_TARGET_DIR=str(runtime/'target'),RUST_TEST_THREADS='2',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
records=[];sources=None
def call(label,args,folder=quality,expected=0):
    if (p/'pause-requested').exists():raise SystemExit(2)
    if sources:assert manifest()==sources
    command=list(map(str,args));started=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (folder/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    row=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-started);records.append(row)
    with (folder/(label+'-terminal.json')).open('x') as f:json.dump(row,f,indent=2)
    print(json.dumps(row),flush=True);assert code==expected
    if sources:assert manifest()==sources
cargo='C:/Users/link_/.cargo/bin/cargo.exe';call('format',[cargo,'fmt','--all'])
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import freeze,manifest
freeze(quality);sources=manifest();base=read(p/'source-snapshot.json')
assert len(sources)==611 and {n for n,h in sources.items() if h!=base[n]}=={'src/replay.rs','src/stages.rs','src/main.rs','tests/replay.rs'}
call('clippy',[cargo,'clippy','--locked','--all-targets','--','-D','warnings'])
call('focused-replay-tests',[cargo,'test','--locked','--release','--test','replay'])
assert '6 passed; 0 failed; 0 ignored;' in (quality/'focused-replay-tests.log').read_text('utf-8')
call('rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,m)) for m in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;', (quality/'rust-tests.log').read_text('utf-8'))]
assert sum(n[0] for n in counts)==277 and all(n[1]==n[2]==0 for n in counts)
call('build',[cargo,'build','--locked','--release'])
solver=p/'solver/soil-fem-3d.exe';solver.parent.mkdir();shutil.copy2(runtime/'target/release/soil-fem-3d.exe',solver);solver_sha=sha(solver)
call('python-tests',[sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
assert 'Ran 171 tests' in (quality/'python-tests.log').read_text('utf-8')
with (quality/'terminal.json').open('x') as f:json.dump(dict(passed=True,rust_tests=277,python_tests=171,source_files=611,only4_private_rust_files_changed=True,solver_sha256=solver_sha,production_not_modified=True,actual_cli_controls_pending=True,records=records,p6_complete=False),f,indent=2)
out=p/'cli-controls';out.mkdir();mesh=read(runtime/'examples/mc-column-mesh.json')
from collections import Counter
faces=Counter(tuple(sorted(t[j] for j in face)) for t in mesh['tetrahedra'] for face in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)])
mesh.setdefault('face_sets',{})['top']=[list(f) for f,n in faces.items() if n==1 and all(mesh['nodes'][i][2]==1. for i in f)]
assert mesh['face_sets']['top']
ground=read(runtime/'examples/ground-geostatic-mc-column.json');case=read(runtime/'examples/stages-geostatic-column.json')
spec=dict(schema='soil-fem-replay/1',material_ids=[1,2],factors=[1.,1.1,1.2],psi_policy='keep',tension_policy='keep')
for n,v in [('mesh.json',mesh),('ground.json',ground),('case.json',case),('spec.json',spec)]:
    with (out/n).open('x') as f:json.dump(v,f,allow_nan=False)
paths=[out/n for n in ['mesh.json','ground.json','case.json','spec.json']]
call('full-cli',[solver,'replay-stages',*paths,out/'full.json'],out)
full=read(out/'full.json');assert all(t['outcome']=='completed' for t in full['trials'])
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
for count in [0,1,len(case['stages'])]:
    dest=out/('gate'+str(count)+'.json');call('gate'+str(count),[solver,'replay-through-stage',*paths,str(count),dest],out)
    gate=read(dest);assert gate['requested_through_construction_stages']==count and gate['stability_classification']=='not_evaluated'
    assert len(gate['original_case']['stages'])==len(case['stages']) and len(gate['replay']['case']['stages'])==count
    for a,b in zip(gate['replay']['trials'],full['trials'],strict=True):
        assert a['factor']==b['factor'] and a['outcome']=='completed'
        assert stable(a['run']['frames'])==stable(b['run']['frames'][:count+1])
        assert stable(a['run']['last_committed'])==stable(b['run']['frames'][count]['state'])
for label,count in [('negative','-1'),('past-end',str(len(case['stages'])+1)),('fraction','1.5')]:
    dest=out/(label+'.json');call(label,[solver,'replay-through-stage',*paths,count,dest],out,expected=1);assert not dest.exists()
sentinel=out/'existing.json';sentinel.write_bytes(b'owned existing artifact');call('no-overwrite',[solver,'replay-through-stage',*paths,'0',sentinel],out,expected=1);assert sentinel.read_bytes()==b'owned existing artifact'
assert manifest()==sources and sha(solver)==solver_sha
with (p/'terminal.json').open('x') as f:json.dump(dict(passed=True,private_candidate=True,rust_tests=277,python_tests=171,legacy_replay_unchanged=True,
    actual_cli_zero_first_last_stage_gates_match_full_factorwise_histories=True,bad_counts_create_no_report=True,existing_output_preserved=True,
    factor_independence_and_partial_failure_restoration_tested=True,physical_stability_or_fs_not_inferred=True,
    strength_reduction_starts_from_original_initialization=True,reduction_only_after_an_unreduced_origin_not_implemented=True,
    solver_sha256=solver_sha,records=records,production_not_modified=True,p6_complete=False),f,indent=2)
print('PASS private stage-gate277/171, actual CLI prefix histories and failure/output guards; no Fs classification',flush=True)
