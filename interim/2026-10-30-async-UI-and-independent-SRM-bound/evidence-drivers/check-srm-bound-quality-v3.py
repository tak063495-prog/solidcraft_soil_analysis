"""Freeze and verify the isolated physical-exclusion candidate and regressions."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
candidate=p/'srm-search-development-v1';runtime=candidate/'runtime';out=candidate/'quality-v3';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
 if (p/'pause-requested').exists():raise SystemExit(2)
active()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
 script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest,freeze
baseline=read(candidate/'baseline635.json');sources=manifest()
added=set(sources)-set(baseline);changed={n for n in baseline if sources[n]!=baseline[n]}
assert added=={'src/srm_bound.rs','tests/srm_bound.rs','docs/SRM_MEAN_STRESS_BOUND.ja.md',
 'examples/srm-mean-stress-mesh.json','examples/srm-mean-stress-problem.json'}
assert changed=={'src/lib.rs','src/main.rs'} and len(sources)==640
assert all(sha(root/n)==h for n,h in baseline.items())
root_solver=root/'target/release/soil-fem-3d.exe';root_solver_sha=sha(root_solver)
assert root_solver_sha=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
freeze(out)
dump(out/'protocol.json',dict(prepared=True,source_files=640,modified=sorted(changed),added=sorted(added),
 old_all_existing_constitutive_and_stage_modules_byte_unchanged=True,
 independent_exclusion_primitive_not_full_SRM_search_or_safety_gate=True,
 required_new_tests=9,all_previous_Rust_and_Python_regressions_required=True,
 actual_new_CLI_and_no_overwrite_and_not_applicable_cases_required=True,
 Root_and_live_Native_frozen_sources_not_modified=True,p6_complete=False))
while True:
 active()
 memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
  'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
  encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
 if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=12*(1<<30) and shutil.disk_usage(out).free>=32*(1<<30):break
 time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,minimum_each_free_bytes=12*(1<<30),
 free_physical_bytes=memory['FreePhysicalMemory']*1024,free_virtual_bytes=memory['FreeVirtualMemory']*1024,
 free_disk_bytes=shutil.disk_usage(out).free))
env=dict(os.environ,CARGO_TARGET_DIR=str(candidate/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',
 PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for name in ('PYTHONPATH','PYTHONHOME'):env.pop(name,None)
records=[]
def call(label,command,expected=0):
 active();assert manifest()==sources
 start=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(list(map(str,command)),cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
   creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(command[0]),arguments=list(map(str,command[1:])),
   created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  code=child.wait()
 unchanged=manifest()==sources and all(sha(root/n)==h for n,h in baseline.items()) and sha(root_solver)==root_solver_sha
 row=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start,
  candidate_and_original_Root_sources_unchanged=unchanged)
 records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
 if code!=expected or not unchanged:
  dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,original_failure_retained=True,p6_complete=False))
  raise SystemExit(1)
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
call('format',[cargo,'fmt','--all','--','--check'])
call('clippy',[cargo,'clippy','--locked','--release','--all-targets','--','-D','warnings'])
call('new-nine-tests',[cargo,'test','--locked','--release','--test','srm_bound','--','--test-threads=2'])
log=(out/'new-nine-tests.log').read_text('utf-8')
assert '9 passed; 0 failed; 0 ignored' in log and 'running 0 tests' not in log
call('all-Rust-tests',[cargo,'test','--locked','--release'])
call('CLI-build',[cargo,'build','--locked','--release','--message-format=json-render-diagnostics'])
executables=[]
for line in (out/'CLI-build.log').read_text('utf-8').splitlines():
 try:item=json.loads(line)
 except json.JSONDecodeError:continue
 if item.get('reason')=='compiler-artifact' and item.get('executable') and item.get('target',{}).get('name')=='soil-fem-3d':executables.append(item['executable'])
assert len(set(executables))==1
binary=Path(executables[-1]);solver=out/'soil-fem-3d.exe'
shutil.copyfile(binary,solver);assert sha(solver)==sha(binary) and sha(solver)!=root_solver_sha
target=runtime/'target/release/soil-fem-3d.exe';target.parent.mkdir(parents=True)
shutil.copyfile(solver,target);assert sha(target)==sha(solver)
call('all-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
mesh=runtime/'examples/srm-mean-stress-mesh.json';problem=runtime/'examples/srm-mean-stress-problem.json'
call('CLI-admissible-mean',[solver,'mc-mean-stress-bound',mesh,problem,out/'mean-not-excluded.json'])
value=read(out/'mean-not-excluded.json');assert value['status']=='evaluated'
assert not value['certificate']['mean_stress_excludes_admissible_discrete_equilibrium']
assert not value['certificate']['lower_bound_or_stability_proven']
modified=read(problem);modified['material']['cohesion_kpa']/=1.5
dump(out/'reduced-problem.json',modified)
call('CLI-exclusion',[solver,'mc-mean-stress-bound',mesh,out/'reduced-problem.json',out/'mean-excluded.json'])
value=read(out/'mean-excluded.json');assert value['status']=='evaluated'
assert value['certificate']['mean_stress_excludes_admissible_discrete_equilibrium']
assert abs(value['certificate']['mean_mc_yield_margin_kpa']-(15.-20./1.5))<1e-10
assert not value['certificate']['continuum_collapse_factor_certified']
modified['boundary']['constraints'][0]['node_set']='all';dump(out/'distributed-support.json',modified)
call('CLI-distributed-not-applicable',[solver,'mc-mean-stress-bound',mesh,out/'distributed-support.json',out/'not-applicable.json'])
assert read(out/'not-applicable.json')['status']=='not_applicable'
before=sha(out/'mean-excluded.json')
call('CLI-refuse-overwrite',[solver,'mc-mean-stress-bound',mesh,problem,out/'mean-excluded.json'],1)
assert sha(out/'mean-excluded.json')==before
dump(out/'terminal.json',dict(passed=True,source_files=640,records=records,actual_new_tests=9,
 all_Rust_and_Python_regressions_passed=True,new_CLI_executed_with_actual_fresh_binary=True,
 CLI_positive_and_negative_and_unavailable_and_no_overwrite_cases_passed=True,
 Root_and_live_Native_sources_and_original_A4_unchanged=True,
 all_existing_model_and_stage_modules_byte_unchanged=True,
 general_SRM_search_Fs_req_and_3d_slope_safety_not_yet_implemented=True,p6_complete=False))
print('PASS isolated640 physical-exclusion primitive and all regressions; full SRM/P6 remain',flush=True)
