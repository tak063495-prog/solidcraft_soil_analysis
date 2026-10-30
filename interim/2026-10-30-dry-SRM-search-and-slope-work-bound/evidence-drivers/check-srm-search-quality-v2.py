"""Freeze and actually exercise the separate645 fresh-history MC search."""
from pathlib import Path
import datetime, hashlib, json, os, re, shutil, subprocess, sys, time

p = Path(__file__).resolve().parent
root = p.parents[1]
candidate = p / 'srm-adaptive-search-v1'
runtime = candidate / 'runtime'
out = candidate / 'quality-v2'
out.mkdir()
read = lambda q: json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
def dump(q, v):
    with q.open('x', encoding='utf-8') as f:
        json.dump(v, f, indent=2, allow_nan=False)
def active():
    if (p / 'pause-requested').exists():
        raise SystemExit(2)
active()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest, freeze
sources = manifest()
base = read(candidate/'baseline640.json')
assert len(sources)==645 and set(sources)-set(base)=={
    'src/srm_search.rs','src/stages/srm_adapter.rs','tests/srm_search.rs',
    'examples/srm-search.json','docs/SRM_SEARCH.ja.md'}
assert {n for n in base if sources[n]!=base[n]}=={
    'src/lib.rs','src/main.rs','src/stages.rs','src/srm_bound.rs','docs/SRM_MEAN_STRESS_BOUND.ja.md'}
root_sources = read(p/'Root-source635-promotion-v1/quality-v1/source-snapshot.json')
prior_sources = read(p/'srm-search-development-v1/quality-v3/source-snapshot.json')
root_solver = root/'target/release/soil-fem-3d.exe'
assert sha(root_solver)=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
def unchanged():
    return (manifest()==sources and all(sha(root/n)==h for n,h in root_sources.items())
        and all(sha(p/'srm-search-development-v1/runtime'/n)==h for n,h in prior_sources.items())
        and sha(root_solver)=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114')
assert unchanged()
freeze(out)
dump(out/'protocol.json',dict(prepared=True,source_files=645,new_search_tests=11,prior_bound_tests=9,
    full320_Rust_and219_Python_required=True,independent_private_target=True,
    original_live_Root635_and640_and_Native_not_modified=True,
    actual_nonlinear_history_and_theory_and_failure_and_CLI_required=True,
    limited_dry_MC_search_not_general_supported_slope_or_safety_certificate=True,p6_complete=False))
while True:
    active()
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=12*(1<<30) and shutil.disk_usage(out).free>=32*(1<<30):break
    time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,free_physical_bytes=memory['FreePhysicalMemory']*1024,
    free_virtual_bytes=memory['FreeVirtualMemory']*1024,free_disk_bytes=shutil.disk_usage(out).free))
env=dict(os.environ,CARGO_TARGET_DIR=str(candidate/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',
    PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for key in ('PYTHONPATH','PYTHONHOME'):env.pop(key,None)
records=[]
def call(label,args,expected=0):
    active();assert unchanged()
    start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(args[0]),
            arguments=list(map(str,args[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    stable=unchanged()
    item=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start,
        all_candidate_and_original_sources_and_A4_unchanged=stable)
    records.append(item);dump(out/(label+'-terminal.json'),item);print(json.dumps(item),flush=True)
    if code!=expected or not stable:
        dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,actual_failure_retained=True,p6_complete=False))
        raise SystemExit(1)
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
call('format',[cargo,'fmt','--all','--','--check'])
call('clippy',[cargo,'clippy','--locked','--release','--all-targets','--','-D','warnings'])
call('new-eleven-search-tests',[cargo,'test','--locked','--release','--test','srm_search','--','--test-threads=2'])
assert '11 passed; 0 failed; 0 ignored' in (out/'new-eleven-search-tests.log').read_text('utf-8')
call('prior-nine-bound-tests',[cargo,'test','--locked','--release','--test','srm_bound','--','--test-threads=2'])
assert '9 passed; 0 failed; 0 ignored' in (out/'prior-nine-bound-tests.log').read_text('utf-8')
call('CLI-build',[cargo,'build','--locked','--release','--message-format=json-render-diagnostics'])
executables=[]
for line in (out/'CLI-build.log').read_text('utf-8').splitlines():
    try:item=json.loads(line)
    except json.JSONDecodeError:continue
    if item.get('reason')=='compiler-artifact' and item.get('target',{}).get('name')=='soil-fem-3d' and item.get('executable'):executables.append(item['executable'])
assert len(set(executables))==1
binary=Path(executables[-1]);solver=out/'soil-fem-3d.exe'
shutil.copyfile(binary,solver);assert sha(binary)==sha(solver) and sha(solver)!=sha(root_solver)
dest=runtime/'target/release/soil-fem-3d.exe';dest.parent.mkdir(parents=True)
shutil.copyfile(solver,dest);assert sha(dest)==sha(solver)
call('CLI-usage',[solver,'invalid-srm-command'],1)
assert 'search-mc-srm MESH' in (out/'CLI-usage.log').read_text('utf-8')
# Construct the same independently specified own cube for an actual CLI run.
problem=read(runtime/'examples/srm-mean-stress-problem.json')
ground=read(runtime/'examples/ground-geostatic-mc-column.json')
seed=ground['initialization']['geostatic']['seed']['states']['LowerClay']
ground['initialization']['geostatic']['seed']['states']={'TestCube':seed}
ground['initialization']['geostatic']['verification_reject_once']=False
ground['materials']=[dict(id=1,name='Own unconfined CLI cube',unit_weight_kn_m3=1e-10,mohr_coulomb=problem['material'])]
ground['bodies']=[dict(body='TestCube',material_id=1)]
case=read(runtime/'examples/stages-geostatic-column.json')
case['boundary']=problem['boundary'];case['monitors']={};case['monitor_sets']={}
case['stages']=[dict(kind='boundary',name='load',target_factor=1.,initial_increment=1.)]
case['verification_reject_once_stages']=[]
dump(out/'ground.json',ground);dump(out/'case.json',case)
mesh=runtime/'examples/srm-mean-stress-mesh.json';spec=read(runtime/'examples/srm-search.json')
call('CLI-bracket',[solver,'search-mc-srm',mesh,out/'ground.json',out/'case.json',runtime/'examples/srm-search.json',out/'search-result.json'])
report=read(out/'search-result.json');bounds=report['interval']
assert report['requirement']=='equilibrium_requirement_met' and report['stop_reason']=='factor_tolerance'
assert bounds['attained_equilibrium_factor']<20./15.<bounds['independently_excluded_factor']
assert bounds['width']<=spec['absolute_factor_tolerance']
assert report['unreduced_source_history_unchanged'] and not report['continuum_stability_or_design_safety_certified']
spec['required_factor']=1.6;dump(out/'excluded-spec.json',spec)
call('CLI-required-excluded',[solver,'search-mc-srm',mesh,out/'ground.json',out/'case.json',out/'excluded-spec.json',out/'excluded-result.json'])
assert read(out/'excluded-result.json')['requirement']=='independently_excluded'
before=sha(out/'search-result.json')
call('CLI-refuse-overwrite',[solver,'search-mc-srm',mesh,out/'ground.json',out/'case.json',runtime/'examples/srm-search.json',out/'search-result.json'],1)
assert sha(out/'search-result.json')==before
ground['initialization']['water_level_m']=0.5;dump(out/'unsupported-wet-ground.json',ground)
call('CLI-reject-unsupported-water',[solver,'search-mc-srm',mesh,out/'unsupported-wet-ground.json',out/'case.json',runtime/'examples/srm-search.json',out/'must-not-exist.json'],1)
assert not (out/'must-not-exist.json').exists()
call('all-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
assert re.search(r'Ran 219 tests in [\d.]+s\s+OK\s*$',(out/'all-Python-tests.log').read_text('utf-8'))
# Different caches allow the focused implementation work above without ever
# writing to the still-live640 target. Do not duplicate the old full suite while
# its own result remains pending.
while not (p/'srm-search-development-v1/quality-v3/terminal.json').exists():
    active();time.sleep(5)
assert read(p/'srm-search-development-v1/quality-v3/terminal.json')['passed']
call('all-Rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',
    (out/'all-Rust-tests.log').read_text('utf-8'))]
assert counts and tuple(map(sum,zip(*counts)))==(320,0,0)
dump(out/'terminal.json',dict(passed=True,source_files=645,records=records,actual_new_search_tests=11,
    actual_previous_bound_tests=9,Rust_passed=320,Python_passed=219,all_ignored_skipped_failed=0,
    actual_closed_form_and_plastic_history_and_cutback_restoration_and_failure_cases_passed=True,
    actual_fresh_CLI_bracket_required_exclusion_water_rejection_and_no_overwrite_passed=True,
    original_Root635_previous640_and_Native_sources_unchanged=True,
    ordinary_distributed_support_physical_upper_water_and_design_safety_and_allP6_not_claimed=True,p6_complete=False))
print('PASS isolated645 limited dry MC search and regressions; general SRM/P6 remain',flush=True)
