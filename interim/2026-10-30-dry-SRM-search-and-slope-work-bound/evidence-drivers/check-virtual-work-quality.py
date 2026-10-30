"""Freeze650, exercise full pointwise work against independent FEM, then regressions."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from qa_log_counts import unittest_summary
import virtual_work_references as reference
p=Path(__file__).resolve().parent;root=p.parents[1]
candidate=p/'srm-general-virtual-work-v1';runtime=candidate/'runtime';out=candidate/'quality-v1';out.mkdir()
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
sources=manifest();base=read(candidate/'baseline645.json')
assert len(sources)==650 and set(sources)-set(base)=={
    'src/srm_virtual_work.rs','tests/srm_virtual_work.rs','docs/MC_VIRTUAL_WORK_BOUND.ja.md',
    'examples/mc-virtual-work-mesh.json','examples/mc-virtual-work-problem.json'}
assert {n for n in base if sources[n]!=base[n]}=={
    'src/lib.rs','src/main.rs','src/srm_search.rs','docs/SRM_SEARCH.ja.md'}
prior=p/'srm-adaptive-search-v1/runtime'
def unchanged():return manifest()==sources and all(sha(prior/n)==h for n,h in base.items())
assert unchanged()
freeze(out)
root635=read(p/'Root-source635-promotion-v1/quality-v1/source-snapshot.json')
root_initial={n:sha(root/n) for n in root635}
dump(out/'source-binding.json',dict(prepared=True,source_files=650,base_frozen645_files=645,
    baseline_private645_source_manifest_sha256=sha(candidate/'baseline645.json'),
    root_at_start_matches_frozen635=root_initial==root635,Root_start_files_sha256=root_initial,
    live_Root_not_a_runtime_dependency=True,
    Root_may_be_intentionally_upgraded_after_its_existing_source_bound_jobs_finish=True,
    private645_and_new650_must_stay_frozen=True,
    driver_sha256=sha(Path(__file__)),reference_helper_sha256=sha(p/'virtual_work_references.py'),
    log_inventory_helper_sha256=sha(p/'qa_log_counts.py'),p6_complete=False))
cases=reference.prepare(runtime,out,candidate/'independent-FEM-v1')
dump(out/'protocol.json',dict(prepared=True,source_files=650,new_virtual_work_tests=12,
    prior_search_tests=11,prior_mean_tests=9,independent_FEM_CLI_cases=6,
    all504_reference_point_records_and_raw_framed_point_hashes_required=True,
    full332_Rust_and219_Python_required=True,
    no_material_kinematics_or_existing_force_yield_tolerance_changes=True,
    actual_slope_and_automatic_candidate_and_fullP6_not_from_cube_tests=True,p6_complete=False))
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
for name in ['PYTHONPATH','PYTHONHOME']:env.pop(name,None)
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
    row=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start,
        both_private650_and645_sources_unchanged=stable,
        current_Root_matches_initial635=all((root/n).exists() and sha(root/n)==h for n,h in root635.items()))
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    if code!=expected or not stable:
        dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,actual_failure_retained=True,p6_complete=False))
        raise SystemExit(1)
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
call('format',[cargo,'fmt','--all','--','--check'])
call('clippy',[cargo,'clippy','--locked','--release','--all-targets','--','-D','warnings'])
call('new-twelve-tests',[cargo,'test','--locked','--release','--test','srm_virtual_work','--','--test-threads=2'])
assert '12 passed; 0 failed; 0 ignored' in (out/'new-twelve-tests.log').read_text('utf-8')
call('prior-search-eleven-tests',[cargo,'test','--locked','--release','--test','srm_search','--','--test-threads=2'])
assert '11 passed; 0 failed; 0 ignored' in (out/'prior-search-eleven-tests.log').read_text('utf-8')
call('prior-mean-nine-tests',[cargo,'test','--locked','--release','--test','srm_bound','--','--test-threads=2'])
assert '9 passed; 0 failed; 0 ignored' in (out/'prior-mean-nine-tests.log').read_text('utf-8')
call('CLI-build',[cargo,'build','--locked','--release','--message-format=json-render-diagnostics'])
executables=[]
for line in (out/'CLI-build.log').read_text('utf-8').splitlines():
    try:item=json.loads(line)
    except json.JSONDecodeError:continue
    if item.get('reason')=='compiler-artifact' and item.get('target',{}).get('name')=='soil-fem-3d' and item.get('executable'):executables.append(item['executable'])
assert len(set(executables))==1
binary=Path(executables[-1]);solver=out/'soil-fem-3d.exe';shutil.copyfile(binary,solver)
assert sha(solver)==sha(binary) and sha(solver)!=sha(p/'srm-adaptive-search-v1/quality-v2/soil-fem-3d.exe')
dest=runtime/'target/release/soil-fem-3d.exe';dest.parent.mkdir(parents=True);shutil.copyfile(solver,dest)
call('CLI-usage',[solver,'invalid-virtual-work-command'],1)
assert 'mc-virtual-work-bound MESH' in (out/'CLI-usage.log').read_text('utf-8')
verified=[]
for case in cases:
    call('CLI-independent-'+str(case['index']),[solver,'mc-virtual-work-bound',case['mesh'],case['problem'],case['output']])
    verified.append(reference.verify(case))
assert sum(r['all_GP_checked'] for r in verified)==504
dump(out/'full-independent-reference-comparison.json',dict(passed=True,all6_actual_CLI_cases=True,
    all504_GP_records_checked=True,all_framed_point_hashes_from_raw_CLI_lexemes_verified=True,
    all_same_reference_data_and_original2e_minus12_numeric_comparison_limits=True,
    cases=verified,p6_complete=False))
before=sha(Path(cases[0]['output']))
call('CLI-refuse-overwrite',[solver,'mc-virtual-work-bound',cases[0]['mesh'],cases[0]['problem'],cases[0]['output']],1)
assert sha(Path(cases[0]['output']))==before
bad=read(cases[0]['problem']);bad['virtual_displacement_m'][0][0]=1e-15
dump(out/'supported-component-violation.json',bad)
call('CLI-supported-component-unavailable',[solver,'mc-virtual-work-bound',cases[0]['mesh'],out/'supported-component-violation.json',out/'not-applicable.json'])
unavailable=read(out/'not-applicable.json')
assert unavailable['assessment']['status']=='not_applicable' and not unavailable['complete_point_audit']
bad=read(cases[0]['problem']);bad['water_policy']='unsupported'
dump(out/'unsupported-water.json',bad)
call('CLI-reject-water',[solver,'mc-virtual-work-bound',cases[0]['mesh'],out/'unsupported-water.json',out/'must-not-exist.json'],1)
assert not (out/'must-not-exist.json').exists()
call('all-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
inventory=unittest_summary((out/'all-Python-tests.log').read_text('utf-8'),219)
dump(out/'Python-inventory.json',inventory)
while not (p/'srm-adaptive-search-v1/quality-v2-continuation-v2/terminal.json').exists():
    active();time.sleep(5)
assert read(p/'srm-adaptive-search-v1/quality-v2-continuation-v2/terminal.json')['passed']
call('all-Rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',
    (out/'all-Rust-tests.log').read_text('utf-8'))]
assert counts and tuple(map(sum,zip(*counts)))==(332,0,0)
dump(out/'terminal.json',dict(passed=True,source_files=650,Rust_passed=332,Python_passed=219,
    ignored_skipped_failed=0,actual_new12_and_prior11_and9_passed=True,
    actual_distributed_support_frictional_search_and_real_plastic_FEM_and_numerical_unknowns_tested=True,
    actual6_fresh_CLI_cases_and_all504_independent_reference_points_and_raw_point_hashes_verified=True,
    actual_CLI_help_no_overwrite_supported_component_unavailable_water_rejection_passed=True,
    private645_and650_sources_unchanged=True,Root_was_not_used_as_a_mutable_runtime_dependency=True,
    automatic_actual_slope_and_full3D_convergence_and_allP6_not_claimed=True,records=records,p6_complete=False))
print('PASS isolated650 pointwise MC work witness and all regressions; actual slope/P6 remain',flush=True)
