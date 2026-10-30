"""Frozen656: real CLI lifecycle/failure recovery and complete regressions."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
from qa_log_counts import unittest_summary
p=Path(__file__).resolve().parent;candidate=p/'srm-durable-job-v1'
runtime=candidate/'runtime';out=candidate/'quality-v3';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active();sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest,freeze
base=read(candidate/'baseline653.json');sources=manifest()
assert len(base)==653 and len(sources)==656
assert set(sources)-set(base)=={'src/srm_search_job.rs','tests/srm_search_job.rs','docs/SRM_SEARCH_JOBS.ja.md'}
assert {n for n in base if sources[n]!=base[n]}=={'src/lib.rs','src/main.rs','src/stages.rs'}
def unchanged():
    return manifest()==sources and all(sha(p/'srm-controlled-search-v1/runtime'/n)==h for n,h in base.items())
assert unchanged();freeze(out)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out/'protocol.json',dict(prepared=True,source_files=656,new_actual_CLI_tests=5,new_stop_boundary_regressions=4,
    frozen_parent653_scientific_restart_tests_required_before_new_job_tests=True,
    all347_Rust_and219_Python_required=True,
    actual_frozen_solver_immutable_inputs_one_shot_cooperative_stop_and_recovery_required=True,
    actual_corrupted_slot_foreign_pointer_wrong_stop_write_failure_and_finished_report_recovery_required=True,
    only_elapsed_seconds_removed_from_scientific_comparisons=True,
    all_original35_requirements_and_original5percent_displacement_FAIL_retained=True,p6_complete=False))
while True:
    active();memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=16*(1<<30) and shutil.disk_usage(out).free>=32*(1<<30):break
    time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,free_physical_bytes=memory['FreePhysicalMemory']*1024,
    free_virtual_bytes=memory['FreeVirtualMemory']*1024,free_disk_bytes=shutil.disk_usage(out).free))
env=dict(os.environ,CARGO_TARGET_DIR=str(candidate/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',
    PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for name in ['PYTHONPATH','PYTHONHOME']:env.pop(name,None)
records=[]
def call(label,args,expected=0):
    active();assert unchanged();start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(args[0]),
            arguments=list(map(str,args[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    stable=unchanged();row=dict(step=label,exit_code=code,expected_exit_code=expected,
        elapsed_seconds=time.monotonic()-start,both656_and653_sources_unchanged=stable)
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    if code!=expected or not stable:
        dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,actual_failure_retained=True,p6_complete=False))
        raise SystemExit(1)
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
call('format',[cargo,'fmt','--all','--','--check'])
call('clippy',[cargo,'clippy','--locked','--release','--all-targets','--','-D','warnings'])
prior=p/'srm-controlled-search-v1/quality-v4/new-six-restart-tests-terminal.json'
while not prior.exists():active();time.sleep(5)
assert read(prior)['exit_code']==0 and read(prior)['both653_and650_sources_unchanged']
call('five-CLI-and-four-boundary-tests',[cargo,'test','--locked','--release','--test','srm_search_job','--','--test-threads=2'])
assert '9 passed; 0 failed; 0 ignored' in (out/'five-CLI-and-four-boundary-tests.log').read_text('utf-8')
call('CLI-build',[cargo,'build','--locked','--release','--message-format=json-render-diagnostics'])
executables=[]
for line in (out/'CLI-build.log').read_text('utf-8').splitlines():
    try:item=json.loads(line)
    except json.JSONDecodeError:continue
    if item.get('reason')=='compiler-artifact' and item.get('target',{}).get('name')=='soil-fem-3d' and item.get('executable'):executables.append(item['executable'])
assert len(set(executables))==1
solver=out/'soil-fem-3d.exe';shutil.copyfile(executables[-1],solver);assert sha(solver)==sha(Path(executables[-1]))
dest=runtime/'target/release/soil-fem-3d.exe';dest.parent.mkdir(parents=True);shutil.copyfile(solver,dest)
call('CLI-command-help',[solver,'unknown-srm-job-command'],1)
helptext=(out/'CLI-command-help.log').read_text('utf-8')
for command in ['prepare-srm-search-job','run-srm-search-job','prepare-srm-search-recovery','stop-srm-search-job','srm-search-job-status','verify-srm-search-job']:
    assert command in helptext
call('all-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
dump(out/'Python-inventory.json',unittest_summary((out/'all-Python-tests.log').read_text('utf-8'),219))
parent=p/'srm-controlled-search-v1/quality-v4/terminal.json'
while not parent.exists():active();time.sleep(5)
assert read(parent)['passed']
call('all-Rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',
    (out/'all-Rust-tests.log').read_text('utf-8'))]
assert counts and tuple(map(sum,zip(*counts)))==(347,0,0)
dump(out/'terminal.json',dict(passed=True,source_files=656,Rust_passed=347,Python_passed=219,
    ignored_skipped_failed=0,actual5_CLI_lifecycle_corruption_write_failure_and_recovery_tests_passed=True,
    all_scientific_fields_match_uninterrupted_except_elapsed_seconds=True,
    completed_search_indeterminate_requirement_not_promoted_to_safety=True,
    constitutive_force_yield_and_trial_search_algorithms_unchanged=True,staged_resume_stop_control_and_atomic_request_publication_fixed=True,
    private653_and656_sources_unchanged=True,
    native_gate_UI_general_slope_exact_Fs_and_fullP6_not_claimed=True,records=records,p6_complete=False))
print('PASS isolated656 durable SRM CLI; fullP6 and native safety gate remain',flush=True)
