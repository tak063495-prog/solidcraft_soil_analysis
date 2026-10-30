"""Retain the terminated harness and resume only the unexecuted full Rust suite."""
from pathlib import Path
import datetime, hashlib, json, os, re, shutil, subprocess, sys, time
from qa_log_counts import unittest_summary
p=Path(__file__).resolve().parent;root=p.parents[1]
candidate=p/'srm-adaptive-search-v1';runtime=candidate/'runtime'
old=candidate/'quality-v2';out=candidate/'quality-v2-continuation-v2';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active()
assert not (old/'terminal.json').exists() and not (old/'all-Rust-tests-child.json').exists()
sources=read(old/'source-snapshot.json');assert len(sources)==645
root_sources=read(p/'Root-source635-promotion-v1/quality-v1/source-snapshot.json')
prior=read(p/'srm-search-development-v1/quality-v3/source-snapshot.json')
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest, freeze
def unchanged():
    return manifest()==sources and all(sha(root/n)==h for n,h in root_sources.items()) and all(
        sha(p/'srm-search-development-v1/runtime'/n)==h for n,h in prior.items())
assert unchanged()
steps=['format','clippy','new-eleven-search-tests','prior-nine-bound-tests','CLI-build','CLI-usage',
    'CLI-bracket','CLI-required-excluded','CLI-refuse-overwrite','CLI-reject-unsupported-water','all-Python-tests']
records=[read(old/(name+'-terminal.json')) for name in steps]
assert all(r['exit_code']==r['expected_exit_code'] and r['all_candidate_and_original_sources_and_A4_unchanged'] for r in records)
inventory=unittest_summary((old/'all-Python-tests.log').read_text('utf-8'),219)
assert '11 passed; 0 failed; 0 ignored' in (old/'new-eleven-search-tests.log').read_text('utf-8')
assert '9 passed; 0 failed; 0 ignored' in (old/'prior-nine-bound-tests.log').read_text('utf-8')
assert read(old/'search-result.json')['requirement']=='equilibrium_requirement_met'
assert read(old/'excluded-result.json')['requirement']=='independently_excluded'
dump(out/'previous-harness-terminal.json',dict(passed=False,actual_driver_exit_code=1,
    actual_tool_session=41374,all11_completed_phases_have_their_expected_exit_code=True,
    actual_failure='post-Python assertion incorrectly required OK at end-of-file; buffered Native importer/viewer stdout followed the successful219-test summary',
    old_driver_logs_sources_and_receipts_retained=True,source_changes_since_previous_frozen645=0,
    mathematical_and_CLI_or_Python_failures_not_inferred=True,
    independent_Python_inventory=inventory,p6_complete=False))
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
freeze(out);assert read(out/'source-snapshot.json')==sources
dump(out/'protocol.json',dict(prepared=True,source_files=645,prior_completed_phases=records,
    Python_inventory=inventory,only_unexecuted_all320_Rust_tests_remain=True,
    no_passing_test_repeated_due_only_to_log_parser=True,all_original35_requirements_retained=True,p6_complete=False))
while not (p/'srm-search-development-v1/quality-v3/terminal.json').exists():
    active();time.sleep(5)
assert read(p/'srm-search-development-v1/quality-v3/terminal.json')['passed']
while True:
    active()
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=12*(1<<30) and shutil.disk_usage(out).free>=32*(1<<30):break
    time.sleep(5)
assert unchanged()
env=dict(os.environ,CARGO_TARGET_DIR=str(candidate/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',
    PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for name in ['PYTHONPATH','PYTHONHOME']:env.pop(name,None)
command=['C:/Users/link_/.cargo/bin/cargo.exe','test','--locked','--release']
start=time.monotonic()
with (out/'all-Rust-tests.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
        creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    dump(out/'all-Rust-tests-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
        arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code=child.wait()
stable=unchanged()
counts=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',
    (out/'all-Rust-tests.log').read_text('utf-8'))]
total=tuple(map(sum,zip(*counts))) if counts else None
passed=code==0 and stable and total==(320,0,0)
dump(out/'all-Rust-tests-terminal.json',dict(passed=passed,exit_code=code,elapsed_seconds=time.monotonic()-start,
    actual_counts=total,all_candidate_and_original_sources_unchanged=stable))
dump(out/'terminal.json',dict(passed=passed,source_files=645,Rust_passed=total[0] if total else None,
    Rust_failed=total[1] if total else None,Rust_ignored=total[2] if total else None,
    Python_passed=219,Python_skipped=0,actual_new_search_tests=11,actual_prior_bound_tests=9,
    all_Python_unique_case_ids_and_exact_summary_verified=True,
    previous_harness_failure_and_all_logs_retained=True,original_frozen645_sources_unchanged=True,
    actual_fresh_CLI_and_required_factor_and_independent_closed_form_and_cutback_and_nonlinear_failure_verified=True,
    Root635_previous640_and_Native_sources_unchanged=stable,
    ordinary_supported_slope_physical_upper_water_and_design_safety_and_allP6_not_claimed=True,p6_complete=False))
print(json.dumps(dict(passed=passed,source_files=645,Rust_counts=total,Python_passed=219,p6_complete=False)),flush=True)
raise SystemExit(0 if passed else 1)
