"""Freeze and type/lint the new responsive-loader candidate without old edits."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent
native=p/'native-async-view'
original=p/'native-owned-v4'
out=p/'native-async-view-quality-v4'
out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,value):
    with q.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False)
def pause_check():
    if (p/'pause-requested').exists():raise SystemExit(2)
pause_check()
assert read(p/'native-async-view-preparation-v1/terminal.json')['passed']
prior=read(p/'native-owned-v4-quality-v1/source-manifest.json')
sources={name:sha(native/name) for name in prior['files_sha256']}
changed={name for name,h in prior['files_sha256'].items() if sources[name]!=h}
assert changed=={'ROADMAP.md','crates/engine/src/cmd/soil_fem.rs','crates/engine/src/soil_results/native.rs',
    'crates/engine/src/soil_results/tests.rs','crates/ui-egui/src/lib.rs','crates/ui-egui/src/soil_results.rs'}
assert len(sources)==73
assert all(sha(original/name)==h for name,h in prior['files_sha256'].items())
snapshot=dict(base_commit=prior['base_commit'],files_sha256=sources)
dump(out/'source-manifest.json',snapshot)
with zipfile.ZipFile(out/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as z:
    for name,h in sorted(sources.items()):
        raw=(native/name).read_bytes();assert hashlib.sha256(raw).hexdigest()==h;z.writestr(name,raw)
with zipfile.ZipFile(out/'source-snapshot.zip') as z:
    assert set(z.namelist())==set(sources) and len(z.namelist())==73 and z.testzip() is None
    assert all(hashlib.sha256(z.read(name)).hexdigest()==h for name,h in sources.items())
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out/'protocol.json',dict(prepared=True,source_files=73,changed=sorted(changed),
    old_owned73_source_Git_fullCI_and_real_workflow_not_modified=True,
    existing_command_path_only_no_new_icons_assets_or_design_serialization=True,
    existing_scientific_parser_and_configure_computation_retained=True,
    background_status_commits_only_after_owned_thread_exit_and_current_input_match=True,
    one_process_wide_reader_permit_and_one_pending_reader_per_session=True,
    cancellation_discards_without_claiming_instant_reader_exit=True,
    no_numerical_FEM_solve_by_this_preflight=True,
    four_new_existing_area_regressions_must_execute=True, explicit_object_only_open_close_and_invalid_input_retains_previous_arc=True, existing_public_null_as_omitted_compatibility_and_direct_leaf_null_rejection=True, default_close_cancels_and_clears_previous_without_resurrection=True,
    required=['format','engine_and_UI_all_targets_clippy'],
    fullCI_wasm_public525_realGUI_still_required=True,
    p6_complete=False))
while True:
    pause_check()
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=12*(1<<30) and shutil.disk_usage(native).free>=32*(1<<30):break
    time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,minimum_each_free_bytes=12*(1<<30),
    free_physical_bytes=memory['FreePhysicalMemory']*1024,free_virtual_bytes=memory['FreeVirtualMemory']*1024,
    free_disk_bytes=shutil.disk_usage(native).free))
env=dict(os.environ,CARGO_TARGET_DIR=str(native/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',PYTHONUTF8='1')
records=[]
for label,args in [('format',['fmt','--all','--','--check']),
    ('engine-and-UI-clippy',['clippy','--locked','--profile','ci','-p','solidcraft-engine','-p','solidcraft-ui-egui','--all-targets','--','-D','warnings']),
    ('result-reader-tests',['test','--locked','--profile','ci','-p','solidcraft-engine','--lib','soil_results::','--','--test-threads=1']),
    ('owned-job-tests',['test','--locked','--profile','ci','-p','solidcraft-engine','--lib','soil_jobs::','--','--test-threads=1']),
    ('UI-library-tests',['test','--locked','--profile','ci','-p','solidcraft-ui-egui','--lib','--','--test-threads=2'])]:
    pause_check()
    assert all(sha(native/name)==h for name,h in sources.items())
    start=time.monotonic()
    command=['C:/Users/link_/.cargo/bin/cargo.exe',*args]
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=native,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],
            created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    unchanged=all(sha(native/name)==h for name,h in sources.items()) and all(sha(original/name)==h for name,h in prior['files_sha256'].items())
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start,candidate_and_original_sources_unchanged=unchanged)
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    assert code==0 and unchanged,label
regressions = [
    'background_public_load_preserves_legacy_mcc_epoch_values_and_integrity_guards',
    'cancelled_background_keeps_previous_result_and_owns_the_reader_until_exit',
    'changed_cad_rejects_delayed_verified_result_and_retains_the_previous_arc',
    'disconnected_or_panicked_background_worker_reports_failure_without_losing_data']
reader_log=(out/'result-reader-tests.log').read_text('utf-8')
assert all(('::'+name+' ... ok') in reader_log for name in regressions)
assert 'running 0 tests' not in reader_log and 'test result: ok.' in reader_log
dump(out/'terminal.json',dict(passed=True,source_files=73,records=records,
    actual_four_regressions_executed=regressions, type_and_lint_preflight_only=False,unit_tests_not_yet_executed=False,
    actual_large_async_GUI_and_fullCI_wasm_and_new_owned_workflow_and_allP6_separate=True,p6_complete=False))
print('PASS frozen async v2 reader/owned/UI regressions; fullCI/wasm/largeGUI remain separate',flush=True)
