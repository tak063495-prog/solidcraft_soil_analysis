"""Make the exercised MCC partial-cut stop regression part of ordinary cargo test."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];adoption=p/'Root-source656-promotion-v1';q=adoption/'quality-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
assert read(q/'terminal.json')['passed'] and not (p/'pause-requested').exists()
sys.path.insert(0,str(root/'tools'));from freeze_sources import manifest,freeze
before=read(q/'source-snapshot.json');assert manifest()==before and len(before)==656
out=adoption/'MCC-regressions-v1';out.mkdir()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
incoming=p/'mcc_stop_boundaries_production.rs';destination=root/'tests/mcc_stop_boundaries.rs'
assert not destination.exists();shutil.copyfile(incoming,destination);assert sha(incoming)==sha(destination)
sources=manifest();assert len(sources)==657 and set(sources)-set(before)=={'tests/mcc_stop_boundaries.rs'}
assert all(sources[n]==h for n,h in before.items());freeze(out)
env=dict(os.environ,CARGO_TARGET_DIR=str(root/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',
    PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
records=[];cargo='C:/Users/link_/.cargo/bin/cargo.exe'
for label,args in [('format',[cargo,'fmt','--all','--','--check']),
    ('new-MCC-clippy',[cargo,'clippy','--locked','--release','--test','mcc_stop_boundaries','--','-D','warnings']),
    ('new-MCC-tests',[cargo,'test','--locked','--release','--test','mcc_stop_boundaries','--','--test-threads=2'])]:
    if (p/'pause-requested').exists():raise SystemExit(2)
    assert manifest()==sources;start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(args,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=args[0],arguments=args[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start,sources_unchanged=manifest()==sources)
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    if code or not row['sources_unchanged']:
        dump(out/'terminal.json',dict(passed=False,records=records,actual_failure_retained=True,p6_complete=False));raise SystemExit(1)
assert '2 passed; 0 failed; 0 ignored' in (out/'new-MCC-tests.log').read_text('utf-8')
assert sha(root/'target/release/soil-fem-3d.exe')==read(q/'terminal.json')['binary_sha256']
dump(out/'terminal.json',dict(passed=True,Root_source_files=657,actual_new_Root_Cargo_tests_passed=2,
    mixed_and_projected_MCC_early_late_partial_local_cut_full_history_byte_equality=True,
    includes_actual_plastic_cutback_and_dormant_GP_and_release_state=True,
    exact_accepted_seal_JSON_roundtrip_and_entire_resumed_scientific_bytes_excluding_elapsed=True,
    only_test_added_all_existing656_sources_and_actual_CLI_bytes_unchanged=True,
    full347_of_identical_numerical_source_previously_passed_and_new_Root2_focused=True,
    expected_combined_CI_Rust_tests=349,full349_Root_suite_not_rerun_or_claimed=True,
    Root219_Python_and_actual_CLI_physics_qualification_reused=True,records=records,p6_complete=False))
print('PASS two ordinary cargo MCC stop regressions; full349 remote CI remains',flush=True)
