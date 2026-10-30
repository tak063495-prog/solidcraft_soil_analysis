"""Freeze the controlled search candidate and verify physical restart equivalence."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
from qa_log_counts import unittest_summary
p=Path(__file__).resolve().parent;candidate=p/'srm-controlled-search-v1'
runtime=candidate/'runtime';out=candidate/'quality-v4';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active()
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest,freeze
base=read(candidate/'baseline650.json');sources=manifest()
assert len(base)==650 and len(sources)==653
assert set(sources)-set(base)=={'src/srm_search_control.rs','tests/srm_search_control.rs','docs/SRM_CONTROLLED_SEARCH.ja.md'}
assert {n for n in base if sources[n]!=base[n]}=={'src/lib.rs','src/srm_search.rs'}
def unchanged():
    return manifest()==sources and all(sha(p/'srm-general-virtual-work-v1/runtime'/n)==h for n,h in base.items())
assert unchanged();freeze(out)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out/'protocol.json',dict(prepared=True,source_files=653,new_scientific_restart_tests=6,
    prior_virtual_tests=12,prior_search_tests=11,prior_mean_tests=9,
    every_saved_snapshot_serialized_and_restored=True,
    all338_Rust_and219_Python_required=True,only_elapsed_seconds_removed_from_history_comparison=True,
    original650_frozen_and_existing_large_jobs_unchanged=True,
    file_job_recovery_and_general_slope_exact_Fs_still_pending=True,p6_complete=False))
while True:
    active()
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
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
        elapsed_seconds=time.monotonic()-start,both653_and650_sources_unchanged=stable)
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    if code!=expected or not stable:
        dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,actual_failure_retained=True,p6_complete=False))
        raise SystemExit(1)
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
call('format',[cargo,'fmt','--all','--','--check'])
call('clippy',[cargo,'clippy','--locked','--release','--all-targets','--','-D','warnings'])
for label,name,count in [('new-six-restart-tests','srm_search_control',6),('prior-twelve-virtual-tests','srm_virtual_work',12),
    ('prior-eleven-search-tests','srm_search',11),('prior-nine-mean-tests','srm_bound',9)]:
    call(label,[cargo,'test','--locked','--release','--test',name,'--','--test-threads=2'])
    assert f'{count} passed; 0 failed; 0 ignored' in (out/(label+'.log')).read_text('utf-8')
call('CLI-build',[cargo,'build','--locked','--release','--message-format=json-render-diagnostics'])
executables=[]
for line in (out/'CLI-build.log').read_text('utf-8').splitlines():
    try:item=json.loads(line)
    except json.JSONDecodeError:continue
    if item.get('reason')=='compiler-artifact' and item.get('target',{}).get('name')=='soil-fem-3d' and item.get('executable'):executables.append(item['executable'])
assert len(set(executables))==1
solver=out/'soil-fem-3d.exe';shutil.copyfile(executables[-1],solver)
assert sha(solver)==sha(Path(executables[-1]))
dest=runtime/'target/release/soil-fem-3d.exe';dest.parent.mkdir(parents=True);shutil.copyfile(solver,dest)
call('all-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
inventory=unittest_summary((out/'all-Python-tests.log').read_text('utf-8'),219)
dump(out/'Python-inventory.json',inventory)
while not (p/'srm-general-virtual-work-v1/quality-v1-continuation-v1/terminal.json').exists():
    active();time.sleep(5)
assert read(p/'srm-general-virtual-work-v1/quality-v1-continuation-v1/terminal.json')['passed']
call('all-Rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',
    (out/'all-Rust-tests.log').read_text('utf-8'))]
assert counts and tuple(map(sum,zip(*counts)))==(338,0,0)
dump(out/'terminal.json',dict(passed=True,source_files=653,Rust_passed=338,Python_passed=219,
    ignored_skipped_failed=0,actual6_new_restart_and_prior12_and11_and9_passed=True,
    all_scientific_fields_match_uninterrupted_except_elapsed_seconds=True,
    all_existing_constitutive_and_staged_FEM_algorithms_unchanged=True,
    private650_and653_sources_unchanged=True,file_jobs_and_general_slope_exact_Fs_and_fullP6_not_claimed=True,
    records=records,p6_complete=False))
print('PASS isolated653 scientific restart equivalence; durable file jobs and fullP6 remain',flush=True)
