"""Format then freeze the new private API; build and run exact old+new checks."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];base=p/'epoch-SRM-job-v1';runtime=base/'runtime';out=base/'quality-v2'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
assert read(base/'preparation.json')['passed'];active();out.mkdir()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
for gate in [p/'epoch-SRM-search-v3/qualification-v2/terminal.json',p/'epoch-SRM-control-v2/quality-v1/terminal.json']:
    while not gate.exists():active();time.sleep(5)
    assert read(gate)['passed'],str(gate)

while True:
    active();mem=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],encoding='utf-8-sig'))
    if min(mem['FreePhysicalMemory'],mem['FreeVirtualMemory'])*1024>=12*(1<<30) and shutil.disk_usage(base).free>=24*(1<<30):break
    time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,minimum_free_physical_virtual_GiB=12,minimum_free_disk_GiB=24,observed=mem))
env=dict(os.environ,CARGO_TARGET_DIR=str(base/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='1',
    OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
env.pop('PYTHONPATH',None);env.pop('PYTHONHOME',None)
cargo='C:/Users/link_/.cargo/bin/cargo.exe';python=p/'source635-offline-setup-v2/verification env with spaces/Scripts/python.exe'
original=read(base/'Root657-source-snapshot.json');records=[];frozen=None
def check():
    active();assert all(sha(root/n)==h for n,h in original.items())
    if frozen is not None:assert all(sha(runtime/n)==h for n,h in frozen.items())
def call(label,command):
    check();start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,command)),cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(command[0]),arguments=list(map(str,command[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    check();row=dict(phase=label,exit_code=code,elapsed_seconds=time.monotonic()-start);records.append(row);dump(out/(label+'-terminal.json'),row)
    print(json.dumps(row),flush=True)
    if code:
        dump(out/'terminal.json',dict(passed=False,failed_phase=label,records=records,
            candidate_sources_frozen_and_Root_and_live_jobs_unchanged=True,p6_complete=False));raise SystemExit(code)
call('format',[cargo,'fmt','--all'])
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest,freeze
frozen=manifest();parent=read(base/'parent662-source-snapshot.json')
assert len(frozen)==664 and set(frozen)-set(parent)=={'src/epoch_srm_job.rs','tests/epoch_srm_job.rs'}
assert {n for n,h in parent.items() if frozen[n]!=h}=={'src/lib.rs','src/main.rs','docs/EPOCH_SRM_SEARCH.ja.md'}
freeze(out)
call('format-check',[cargo,'fmt','--all','--check'])
call('clippy',[cargo,'clippy','--locked','--all-targets','--','-D','warnings'])
call('focused-job-tests',[cargo,'test','--locked','--release','--test','epoch_srm_job','--','--nocapture'])
focused=(out/'focused-job-tests.log').read_text('utf-8');assert re.search(r'test result: ok\. 10 passed; 0 failed; 0 ignored;',focused)
call('actual-CLI-build',[cargo,'build','--locked','--release','--message-format=json'])
artifacts=[]
for line in (out/'actual-CLI-build.log').read_text('utf-8').splitlines():
    try:item=json.loads(line)
    except json.JSONDecodeError:continue
    if item.get('reason')=='compiler-artifact' and item['target']['name']=='soil-fem-3d' and item.get('executable'):artifacts.append(Path(item['executable']))
assert len(set(artifacts))==1
cli=artifacts[0];local=runtime/'target/release/soil-fem-3d.exe'
assert cli.resolve()==(base/'target/release/soil-fem-3d.exe').resolve() and not local.exists()
local.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(cli,local)
assert sha(cli)==sha(local);check();assert manifest()==frozen
dump(out/'Python-CLI-placement.json',dict(passed=True,source=str(cli),existing_test_location=str(local),binary_sha256=sha(cli),sources_and_tests_unchanged=True))

call('all-Rust-tests',[cargo,'test','--locked','--release'])
log=(out/'all-Rust-tests.log').read_text('utf-8');rows=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',log)
assert rows and sum(int(n) for n,_,_ in rows)==380 and all(int(f)==0 and int(i)==0 for _,f,i in rows)
call('all-Python-tests',[python,'-X','utf8','-W','error::SyntaxWarning','-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
log=(out/'all-Python-tests.log').read_text('utf-8');assert re.search(r'Ran 219 tests in ',log) and re.search(r'(?m)^OK\s*$',log)
check();assert manifest()==frozen
dump(out/'terminal.json',dict(passed=True,source_files=664,Rust_tests=380,Python_tests=219,
    old349_Rust_and9_epoch_and12_control_and10_new_job_tests_retained=True,all_sources_frozen_and_Root_and_live_jobs_unchanged=True,
    existing_actual_CLI_covered_by_integration_test=True,new_job_CLI_and_actual_owned_kill_recovery_covered=True,production_large_job_stop_guard_and_Native_NOT_yet_implemented=True,records=records,
    production_Root_adoption_durable_epoch_search_Native_gate_and_allP6_still_pending=True,p6_complete=False))
print('PASS private664 epoch SRM disk job,380Rust219Python; production adoption/large job stop guard/Native still pending',flush=True)
