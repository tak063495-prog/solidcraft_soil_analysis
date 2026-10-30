"""Link the same independent diagnostic against the frozen old/new Rust libs."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;out=p/'mcc-local-cut-stop-verification-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
source=p/'mcc-local-cut-stop-probe.rs';binding={str(source):sha(source),str(Path(__file__)):sha(Path(__file__))}
records=[]
def call(label,args,cwd):
    active();start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=cwd,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(args[0]),
            arguments=list(map(str,args[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(label=label,exit_code=code,elapsed_seconds=time.monotonic()-start);records.append(row);dump(out/(label+'-terminal.json'),row)
    print(json.dumps(row),flush=True)
    if code:
        dump(out/'terminal.json',dict(passed=False,failed_phase=label,records=records,p6_complete=False));raise SystemExit(1)
for version,name,quality in [(653,'srm-controlled-search-v1','quality-v4'),(656,'srm-durable-job-v1','quality-v3')]:
    base=p/name;runtime=base/'runtime';manifest=read(base/quality/'source-snapshot.json')
    assert len(manifest)==version and all(sha(runtime/n)==h for n,h in manifest.items())
    artifacts=[]
    for line in (base/quality/'CLI-build.log').read_text('utf-8').splitlines():
        try:item=json.loads(line)
        except json.JSONDecodeError:continue
        if item.get('reason')=='compiler-artifact':artifacts.append(item)
    libs={}
    for target in ['soil_fem_3d','serde_json']:
        values={Path(q) for item in artifacts if item['target']['name']==target for q in item['filenames'] if q.endswith('.rlib')}
        assert len(values)==1,(target,values)
        libs[target]=values.pop();binding[str(libs[target])]=sha(libs[target])
    exe=out/f'probe-{version}.exe'
    call(f'link-{version}',['C:/Users/link_/.cargo/bin/rustc.exe','--edition=2024','-C','opt-level=2',
        '-L',f'dependency={base/"target/release/deps"}',
        '--extern',f'soil_fem_3d={libs["soil_fem_3d"]}','--extern',f'serde_json={libs["serde_json"]}',source,'-o',exe],runtime)
    binding[str(exe)]=sha(exe)
    call(f'actual-MCC-{version}',[exe,runtime,out/f'results-{version}','fixed' if version==656 else 'old'],runtime)
    assert all(sha(runtime/n)==h for n,h in manifest.items())
    assert all(sha(Path(n))==h for n,h in binding.items())
dump(out/'source-binding.json',binding)
dump(out/'terminal.json',dict(passed=True,actual_old_failure_reproduced_and_fixed_candidate_passed=True,
    formulations=['mixed-MCC-Tet10','projected-MCC-Tet10'],boundaries=['partial-local-cut-early-stop','partial-local-cut-late-stop'],
    all_4_new_resume_histories_match_uninterrupted_except_elapsed_seconds=True,
    no_frozen_runtime_source_or_existing_job_changed=True,records=records,p6_complete=False))
