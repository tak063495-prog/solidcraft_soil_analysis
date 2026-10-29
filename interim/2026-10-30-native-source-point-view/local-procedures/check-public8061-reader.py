"""Read real version4 data through Native and compare direct values with independent VTUs."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
import meshio,numpy as np
p=Path(__file__).resolve().parent;root=p.parents[1];native=root/'runs/native-typed-source-integration-20261029-v1/native';out=p/'public8061-reader-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'native-quality-v2/terminal.json')['passed']
assert read(p/'public-legacy-reader-v1/terminal.json')['passed']
out.mkdir();dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sources=read(p/'native-quality-v2/source73.json')['files_sha256']
def check():
 assert all(sha(native/name)==h for name,h in sources.items())
 if (p/'pause-requested').exists():raise SystemExit(2)
bindings=read(p/'public-legacy-reader-v1/library-bindings.json');assert all(sha(Path(q))==h for q,h in bindings.items())
engine=next(Path(q) for q in bindings if Path(q).name=='libsolidcraft_engine.rlib')
serde=next(Path(q) for q in bindings if Path(q).name.startswith('libserde_json') and q.endswith('.rlib'))
probe=p/'source-point-reader-probe.rs';frozen=out/probe.name;frozen.write_bytes(probe.read_bytes());binary=out/'public-source-point-reader.exe'
def call(label,command):
 check();start=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=native,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=list(map(str,command[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 check();row=dict(exit_code=code,elapsed_seconds=time.monotonic()-start);dump(out/(label+'-terminal.json'),row);print(json.dumps(dict(step=label,**row)),flush=True);assert code==0
call('compile-probe',['C:/Users/link_/.cargo/bin/rustc.exe','--edition','2024','-L','dependency='+str(serde.parent),'--extern','solidcraft_engine='+str(engine),'--extern','serde_json='+str(serde),str(frozen),'-o',str(binary)])
dump(out/'compiled-probe-binding.json',dict(passed=True,engine=str(engine),engine_sha256=sha(engine),serde=str(serde),serde_sha256=sha(serde),
 probe_sha256=sha(frozen),binary_sha256=sha(binary),source73_manifest_sha256=sha(p/'native-quality-v2/source73.json')))
if sys.argv[1:]==['compile-only']:
 print('PASS public version4 probe compiled and bound; actual data phase must be invoked separately',flush=True);raise SystemExit(0)
raise SystemExit('Run compiled probe through the separate read-real8061.py data-phase driver once the existing export finishes')
