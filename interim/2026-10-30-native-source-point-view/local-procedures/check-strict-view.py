"""Reverify the unchanged accepted packet with strict JSON and relative-file admission."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=root/'runs/native-typed-source-integration-20261029-v1/runtime';out=p/'strict-view-verification-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'actual8061-v2/terminal.json')['passed'] and read(p/'public8061-reader-v1/terminal.json')['passed']
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest
sources=manifest();dump(out/'source628-binding.json',sources)
workflow=p/'actual8061-v2/workflow';files=[workflow/'view/view-manifest.json',workflow/'view-verification.json',workflow/'analysis/stage-results.json']
binding={str(q):sha(q) for q in files};dump(out/'original-packet-binding.json',binding)
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONPATH=str(root/'runs/native-typed-source-development-20261029-v1/deps'),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
command=[sys.executable,'-X','utf8',str(runtime/'tools/verify_solidcraft_view.py'),str(workflow),str(workflow/'view'),str(out/'view-verification.json')];start=time.monotonic()
with (out/'verification.log').open('x',encoding='utf-8') as log:
 child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
 dump(out/'worker-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
assert manifest()==sources and binding=={q:sha(Path(q)) for q in binding}
dump(out/'terminal.json',dict(passed=code==0,exit_code=code,source628_unchanged=True,original_packet_and_old_receipt_unchanged=True,
 strict_duplicate_JSON_keys_and_inventory_and_units_admission_used=True,elapsed_seconds=time.monotonic()-start,p6_complete=False))
print(json.dumps(dict(passed=code==0)),flush=True);assert code==0
