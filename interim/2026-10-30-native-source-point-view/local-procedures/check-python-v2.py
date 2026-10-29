from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];n=root/'runs/native-typed-source-integration-20261029-v1';runtime=n/'runtime'
out=p/'python-quality-v2';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest,freeze
sources=manifest();prior=read(n/'producer-python-quality-v2/source-snapshot.json')
assert len(sources)==628 and set(sources)==set(prior)
assert {name for name,h in prior.items() if sources[name]!=h}=={'tools/export_solidcraft_view.py','tools/verify_solidcraft_view.py','tests/test_source_point_descriptor.py','tools/source_point_descriptor.py'}
freeze(out);solver=runtime/'target/release/soil-fem-3d.exe'
assert sha(solver)=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONPATH=str(root/'runs/native-typed-source-development-20261029-v1/deps'),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
command=[sys.executable,'-X','utf8','-W','error::SyntaxWarning','-m','unittest','discover','-s','tests','-p','test_*.py','-v'];start=time.monotonic()
with (out/'tests.log').open('x',encoding='utf-8') as log:
 child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
 dump(out/'tests-child.json',dict(pid=child.pid,executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
assert manifest()==sources
text=(out/'tests.log').read_text('utf-8');passed=code==0 and 'Ran 208 tests' in text and '\nOK\n' in text
dump(out/'terminal.json',dict(passed=passed,exit_code=code,tests=208,all208_count_observed='Ran 208 tests' in text,
 source628_unchanged=True,Rust_and_A4_unchanged=True,elapsed_seconds=time.monotonic()-start,strict_view_reverification_pending=True,p6_complete=False))
print(json.dumps(dict(passed=passed,tests=208)));assert passed
