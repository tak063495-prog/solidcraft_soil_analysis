"""Check the producer changes with the unchanged accepted A4 FEM executable."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime';out=p/'producer-python-quality-v2';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest,freeze
sources=manifest();assert len(sources)==628,len(sources)
old=read(p.parent/'native-mc-epoch-view-development-20261029-v1/python-quality-v1/source-snapshot.json')
assert {name for name,h in old.items() if sources[name]!=h}=={'tools/export_solidcraft_view.py'}
assert all(sources[name]==h for name,h in old.items() if name.endswith('.rs'))
freeze(out)
solver=root/'target/release/soil-fem-3d.exe';assert sha(solver)=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
(runtime/'target/release').mkdir(parents=True,exist_ok=True);shutil.copyfile(solver,runtime/'target/release/soil-fem-3d.exe')
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONPATH=str(p.parent/'native-typed-source-development-20261029-v1/deps'),
 OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
command=[sys.executable,'-X','utf8','-W','error::SyntaxWarning','-m','unittest','discover','-s','tests','-p','test_*.py','-v'];start=time.monotonic()
if (p/'pause-requested').exists():raise SystemExit(2)
with (out/'tests.log').open('x',encoding='utf-8') as log:
 child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
 dump(out/'tests-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
assert manifest()==sources and sha(runtime/'target/release/soil-fem-3d.exe')==sha(solver)
text=(out/'tests.log').read_text('utf-8')
dump(out/'terminal.json',dict(passed=code==0,exit_code=code,expected_tests=202,
 all202_count_observed='Ran 202 tests' in text,source_files=628,source_unchanged=True,Rust_and_A4_unchanged=True,
 elapsed_seconds=time.monotonic()-start,new_format_actual_export_and_Native_consumer_not_yet_accepted=True,p6_complete=False))
print(json.dumps(dict(exit_code=code,expected202_observed='Ran 202 tests' in text)),flush=True)
assert code==0 and 'Ran 202 tests' in text and '\nOK\n' in text
