"""Full repaired integration in a fresh offline hash-locked CP313 environment."""
from pathlib import Path
import datetime,hashlib,json,os,re,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'python-blocked-integration-v1';runtime=out/'runtime';quality=out/'quality-v2';quality.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(out/'quality-v1/terminal.json')['passed'] is False
assert read(out/'quality-v1/terminal.json')['tests']==218
assert read(out/'wheel-preparation.json')['passed']
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest,freeze
sources=manifest();prior=read(p/'python-quality-v2/source-snapshot.json')
assert len(sources)==633
assert set(sources)-set(prior)=={'tools/streamed_source_v2.py','tests/test_streamed_source_v2.py','requirements-verification-win-cp313.lock.txt','tools/verification_wheels_windows.json','docs/STORAGE_BLOCKS.ja.md'}
assert {name for name,h in prior.items() if sources[name]!=h}=={'tools/export_solidcraft_view.py','tools/verify_solidcraft_view.py','tests/test_streamed_export_source.py','tools/freeze_sources.py','requirements.lock.txt'}
freeze(quality);assert sha(runtime/'target/release/soil-fem-3d.exe')=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
dump(quality/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(quality/'retained-first-failure.json',dict(first_quality_failed=True,failure='existing plate barycentric tuple was changed to list by overly broad Sequence cleaning',
 repair='Clean only explicit stored block sequences; retain every existing tuple API and original assertions',
 original_source630_zip_sha256=sha(out/'quality-v1/source-snapshot.zip'),original_failure_log_sha256=sha(out/'quality-v1/tests.log'),
 numerical_Rust_source_A4_and_original628_unchanged=True))
env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
env.pop('PYTHONPATH',None);env.pop('PYTHONHOME',None)
records=[]
def call(label,command):
 if (p/'pause-requested').exists():raise SystemExit(2)
 command=list(map(str,command));start=time.monotonic()
 with (quality/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(quality/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 assert manifest()==sources
 row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start,sources_unchanged=True);records.append(row);dump(quality/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
 if code:raise SystemExit(code)
cold=out/'clean-env';call('create-env',[sys.executable,'-m','venv',cold]);python=cold/'Scripts/python.exe'
call('offline-install',[python,'-m','pip','install','--disable-pip-version-check','--no-index','--find-links',out/'wheelhouse','--require-hashes','-r',runtime/'requirements-verification-win-cp313.lock.txt'])
call('pip-check',[python,'-m','pip','check'])
call('isolation',[python,'-I',p/'check-blocked-isolation.py',quality/'isolation.json',runtime/'tools/verification_wheels_windows.json'])
assert read(quality/'isolation.json')['passed']
call('tests',[python,'-X','utf8','-W','error::SyntaxWarning','-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
text=(quality/'tests.log').read_text('utf-8');count=re.search(r'Ran (\d+) tests',text);count=int(count[1]) if count else None
assert count==218 and '\nOK\n' in text
base=root/'runs/native-typed-source-integration-20261029-v1/runtime';assert all(sha(base/name)==h for name,h in prior.items())
dump(quality/'terminal.json',dict(passed=True,tests=218,all218_count_observed=True,source_files=633,sources_unchanged=True,
 fresh_same_host_environment=True,all9_wheels_SHA256_pinned_and_installed_offline=True,all9_exact_package_imports_isolated=True,
 prior218_tuple_failure_fixed_without_changing_original_assertions=True,numerical_Rust_source_A4_and_original628_unchanged=True,
 records=records,new_PC_other_GPU_and_actual633_producer_and_formal_adoption_pending=True,p6_complete=False))
print('PASS source633 all218 in new offline hash-locked CP313 env; formal producer/Native/newPC gates pending',flush=True)
