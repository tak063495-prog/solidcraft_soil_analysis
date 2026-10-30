"""Compare every scientific field in the storage-backed full8061 diagnostic."""
from pathlib import Path,PurePosixPath
import datetime,hashlib,json,os,shutil,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1]
out=p/'storage-backed-fields8061-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
 if (p/'pause-requested').exists():raise SystemExit(2)
active()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
 script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
adapter=p/'storage-backed-field-audit-v1'
assert read(adapter/'preparation.json')['prepared']
assert sha(adapter/'audit-fields.py')==read(adapter/'preparation.json')['storage_driver_sha256']
original=root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061'
assert read(original/'study-terminal.json')['passed']
assert read(p/'geometry-audit-real8061-v1/worker-passed.json')['passed']
assert read(original/'field-verification.json')['seed_load_and_recorded_cutback_gates_passed']
source_binding={str(q):sha(q) for q in [original/'workflow/analysis/stage-results.json',
 original/'field-verification.json',original/'preparation.json',original/'terminal.json',
 adapter/'audit-fields.py',adapter/'preparation.json',Path(__file__)]}
dump(out/'protocol.json',dict(prepared=True,elements=8061,
 full_original_scientific_field_report_must_match_exactly=True,
 only_binding_paths_and_driver_identity_excluded_from_report_equality=True,
 immutable_originals_must_remain_unchanged=True,source_binding=source_binding,
 readonly_diagnostic_not_new_numerical_solve=True,p6_complete=False))
while True:
 active()
 memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
  'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
  encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
 if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=28*(1<<30) and shutil.disk_usage(out).free>=32*(1<<30):break
 time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,minimum_each_free_bytes=28*(1<<30),
 free_physical_bytes=memory['FreePhysicalMemory']*1024,free_virtual_bytes=memory['FreeVirtualMemory']*1024,
 free_disk_bytes=shutil.disk_usage(out).free))
runtime=out/'runtime';runtime.mkdir()
quality=p/'Root-source635-promotion-v1/quality-v1'
sources=read(quality/'source-snapshot.json');assert len(sources)==635
with zipfile.ZipFile(quality/'source-snapshot.zip') as z:
 assert z.testzip() is None and set(z.namelist())==set(sources) and len(z.namelist())==635
 for name,h in sources.items():
  q=PurePosixPath(name);assert not q.is_absolute() and '..' not in q.parts
  raw=z.read(name);assert hashlib.sha256(raw).hexdigest()==h
  dest=runtime/name;dest.parent.mkdir(parents=True,exist_ok=True)
  with dest.open('xb') as f:f.write(raw)
for name in ['fem_reference.py','mcc_rate_fem.py']:
 assert sha(runtime/'tools'/name)==sha(original/'runtime/tools'/name),'original numerical reference code differs'
dump(out/'source635.json',sources)
analysis=out/'workflow/analysis';(analysis/'ground').mkdir(parents=True)
for name in ['ground/mesh.json','ground/ground.json','stages.json','stage-results.json']:
 a=original/'workflow/analysis'/name;b=analysis/name
 shutil.copyfile(a,b);assert sha(a)==sha(b)
for name in ['preparation.json','terminal.json']:
 shutil.copyfile(original/name,out/name);assert sha(original/name)==sha(out/name)
# terminal.json above is the legacy diagnostic's required mechanical-package
# input. This manager writes its own outcome under adapter-terminal.json.
shutil.copyfile(adapter/'audit-fields.py',out/'audit-fields.py')
active()
env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',
 OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for key in ('PYTHONPATH','PYTHONHOME'):env.pop(key,None)
command=[sys.executable,'-I','-X','utf8',str(out/'audit-fields.py')];start=time.monotonic()
with (out/'audit.log').open('x',encoding='utf-8') as log:
 child=subprocess.Popen(command,cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
  creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
 dump(out/'audit-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],
  created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 code=child.wait()
dump(out/'audit-process-terminal.json',dict(passed=code==0,exit_code=code,elapsed_seconds=time.monotonic()-start))
assert code==0
actual=read(out/'field-verification.json');expected=read(original/'field-verification.json')
assert actual.pop('sha256') and expected.pop('sha256')
assert actual==expected,'original full seed/load/history/surface diagnostic differs'
assert all(sha(Path(q))==h for q,h in source_binding.items())
assert all(sha(runtime/n)==h for n,h in sources.items())
dump(out/'adapter-terminal.json',dict(passed=True,elements=8061,
 all_original_scientific_seed_load_history_and_surface_report_values_exact=True,
 original_result_and_receipts_and_sources_unchanged=True,
 only_report_binding_paths_and_driver_identity_excluded=True,
 actual_Result_IO_storage_backed=True,all635_runtime_sources_unchanged=True,
 full_original_numerical_tail_byte_identical=True,
 actual_adapter_test_not_a_new_FEM_solution=True,
 general34350_full_solve_and_allP6_still_pending=True,p6_complete=False))
print('PASS storage-backed8061 full seed/load/history/surface diagnostic exact',flush=True)
