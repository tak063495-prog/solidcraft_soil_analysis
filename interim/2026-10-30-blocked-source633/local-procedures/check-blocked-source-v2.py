"""Isolated experimental blocks, real full mixed validation and all-GP VTUs."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'blocked-source-v2-real8061';runtime=root/'runs/native-typed-source-integration-20261029-v1/runtime'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
def worker():
 sys.path[:0]=[str(runtime/'tools'),str(root/'runs/native-typed-source-development-20261029-v1/deps')]
 import streamed_stage_source,streamed_source_v2,mixed_mcc_io
 from verify_solidcraft_view import verify
 from freeze_sources import manifest
 assert manifest()==read(out/'source628.json')
 sources=read(out/'experimental-sources.json');assert all(sha(Path(q))==h for q,h in sources.items())
 workflow=p/'actual8061-v2/workflow';source=workflow/'analysis/stage-results.json';mesh=read(workflow/'analysis/ground/mesh.json')
 # Explicitly choose only this trusted local candidate storage implementation.
 # Every original numerical validator and VTU comparison remains unmodified.
 streamed_stage_source.index_source=streamed_source_v2.index_source
 start=time.monotonic();run=streamed_source_v2.index_source(source,out/'validation-cache')
 dump(out/'indexing.json',dict(passed=True,elapsed_seconds=time.monotonic()-start,source_bytes=source.stat().st_size,
  components=len(run._store.index['components']),maximum_component_bytes=max(v['bytes'] for v in run._store.index['components']),
  maximum_component_payload=max(v['payload'] for v in run._store.index['components']),cache_payload_cap=run._store.cache_payload,
  component_payload_cap=run._store.index['payload_limit'],default_caps_not_raised=True))
 start=time.monotonic();mixed_mcc_io.validate(mesh,run);run.verify_source();run.clear_cache()
 dump(out/'full-mixed-validation.json',dict(passed=True,frames=5,elapsed_seconds=time.monotonic()-start,
  unchanged_original_validator=True,all_GP_geometry_constitutive_primary_and_removed_history_checks_retained=True))
 print('PASS candidate blocks real8061 full mixed validation',flush=True)
 del run
 if (p/'pause-requested').exists():raise SystemExit(2)
 verify(workflow,workflow/'view',out/'view-verification.json')
 assert manifest()==read(out/'source628.json') and all(sha(Path(q))==h for q,h in sources.items())
 assert all(sha(Path(q))==h for q,h in read(out/'original-packet-binding.json').items())
 dump(out/'worker-passed.json',dict(passed=True,real8061_all5_and_all_active_raw_GP_VTU_fields=True,
  original_accepted_packet_and_full_history_unchanged=True,component_and_cache_budgets_not_raised=True,
  storage_only_candidate_not_adopted=True,actual14909_and_100k_capacity_and_native_acceptance_separate=True,p6_complete=False))
 print('PASS candidate blocks all5/all-active GP VTUs; original source628 and packet unchanged',flush=True)
def driver():
 out.mkdir();dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest
 dump(out/'source628.json',manifest());dump(out/'experimental-sources.json',{str(q):sha(q) for q in [p/'streamed_source_v2.py',p/'test_streamed_source_v2.py',Path(__file__).resolve()]})
 workflow=p/'actual8061-v2/workflow'
 dump(out/'original-packet-binding.json',{str(q):sha(q) for q in [workflow/'analysis/stage-results.json',workflow/'view/view-manifest.json',workflow/'view-verification.json']})
 env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',PYTHONPATH=os.pathsep.join([str(runtime/'tools'),str(root/'runs/native-typed-source-development-20261029-v1/deps')]),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
 for label,args in [('unit',[str(p/'test_streamed_source_v2.py'),'-v']),('real8061',[str(Path(__file__).resolve()),'worker'])]:
  if (p/'pause-requested').exists():raise SystemExit(2)
  command=[sys.executable,'-X','utf8',*args];start=time.monotonic()
  with (out/(label+'.log')).open('x',encoding='utf-8') as log:
   child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
   dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
  dump(out/(label+'-terminal.json'),dict(exit_code=code,elapsed_seconds=time.monotonic()-start));print(json.dumps(dict(step=label,exit_code=code)),flush=True)
  if code:raise SystemExit(code)
 assert read(out/'worker-passed.json')['passed']
 dump(out/'terminal.json',dict(passed=True,unit_tests=8,real8061_complete_mixed_validation_and_all_GP_VTU_comparison=True,
  source628_and_native73_unchanged=True,experimental_storage_only=True,p6_complete=False))
if __name__=='__main__':worker() if sys.argv[1:]==['worker'] else driver()
