"""Real14909 export using integrated source633 and the verified clean dependency environment."""
from pathlib import Path
import datetime,hashlib,importlib.util,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'actual14909-blocked-v2';runtime=p/'python-blocked-integration-v1/runtime'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
def worker():
 sys.path[:0]=[str(runtime/'tools')]
 from freeze_sources import manifest
 assert manifest()==read(out/'source633.json')
 binding=read(out/'experimental-sources.json');assert all(sha(Path(q))==h for q,h in binding.items())
 legacy=root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061/runtime/tools/plate_surfaces.py'
 expected='9b43db5fdff4a2d80c6946c1ade67db5e71e4f68d63c37829b739240c82f17c9'
 assert sha(legacy)==expected
 spec=importlib.util.spec_from_file_location('plate_surfaces',legacy);module=importlib.util.module_from_spec(spec);sys.modules['plate_surfaces']=module;spec.loader.exec_module(module)
 from export_solidcraft_view import export
 from verify_solidcraft_view import verify
 original=p/'mesh14909';assert read(original/'study-terminal.json')['passed']
 prepared=root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v5/workflow'
 assert Path(read(original/'preparation.json')['original_preparation'])==prepared.parent
 for name in ['ground/mesh.json','ground/ground.json','stages.json']:assert sha(original/'workflow/analysis'/name)==sha(prepared/'analysis'/name)
 workflow=out/'workflow';workflow.mkdir();shutil.copytree(prepared/'conformal',workflow/'conformal')
 analysis=workflow/'analysis';analysis.mkdir();shutil.copytree(prepared/'analysis/ground',analysis/'ground')
 for name in ['stages.json','stage-results.json','accepted-stage-run.json','stage-verification.json']:
  shutil.copyfile(original/'workflow/analysis'/name,analysis/name);assert sha(analysis/name)==sha(original/'workflow/analysis'/name)
 shutil.copytree(original/'workflow/analysis/results',analysis/'results')
 solver=original/'solver/soil-fem-3d.exe';assert sha(solver)=='f2ab7b5a1ea8684d2f05f6edc1e35e51258b9f9890ac3bc9b0d849a370e00f27'
 (analysis/'solver').mkdir();shutil.copyfile(solver,analysis/'solver/soil-fem-3d.exe')
 paths=[analysis/name for name in ['ground/mesh.json','ground/ground.json','stages.json','stage-results.json','accepted-stage-run.json','stage-verification.json']]
 originals={str(q):sha(q) for q in paths}
 dump(out/'preparation.json',dict(passed=True,original_three_inputs_and_result_receipts_byte_exact=True,
  original_result_bytes=(analysis/'stage-results.json').stat().st_size,original_result_sha256=sha(analysis/'stage-results.json'),
  historical_plate_module_sha256=expected,original_receipts_unmodified=True,all_original_mechanical_and_CAD_audits_retained=True,sha256=originals,p6_complete=False))
 if (p/'pause-requested').exists():raise SystemExit(2)
 start=time.monotonic();export(workflow,workflow/'view',source_points=True)
 view=read(workflow/'view/view-manifest.json');assert view['schema']=='soil-fem-solidcraft-view/4' and len(view['frames'])==5
 dump(out/'export-passed.json',dict(passed=True,frames=5,elements=14909,GP_count=402543,
  packet_bytes=sum((workflow/'view'/name).stat().st_size for name in view['files_sha256']),
  elapsed_seconds=time.monotonic()-start,full_original_mechanical_and_input_linkage_audits_retained=True,
  block_storage_candidate=True,Native_and_GUI_acceptance_not_yet_claimed=True,p6_complete=False))
 print('PASS real14909 full mechanical source-point export with bounded block storage',flush=True)
 if (p/'pause-requested').exists():raise SystemExit(2)
 verify(workflow,workflow/'view',out/'view-verification.json')
 assert all(sha(Path(q))==h for q,h in originals.items()) and all(sha(Path(q))==h for q,h in binding.items()) and manifest()==read(out/'source633.json')
 dump(out/'worker-passed.json',dict(passed=True,elements=14909,GP_count=402543,all5_all_active_GP_VTU_comparison=True,
  full_original_mechanical_audits_and_CAD_provenance_retained=True,source633_and_original_results_unchanged=True,
  explicit_candidate_block_storage_only=True,Native_GUI_and_general100k_and_spatial_acceptance_separate=True,p6_complete=False))
 print('PASS real14909 all5/all-active raw GP VTUs; Native and general100k acceptance pending',flush=True)
def driver():
 out.mkdir();dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 gates=[p/'mesh14909/study-terminal.json',p/'python-blocked-integration-v1/quality-v2/terminal.json']
 while not all(q.exists() for q in gates):
  if (p/'pause-requested').exists():raise SystemExit(2)
  time.sleep(5)
 if (p/'pause-requested').exists():raise SystemExit(2)
 assert all(read(q)['passed'] for q in gates)
 sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest
 assert manifest()==read(p/'python-blocked-integration-v1/quality-v2/source-snapshot.json') and len(manifest())==633
 assert Path(sys.prefix).resolve()==(p/'python-blocked-integration-v1/clean-env').resolve()
 dump(out/'source633.json',manifest());dump(out/'experimental-sources.json',{str(q):sha(q) for q in [runtime/'tools/streamed_source_v2.py',Path(__file__).resolve()]})
 env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
 env.pop('PYTHONPATH',None);env.pop('PYTHONHOME',None)
 command=[sys.executable,'-I','-X','utf8',str(Path(__file__).resolve()),'worker'];start=time.monotonic()
 with (out/'export-and-verify.log').open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/'worker-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 dump(out/'terminal.json',dict(passed=code==0,exit_code=code,elapsed_seconds=time.monotonic()-start,p6_complete=False))
 print(json.dumps(dict(exit_code=code)),flush=True)
 if code:raise SystemExit(code)
 assert read(out/'worker-passed.json')['passed']
if __name__=='__main__':worker() if sys.argv[1:]==['worker'] else driver()
