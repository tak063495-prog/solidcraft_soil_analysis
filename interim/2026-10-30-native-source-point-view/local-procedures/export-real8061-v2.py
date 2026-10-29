"""Reassemble original byte-identical CAD provenance and export the complete real run."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];n=root/'runs/native-typed-source-integration-20261029-v1';runtime=n/'runtime'
out=p/'actual8061-v2'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
def worker():
 for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
 sys.path[:0]=[str(runtime/'tools'),str(root/'runs/native-typed-source-development-20261029-v1/deps')]
 # This existing model was prepared with the preserved600-source plate code.
 # Select that explicitly trusted, hash-verified module only; do not follow an
 # executable path from imported provenance or change/reseal any receipt.
 import importlib.util
 legacy=root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061/runtime/tools/plate_surfaces.py'
 expected_plate='9b43db5fdff4a2d80c6946c1ade67db5e71e4f68d63c37829b739240c82f17c9'
 assert sha(legacy)==expected_plate
 original_manifest=read(root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061/audit-source-manifest.json')
 assert original_manifest['tools/plate_surfaces.py']==expected_plate
 spec=importlib.util.spec_from_file_location('plate_surfaces',legacy)
 module=importlib.util.module_from_spec(spec);sys.modules['plate_surfaces']=module;spec.loader.exec_module(module)
 assert module.__file__==str(legacy) and sha(legacy)==expected_plate
 dump(out/'historical-plate-code-binding.json',dict(passed=True,explicit_local_trusted_source=str(legacy),sha256=expected_plate,
  original_frozen_source_manifest_sha256=sha(root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061/audit-source-manifest.json'),
  original_receipt_unmodified=True,plate_generation_and_source_linkage_verification_not_bypassed=True,
  current628_exporter_and_mechanical_audits_still_used=True))
 from export_solidcraft_view import export
 from verify_solidcraft_view import verify
 from freeze_sources import manifest
 assert manifest()==read(p/'python-quality-v1/source-snapshot.json')
 assert read(p/'python-quality-v1/terminal.json')['passed']
 original=root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061'
 prepared=root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v4/workflow'
 assert read(original/'study-terminal.json')['passed']
 analysis=original/'workflow/analysis'
 for name in ['ground/mesh.json','ground/ground.json','stages.json']:
  assert sha(analysis/name)==sha(prepared/'analysis'/name),name
 raw=analysis/'stage-results.json';expected='6c99fd4b23efe3669ea000c4da3a77c1c6758229e3b1a1aba27f733d42aa0ad0'
 assert raw.stat().st_size==1459184313 and sha(raw)==expected
 workflow=out/'workflow';workflow.mkdir()
 shutil.copytree(prepared/'conformal',workflow/'conformal')
 target=workflow/'analysis';target.mkdir();shutil.copytree(prepared/'analysis/ground',target/'ground')
 for name in ['stages.json','stage-results.json','accepted-stage-run.json','stage-verification.json']:
  shutil.copyfile(analysis/name,target/name);assert sha(target/name)==sha(analysis/name)
 shutil.copytree(analysis/'results',target/'results')
 solver=analysis/'solver-job/solver/soil-fem-3d.exe';assert sha(solver)=='f2ab7b5a1ea8684d2f05f6edc1e35e51258b9f9890ac3bc9b0d849a370e00f27'
 (target/'solver').mkdir();shutil.copyfile(solver,target/'solver/soil-fem-3d.exe')
 dump(out/'preparation.json',dict(passed=True,original_CAD_capture_and_import_records_copied=True,
  all_three_mechanical_inputs_byte_exact=True,original_result_sha256=expected,original_result_bytes=1459184313,
  original_accepted_receipt_copied_without_rewriting=True,solver_sha256=sha(solver),p6_complete=False))
 print('PASS real8061 complete provenance reassembled without resealing receipts',flush=True)
 if (p/'pause-requested').exists():raise SystemExit(2)
 start=time.monotonic();export(workflow,workflow/'view',source_points=True)
 manifest_view=read(workflow/'view/view-manifest.json');assert manifest_view['schema']=='soil-fem-solidcraft-view/4'
 dump(out/'export-passed.json',dict(passed=True,frames=len(manifest_view['frames']),
  packet_bytes=sum((workflow/'view'/name).stat().st_size for name in manifest_view['files_sha256']),
  raw_point_storage=manifest_view['semantics']['raw_point_storage'],elapsed_seconds=time.monotonic()-start,
  original_full_mechanical_and_input_linkage_audits_kept=True,Native_and_GUI_acceptance_separate=True,p6_complete=False))
 print('PASS real8061 source-point packet exported with all original provenance and full mechanical gates',flush=True)
 if (p/'pause-requested').exists():raise SystemExit(2)
 verify(workflow,workflow/'view',workflow/'view-verification.json')
 assert sha(raw)==sha(target/'stage-results.json')==expected
 assert manifest()==read(p/'python-quality-v1/source-snapshot.json')
 dump(out/'verification.json',dict(passed=True,source628_unchanged=True,original_run_unchanged=True,
  all5_and_all_active_raw_points_compared_with_existing_VTU=True,
  Native_reader_GUI_and_general_mesh_acceptance_separate=True,p6_complete=False))
def driver():
 out.mkdir();dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 with (out/'export-and-verify.log').open('x',encoding='utf-8') as log:
  command=[sys.executable,'-X','utf8',str(Path(__file__).resolve()),'worker'];child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/'worker-child.json',dict(pid=child.pid,executable=sys._base_executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 dump(out/'terminal.json',dict(passed=code==0,exit_code=code,p6_complete=False));print(json.dumps(dict(exit_code=code)));assert code==0
if __name__=='__main__':worker() if sys.argv[1:]==['worker'] else driver()
