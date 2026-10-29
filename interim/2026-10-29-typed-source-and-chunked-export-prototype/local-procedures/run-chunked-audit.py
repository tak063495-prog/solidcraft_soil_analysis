"""Run the unchanged full projected mechanical audit through the chunked source."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'chunked-mechanical-audit-v1'
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)

def worker():
 for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
 sys.path[:0]=[str(p/'runtime/tools'),str(p.parent/'native-typed-source-development-20261029-v1/deps')]
 from streamed_stage_source import index_source
 from verify_projected_stage_package import audit
 import export_solidcraft_view as view
 source_case=p.parent/'mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061'
 analysis=source_case/'workflow/analysis';raw=analysis/'stage-results.json'
 files=[raw,analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stages.json']
 binding={str(q):sha(q) for q in files}
 expected=next(v for k,v in read(source_case/'field-verification.json')['sha256'].items() if k.endswith('analysis\\stage-results.json'))
 assert binding[str(raw)]==expected
 start=time.monotonic();run=index_source(raw,out/'components')
 print(json.dumps(dict(step='fresh_index_complete',elapsed_seconds=time.monotonic()-start)),flush=True)
 mesh,ground,case=map(read,files[1:])
 begin=time.monotonic();audit(mesh,ground,case,run)
 run.verify_source()
 dump(out/'mechanical-audit-passed.json',dict(passed=True,all_original_projected_audit_gates_retained=True,
    frames=len(run['frames']),elapsed_seconds=time.monotonic()-begin,
    source_sha256=expected,independent_nonlinear_solver_comparison_not_added=True,
    scope='Existing independent geometry/load/equilibrium/projected equation/history audit through the chunked source',
    accepted_Native_packet_not_yet_created=True,p6_complete=False))
 print('PASS unchanged full projected audit through chunked source',flush=True)
 run.clear_cache();sizes=[]
 def encoded_size(value):
  return sum(len(part.encode('utf-8')) for part in json.JSONEncoder(ensure_ascii=False,allow_nan=False,separators=(',',':')).iterencode(value))
 geometry=dict(schema='soil-fem-view-geometry/1',length_unit='m',axis='Z-up',nodes_m=mesh['nodes'],
    cells=mesh.get('tetrahedra10') or mesh['tetrahedra'],cell_kind='tetra10' if mesh.get('tetrahedra10') else 'tetra4',
    regions=mesh['regions'],node_sets=mesh.get('node_sets',{}),crs=mesh.get('crs'),origin_m=mesh.get('origin_m',[0.,0.,0.]))
 geometry_bytes=encoded_size(geometry)
 for i in range(len(run['frames'])):
  if (p/'pause-requested').exists():raise SystemExit(2)
  frame,points=view._scene_frame(mesh,run,i,validated_mixed=True)
  row=dict(frame=i,frame_bytes=encoded_size(frame),points_bytes=encoded_size(points));sizes.append(row)
  print(json.dumps(row),flush=True);del frame,points;run.clear_cache()
 total=geometry_bytes+sum(r['frame_bytes']+r['points_bytes'] for r in sizes)
 run.verify_source();assert all(sha(Path(q))==h for q,h in binding.items())
 dump(out/'packet-size-diagnostic.json',dict(passed=True,actual_all5_encoded_size_measured=True,
    geometry_bytes=geometry_bytes,frames=sizes,total_bytes=total,
    existing_total_budget_bytes=view.MAX_TOTAL,existing_file_budget_bytes=view.MAX_FILE,
    fits_existing_total_budget=total<=view.MAX_TOTAL,
    fits_existing_file_budget=max([geometry_bytes,*[r[k] for r in sizes for k in ['frame_bytes','points_bytes']]])<=view.MAX_FILE,
    no_accepted_manifest_emitted=True,whole_source_immutable=True,elapsed_seconds=time.monotonic()-start,p6_complete=False))
 print('PASS size diagnostic; this is not an accepted Native export',flush=True)

def driver():
 out.mkdir();dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
  script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 prior=p/'producer-source-index-v2'
 while not (prior/'terminal.json').exists():
  if (p/'pause-requested').exists():raise SystemExit(2)
  time.sleep(5)
 assert read(prior/'terminal.json')['passed'] and read(prior/'verification.json')['passed']
 assert all(sha(Path(q))==h for q,h in read(prior/'source-binding.json').items())
 original=read(p.parent/'native-mc-epoch-view-development-20261029-v1/python-quality-v1/source-snapshot.json')
 names=set(original)|{'tools/streamed_stage_source.py','tests/test_streamed_stage_source.py'}
 binding={name:sha(p/'runtime'/name) for name in sorted(names)}
 dump(out/'runtime625-source-binding.json',binding)
 env=dict(os.environ,PYTHONPATH=os.pathsep.join([str(root/'.venv/Lib/site-packages'),str(p/'runtime/tools'),str(p.parent/'native-typed-source-development-20261029-v1/deps')]),
    OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
 command=[sys._base_executable,'-X','utf8',str(Path(__file__).resolve()),'--worker'];start=time.monotonic()
 if (p/'pause-requested').exists():raise SystemExit(2)
 with (out/'audit.log').open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/'worker-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  code=child.wait()
 assert all(sha(p/'runtime'/name)==h for name,h in binding.items())
 dump(out/'terminal.json',dict(exit_code=code,passed=code==0,elapsed_seconds=time.monotonic()-start,
    all_runtime_sources_unchanged=True,no_Native_export_acceptance_claimed=True,p6_complete=False))
 print(json.dumps(dict(exit_code=code,elapsed_seconds=time.monotonic()-start)),flush=True)
 raise SystemExit(code)

if '--worker' in sys.argv:worker()
else:driver()
