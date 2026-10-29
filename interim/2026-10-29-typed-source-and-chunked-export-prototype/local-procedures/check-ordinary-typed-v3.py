"""Exercise the actual public Native reader with the unchanged535MB MCC model."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];native=p/'native';out=p/'ordinary-typed-large-v3';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
quality=p/'quality-v5'
while not (quality/'terminal.json').exists():
 if (p/'pause-requested').exists():raise SystemExit(2)
 for q in quality.glob('*-terminal.json'):
  if read(q)['exit_code']!=0:raise SystemExit('Native quality failed: '+q.name)
 time.sleep(3)
assert read(quality/'terminal.json')['passed']
sources=read(quality/'source72.json')['files_sha256']
def check():
 assert all(sha(native/name)==h for name,h in sources.items())
 if (p/'pause-requested').exists():raise SystemExit(2)
probe=p.parent/'native-mc-epoch-view-development-20261029-v1/ordinary-reader-probe.rs'
frozen=out/'ordinary-reader-probe.rs';shutil.copyfile(probe,frozen);assert sha(probe)==sha(frozen)
env=dict(os.environ,CARGO_TARGET_DIR=str(native/'target'),CARGO_BUILD_JOBS='2')
def call(label,command):
 check();start=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=native,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=list(map(str,command[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  code=child.wait()
 check();dump(out/(label+'-terminal.json'),dict(exit_code=code,elapsed_seconds=time.monotonic()-start));assert code==0
check()
prior_build=p/'ordinary-typed-large-v2'
assert read(prior_build/'engine-build-terminal.json')['exit_code']==0
shutil.copyfile(prior_build/'engine-build.log',out/'engine-build.log')
assert all(sha(Path(q))==h for q,h in read(prior_build/'library-bindings.json').items())
artifacts=[]
for line in (out/'engine-build.log').read_text('utf-8').splitlines():
 try:r=json.loads(line)
 except json.JSONDecodeError:continue
 if r.get('reason')=='compiler-artifact':artifacts.append(r)
def artifact(name,required_feature=None):
 rows=[r for r in artifacts if r['target']['name']==name and (required_feature is None or required_feature in r['features'])]
 paths={f for r in rows for f in r['filenames'] if f.endswith('.rlib')};assert len(paths)==1,(name,paths)
 return Path(paths.pop())
engine=artifact('solidcraft_engine');serde=artifact('serde_json','float_roundtrip')
binary=out/'ordinary-reader.exe';bindings={str(q):sha(q) for q in [engine,serde,frozen]}
dump(out/'library-bindings.json',bindings)
call('compile-probe',['C:/Users/link_/.cargo/bin/rustc.exe','--edition','2024','-L','dependency='+str(serde.parent),
 '--extern','solidcraft_engine='+str(engine),'--extern','serde_json='+str(serde),str(frozen),'-o',str(binary)])
case=root/'runs/mcc-large-structures-development-20261028-v1/large-case'
document=case/'source.solidcraft.json';manifest=case/'workflow/view/view-manifest.json';raw=case/'workflow/analysis/stage-results.json'
original={str(q):sha(q) for q in [document,manifest,raw,binary]}
prior=read(p.parent/'native-mc-epoch-view-development-20261029-v1/ordinary-large-v4/terminal.json')
assert prior['passed'] and raw.stat().st_size==prior['raw_source_bytes']==535711781
call('actual-reader',[str(binary),str(document),str(manifest),str(out/'reader-receipt.json')])
receipt=read(out/'reader-receipt.json');source=read(raw);checks=0
assert len(receipt['frames'])==len(source['frames'])==5
for frame,actual in zip(source['frames'],receipt['frames'],strict=True):
 assert actual['configured']['frame_name']==frame['name'] and actual['configured']['triangles']>0
 storage=actual['configured']['source_storage'];assert storage['representation']=='typed_projected_MCC_selected_frame'
 assert storage['material_history_count']==storage['integration_layout_count']==75789
 assert storage['conservative_retained_accounted_bytes']<=512<<20
 assert 'mc_epoch' not in actual['configured']
 for entry in actual['points']:
  expected=frame['state']['material_states'][entry['id']];value=entry['value']
  for key in ['effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio','preconsolidation_pressure_kpa']:
   assert value[key]==expected[key],(frame['name'],entry['id'],key);checks+=1
 for entry in actual['nodes']:
  assert entry['value']['stored_node_fields']['total_displacement_from_seed_m']==frame['state']['displacement_m'][entry['id']];checks+=1
assert original=={q:sha(Path(q)) for q in original}
assert bindings=={q:sha(Path(q)) for q in bindings}
dump(out/'terminal.json',dict(passed=True,source_files=72,frames=5,original_raw_bytes=raw.stat().st_size,
 actual_public_reader_and_original_CAD_input_identity_used=True,typed_representation_confirmed_for_every_frame=True,
 all_active_cell_and_projected_source_checks_applied=True,representative_raw_value_checks=checks,
 original_raw_CAD_manifest_unchanged=True,library_bindings_unchanged=True,reader_binary_sha256=sha(binary),
 additional8061_packet_and_GUI_and_full_CI_still_pending=True,reference_verifier_Python_DOM_memory_not_native_reader_memory=True,p6_complete=False))
print('PASS actual535MB all5 through typed Native public reader; raw fields and CAD identity verified',flush=True)
