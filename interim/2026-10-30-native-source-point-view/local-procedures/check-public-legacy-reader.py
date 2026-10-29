"""Keep the real existing535MB public reader regression under the new consumer source."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];native=root/'runs/native-typed-source-integration-20261029-v1/native';out=p/'public-legacy-reader-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'native-quality-v2/terminal.json')['passed']
sources=read(p/'native-quality-v2/source73.json')['files_sha256']
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
def check():
 assert all(sha(native/name)==h for name,h in sources.items())
 if (p/'pause-requested').exists():raise SystemExit(2)
env=dict(os.environ,CARGO_TARGET_DIR=str(native/'target'),CARGO_BUILD_JOBS='2')
def call(label,command):
 check();start=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=native,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=list(map(str,command[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 check();row=dict(exit_code=code,elapsed_seconds=time.monotonic()-start);dump(out/(label+'-terminal.json'),row);print(json.dumps(dict(step=label,**row)),flush=True);assert code==0
call('engine-build',['C:/Users/link_/.cargo/bin/cargo.exe','build','--locked','--profile','ci','-p','solidcraft-engine','--lib','--message-format=json-render-diagnostics'])
artifacts=[]
for line in (out/'engine-build.log').read_text('utf-8').splitlines():
 try:item=json.loads(line)
 except json.JSONDecodeError:continue
 if item.get('reason')=='compiler-artifact':artifacts.append(item)
def artifact(name,feature=None):
 paths={f for item in artifacts if item['target']['name']==name and (feature is None or feature in item['features']) for f in item['filenames'] if f.endswith('.rlib')}
 assert len(paths)==1,(name,paths)
 return Path(paths.pop())
engine=artifact('solidcraft_engine');serde=artifact('serde_json','float_roundtrip')
source=root/'runs/native-mc-epoch-view-development-20261029-v1/ordinary-reader-probe.rs'
probe=out/'public-reader-probe.rs';probe.write_text(source.read_text('utf-8').replace('use serde_json::{Value,json};','use serde_json::json;'),'utf-8')
binary=out/'public-reader.exe';bindings={str(q):sha(q) for q in [engine,serde,probe]};dump(out/'library-bindings.json',bindings)
call('compile-probe',['C:/Users/link_/.cargo/bin/rustc.exe','--edition','2024','-L','dependency='+str(serde.parent),'--extern','solidcraft_engine='+str(engine),'--extern','serde_json='+str(serde),str(probe),'-o',str(binary)])
case=root/'runs/mcc-large-structures-development-20261028-v1/large-case'
document=case/'source.solidcraft.json';manifest=case/'workflow/view/view-manifest.json';raw=case/'workflow/analysis/stage-results.json'
original={str(q):sha(q) for q in [document,manifest,raw,binary]};assert raw.stat().st_size==535711781
call('actual-reader',[str(binary),str(document),str(manifest),str(out/'reader-receipt.json')])
receipt=read(out/'reader-receipt.json');run=read(raw);checks=0
assert len(receipt['frames'])==len(run['frames'])==5
for frame,actual in zip(run['frames'],receipt['frames'],strict=True):
 assert actual['configured']['frame_name']==frame['name'] and actual['configured']['triangles']>0
 assert actual['configured']['source_storage']['representation']=='typed_projected_MCC_selected_frame'
 assert actual['configured']['raw_point_storage']=='separate_point_file'
 for entry in actual['points']:
  expected=frame['state']['material_states'][entry['id']]
  for key in ['effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio','preconsolidation_pressure_kpa']:
   assert entry['value'][key]==expected[key];checks+=1
 for entry in actual['nodes']:
  assert entry['value']['stored_node_fields']['total_displacement_from_seed_m']==frame['state']['displacement_m'][entry['id']];checks+=1
assert original=={q:sha(Path(q)) for q in original} and bindings=={q:sha(Path(q)) for q in bindings}
dump(out/'terminal.json',dict(passed=True,source_files=73,frames=5,representative_raw_and_node_checks=checks,
 original535MB_source_CAD_and_manifest_unchanged=True,actual_public_reader_used=True,
 old_view_packet_compatibility_confirmed=True,new8061_version4_and_GUI_acceptance_still_separate=True,p6_complete=False))
print('PASS real535MB legacy packet under source73; new8061 packet and GUI remain separate',flush=True)
