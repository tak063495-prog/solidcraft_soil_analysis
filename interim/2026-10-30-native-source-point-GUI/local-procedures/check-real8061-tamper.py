"""Deliberately forged negative controls; none is an accepted analysis packet."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'real8061-negative-controls-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 q.parent.mkdir(parents=True,exist_ok=True)
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'public8061-reader-v1/terminal.json')['passed']
binding=read(p/'public8061-reader-v1/compiled-probe-binding.json');binary=p/'public8061-reader-v1/public-source-point-reader.exe'
assert sha(binary)==binding['binary_sha256']
source=p/'actual8061-v2/workflow';m=read(source/'view/view-manifest.json')
baseline={str(q):sha(q) for q in [source/'view/view-manifest.json',source/'view-verification.json',source/'analysis/stage-results.json',source/'view/points-000.json',binary]}
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
records=[]
for label in ['resealed-active-count','duplicate-last-value-wins','descriptor-over8KiB','valid-JSON-source-byte-change']:
 if (p/'pause-requested').exists():raise SystemExit(2)
 case=out/label;workflow=case/'workflow';workflow.mkdir(parents=True);start=time.monotonic()
 for name,digest in m['job_sources_sha256'].items():
  target=workflow/name;target.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(source/name,target);assert sha(target)==digest
 for name,digest in m['files_sha256'].items():
  target=workflow/'view'/name;target.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(source/'view'/name,target);assert sha(target)==digest
 doc=workflow/'source.solidcraft.json';shutil.copyfile(source/'conformal/source.solidcraft.json',doc)
 manifest=json.loads(json.dumps(m));descriptor=workflow/'view/points-000.json'
 if label=='resealed-active-count':
  value=read(descriptor);value['active_integration_point_count']+=1;descriptor.write_text(json.dumps(value),'utf-8')
 elif label=='duplicate-last-value-wins':
  value=read(descriptor);encoded=json.dumps(value).replace('"stress_unit": "kPa"','"stress_unit": "Pa", "stress_unit": "kPa"')
  assert encoded!=json.dumps(value) and json.loads(encoded)==value;descriptor.write_text(encoded,'utf-8')
 elif label=='descriptor-over8KiB':
  descriptor.write_bytes(descriptor.read_bytes()+b' '*8193)
 else:
  raw=workflow/'analysis/stage-results.json'
  with raw.open('r+b') as f:
   header=f.read(4096);index=header.index(b' ');f.seek(index);f.write(b'\t')
  assert sha(raw)!=manifest['job_sources_sha256']['analysis/stage-results.json']
 manifest['files_sha256']['points-000.json']=sha(descriptor)
 manifest_path=workflow/'view/view-manifest.json';manifest_path.write_text(json.dumps(manifest),'utf-8')
 # Deliberately forge only this test copy's view receipt. This challenges the
 # Native reader's independent checks instead of testing a stale outer hash.
 audit=read(source/'view-verification.json');audit['view_manifest_sha256']=sha(manifest_path)
 audit['scope']='Deliberately forged negative control only; never an accepted result'
 dump(workflow/'view-verification.json',audit)
 command=[str(binary),str(doc),str(manifest_path),str(p/'public8061-reader-v1/query-plan.json'),str(case/'must-not-exist.json')]
 with (case/'reader.log').open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(case/'reader-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 message=(case/'reader.log').read_text('utf-8')
 expected={'resealed-active-count':'source-point descriptor differs','duplicate-last-value-wins':'duplicate field',
  'descriptor-over8KiB':'limit','valid-JSON-source-byte-change':'source changed'}[label]
 okay=code!=0 and expected in message and not (case/'must-not-exist.json').exists()
 row=dict(case=label,expected_rejection=okay,exit_code=code,expected_message_fragment=expected,elapsed_seconds=time.monotonic()-start,
  deliberately_forged_receipt_for_negative_control=True,accepted_analysis_packet=False)
 dump(case/'negative-control-result.json',row);records.append(row);print(json.dumps(row),flush=True);assert okay
assert baseline=={name:sha(Path(name)) for name in baseline}
dump(out/'terminal.json',dict(passed=True,all4_deliberately_forged_copies_rejected=True,records=records,
 original_accepted_packet_and_source_and_binary_unchanged=True,positive_public_all5_control_preserved=True,p6_complete=False))
print('PASS all four real8061 Native negative controls; original accepted workflow unchanged',flush=True)
