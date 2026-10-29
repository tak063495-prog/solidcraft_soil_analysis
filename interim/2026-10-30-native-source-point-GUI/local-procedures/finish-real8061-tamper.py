"""Preserve the original harness failure; correct only its error wording and finish case4."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'real8061-negative-controls-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 q.parent.mkdir(parents=True,exist_ok=True)
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
records=[read(out/label/'negative-control-result.json') for label in ['resealed-active-count','duplicate-last-value-wins']]
assert all(v['expected_rejection'] for v in records)
failed=read(out/'descriptor-over8KiB/negative-control-result.json');message=(out/'descriptor-over8KiB/reader.log').read_text('utf-8')
assert failed['exit_code']==1 and not failed['expected_rejection'] and 'file too large' in message
assert (out/'descriptor-over8KiB/workflow/view/points-000.json').stat().st_size>8192
assert not (out/'descriptor-over8KiB/must-not-exist.json').exists()
correction=dict(case='descriptor-over8KiB',expected_rejection=True,exit_code=1,expected_message_fragment='file too large',
 original_harness_only_expected_the_wrong_word_limit=True,original_failed_result_and_log_unchanged=True,
 no_product_code_or_size_budget_changed=True,accepted_analysis_packet=False)
dump(out/'descriptor-over8KiB/harness-message-correction.json',correction);records.append(correction)
source=p/'actual8061-v2/workflow';manifest=read(source/'view/view-manifest.json');binding=read(p/'public8061-reader-v1/compiled-probe-binding.json')
binary=p/'public8061-reader-v1/public-source-point-reader.exe';assert sha(binary)==binding['binary_sha256']
baseline={str(q):sha(q) for q in [source/'view/view-manifest.json',source/'view-verification.json',source/'analysis/stage-results.json',source/'view/points-000.json',binary]}
assert baseline[str(source/'analysis/stage-results.json')]=='6c99fd4b23efe3669ea000c4da3a77c1c6758229e3b1a1aba27f733d42aa0ad0'
dump(out/'continuation-driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
if (p/'pause-requested').exists():raise SystemExit(2)
case=out/'valid-JSON-source-byte-change';workflow=case/'workflow';workflow.mkdir(parents=True);start=time.monotonic()
for name,digest in manifest['job_sources_sha256'].items():
 target=workflow/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/name,target);assert sha(target)==digest
for name,digest in manifest['files_sha256'].items():
 target=workflow/'view'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/'view'/name,target);assert sha(target)==digest
doc=workflow/'source.solidcraft.json';shutil.copyfile(source/'conformal/source.solidcraft.json',doc)
raw=workflow/'analysis/stage-results.json'
with raw.open('r+b') as f:
 header=f.read(4096);index=header.index(b'\n');assert index<header.index(b'"');f.seek(index);f.write(b'\r')
assert sha(raw)!=manifest['job_sources_sha256']['analysis/stage-results.json']
manifest_path=workflow/'view/view-manifest.json';manifest_path.write_text(json.dumps(manifest),'utf-8')
audit=read(source/'view-verification.json');audit['view_manifest_sha256']=sha(manifest_path);audit['scope']='Deliberately forged negative control only; never an accepted result'
dump(workflow/'view-verification.json',audit)
command=[str(binary),str(doc),str(manifest_path),str(p/'public8061-reader-v1/query-plan.json'),str(case/'must-not-exist.json')]
with (case/'reader.log').open('x',encoding='utf-8') as log:
 child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
 dump(case/'reader-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
message=(case/'reader.log').read_text('utf-8');okay=code!=0 and 'source changed' in message and not (case/'must-not-exist.json').exists()
row=dict(case=case.name,expected_rejection=okay,exit_code=code,JSON_whitespace_only_changed=True,
 expected_message_fragment='source changed',elapsed_seconds=time.monotonic()-start,accepted_analysis_packet=False)
dump(case/'negative-control-result.json',row);records.append(row);assert okay
assert baseline=={q:sha(Path(q)) for q in baseline}
dump(out/'terminal.json',dict(passed=True,all4_deliberately_forged_copies_rejected=True,records=records,
 original_accepted_packet_and_source_and_binary_unchanged=True,positive_public_all5_control_preserved=True,
 initial_harness_message_failure_preserved_and_explained=True,p6_complete=False))
print('PASS all four real8061 Native negative controls; original accepted source unchanged',flush=True)
