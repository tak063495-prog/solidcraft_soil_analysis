"""Strengthen staged-job matching; retain guard101's safe refusal as evidence."""
from pathlib import Path
from collections import Counter
import ast,copy,datetime,hashlib,json,subprocess,sys
p=Path(__file__).resolve().parent;root=p.parents[1];j=root/'runs/goal-progress-20261027-v1';out=j/'pause-remaining20-20261030-183721/guard102-verification-v1';out.mkdir()
sys.path.insert(0,str(j));from match_stage_job_process import matches_stage_job_process
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
policyfile=j/'active-pause-policy.json';policy=read(policyfile);old=Path(policy['script'])
assert policy['pause_trigger_observed'] and policy['latest_observed_weekly_usage_percent']==81 and (p/'pause-requested').exists()
assert sha(old)==policy['script_sha256']=='a6a1893a66f7add60a449ceff3c8de0f524b56a2cc24376db7b79103d8ec022c'
dump(out/'retained-guard101-refusal.json',dict(observed_failure='AssertionError: no verified accepted snapshot; do not kill the pending solver',
    source=str(old),source_sha256=sha(old),failure_line=397,forced_held_process_stop_loop_not_reached=True,
    final9600_checkpointed_source_preserved_and_under_separate_deep_verification=True,
    exact_failing_stale_job_not_logged_so_not_claimed_proven=True))
job=p/'selector-fixture/job';exe=job/'solver/soil-fem-3d.exe';status=dict(state='running',process_id=123)
record=dict(ExecutablePath=str(exe),CommandLine=json.dumps([str(exe),'run-stages-job',str(job)]))
cases=[('exact',job,status,{123},{123:record},True)]
def case(tag,**kw):
    row=dict(job=job,status=copy.deepcopy(status),owned={123},records={123:copy.deepcopy(record)},expected=False);row.update(kw)
    cases.append((tag,row['job'],row['status'],row['owned'],row['records'],row['expected']))
case('PID not held',owned=set());case('missing record',records={});case('status not running',status=dict(state='checkpointed',process_id=123))
case('bool PID',status=dict(state='running',process_id=True))
foreign=copy.deepcopy(record);foreign['ExecutablePath']=str(p/'different-job/solver/soil-fem-3d.exe');case('same saved PID but different actual executable',records={123:foreign})
for tag,args in [('different job',[str(exe),'run-stages-job',str(p/'different-job')]),('wrong verb',[str(exe),'run-srm-search-job',str(job)]),
    ('relative job',[str(exe),'run-stages-job','job']),('spoofed argv0',[str(p/'foreign.exe'),'run-stages-job',str(job)]),
    ('extra argument',[str(exe),'run-stages-job',str(job),'ignored'])]:
    altered=copy.deepcopy(record);altered['CommandLine']=json.dumps(args);case(tag,records={123:altered})
altered=copy.deepcopy(record);altered['ExecutablePath']='soil-fem-3d.exe';case('relative image',records={123:altered})
altered=copy.deepcopy(record);altered['CommandLine']='invalid';case('unparseable command',records={123:altered})
rows=[]
for tag,q,s,owned,records,expected in cases:
    actual=matches_stage_job_process(q,s,owned,records,json.loads);assert actual is expected,tag
    rows.append(dict(tag=tag,matched=actual,expected=expected))
assert len(rows)==13
source=old.read_text('utf-8');start=source.index('    def windows_argv(command):');end=source.index('    for manifest_path in (q for scope in mc_scopes',start)
parser=source[start:end];source=source[:start]+source[end:]
anchor="    for job in jobs:\n        status=job/'status.json'\n";assert source.count(anchor)==2
source=source.replace(anchor,'    from match_stage_job_process import matches_stage_job_process\n'+parser+anchor,1)
anchor2="        if value.get('process_id') not in owned:continue\n        solver=job/'solver/soil-fem-3d.exe'"
assert source.count(anchor2)==1
source=source.replace(anchor2,"        if value.get('process_id') not in owned:continue\n        if not matches_stage_job_process(job,value,owned,records,windows_argv):continue\n        solver=job/'solver/soil-fem-3d.exe'")
new=j/'safe-pause-latest102.py';ast.parse(source)
with new.open('x',encoding='utf-8',newline='\n') as f:f.write(source)
results={}
for label,script,args in [('old-inspect',old,[]),('new-inspect',new,[]),('strict-boundary',new,['--execute','--usage-percent','80'])]:
    result=subprocess.run([sys.executable,str(script),*args],capture_output=True,text=True,encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW,timeout=120)
    (out/(label+'.stdout')).write_text(result.stdout,'utf-8');(out/(label+'.stderr')).write_text(result.stderr,'utf-8');results[label]=result.returncode
    assert result.returncode==(2 if label=='strict-boundary' else 0)
before=read(out/'old-inspect.stdout');after=read(out/'new-inspect.stdout');canonical=lambda rows:Counter(json.dumps(r,sort_keys=True) for r in rows)
assert before['inspect_only'] and after['inspect_only'] and before['owned_processes']==after['owned_processes']
assert canonical(before['ownership_seed_receipts'])==canonical(after['ownership_seed_receipts'])
receipt=dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),new_guard=str(new),new_guard_sha256=sha(new),
    exact_stage_selector=str(j/'match_stage_job_process.py'),exact_stage_selector_sha256=sha(j/'match_stage_job_process.py'),
    positive=1,negative=12,cases=rows,read_only_old_new_ownership_identical=True,owned_processes_observed=after['owned_processes'],
    strict80_boundary_rejected=True,all_previous_scope_and_SRM_registrations_retained=True,
    no_force_stop_or_new_scientific_phase_executed_by_this_verifier=True,p6_complete=False)
dump(out/'terminal.json',receipt)
policy.update(updated_utc=receipt['utc'],script=str(new),script_sha256=receipt['new_guard_sha256'],
    exact_stage_selector=str(j/'match_stage_job_process.py'),exact_stage_selector_sha256=receipt['exact_stage_selector_sha256'],
    pause_guard101_refusal_retained=str(out/'retained-guard101-refusal.json'),pause_guard102_verification_receipt=str(out/'terminal.json'))
temp=out/'active-policy.tmp';dump(temp,policy);temp.replace(policyfile);print(json.dumps({k:v for k,v in receipt.items() if k!='cases'}),flush=True)
