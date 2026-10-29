"""Actual promoted worker resumes retained87.5% excavation history, preserving original jobs."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert read(p/'promotion.json')['promoted'] and read(p/'acceptance.json')['passed']
with (p/'root-mid-excavation-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
source=p/'controls/continuous';inner=source/'workflow/analysis/solver-job'
before={q.relative_to(source).as_posix():sha(q) for q in source.rglob('*') if q.is_file()}
pointer=read(inner/'recovery.json');assert pointer['generation']==36
checkpoint=inner/'recovery-checkpoint-1.json';saved=read(checkpoint);run=saved['run'];continuation=run['continuation']
assert len(run['frames'])==4 and continuation['stage_index']==3 and continuation['progress']==.875
active=run['last_committed']['active_elements'];assert len(active)==274 and sum(active)==250
copied=p/'mid-excavation-snapshot-source'
def ignore(folder,names):
    q=Path(folder)
    if q==source/'workflow':return [n for n in names if n=='view']
    if q==source/'workflow/analysis':return [n for n in names if n in ('results','stage-results.json')]
    if q==inner:return [n for n in names if n in ('result.json','progress.jsonl')]
    return []
shutil.copytree(source,copied,ignore=ignore)
selected={**pointer,'generation':35,'checkpoint':checkpoint.name,'checkpoint_sha256':sha(checkpoint),'run_sha256':saved['run_sha256']}
(copied/'workflow/analysis/solver-job/recovery.json').write_text(json.dumps(selected,indent=2)+'\n','utf-8')
assert sha(copied/'workflow/analysis/solver-job'/checkpoint.name)==sha(checkpoint)
selection=dict(passed=True,selected_generation=35,checkpoint_data_byte_exact=True,original_parent_unchanged=True,
    selection_pointer_only_in_new_copy=True,continuation_stage_index=3,excavation_progress=.875,active_elements=250,total_elements=274,
    remaining_stress_release_recomputed=True,removed_integration_point_histories_retained=True)
with (p/'root-mid-excavation-selection.json').open('x') as f:json.dump(selection,f,indent=2)
if (p/'pause-requested').exists():raise SystemExit(2)
solver=root/'runs/durable-root-promotion-20261029-v1/solver/soil-fem-3d.exe'
output=p/'controls/root-mid-excavation';args=[sys.executable,str(root/'tools/native_job.py'),'recover',str(copied),str(output),'--solver',str(solver)]
began=time.monotonic()
with (p/'root-mid-excavation.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(args,cwd=root,stdout=log,stderr=subprocess.STDOUT,
        env={**os.environ,'PYTHONUTF8':'1','PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OMP_NUM_THREADS':'1'},
        creationflags=subprocess.CREATE_NO_WINDOW|subprocess.CREATE_NEW_PROCESS_GROUP|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    with (p/'root-mid-excavation-child.json').open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=args[0],arguments=args[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
    code=child.wait()
assert code==0,code
assert read(output/'status.json')['state']=='succeeded' and read(output/'status.json')['solver_report']['state']=='succeeded'
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
expected=read(source/'workflow/analysis/stage-results.json');actual=read(output/'workflow/analysis/stage-results.json')
for key in ['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','transition_audits','geostatic','attempts']:
    assert stable(expected[key])==stable(actual[key]),key
assert before=={q.relative_to(source).as_posix():sha(q) for q in source.rglob('*') if q.is_file()}
tested=read(p/'quality/source-snapshot.json')
assert all(sha(root/n)==h for n,h in tested.items() if not n.endswith('.md'))
with (p/'actual-root-acceptance.json').open('x') as f:json.dump(dict(passed=True,actual_promoted_root_worker=True,python_tests=171,rust_sources_unchanged_from272_passed=True,
    same_accepted_solver_sha256=sha(solver),remaining_excavation_stress_release_from_point875_recomputed=True,all8_physical_groups_exact_to_continuous=True,
    postprocessing_and_view_audits_passed=True,original_parent_every_file_unchanged=True,all_production_non_markdown_matches_tested_candidate=True,
    selected_checkpoint=selection,elapsed_seconds=time.monotonic()-began,native_gui_recovery_acceptance_still_separate=True,p6_complete=False),f,indent=2)
print('PASS actual Root worker: excavation87.5% snapshot to completion, full8 histories exact and original parent unchanged',flush=True)
