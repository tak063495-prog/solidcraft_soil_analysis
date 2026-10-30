"""Use closed real CLI fault-test jobs, without changing inputs or running a solve."""
from pathlib import Path
import datetime,hashlib,json,os,struct,sys
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'epoch-SRM-preserver-v1';out.mkdir()
sys.path.insert(0,str(root/'runs/goal-progress-20261027-v1'))
from preserve_epoch_srm_search_job import preserve_epoch_search_job,read,sha
dump=lambda q,v:q.write_text(json.dumps(v,indent=2,allow_nan=False),encoding='utf-8')
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
base=p/'epoch-SRM-job-v1';quality=base/'quality-v2';frozen=read(quality/'source-snapshot.json')
assert read(quality/'focused-job-tests-terminal.json')['exit_code']==0
assert all(sha(base/'runtime'/n)==h for n,h in frozen.items())
fixtures=base/'runtime/runs'
def job(tag,sub='first'):
    matches=list(fixtures.glob('epoch SRM job '+tag+' 17136-*'));assert len(matches)==1
    q=matches[0]/sub;assert (q/'job.json').is_file();return q
def inventory(q):return {n.relative_to(q).as_posix():sha(n) for n in q.rglob('*') if n.is_file()}
def signature(v):
    if type(v) is float:return ('f64',struct.pack('>d',v).hex())
    if isinstance(v,dict):return ('object',tuple((k,signature(x)) for k,x in sorted(v.items())))
    if isinstance(v,list):return ('array',tuple(map(signature,v)))
    return (type(v).__name__,v)
rows=[]
def preserve(q,dest,**overrides):
    manifest=read(q/'job.json');solver=q/'solver/soil-fem-3d.exe'
    kwargs=dict(expected_solver_sha256=sha(solver),expected_job_id=manifest['job_id'],expected_manifest_sha256=sha(q/'job.json'))
    kwargs.update(overrides)
    return preserve_epoch_search_job(q,dest,solver,**kwargs)
for tag,sub,epoch,complete in [('early stop','first',False,False),('roundtrip','first',True,False),
    ('nonlinear cooperative stop','first',True,False),('actual process kill','first',True,False),
    ('finished fork','first',True,True),('unknown','first',True,True),('roundtrip','recovered',True,True)]:
    q=job(tag,sub);before=inventory(q);pointer=read(q/'recovery.json');original=read(q/pointer['checkpoint'],32*1024*1024)
    receipt=preserve(q,out/f'positive-{len(rows)}');child=Path(receipt['frozen_job'])
    copied=read(child/'inputs/checkpoint.json',32*1024*1024)
    assert signature(original['checkpoint'])==signature(copied),'all fields and float bits must be identical'
    assert receipt['saved_material_epoch']==epoch and receipt['search_complete']==complete
    assert receipt['selected_generation']==pointer['generation']
    assert receipt['selected_source_snapshot_sha256']==pointer['checkpoint_sha256']
    assert receipt['selected_canonical_snapshot_sha256']==pointer['payload_sha256']
    assert before==inventory(q),'preserver altered original job'
    rows.append(dict(kind='positive',tag=tag,sub=sub,passed=True,full_checkpoint_fields_and_f64_bits_equal=True,receipt=receipt))
    print(f'PASS preserve {tag}/{sub}',flush=True)
normal=job('roundtrip');negative=[('foreign selected pointer',job('corruption','foreign'),{}),
    ('truncated selected slot',job('corruption','truncated'),{}),('resealed semantic change',job('corruption','resealed'),{}),
    ('wrong registered manifest',normal,dict(expected_manifest_sha256='0'*64)),
    ('wrong registered identity',normal,dict(expected_job_id='0'*64)),
    ('wrong registered binary',normal,dict(expected_solver_sha256='0'*64))]
for tag,q,overrides in negative:
    dest=out/f'negative-{len(rows)}';before=inventory(q)
    try:preserve(q,dest,**overrides)
    except (ValueError,OSError) as error:reason=str(error)
    else:raise AssertionError('invalid source accepted: '+tag)
    assert not (dest/'receipt.json').exists() and inventory(q)==before
    rows.append(dict(kind='negative',tag=tag,passed=True,rejected=True,reason=reason,original_job_unchanged=True))
    print('PASS reject '+tag,flush=True)
assert all(sha(base/'runtime'/n)==h for n,h in frozen.items())
baseline=read(base/'Root657-source-snapshot.json');assert all(sha(root/n)==h for n,h in baseline.items())
result=dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),positive=7,negative=6,
    selected_complete_checkpoint_fields_and_float_bits_identical=True,source_jobs_and_Root657_and_immutable664_sources_unchanged=True,
    no_scientific_solver_or_process_signal_started=True,actual_job_safety_guard_registration_still_pending=True,
    source_sha256=sha(root/'runs/goal-progress-20261027-v1/preserve_epoch_srm_search_job.py'),
    driver_sha256=sha(Path(__file__)),rows=rows,p6_complete=False)
dump(out/'terminal.json',result);print('PASS epoch SRM preserver:7 real checkpoints,6 fail-closed counterexamples',flush=True)
