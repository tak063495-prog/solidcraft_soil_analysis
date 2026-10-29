"""Fresh-process period-boundary restarts of the exact audited82944 Tet4 case."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'rain-refined-restart-v2'
original=root/'runs/richards-local-rain-3d-refinement-20261029-v1';case_source=original/'rt0-24-dt100'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def pause_check():
    if (p/'pause-requested').exists():raise SystemExit(2)
out.mkdir()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
 script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
pause_check()
assert read(p/'rain-restart-comparator-v1/terminal.json')['passed']
assert read(original/'study-terminal.json')['passed'] and read(original/'comparison.json')['passed']
audit=read(case_source/'audit.json');assert audit['passed'] and len(audit['frames'])==108
assert [audit['frames'][i-1]['time_s'] for i in (36,72,108)]==[3600.,7200.,10800.]
expected_solver='d55ffe42c3b61e0891fec513489c9f5881830733ac9e003a01a0cf69f495905d'
expected_result='4e0559856cb3124b72a05a6d818bc7867403db3764c6ca195135872b248f3379'
baseline=case_source/'result.json';assert sha(baseline)==expected_result
assert expected_result in audit['sha256'].values()
solver=out/'soil-fem-3d.exe';shutil.copyfile(original/'solver/soil-fem-3d.exe',solver);assert sha(solver)==expected_solver
for name in ['mesh.json','case.json','seed-checkpoint.json']:shutil.copyfile(case_source/name,out/name)
assert sha(out/'seed-checkpoint.json')==read(case_source/'seed-identity.json')['checkpoint_sha256']
assert all(sha(out/name)==sha(case_source/name) for name in ['mesh.json','case.json','seed-checkpoint.json'])
binding={str(q.resolve()):sha(q) for q in [baseline,case_source/'audit.json',case_source/'seed-checkpoint.json',
    case_source/'mesh.json',case_source/'case.json',solver,out/'mesh.json',out/'case.json',out/'seed-checkpoint.json',
    p/'compare_richards_stream.py',Path(__file__).resolve()]}
runtime=p/'python-blocked-integration-v1/runtime'
sources=read(p/'python-blocked-integration-v1/quality-v2/source-snapshot.json')
assert all(sha(runtime/name)==h for name,h in sources.items())
binding.update({str(runtime/name):h for name,h in sources.items()})
binding[str(p/'rain-refined-restart-v1/retained-harness-failure.json')]=sha(p/'rain-refined-restart-v1/retained-harness-failure.json')
dump(out/'source-binding.json',binding)
normalizer=(original/'runtime/src/richards.rs').read_text('utf-8').split('fn normalized(run: &Report)',1)[1].split('\n}\n',1)[0]
assert normalizer.count('linear.elapsed_seconds = 0.')==1
assert normalizer.count('= 0.')==1
dump(out/'protocol.json',dict(elements=82944,period_boundaries_s=[3600,7200],accepted_counts=[36,72],
    original_solver_sha256=expected_solver,original_audited_result_sha256=expected_result,
    original_Rust_normalizer_sha256=hashlib.sha256(normalizer.encode()).hexdigest(),
    timer_exclusion='attempts.item.iterations.item.linear.elapsed_seconds',
    all109_frames_including_initial_and_all_remaining_fields_must_match=True,
    original_inputs_periods_materials_boundaries_tolerances_and_binary_unchanged=True,
    baseline_first_accepted_checkpoint_preserved_before_new_solve=True,
    no_original_audit_or_result_or_checkpoint_will_be_overwritten=True,p6_complete=False))
sys.path[:0]=[str(p),str(runtime/'tools')]
from compare_richards_stream import compare
from streamed_stage_source import CheckedReader
import ijson
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',
    PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
def unchanged():
    assert all(sha(Path(q))==h for q,h in binding.items())
def call(label,command,folder,*,expected_exit=0,error_text=None):
    pause_check();command=list(map(str,command));start=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(folder/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(exit_code=code,elapsed_seconds=time.monotonic()-start,expected_exit=expected_exit)
    dump(folder/(label+'-terminal.json'),row);print(json.dumps(dict(step=label,**row)),flush=True)
    assert code==expected_exit,label
    if error_text:assert error_text in (folder/(label+'.log')).read_text('utf-8')
    return row
def resource(label,folder):
    while True:
        pause_check()
        raw=subprocess.check_output(['pwsh','-NoProfile','-Command',
            '(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory'],encoding='utf-8-sig',
            creationflags=subprocess.CREATE_NO_WINDOW)
        free=int(raw.strip())*1024
        if free>=24*(1<<30):break
        time.sleep(5)
    dump(folder/(label+'-resource.json'),dict(free_physical_bytes=free,required_free_bytes=24*(1<<30),
        exclusive_resources_or_speed_not_claimed=True))
# Use the preserved one-step state for fast negative controls on the same82944 mesh.
negative=out/'negative-controls';negative.mkdir()
bad=negative/'changed-payload.json'
with (out/'seed-checkpoint.json').open('rb') as source,bad.open('xb') as target:
    first=source.read(4096);match=re.search(rb'("payload_sha256"\s*:\s*")([0-9a-f]{64})(")',first)
    assert match
    value=match.group(2);replacement=(b'0' if value[:1]!=b'0' else b'1')+value[1:]
    target.write(first[:match.start(2)]+replacement+first[match.end(2):]);shutil.copyfileobj(source,target)
call('payload-refusal',[solver,'resume-richards',out/'mesh.json',out/'case.json',bad,negative/'forbidden-payload-result.json'],negative,
    expected_exit=1,error_text='hydraulic checkpoint payload changed')
assert not (negative/'forbidden-payload-result.json').exists()
changed=read(out/'case.json');rain=next(b['condition'] for b in changed['analysis']['periods'][1]['boundaries'] if b['face_set']=='rain')
assert rain['rainfall_vertical_m_s']==0.;rain['rainfall_vertical_m_s']=1e-8
dump(negative/'changed-future-period.json',changed)
call('future-input-refusal',[solver,'resume-richards',out/'mesh.json',negative/'changed-future-period.json',out/'seed-checkpoint.json',negative/'forbidden-input-result.json'],negative,
    expected_exit=1,error_text='checkpoint inputs, executable, runtime, frames or hydraulic history changed')
assert not (negative/'forbidden-input-result.json').exists();unchanged()
dump(negative/'terminal.json',dict(passed=True,refined82944_mesh=True,
    payload_mutation_and_future_inactive_period_input_change_rejected=True,
    no_rejected_result_written=True,original_inputs_and_seed_preserved=True,p6_complete=False))
def checkpoint_values(path,prefix):
    with path.open('rb') as stream:
        yield from ijson.items(CheckedReader(stream,32*1024**3),prefix,use_float=True)
summaries=[]
for count,time_s in [(36,3600.),(72,7200.)]:
    pause_check();folder=out/('accepted-%03d'%count);folder.mkdir();resource('checkpoint',folder)
    checkpoint=folder/'checkpoint.json'
    if count==36:
        old=p/'rain-refined-restart-v1';retained=read(old/'retained-harness-failure.json')
        saved=old/'accepted-036/checkpoint.json'
        assert retained['passed'] and sha(saved)==retained['checkpoint_sha256']
        assert retained['last_committed_time_s']==3600. and retained['accepted_steps']==36
        assert all(sha(out/name)==sha(old/name) for name in ['mesh.json','case.json','soil-fem-3d.exe'])
        shutil.copyfile(saved,checkpoint);assert sha(checkpoint)==retained['checkpoint_sha256']
        dump(folder/'checkpoint-reused.json',dict(passed=True,original_valid_checkpoint_byte_exact=True,
            source=str(saved),sha256=sha(checkpoint),original_failed_harness_preserved=True,
            no_physical_inputs_or_tolerances_or_time_counts_changed=True,p6_complete=False))
    else:
        call('checkpoint',[solver,'checkpoint-richards',out/'mesh.json',out/'case.json',count,checkpoint],folder)
    digest=sha(checkpoint)
    frame_times=list(checkpoint_values(checkpoint,'run.frames.item.state.time_s'))
    accepted=sum(checkpoint_values(checkpoint,'run.attempts.item.accepted'));frames=len(frame_times)
    header={prefix:next(checkpoint_values(checkpoint,prefix)) for prefix in
        ['schema','run.schema','run.complete','run.last_committed.time_s']}
    checkpoint_runtime=next(checkpoint_values(checkpoint,'run.runtime'))
    header.update({'run.runtime.'+key:value for key,value in checkpoint_runtime.items()})
    assert header['schema']=='soil-fem-richards-checkpoint/1' and header['run.schema']=='soil-fem-richards-result/1'
    assert not header['run.complete'] and header['run.last_committed.time_s']==time_s
    assert accepted==count and frames==count+1 and frame_times[0]==0. and frame_times[-1]==time_s
    dump(folder/'checkpoint-identity.json',dict(passed=True,sha256=digest,accepted_steps=accepted,frames=frames,header=header,
        own_frozen_solver_produced_checkpoint=True,will_be_validated_by_own_resume_API=True,p6_complete=False))
    resource('resume',folder)
    result=folder/'result.json'
    call('resume',[solver,'resume-richards',out/'mesh.json',out/'case.json',checkpoint,result],folder)
    assert sha(checkpoint)==digest
    pause_check();start=time.monotonic()
    report=compare(baseline,result,expected_left_sha256=expected_result)
    assert report['complete'] and report['schema']=='soil-fem-richards-result/1'
    assert report['frames']==109 and report['last_committed_time_s']==10800.
    report.update(restart_time_s=time_s,checkpoint_sha256=digest,original_full_discrete_audit_sha256=sha(case_source/'audit.json'),
        original_audit_inherited_by_all_non_time_fields_identity_not_rerun=True,elapsed_seconds=time.monotonic()-start,p6_complete=False)
    dump(folder/'comparison.json',report);unchanged();summaries.append(dict(restart_time_s=time_s,accepted_steps=count,
        comparison_sha256=sha(folder/'comparison.json'),all_non_time_fields_equal=True))
    print(json.dumps(dict(restart_time_s=time_s,all109_frames_including_initial_and_non_time_fields_equal=True)),flush=True)
unchanged()
dump(out/'terminal.json',dict(passed=True,elements=82944,restarts=summaries,negative_controls_passed=True,
    exact_original_solver_and_inputs_and_baseline_and_checkpoints_preserved=True,
    all109_frames_including_initial_heads_fluxes_storage_surface_decisions_attempts_runtime_and_bindings_equal=True,
    only_original_Rust_solver_timer_prefix_excluded=True,whole_result_JSON_not_materialized=True,
    full_independent_discrete_audit_inherited_from_identity_not_reexecuted=True,
    general_ground_ponding_mechanics_stability_wet_contact_and_P6_separate=True,p6_complete=False))
print('PASS refined82944 rain period-boundary restarts; all109 frames including initial and every non-time field exact',flush=True)
