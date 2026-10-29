"""Verify every endpoint of the real CAD BIRTH/RESET/MATSET history independently.

The retained CAD example is elastic for factors1/2/3. The new Rust tests cover
plastic factor histories separately. Neither check creates an Fs bound.
"""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'cad-controls';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
with (p/'cad-control-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
print('Waiting for full277/171 and actual column CLI guards before CAD endpoint controls',flush=True)
while not (p/'terminal.json').exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    for q in (p/'quality-v2').glob('*-terminal.json'):
        if read(q).get('exit_code')!=read(q).get('expected_exit_code',0):raise RuntimeError('Stage-gate QA failed before CAD controls: '+q.name)
    time.sleep(5)
assert read(p/'terminal.json')['passed']
sources=read(p/'quality-v2/source-snapshot.json');assert all(sha(p/'runtime'/n)==h for n,h in sources.items())
solver=p/'solver/soil-fem-3d.exe';assert sha(solver)==read(p/'terminal.json')['solver_sha256']
baseline=root/'runs/reset-replay-test-20261024-012858-00a4e04d';assert read(baseline/'verification.json')['passed']
for old,new in [('cad/ground/mesh.json','mesh.json'),('cad/ground/ground.json','ground.json'),('cad/stages.json','case.json'),('replay/specification.json','spec.json')]:
    shutil.copy2(baseline/old,out/new);assert sha(baseline/old)==sha(out/new)
binding={n:sha(out/n) for n in ['mesh.json','ground.json','case.json','spec.json']};case=read(out/'case.json');spec=read(out/'spec.json')
assert spec['factors']==[1.,2.,3.] and len(case['stages'])==10
records=[];reports=[];env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
def call(label,args,folder=out):
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));started=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        with (folder/(label+'-child.json')).open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-started);records.append(row)
    with (folder/(label+'-terminal.json')).open('x') as f:json.dump(row,f,indent=2)
    print(json.dumps(row),flush=True);assert code==0
paths=[out/n for n in ['mesh.json','ground.json','case.json','spec.json']]
call('full',[solver,'replay-stages',*paths,out/'full.json']);full=read(out/'full.json');assert all(t['outcome']=='completed' for t in full['trials'])
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
for count in range(len(case['stages'])+1):
    target=out/('gate'+str(count));target.mkdir();gatepath=target/'report.json'
    call('gate'+str(count),[solver,'replay-through-stage',*paths,str(count),gatepath]);gate=read(gatepath)
    assert gate['requested_through_construction_stages']==count and gate['stability_classification']=='not_evaluated'
    assert gate['strength_applied_from']=='original_inputs_before_initial_equilibration'
    for index,(a,b) in enumerate(zip(gate['replay']['trials'],full['trials'],strict=True)):
        assert a['factor']==b['factor'] and a['outcome']=='completed'
        assert stable(a['run']['frames'])==stable(b['run']['frames'][:count+1])
        assert stable(a['run']['last_committed'])==stable(b['run']['frames'][count]['state'])
        trial=target/('trial'+str(index));(trial/'ground').mkdir(parents=True);shutil.copy2(out/'mesh.json',trial/'ground/mesh.json')
        with (trial/'ground/ground.json').open('x') as f:json.dump(a['adjusted_ground'],f)
        with (trial/'stage-results.json').open('x') as f:json.dump(a['run'],f)
        audit=trial/'independent.json'
        call('independent',[sys.executable,p/'runtime/tools/verify_construction_linear.py',trial,audit],trial)
        report=read(audit)
        reports.append(dict(gate=count,factor=a['factor'],frames=len(report['frames']),maximum_displacement_relative_error=max(r['relative_3d_displacement_error'] for r in report['frames']),maximum_stress_relative_error=max(r['relative_all_point_stress_error'] for r in report['frames']),source_sha256=sha(gatepath),independent_audit_sha256=sha(audit)))
assert all(sha(out/n)==h for n,h in binding.items()) and all(sha(p/'runtime'/n)==h for n,h in sources.items())
with (out/'terminal.json').open('x') as f:json.dump(dict(passed=True,actual_retained_cad_case=True,all11_endpoints_and3_factors=True,all_prefix_frames_and_end_histories_exact_to_full_replay=True,
    independent_3d_elastic_history_and_duffy_integration_and_svd_for_all33_cases=True,independent_frame_comparisons=sum(r['frames'] for r in reports),reports=reports,
    birth_death_material_set_display_reset_and_explicit_state_reset_covered=True,plastic_prefix_gate_coverage_in_separate_rust_tests=True,
    fs_search_and_physical_stability_classification_not_implemented=True,post_unreduced_origin_reduction_not_implemented=True,
    original_inputs_sha256=binding,solver_sha256=sha(solver),records=records,production_not_modified=True,p6_complete=False),f,indent=2)
print('PASS real CAD11 endpoints/3 factors, exact prefixes and independent elastic construction histories; no Fs bound',flush=True)
