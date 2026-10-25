"""Observe a live actual-CAD mechanics job, pause it, and resume in another process."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from verify_default_solver import deterministic


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,value):
    with p.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2,allow_nan=False)


def verify(solver,out):
    root=Path(__file__).resolve().parents[1];base=root/'runs/mcc-tet10-cycle-acceptance-20261025-v1/inputs'
    expected=root/'runs/mixed-mcc-ordering-legacy-20261025-v1/legacy-mcc.json'
    sources=[solver,Path(__file__),root/'tools/verify_default_solver.py',root/'tools/run_logged.py',expected,
             base/'mesh.json',base/'ground.json',base/'slope-cycle-case.json']
    fingerprints={str(p):sha(p) for p in sources}
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    binary=out/'solver/soil-fem-3d.exe'
    def invoke(tag,exe,args,code=0):
        command=[sys.executable,'-X','utf8',str(root/'tools/run_logged.py'),str(out/(tag+'.log')),str(exe),*map(str,args)]
        completed=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',cwd=root)
        if completed.returncode!=code:raise AssertionError(tag+': '+completed.stdout[-2000:])
        print(tag+' completed with exit '+str(code),flush=True)
    job=out/'job'
    invoke('prepare',binary,['prepare-stages-job',base/'mesh.json',base/'ground.json',base/'slope-cycle-case.json',job])
    frozen=job/'solver/soil-fem-3d.exe';logpath=out/'live-job.log'
    with logpath.open('x',encoding='utf-8') as log:
        command=[str(frozen),'run-stages-job',str(job)]
        log.write('Command: '+repr(command)+'\n');log.flush()
        child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,cwd=root)
        start=time.monotonic();observed=None
        try:
            while time.monotonic()-start<60:
                if child.poll() is not None:raise AssertionError('job exited before a live Newton observation')
                path=job/'progress.jsonl'
                if path.exists():
                    with path.open(encoding='utf-8') as stream:
                        for line in stream:
                            try:row=json.loads(line)
                            except json.JSONDecodeError:continue
                            if row['progress']['kind']=='newton' and row['progress']['stage_index']>=1:observed=row;break
                if observed is not None:break
                time.sleep(.02)
            assert observed is not None and child.poll() is None
            write(out/'observed-live-newton.json',dict(process_id=child.pid,owned_process_poll=None,event=observed))
            invoke('request-stop',frozen,['stop-stages-job',job])
            while child.poll() is None and time.monotonic()-start<120:time.sleep(.05)
            assert child.poll()==2,'cooperative pause did not finish with checkpoint exit 2'
            log.write('\nExit code: 2\n')
        finally:
            if child.poll() is None:
                # This test owns the Popen handle. Never terminate a recorded
                # PID from a status file or any unrelated application.
                child.terminate();child.wait(timeout=20)
    state=read(job/'status.json');saved=read(job/'checkpoint.json')
    assert state['state']=='checkpointed' and not (job/'result.json').exists()
    assert saved['run']['failure_reason']=='checkpoint_requested'
    assert saved['run']['geostatic']['complete'] and saved['run']['attempts'][-1]['accepted']
    assert saved['run']['last_committed']['stage_index']>=1 and saved['run']['continuation'] is not None
    assert any(any(abs(x)>1e-10 for x in s['plastic_strain_compression_positive']) for s in saved['run']['last_committed']['material_states'])
    resumed=out/'resumed'
    invoke('prepare-resume',frozen,['prepare-stages-resume',job,resumed])
    invoke('resume',resumed/'solver/soil-fem-3d.exe',['run-stages-job',resumed])
    actual=read(resumed/'result.json');reference=read(expected)
    assert actual['complete'] and read(resumed/'status.json')['state']=='succeeded'
    assert actual['provenance']['solver_binary_sha256']==sha(binary)
    assert json.dumps(deterministic(actual),sort_keys=True,allow_nan=False)==json.dumps(deterministic(reference),sort_keys=True,allow_nan=False),'non-time/non-binary result differs from the accepted legacy CAD cycle'
    for folder in (job,resumed):
        status=read(folder/'status.json');rows=[json.loads(s) for s in (folder/'progress.jsonl').read_text(encoding='utf-8').splitlines()]
        assert len(rows)==status['progress_events']
        assert all(r['job_id']==status['job_id'] and r['sequence']==i+1 for i,r in enumerate(rows))
        assert sha(folder/status['artifact'])==status['artifact_sha256']
    assert fingerprints=={str(p):sha(p) for p in sources},'source/input changed during the verification'
    report=dict(schema='soil-fem-stage-job-cad-verification/1',passed=True,goal_complete=False,p6_complete=False,
        scope='External controlled CLI, original SolidCraft slope/5 kPa cycle/local pit, existing displacement MCC formulation',
        live_newton_observed_before_stop_request=True,accepted_construction_boundary_checkpoint=True,
        plastic_history_present=True,separate_process_resume=True,all_nonruntime_fields_exact=True,
        excluded_fields=['elapsed_seconds','solver_binary_sha256'],solver_sha256=sha(binary),
        native_ui_stop_and_resume_connected=False,input_and_source_sha256=fingerprints,
        evidence_sha256={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    write(out/'verification.json',report)
    print('Actual CAD controlled job: live Newton, cooperative checkpoint, separate-process resume, all nonruntime fields exact.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('solver',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.solver.resolve(),a.output.resolve())
