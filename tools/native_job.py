"""External FEM job worker for the native SolidCraft host.

The host starts this worker as a child and reads status.json. Imported CAD never
selects a program. Cancellation is a request to this worker, which owns the
PowerShell child and its descendants; observers never kill a recorded PID.
"""
import argparse
from datetime import datetime,timezone
import json
import math
import os
from pathlib import Path
import queue
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from import_solidcraft_capture import read,validate
from verify_acceptance_packet import sha
import native_workflow_snapshot

ROOT=Path(__file__).resolve().parents[1]
TERMINAL={'succeeded','failed','cancelled','checkpointed'}


def watch_host(stream):
    """EOF on the dedicated pipe requests cleanup, including after a host crash.

    Only the native host opts into this protocol. No command is accepted on the
    pipe; data cannot select programs or modify the immutable analysis input.
    """
    disconnected=threading.Event()
    def observe():
        try:
            # A daemon blocked in BufferedReader.read holds Python's buffered-I/O
            # lock during interpreter shutdown. Use the raw descriptor instead.
            while os.read(stream.fileno(),1):pass
        finally:disconnected.set()
    threading.Thread(target=observe,name='fem-host-lifetime',daemon=True).start()
    return disconnected


def utc():return datetime.now(timezone.utc).isoformat()


def atomic(path,value):
    temporary=path.with_name(path.name+'.'+uuid.uuid4().hex+'.part')
    with temporary.open('x',encoding='utf-8') as stream:json.dump(value,stream,ensure_ascii=False,allow_nan=False)
    os.replace(temporary,path)


def copy_checked(source,destination,digest):
    with source.open('rb') as incoming,destination.open('xb') as outgoing:shutil.copyfileobj(incoming,outgoing,1024*1024)
    if sha(source)!=digest or sha(destination)!=digest:raise ValueError('input changed during job snapshot: '+source.name)


def status(folder):
    value=read(folder.resolve()/'status.json',256*1024)
    if value.get('schema') not in ('soil-fem-native-job-status/1','soil-fem-native-job-status/2'):raise ValueError('unsupported job status')
    # This is a last-reported snapshot, not proof that a process is still alive.
    value['observation']='last worker report; the host must also observe its owned child process'
    return value


def request_cancel(folder):
    folder=folder.resolve();owner=read(folder/'job.json',256*1024);current=status(folder)
    if owner.get('schema') not in ('soil-fem-native-job/1','soil-fem-native-job/2') or owner['job_id']!=current['job_id']:raise ValueError('job identity mismatch')
    if current['state'] in TERMINAL:return dict(requested=False,state=current['state'],job_id=owner['job_id'])
    request=folder/'cancel-request.json';value=dict(job_id=owner['job_id'],requested_utc=utc())
    try:
        with request.open('x',encoding='utf-8') as stream:json.dump(value,stream)
    except FileExistsError:
        if read(request,4096).get('job_id')!=owner['job_id']:raise ValueError('cancellation request identity mismatch')
    return dict(requested=True,job_id=owner['job_id'])


def stop_child(child):
    """Stop only the still-owned process tree; retain the Popen process handle."""
    if child.poll() is not None:return
    if os.name=='nt':
        result=subprocess.run(['taskkill.exe','/PID',str(child.pid),'/T','/F'],capture_output=True,
            creationflags=subprocess.CREATE_NO_WINDOW,timeout=20)
        if result.returncode and child.poll() is None:raise RuntimeError('could not stop the owned job process tree')
    else:
        try:os.killpg(child.pid,signal.SIGTERM)
        except ProcessLookupError:return
        try:child.wait(timeout=3)
        except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL)
    child.wait(timeout=20)


def phase(folder):
    workflow=folder/'workflow'
    for name,path in [('checking_view',workflow/'view/view-manifest.json'),('preparing_view',workflow/'analysis/accepted-stage-run.json'),
        ('checking_result',workflow/'analysis/stage-results.json'),('solving',workflow/'analysis/ground/mesh.json'),
        ('meshing',workflow/'conformal/native-capture-import.json')]:
        if path.is_file():return name
    return 'importing_geometry'


def run(capture,folder,*,solver=None,mesh_size=2.,order=2,quadrature=3,split=None,expected_input=None,host_disconnected=None,host_token=None,cooperative=False,resume_from=None):
    capture=capture.resolve();folder=folder.resolve()
    if folder.is_relative_to(capture):raise ValueError('job output must be outside the source capture')
    if not math.isfinite(mesh_size) or mesh_size<=0 or order not in (1,2) or quadrature not in (3,4,5):raise ValueError('invalid mesh settings')
    if host_token is not None and (not isinstance(host_token,str) or not host_token or len(host_token)>128):raise ValueError('invalid native host token')
    recorded,_,_,_=validate(capture)
    if expected_input is not None and recorded['input_sha256']!=expected_input:raise ValueError('capture differs from the host input snapshot')
    solver=(solver or ROOT/'target/release/soil-fem-3d.exe').resolve();solver_digest=sha(solver)
    shell=shutil.which('pwsh')
    if shell is None:raise ValueError('PowerShell 7 is required for the external FEM workflow')
    if split is not None:split=split.resolve();read(split)
    sources=[*sorted((ROOT/'tools').glob('*.py')),ROOT/'scripts/Invoke-SolidCraftCapture.ps1',ROOT/'scripts/Invoke-Stages.ps1',ROOT/'requirements.lock.txt']
    source_hashes={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sources}
    previous_owner=None;receipt=None
    if resume_from is not None:
        import native_stage_profile
        resume_from=resume_from.resolve();cooperative=True
        previous_owner,receipt=native_stage_profile.previous(resume_from,solver_digest,source_hashes,expected_input)
        if recorded['input_sha256']!=previous_owner['input_sha256']:raise ValueError('resume capture differs')
        mesh_size,order,quadrature=(previous_owner[k] for k in ('mesh_size_m','element_order','quadrature_order'))
    folder.mkdir(parents=True,exist_ok=False)
    job_id=uuid.uuid4().hex;started=time.monotonic()
    owner=dict(schema='soil-fem-native-job/2' if cooperative else 'soil-fem-native-job/1',job_id=job_id,created_utc=utc(),input_sha256=recorded['input_sha256'],
        host_token=host_token,
        native_exporter_sha256=recorded['exporter_binary_sha256'],solver_sha256=solver_digest,
        mesh_size_m=mesh_size,element_order=order,quadrature_order=quadrature,tool_sources_sha256=source_hashes,
        workflow_sources='runtime',profile='cooperative_staged' if cooperative else 'legacy_process',
        parent_job_id=previous_owner['job_id'] if previous_owner else None)
    atomic(folder/'job.json',owner)
    current=dict(schema='soil-fem-native-job-status/2' if cooperative else 'soil-fem-native-job-status/1',job_id=job_id,input_sha256=owner['input_sha256'],
        host_token=host_token,
        state='preparing',phase='snapshot',worker_pid=os.getpid(),child_pid=None,created_utc=owner['created_utc'],
        updated_utc=utc(),elapsed_seconds=0.,exit_code=None,result=None,log='worker.log',tail=[],
        stop_policy='cooperative_staged_checkpoint' if cooperative else 'process_stop_without_new_checkpoint',resume_available=False)
    def publish(**changes):
        current.update(changes);current.update(updated_utc=utc(),elapsed_seconds=time.monotonic()-started)
        atomic(folder/'status.json',current)
    def cancelled():
        path=folder/'cancel-request.json'
        return (host_disconnected is not None and host_disconnected.is_set()) or (path.exists() and read(path,4096).get('job_id')==job_id)
    publish();child=None;reader=None
    try:
        runtime=native_workflow_snapshot.capture(ROOT,folder/'runtime',source_hashes)
        owned=folder/'capture';owned.mkdir()
        for name,digest in recorded['files_sha256'].items():
            if cancelled():publish(state='cancelled',phase='snapshot');return 2
            copy_checked(capture/name,owned/name,digest)
        copy_checked(capture/'capture.json',owned/'capture.json',sha(capture/'capture.json'))
        copied,_,_,_=validate(owned)
        if copied!=recorded:raise ValueError('capture changed during job preparation')
        (folder/'solver').mkdir();pinned=folder/'solver/soil-fem-3d.exe';copy_checked(solver,pinned,solver_digest)
        if resume_from is not None:
            import native_stage_profile
            native_stage_profile.clone_context(resume_from,folder,pinned,owner)
            return native_stage_profile.execute(folder,runtime,pinned,owner,publish,cancelled,stop_child,previous_folder=resume_from,receipt=receipt)
        command=[shell,'-NoProfile','-File',str(runtime/'scripts/Invoke-SolidCraftCapture.ps1'),'-CaptureFolder',str(owned),
            '-MeshSizeM',str(mesh_size),'-ElementOrder',str(order),'-QuadratureOrder',str(quadrature),
            '-SolverPath',str(pinned),'-OutputDirectory',str(folder/'workflow'),'-PythonPath',sys.executable]
        if split is not None:
            spec=folder/'joint-split.json';copy_checked(split,spec,sha(split));command+=['-SplitFile',str(spec)]
        if cooperative:command+=['-PrepareOnly']
        if cancelled():publish(state='cancelled',phase='snapshot');return 2
        lines=queue.Queue(maxsize=256)
        options=dict(creationflags=subprocess.CREATE_NEW_PROCESS_GROUP|subprocess.CREATE_NO_WINDOW) if os.name=='nt' else dict(start_new_session=True)
        with (folder/'worker.log').open('x',encoding='utf-8') as log:
            child=subprocess.Popen(command,cwd=runtime,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                text=True,encoding='utf-8',errors='replace',env={**os.environ,'PYTHONUTF8':'1'},**options)
            def collect():
                try:
                    for line in child.stdout:log.write(line);log.flush();lines.put(line)
                except Exception as error:lines.put(('reader_error',str(error)))
                finally:lines.put(None)
            reader=threading.Thread(target=collect,name='fem-job-output',daemon=True);reader.start()
            publish(state='running',phase='importing_geometry',child_pid=child.pid)
            last_write=0.;eof=False;tail=[];cancelled_run=False
            while child.poll() is None or not eof:
                if cancelled() and not cancelled_run:
                    publish(state='cancelling',phase=phase(folder));stop_child(child);cancelled_run=True
                try:
                    line=lines.get(timeout=.2)
                    if line is None:eof=True
                    elif isinstance(line,tuple):raise RuntimeError('job output reader failed: '+line[1])
                    else:tail=(tail+[line.rstrip()[-1024:]])[-30:]
                except queue.Empty:pass
                if time.monotonic()-last_write>.5:
                    publish(phase=phase(folder),tail=tail);last_write=time.monotonic()
            reader.join(timeout=5);code=child.wait();child.stdout.close()
            if cancelled_run or cancelled():publish(state='cancelled',exit_code=code,phase=phase(folder),tail=tail);return 2
        if code!=0:publish(state='failed',exit_code=code,phase=phase(folder),tail=tail);return 1
        if cooperative:
            import native_stage_profile
            return native_stage_profile.execute(folder,runtime,pinned,owner,publish,cancelled,stop_child)
        native_workflow_snapshot.verify(runtime,source_hashes)
        manifest_path=folder/'workflow/view/view-manifest.json';manifest=read(manifest_path)
        audit=read(folder/'workflow/view-verification.json')
        if manifest.get('complete') is not True or manifest.get('input_sha256')!=owner['input_sha256'] or manifest.get('solver_binary_sha256')!=solver_digest:
            raise ValueError('completed result differs from the owned job snapshot')
        if audit.get('passed') is not True or audit.get('view_manifest_sha256')!=sha(manifest_path):raise ValueError('view audit does not match the completed result')
        publish(state='succeeded',phase='complete',exit_code=0,result='workflow/view/view-manifest.json',
            result_sha256=sha(manifest_path),view_audit_sha256=sha(folder/'workflow/view-verification.json'))
        return 0
    except BaseException as error:
        if child is not None and child.poll() is None:stop_child(child)
        publish(state='cancelled' if isinstance(error,KeyboardInterrupt) else 'failed',error=str(error),phase=phase(folder))
        if isinstance(error,KeyboardInterrupt):return 2
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();commands=parser.add_subparsers(dest='command',required=True)
    runner=commands.add_parser('run');runner.add_argument('capture',type=Path);runner.add_argument('job',type=Path)
    runner.add_argument('--solver',type=Path);runner.add_argument('--mesh-size-m',type=float,default=2.);runner.add_argument('--order',type=int,default=2)
    runner.add_argument('--quadrature',type=int,default=3);runner.add_argument('--split',type=Path);runner.add_argument('--expected-input-sha256')
    runner.add_argument('--watch-stdin',action='store_true',help='native host pipe: EOF requests owned-process cleanup')
    runner.add_argument('--host-token',help='opaque identity supplied by the native process owner, never a PID to signal')
    runner.add_argument('--cooperative',action='store_true',help='save solver checkpoints at accepted increments; requires staged-job CLI support')
    resume=commands.add_parser('resume');resume.add_argument('previous',type=Path);resume.add_argument('job',type=Path)
    resume.add_argument('--solver',type=Path);resume.add_argument('--expected-input-sha256');resume.add_argument('--host-token')
    resume.add_argument('--watch-stdin',action='store_true')
    for name in ('status','cancel'):commands.add_parser(name).add_argument('job',type=Path)
    args=parser.parse_args()
    if args.command=='run':sys.exit(run(args.capture,args.job,solver=args.solver,mesh_size=args.mesh_size_m,order=args.order,quadrature=args.quadrature,split=args.split,expected_input=args.expected_input_sha256,host_disconnected=watch_host(sys.stdin.buffer) if args.watch_stdin else None,host_token=args.host_token,cooperative=args.cooperative))
    if args.command=='resume':sys.exit(run(args.previous/'capture',args.job,solver=args.solver,expected_input=args.expected_input_sha256,host_disconnected=watch_host(sys.stdin.buffer) if args.watch_stdin else None,host_token=args.host_token,cooperative=True,resume_from=args.previous))
    print(json.dumps(request_cancel(args.job) if args.command=='cancel' else status(args.job),ensure_ascii=False))
