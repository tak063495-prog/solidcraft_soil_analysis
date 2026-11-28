"""Resume complete phases of an interrupted frozen MCC reference study.

This does not recover any increment inside an unfinished phase. Original
source copies, input files and saved results stay in place without replacement.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import signal
import threading
from contextlib import contextmanager

INPUTS=('mesh.json','ground.json','case.json','mesh-identity.json')
BACKENDS={
 'mcc_projected_sparse_reference.SparseProjectedSystem':(),
 'mcc_projected_factor_reference.ReusedFactorProjectedSystem':('--factor-reuse',),
 'mcc_projected_matrixfree_reference.MatrixFreeFactorProjectedSystem':('--matrix-free-factor-reuse',),
 'mcc_projected_coercive_reference.CoerciveSharedFactorProjectedSystem':('--coercive-factor-reuse',),
}

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def path_identity(path):
    # Windows TEMP may use an 8.3 alias or a different drive/name case from
    # Python's __file__. Compare canonical file identity, retaining raw hashes.
    return os.path.normcase(str(Path(path).resolve()))
def atomic(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(path)

def preflight(folder):
    folder=folder.resolve()
    progress=read(folder/'progress.json');steps=progress['steps']
    if not isinstance(steps,int) or isinstance(steps,bool) or steps<2:
        raise ValueError('invalid frozen reference step count')
    snapshot=read(folder/'source-snapshot.json');frozen=folder/'reference-sources'
    required=('mcc_projected_sparse_run.py','mcc_projected_sparse_reference.py',
      'mcc_projected_volume_diagnostic.py','mcc_projected_cycle_diagnostic.py',
      'mcc_rate_fem.py','mcc_rate_midpoint.py','fem_reference.py',
      'verify_high_order.py','high_order_io.py','requirements.lock.txt','requirements-reference.lock.txt')
    if not all((frozen/name).is_file() for name in required):
        raise ValueError('a required frozen reference source or dependency lock is missing')
    fingerprints={}
    # Some refinement studies deliberately keep the completed coarser study's
    # program while the root snapshot contains later development. Require that
    # explicit historical binding AND the completed earlier cycle, not a silent
    # exception based only on today's file contents.
    overrides={};version_file=folder/'reference-version-binding.json'
    if version_file.exists():
        version=read(version_file);prior_files=version.get('files_sha256',{})
        if version.get('same_programs_as_completed_200') is not True or version.get('steps')!=steps or not prior_files:
            raise ValueError('invalid historical reference-version binding')
        parents={Path(n).resolve().parent for n in prior_files}
        if len(parents)!=1:raise ValueError('historical sources do not belong to one frozen reference')
        prior_sources=parents.pop();prior=prior_sources.parent
        if prior_sources.name!='reference-sources' or not prior.is_relative_to(folder.parent):
            raise ValueError('historical reference is outside the study collection')
        prior_cycle=prior/'cycle.json'
        if sha(prior_cycle)!=version.get('previous_completed_cycle_sha256'):
            raise ValueError('historical completed cycle changed')
        prior_bindings=read(prior_cycle)['input_and_source_sha256']
        by_identity={}
        for name,digest in prior_bindings.items():
            key=path_identity(name)
            if key in by_identity:raise ValueError('ambiguous historical source identity')
            by_identity[key]=digest
        for name,digest in prior_files.items():
            if sha(Path(name))!=digest or by_identity.get(path_identity(name))!=digest:
                raise ValueError('historical completed source binding changed')
            overrides[Path(name).name]=digest
            fingerprints[name]=digest
        fingerprints[str(version_file)]=sha(version_file)
        fingerprints[str(prior_cycle)]=sha(prior_cycle)
    for path in sorted(frozen.iterdir()):
        if not path.is_file() or path.suffix not in ('.py','.txt'):continue
        key='tools/'+path.name if path.suffix=='.py' else path.name
        if key not in snapshot or sha(path)!=overrides.get(path.name,snapshot[key]):
            raise ValueError('frozen reference source changed: '+path.name)
        fingerprints[str(path)]=sha(path)
    for name in INPUTS:
        path=folder/'inputs'/name
        read(path)
        fingerprints[str(path)]=sha(path)
    backend='mcc_projected_sparse_reference.SparseProjectedSystem'
    if (frozen/'mcc_projected_factor_reference.py').exists():
        backend='mcc_projected_factor_reference.ReusedFactorProjectedSystem'
    if (frozen/'mcc_projected_matrixfree_reference.py').exists():
        backend='mcc_projected_matrixfree_reference.MatrixFreeFactorProjectedSystem'
    if (frozen/'mcc_projected_coercive_reference.py').exists():
        backend='mcc_projected_coercive_reference.CoerciveSharedFactorProjectedSystem'
    saved=[]
    for phase,schema in [('initial','soil-fem-projected-volume-mcc-diagnostic/1'),
                         ('cycle','soil-fem-projected-mcc-cycle-diagnostic/1')]:
        path=folder/(phase+'.json')
        if not path.exists():continue
        data=read(path)
        if data.get('schema')!=schema or data.get('reference_steps')!=steps or data.get('reference_backend')!=backend:
            raise ValueError('saved phase schema, steps or backend differs: '+phase)
        if phase=='initial' and data.get('load_fraction')!=1.:
            raise ValueError('saved initialization is incomplete')
        if phase=='cycle' and ('initial' not in saved or data.get('initial_reference_sha256')!=sha(folder/'initial.json')):
            raise ValueError('cycle is detached from its saved initialization')
        bindings=data.get('input_and_source_sha256',{})
        for name in INPUTS:
            if bindings.get(str(folder/'inputs'/name))!=fingerprints[str(folder/'inputs'/name)]:
                raise ValueError('saved phase input binding differs: '+name)
        for name,digest in bindings.items():
            bound=Path(name).resolve()
            if not bound.is_relative_to(folder) or not bound.is_file() or sha(bound)!=digest:
                raise ValueError('saved phase source/input changed: '+name)
        marker=folder/'user-usage-pause.json'
        if marker.exists():
            for prior in read(marker).get('completed_phases',[]):
                if prior['name']==path.name and prior['sha256']!=sha(path):
                    raise ValueError('saved phase differs from the recorded pause: '+phase)
        saved.append(phase)
        fingerprints[str(path)]=sha(path)
    return dict(folder=folder,steps=steps,backend=backend,saved=saved,fingerprints=fingerprints)

@contextmanager
def retain_launch_handle():
    """Defer cooperative interrupts until Popen's handle has been assigned."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous={};pending=[]
    def defer(signum, frame):
        pending.append(signum)
    try:
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous[signum]=signal.signal(signum,defer)
        yield
    finally:
        for signum,handler in previous.items():signal.signal(signum,handler)
        if pending:raise KeyboardInterrupt('interrupted while launching owned reference child')


def stop_owned(process):
    """Use only our retained Popen handle; unlock requires a successful wait."""
    try:
        if process.poll() is not None:
            process.wait(timeout=5)
            return True
    except BaseException:
        pass
    for action in (process.terminate, process.kill):
        try:
            action()
        except BaseException:
            pass
        try:
            process.wait(timeout=5)
            return True
        except BaseException:
            pass
    return False


def run(folder):
    folder=folder.resolve();lock=folder/'.reference-resume.lock'
    # A live or stale lock is never silently discarded. Check its owner before
    # removing a stale lock after a hard interruption.
    try:handle=lock.open('x',encoding='utf-8')
    except FileExistsError:raise ValueError('reference resume lock exists; verify its owner before retrying')
    with handle:json.dump(dict(pid=os.getpid(),parent_pid=os.getpid(),child_pid=None,interrupted=False),handle)
    process=None;progress=None;attempt=None;interrupted=False
    try:
        cfg=preflight(folder)
        attempt=folder/'resume-attempts'/datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        attempt.mkdir(parents=True,exist_ok=False)
        progress=dict(schema='soil-fem-frozen-reference-resume/1',status='prepared',
          original_progress_sha256=sha(folder/'progress.json'),resume_tool_sha256=sha(Path(__file__)),
          steps=cfg['steps'],backend=cfg['backend'],reused_complete_phases=cfg['saved'],
          initial_source_and_input_sha256=cfg['fingerprints'],phases=[],
          no_midphase_recovery=True,goal_complete=False,p6_complete=False)
        def save():atomic(attempt/'progress.json',progress)
        save();atomic(folder/'latest-resume.json',dict(attempt=str(attempt)))
        for phase in ('initial','cycle'):
            if phase in cfg['saved']:continue
            # Recheck frozen files and complete earlier phases before launch.
            preflight(folder)
            if any(sha(Path(n))!=h for n,h in cfg['fingerprints'].items()):
                raise ValueError('frozen study changed after resume preparation')
            command=[sys.executable,'-X','utf8',str(folder/'reference-sources/mcc_projected_sparse_run.py'),
              phase,*[str(folder/'inputs'/n) for n in INPUTS],str(folder/(phase+'.json')),
              '--steps',str(cfg['steps']),'--dependency-lock',str(folder/'reference-sources/requirements-reference.lock.txt'),
              *BACKENDS[cfg['backend']]]
            if phase=='cycle':command.extend(['--initial',str(folder/'initial.json')])
            progress.update(status=phase+'_running',command=command);save()
            with retain_launch_handle():
                process=subprocess.Popen(command)
            progress.update(child_pid=process.pid,parent_pid=os.getpid());save()
            code=process.wait()
            process=None  # wait confirmed the owned child's termination
            progress['phases'].append(dict(phase=phase,exit_code=code))
            progress.pop('child_pid',None)
            if code:
                progress['status']='interrupted_or_failed';save();return code
            verified=preflight(folder)
            if phase not in verified['saved']:raise ValueError('reference did not save its completed phase')
            progress['status']=phase+'_saved';save()
        progress['status']='completed';save()
        return 0
    except BaseException:
        interrupted=True
        raise
    finally:
        child_pid=None if process is None else process.pid
        confirmed=process is None or stop_owned(process)
        record=dict(pid=os.getpid(),parent_pid=os.getpid(),child_pid=child_pid,
            interrupted=interrupted,child_exit_confirmed=confirmed,
            status='interrupted' if interrupted else 'finished',
            manual_check='Verify both recorded PIDs and child termination before removing this lock.' if not confirmed else None)
        # Write the ownership fallback independently of progress/status writes.
        # If any record write fails, the existing exclusive lock stays in place.
        if not confirmed:
            with lock.open('w',encoding='utf-8') as stream:
                json.dump(record,stream,indent=2);stream.flush();os.fsync(stream.fileno())
        if progress is not None and interrupted:
            progress.update(record)
            try:atomic(attempt/'progress.json',progress)
            except BaseException:
                if not confirmed:
                    pass  # the lock itself carries the recovery instructions
        if confirmed:lock.unlink()

if __name__=='__main__':
    def terminate_parent(signum, frame):
        raise KeyboardInterrupt('reference resume termination requested')
    signal.signal(signal.SIGTERM,terminate_parent)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    args=parser.parse_args()
    try:sys.exit(run(args.study))
    except (ValueError,KeyError,OSError) as error:parser.exit(2,str(error)+'\n')
