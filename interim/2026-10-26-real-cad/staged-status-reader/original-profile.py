"""Owned staged solver/checkpoint profile for the external SolidCraft worker.

Executable selection remains a trusted host argument. Saved jobs only provide
bound data; they never select programs or a different workflow implementation.
"""
import json
import os
from pathlib import Path
import queue
import subprocess
import sys
import threading
import time
import uuid
from import_solidcraft_capture import read, validate
from native_workflow_snapshot import sha, verify, atomic_json

INPUT_FILES = ('workflow/analysis/ground/mesh.json', 'workflow/analysis/ground/ground.json',
               'workflow/analysis/ground/model-input-package.json', 'workflow/analysis/stages.json',
               'workflow/analysis/prepared-stage-input.json', 'workflow/conformal/prepared-model.json')


def atomic(path, value):
    atomic_json(path, value)


def copy_file(source, target, digest):
    if sha(source) != digest:
        raise ValueError('saved input changed: '+source.name)
    with source.open('rb') as incoming, target.open('xb') as outgoing:
        for block in iter(lambda: incoming.read(1024*1024), b''):
            outgoing.write(block)
    if sha(source) != digest or sha(target) != digest:
        raise ValueError('saved input changed during copy: '+source.name)


def plain(name):
    if not isinstance(name, str) or not name or Path(name).name != name or name in ('.', '..') or any(c in name for c in ('/', '\\', ':')):
        raise ValueError('saved inventory requires plain filenames')
    return name


def previous(folder, solver_digest, sources, expected_input=None):
    owner, status = read(folder/'job.json', 256*1024), read(folder/'status.json', 256*1024)
    if owner.get('schema') != 'soil-fem-native-job/2' or status.get('schema') != 'soil-fem-native-job-status/2':
        raise ValueError('cooperative native job required for resume')
    if status.get('state') != 'checkpointed' or status.get('job_id') != owner.get('job_id'):
        raise ValueError('saved native job identity/state differs')
    if owner.get('solver_sha256') != solver_digest or owner.get('tool_sources_sha256') != sources:
        raise ValueError('resume requires the same installed solver and workflow version')
    captured, _, _, _ = validate(folder/'capture')
    if captured['input_sha256'] != owner.get('input_sha256') or (expected_input is not None and expected_input != captured['input_sha256']):
        raise ValueError('saved job differs from the current input snapshot')
    receipt_path = folder/'checkpoint-receipt.json'
    if status.get('checkpoint_receipt_sha256') != sha(receipt_path):
        raise ValueError('saved checkpoint receipt changed')
    receipt = read(receipt_path, 256*1024)
    if receipt.get('schema') != 'soil-fem-native-checkpoint/1' or receipt.get('job_id') != owner['job_id'] or receipt.get('input_sha256') != captured['input_sha256'] or receipt.get('solver_sha256') != solver_digest:
        raise ValueError('saved checkpoint receipt identity differs')
    kind = receipt.get('kind')
    artifact = 'workflow/analysis/solver-job/checkpoint.json' if kind == 'solver_checkpoint' else 'workflow/analysis/stage-results.json' if kind == 'completed_solver_result' else None
    if artifact is None or receipt.get('artifact') != artifact or sha(folder/artifact) != receipt.get('artifact_sha256'):
        raise ValueError('saved checkpoint/result artifact changed')
    if set(receipt.get('files_sha256', {})) != set(INPUT_FILES):
        raise ValueError('saved input inventory differs')
    for name in INPUT_FILES:
        if sha(folder/name) != receipt['files_sha256'][name]:
            raise ValueError('saved analysis input changed: '+name)
    return owner, receipt


def clone_context(previous_folder, folder, solver, owner):
    workflow = folder/'workflow'
    conformal = workflow/'conformal'
    conformal.mkdir(parents=True)
    old = previous_folder/'workflow/conformal'
    inventory = read(old/'prepared-model.json', 256*1024)
    if inventory.get('schema') != 'soil-fem-prepared-model/1':
        raise ValueError('saved geometry inventory missing')
    for name, digest in inventory['files_sha256'].items():
        copy_file(old/plain(name), conformal/name, digest)
    copy_file(old/'prepared-model.json', conformal/'prepared-model.json', sha(old/'prepared-model.json'))
    analysis = workflow/'analysis'
    ground = analysis/'ground'
    ground.mkdir(parents=True)
    old_ground = previous_folder/'workflow/analysis/ground'
    imported = read(old_ground/'model-input-package.json', 256*1024)
    if imported.get('schema') != 'soil-fem-prepared-model-import/1':
        raise ValueError('saved prepared-input receipt missing')
    for name, digest in imported['original_inventory']['files_sha256'].items():
        actual = 'prepared-ground.json' if plain(name) == 'ground.json' else name
        copy_file(old_ground/actual, ground/actual, digest)
    copy_file(old_ground/'ground.json', ground/'ground.json', imported['ground_sha256'])
    copy_file(old_ground/'model-input-package.json', ground/'model-input-package.json', sha(old_ground/'model-input-package.json'))
    for name in ('stages.json', 'prepared-stage-input.json'):
        source = previous_folder/'workflow/analysis'/name
        copy_file(source, analysis/name, sha(source))
    (analysis/'solver').mkdir()
    copy_file(solver, analysis/'solver/soil-fem-3d.exe', owner['solver_sha256'])


def command(args, runtime, log_path, publish, cancelled, stop_child, *, phase, solver_job=None, solver=None):
    lines = queue.Queue(maxsize=256)
    options = dict(creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW) if os.name == 'nt' else dict(start_new_session=True)
    child = None
    with log_path.open('x', encoding='utf-8') as log:
        try:
            child = subprocess.Popen(args, cwd=runtime, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, encoding='utf-8', errors='replace', env={**os.environ, 'PYTHONUTF8':'1'}, **options)
            def collect():
                try:
                    for line in child.stdout:
                        log.write(line)
                        log.flush()
                        lines.put(line)
                except Exception as error:
                    lines.put(('reader_error', str(error)))
                finally:
                    lines.put(None)
            reader = threading.Thread(target=collect, daemon=True)
            reader.start()
            publish(state='running', phase=phase, child_pid=child.pid)
            eof, stopped, last_write, tail = False, False, 0., []
            while child.poll() is None or not eof:
                if cancelled() and not stopped:
                    stopped = True
                    publish(state='cancelling', phase='checkpointing' if solver_job else phase)
                    if solver_job:
                        request = subprocess.run([str(solver), 'stop-stages-job', str(solver_job)], capture_output=True,
                                                 text=True, encoding='utf-8', **(dict(creationflags=subprocess.CREATE_NO_WINDOW) if os.name == 'nt' else {}))
                        if request.returncode:
                            raise ValueError('cooperative solver stop request failed: '+request.stderr[-1024:])
                    else:
                        stop_child(child)
                try:
                    line = lines.get(timeout=.1)
                    if line is None:
                        eof = True
                    elif isinstance(line, tuple):
                        raise RuntimeError(line[1])
                    else:
                        tail = (tail+[line.rstrip()[-1024:]])[-30:]
                except queue.Empty:
                    pass
                if time.monotonic()-last_write > .3:
                    extra = {}
                    if solver_job and (solver_job/'status.json').is_file():
                        status = read(solver_job/'status.json', 256*1024)
                        if status.get('schema') != 'soil-fem-stages-job-status/1':
                            raise ValueError('unsupported solver progress schema')
                        extra['solver_report'] = status
                    publish(phase='checkpointing' if stopped and solver_job else phase, tail=tail, **extra)
                    last_write = time.monotonic()
            reader.join(timeout=5)
            code = child.wait()
            child.stdout.close()
            return code, stopped
        except BaseException:
            if child is not None and child.poll() is None:
                stop_child(child)
            raise


def saved(folder, owner, kind, artifact, publish, solver_report=None):
    relative = artifact.relative_to(folder).as_posix()
    receipt = dict(schema='soil-fem-native-checkpoint/1', job_id=owner['job_id'], input_sha256=owner['input_sha256'],
                   solver_sha256=owner['solver_sha256'], kind=kind, artifact=relative,
                   artifact_sha256=sha(artifact), artifact_bytes=artifact.stat().st_size,
                   files_sha256={name: sha(folder/name) for name in INPUT_FILES})
    atomic(folder/'checkpoint-receipt.json', receipt)
    publish(state='checkpointed', phase='saved', exit_code=3, resume_available=True,
            checkpoint_receipt='checkpoint-receipt.json', checkpoint_receipt_sha256=sha(folder/'checkpoint-receipt.json'),
            solver_report=solver_report, stop_policy='cooperative_staged_checkpoint')
    return 3


def execute(folder, runtime, solver, owner, publish, cancelled, stop_child, *, previous_folder=None, receipt=None):
    analysis = folder/'workflow/analysis'
    mesh, ground, case = analysis/'ground/mesh.json', analysis/'ground/ground.json', analysis/'stages.json'
    solver_job = analysis/'solver-job'
    complete = analysis/'stage-results.json'
    if previous_folder is not None and receipt['kind'] == 'completed_solver_result':
        copy_file(previous_folder/'workflow/analysis/stage-results.json', complete, receipt['artifact_sha256'])
    else:
        args = ([str(solver), 'prepare-stages-resume', str(previous_folder/'workflow/analysis/solver-job'), str(solver_job)]
                if previous_folder is not None else [str(solver), 'prepare-stages-job', str(mesh), str(ground), str(case), str(solver_job)])
        code, stopped = command(args, runtime, folder/'solver-prepare.log', publish, cancelled, stop_child, phase='preparing_solver')
        if stopped:
            publish(state='cancelled', phase='preparing_solver', exit_code=2)
            return 2
        if code:
            raise ValueError('installed solver could not prepare a cooperative staged job; inspect solver-prepare.log')
        pinned = solver_job/'solver/soil-fem-3d.exe'
        if sha(pinned) != owner['solver_sha256']:
            raise ValueError('staged job solver fingerprint differs')
        code, stopped = command([str(pinned), 'run-stages-job', str(solver_job)], runtime, folder/'solver-run.log',
                                publish, cancelled, stop_child, phase='solving', solver_job=solver_job, solver=pinned)
        status = read(solver_job/'status.json', 256*1024)
        manifest = read(solver_job/'job.json', 256*1024)
        if status.get('job_id') != manifest.get('job_id'):
            raise ValueError('solver status and job identity differ')
        if not ((status.get('state') == 'checkpointed' and code == 2) or (status.get('state') == 'succeeded' and code == 0)):
            raise ValueError('staged solver did not complete; inspect its retained status/incomplete result')
        artifact = solver_job/('checkpoint.json' if status.get('state') == 'checkpointed' else 'result.json')
        if status.get('artifact_sha256') != sha(artifact):
            raise ValueError('staged solver artifact changed')
        if status.get('state') == 'checkpointed' and code == 2:
            verify(runtime, owner['tool_sources_sha256'])
            return saved(folder, owner, 'solver_checkpoint', artifact, publish, status)
        if code or status.get('state') != 'succeeded':
            raise ValueError('staged solver did not complete; inspect its retained status/incomplete result')
        copy_file(artifact, complete, status['artifact_sha256'])
    if cancelled():
        return saved(folder, owner, 'completed_solver_result', complete, publish)
    commands = [
        ('exporting', [sys.executable, str(runtime/'tools/export_stages.py'), str(mesh), str(complete), str(analysis/'results')]),
        ('checking_result', [sys.executable, str(runtime/'tools/verify_stage_package.py'), str(mesh), str(ground), str(complete), str(analysis/'results'), '--solver', str(analysis/'solver/soil-fem-3d.exe')]),
        ('preparing_view', [sys.executable, str(runtime/'tools/export_solidcraft_view.py'), str(folder/'workflow'), str(folder/'workflow/view')]),
        ('checking_view', [sys.executable, str(runtime/'tools/verify_solidcraft_view.py'), str(folder/'workflow'), str(folder/'workflow/view'), str(folder/'workflow/view-verification.json')])]
    for phase, args in commands:
        code, stopped = command(args, runtime, folder/(phase+'.log'), publish, cancelled, stop_child, phase=phase)
        if stopped:
            return saved(folder, owner, 'completed_solver_result', complete, publish)
        if code:
            raise ValueError('native staged postprocessing failed: '+phase)
    verify(runtime, owner['tool_sources_sha256'])
    view = folder/'workflow/view/view-manifest.json'
    audit = folder/'workflow/view-verification.json'
    value, checked = read(view), read(audit)
    if not value.get('complete') or value.get('input_sha256') != owner['input_sha256'] or value.get('solver_binary_sha256') != owner['solver_sha256'] or checked.get('passed') is not True or checked.get('view_manifest_sha256') != sha(view):
        raise ValueError('completed native view differs from captured input or solver')
    publish(state='succeeded', phase='complete', exit_code=0, resume_available=False,
            result='workflow/view/view-manifest.json', result_sha256=sha(view), view_audit_sha256=sha(audit))
    return 0
