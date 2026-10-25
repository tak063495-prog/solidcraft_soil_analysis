"""Real CAD worker pause/resume comparison using owned process handles."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from freeze_sources import freeze, manifest
from import_solidcraft_capture import read, validate
from native_job import request_cancel, stop_child
from native_workflow_snapshot import sha
from verify_default_solver import deterministic


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def verify(capture, solver, output):
    root = Path(__file__).resolve().parents[1]
    captured, _, _, _ = validate(capture)
    output.mkdir(parents=True, exist_ok=False)
    freeze(output)
    sources = manifest()
    (output/'solver').mkdir()
    pinned = output/'solver/soil-fem-3d.exe'
    shutil.copyfile(solver, pinned)
    digest = sha(pinned)
    assert digest == sha(solver)
    def launch(name, previous=None, watch=False):
        folder = output/name
        args = ([sys.executable, '-X', 'utf8', str(root/'tools/native_job.py'), 'resume', str(previous), str(folder)]
                if previous else [sys.executable, '-X', 'utf8', str(root/'tools/native_job.py'), 'run', str(capture), str(folder),
                                  '--mesh-size-m', '4', '--order', '2', '--quadrature', '3', '--cooperative'])
        args += ['--solver', str(pinned), '--expected-input-sha256', captured['input_sha256'], '--host-token', name]
        if watch:
            args.append('--watch-stdin')
        log = (output/(name+'.log')).open('x', encoding='utf-8')
        child = subprocess.Popen(args, cwd=root, stdin=subprocess.PIPE if watch else subprocess.DEVNULL,
                                 stdout=log, stderr=subprocess.STDOUT)
        return folder, child, log
    def finish(child, log, expected, deadline=180):
        try:
            code = child.wait(timeout=deadline)
            assert code == expected, f'worker exit {code}, expected {expected}; inspect {log.name}'
            log.write('\nExit code: '+str(code)+'\n')
            log.flush()
        finally:
            if child.poll() is None:
                stop_child(child)
            if child.stdin is not None:
                child.stdin.close()
            log.close()
    continuous, child, log = launch('continuous', watch=True)
    finish(child, log, 0)
    result = read(continuous/'workflow/analysis/stage-results.json', 256*1024*1024)
    assert result['complete'] and read(continuous/'status.json')['state'] == 'succeeded'
    print('Continuous cooperative profile completed with audited native view', flush=True)
    observations = []
    for label, host_eof, postprocessing in (('requested-stop', False, False), ('host-disconnect', True, False), ('postprocessing-stop', False, True)):
        saved, child, log = launch(label, watch=host_eof)
        started = time.monotonic()
        observed = None
        try:
            while time.monotonic()-started < 120:
                if child.poll() is not None:
                    raise AssertionError('worker ended before a live solver observation')
                path = saved/('status.json' if postprocessing else 'workflow/analysis/solver-job/status.json')
                if path.exists():
                    report = read(path, 256*1024)
                    progress = None
                    events = path.parent/'progress.jsonl'
                    if not postprocessing and events.exists():
                        for line in events.read_text(encoding='utf-8').splitlines():
                            try:
                                event = json.loads(line)
                            except json.JSONDecodeError:
                                continue
                            point = event.get('progress') or {}
                            if event.get('job_id') == report.get('job_id') and point.get('kind') == 'newton' and point.get('stage_index', 0) >= 1:
                                progress = event
                                break
                    if report.get('state') == 'running' and ((postprocessing and report.get('phase') == 'exporting') or
                        (not postprocessing and progress is not None)):
                        observed = dict(owned_child_poll=child.poll(), solver_report=report, newton_event=progress)
                        break
                time.sleep(.02)
            assert observed is not None and observed['owned_child_poll'] is None
            if host_eof:
                child.stdin.close()
            else:
                assert request_cancel(saved)['requested']
            finish(child, log, 3)
        finally:
            if child.poll() is None:
                stop_child(child)
            if not log.closed:
                log.close()
        status = read(saved/'status.json')
        receipt = read(saved/'checkpoint-receipt.json')
        assert status['state'] == 'checkpointed' and status['resume_available'] is True
        assert receipt['kind'] == ('completed_solver_result' if postprocessing else 'solver_checkpoint')
        assert (saved/'workflow/analysis/stage-results.json').exists() is postprocessing
        assert not request_cancel(saved)['requested']
        write(output/(label+'-observation.json'), observed)
        resumed, child, log = launch(label+'-resumed', previous=saved, watch=True)
        finish(child, log, 0)
        owner = read(resumed/'job.json')
        assert owner['parent_job_id'] == status['job_id'] and owner['job_id'] != status['job_id']
        actual = read(resumed/'workflow/analysis/stage-results.json', 256*1024*1024)
        assert json.dumps(deterministic(actual), sort_keys=True, allow_nan=False) == json.dumps(deterministic(result), sort_keys=True, allow_nan=False)
        assert read(resumed/'status.json')['state'] == 'succeeded'
        observations.append(dict(trigger=label, saved_job=status['job_id'], resumed_job=owner['job_id'],
                                 nonruntime_result_exact=True, checkpoint_kind=receipt['kind']))
        print(label+': saved '+receipt['kind']+', separate-process resume and complete result exact', flush=True)
    assert manifest() == sources
    write(output/'verification.json', dict(schema='soil-fem-native-cooperative-verification/1', passed=True,
        actual_native_capture=True, input_sha256=captured['input_sha256'], solver_sha256=digest,
        continuous_completion=True, checkpoint_trials=observations, native_ui_connected=False,
        plastic_history_required=False, goal_complete=False, p6_complete=False,
        source_sha256=sources,
        evidence_sha256={p.relative_to(output).as_posix(): sha(p) for p in output.rglob('*') if p.is_file()}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('capture', 'solver', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    verify(args.capture.resolve(), args.solver.resolve(), args.output.resolve())
