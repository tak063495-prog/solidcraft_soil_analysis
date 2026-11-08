"""Bounded observation of the existing 5 kPa NC free-contact case.

This is a diagnostic driver, not a plastic 3D reference or an acceptance tool.
It never changes the production solver or the generated numerical settings.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from make_free_contact import make
import fem_reference as fem
from joint_reference import System
from verify_high_order import internal, physical, prepare


NAMES = ('mesh', 'ground', 'stages')


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)


def expected_inputs():
    return make(mcc=True, pc=100., preload_kpa=5., force=(1., .5, 0.))


def check_inputs(folder, job=False):
    paths = [folder/('case.json' if job and n == 'stages' else n+'.json') for n in NAMES]
    for name, path, expected in zip(NAMES, paths, expected_inputs()):
        if read(path) != expected:
            raise ValueError('5 kPa NC input mismatch: '+name)
    return {p.name: sha(p) for p in paths}


def accepted_count(run):
    return sum(a['accepted'] for a in run.get('geostatic', {}).get('attempts', [])) + sum(
        a['accepted'] for a in run['attempts'])


def deterministic(value):
    """Measured linear-solve time is the only recursively excluded field."""
    if isinstance(value, dict):
        return {k: deterministic(v) for k, v in value.items() if k != 'elapsed_seconds'}
    if isinstance(value, list):
        return [deterministic(v) for v in value]
    return value


def compare_prefix(a, b):
    # The stop mechanism's top-level reason/continuation is intentionally not
    # compared. Every physical state, frame, audit and attempt is compared.
    keys = ('last_committed', 'frames', 'excavation_audits', 'transition_audits',
            'reference_displacement_m', 'attempts', 'geostatic')
    for key in keys:
        if deterministic(a.get(key)) != deterministic(b.get(key)):
            raise ValueError('restart/prefix differs: '+key)
    return dict(exact=True, excluded_recursive_fields=['elapsed_seconds'], compared_fields=list(keys))


def snapshot(job):
    status = read(job/'status.json')
    if status.get('artifact'):
        path = job/status['artifact']
        if sha(path) != status['artifact_sha256']:
            raise ValueError('job artifact checksum mismatch')
    elif (job/'recovery.json').exists():
        pointer = read(job/'recovery.json')
        name = pointer['checkpoint']
        if name not in ('recovery-checkpoint-0.json', 'recovery-checkpoint-1.json'):
            raise ValueError('unexpected recovery path')
        path = job/name
        if sha(path) != pointer['checkpoint_sha256']:
            raise ValueError('recovery checksum mismatch')
    elif (job/'inputs/checkpoint.json').exists():
        # A resumed run can spend its whole budget in one unaccepted trial.
        # Its immutable input checkpoint is still the latest known state.
        path = job/'inputs/checkpoint.json'
        expected = read(job/'job.json')['binding']['files_sha256']['inputs/checkpoint.json']
        if sha(path) != expected:
            raise ValueError('input checkpoint checksum mismatch')
    else:
        return None, None
    data = read(path)
    return path.relative_to(job).as_posix(), data.get('run', data)


def force_check(mesh, ground, case, state, q, offsets, system):
    active = np.array(state['active_elements'])
    bulk, yf = internal(mesh, ground, state['material_states'], active, q, offsets)
    force = bulk + system.audit(state)
    load = np.asarray(state['applied_nodal_load_kn']).reshape(-1, 3)
    basis, error = fem.staged_constraint_space(mesh, case, state)
    residual = np.concatenate([b.T@(force[n]-load[n]) for n, b in basis.items()])
    relative = float(np.linalg.norm(residual)/max(np.linalg.norm(load), 1.))
    if relative >= case['settings']['relative_force_tolerance'] or error >= 1e-10:
        raise ValueError('accepted-state force/constraint check failed')
    return dict(relative_force_residual=relative, constraint_error_m=error,
                maximum_normalized_mcc_yield=yf)


def inspect_run(mesh, ground, case, run):
    """Reassemble stored-stress forces; do not call this a displacement reference."""
    q, offsets = prepare(mesh)
    system = System(mesh, case['joint_system'])
    state = run['last_committed']
    attempts = run['geostatic']['attempts'] + run['attempts']
    if any(not a['accepted'] and not a['committed_state_unchanged_on_rejection'] for a in attempts):
        raise ValueError('rejection changed committed state')
    history = state['material_states']
    rows = []
    for frame in run['frames']:
        rows.append(dict(frame=frame['name'], **force_check(mesh, ground, case, frame['state'], q, offsets, system)))
    if run['frames']:
        current = force_check(mesh, ground, case, state, q, offsets, system)
    else:
        current = None  # geostatic seed/partial holding is not a physical endpoint
    report = dict(accepted_increments=accepted_count(run),
        completed_frames=[f['name'] for f in run['frames']],
        stage_index=state['stage_index'], stage_progress=state['stage_progress'],
        geostatic_progress=run['geostatic']['progress'],
        pc_range_kpa=[min(p['preconsolidation_pressure_kpa'] for p in history),
                      max(p['preconsolidation_pressure_kpa'] for p in history)],
        plastic_points=sum(any(x != 0. for x in p['plastic_strain_compression_positive']) for p in history),
        rejections=sum(not a['accepted'] for a in attempts),
        rejected_plastic_trial_changed=any(a.get('rejected_trial_plastic_changed', False) for a in attempts if not a['accepted']),
        rejected_contact_trial_changed=any(a.get('rejected_trial_structure_changed', False) for a in attempts if not a['accepted']),
        all_rejections_preserve_committed_state=True,
        stored_stress_force_checks=rows, last_committed_force_check=current,
        independent_plastic_3d_displacement=False,
        zero_release_check=None, fully_released_and_frozen=False)
    if run['excavation_audits']:
        origin = run['frames'][2]['state']
        active = np.array(state['active_elements'])
        removed = np.array(origin['active_elements']) & ~active
        removed_force, _ = internal(mesh, ground, origin['material_states'], removed, q, offsets)
        endload, _ = physical(mesh, ground, active, 0., True, joint_system=system)
        retained = {**origin, 'active_elements': active.tolist()}
        holding = (np.asarray(origin['applied_nodal_load_kn']).reshape(-1, 3)-endload-removed_force
                   -system.audit(origin)+system.audit(retained))
        used = np.unique(fem.cells(mesh)[active])
        holding[np.setdiff1d(np.arange(len(mesh['nodes'])), used)] = 0.
        audit = run['excavation_audits'][0]
        if not np.allclose(holding.ravel(), audit['initial_holding_load_kn'], rtol=1e-9, atol=1e-7):
            raise ValueError('lambda=0 holding load mismatch')
        retained = {**retained, 'applied_nodal_load_kn': (endload+holding).ravel().tolist()}
        report['zero_release_check'] = force_check(mesh, ground, case, retained, q, offsets, system)
        if audit['zero_release_relative_force_residual'] >= case['settings']['relative_force_tolerance']:
            raise ValueError('native lambda=0 residual exceeds unchanged tolerance')
        if state['joints'] != origin['joints']:
            raise ValueError('removed joint history changed')
        for ei in np.flatnonzero(removed):
            if history[offsets[ei]:offsets[ei+1]] != origin['material_states'][offsets[ei]:offsets[ei+1]]:
                raise ValueError('removed material history changed')
        report['removed_joint_and_material_history_frozen'] = True
        if run['complete']:
            if np.linalg.norm(state['excavation_release_load_kn']) >= 1e-9:
                raise ValueError('final holding load not released')
            report['fully_released_and_frozen'] = True
    return report


def invoke(command, log):
    with log.open('x', encoding='utf-8') as stream:
        result = subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT,
                                timeout=30, check=False)
    if result.returncode:
        raise ValueError('control command failed: '+str(log))


def probe(solver, output, budget_seconds=60., grace_seconds=15., resume_from=None, stop_after_accepted=None):
    if not math.isfinite(budget_seconds) or not 0 < budget_seconds <= 600:
        raise ValueError('budget must be finite, positive and <= 600 seconds')
    if not math.isfinite(grace_seconds) or not 0 < grace_seconds <= 60:
        raise ValueError('grace must be finite, positive and <= 60 seconds')
    if stop_after_accepted is not None and stop_after_accepted < 1:
        raise ValueError('accepted count must be positive')
    if resume_from:
        check_inputs(resume_from/'inputs', job=True)
    output.mkdir(parents=True, exist_ok=False)
    inputs = output/'input'; inputs.mkdir()
    for name, value in zip(NAMES, expected_inputs()):
        write(inputs/(name+'.json'), value)
    fingerprints = check_inputs(inputs)
    solver_sha = sha(solver)
    job = output/'job'
    resume_command = None
    if resume_from:
        # A cooperative checkpoint can be newer than the rotating recovery
        # pointer. Use that final sealed artifact whenever it is available.
        resume_command = ('prepare-stages-resume' if read(resume_from/'status.json')['state'] == 'checkpointed'
                          else 'prepare-stages-recovery')
        command = [solver, resume_command, resume_from, job]
    else:
        command = [solver, 'prepare-durable-stages-job', *[inputs/(n+'.json') for n in NAMES], 1, job]
    invoke(command, output/'prepare.log')
    frozen = job/'solver/soil-fem-3d.exe'
    if sha(frozen) != solver_sha:
        raise ValueError('resume requires the identical solver binary')
    start = time.monotonic(); reason = None; requested = None; terminated = False
    events = []; position = 0; pending = ''; event_file = job/'progress.jsonl'
    with (output/'solve.log').open('x', encoding='utf-8') as log:
        child = subprocess.Popen([str(frozen), 'run-stages-job', str(job)], stdout=log, stderr=subprocess.STDOUT)
        try:
            while child.poll() is None:
                if event_file.exists():
                    with event_file.open(encoding='utf-8') as stream:
                        stream.seek(position); pending += stream.read(); position = stream.tell()
                    lines = pending.split('\n'); pending = lines.pop()
                    events.extend(json.loads(line) for line in lines if line)
                accepted = sum(e['progress']['kind'] == 'accepted' for e in events)
                elapsed = time.monotonic()-start
                if requested is None and (elapsed >= budget_seconds or
                        (stop_after_accepted is not None and accepted >= stop_after_accepted)):
                    reason = 'wall_time_budget' if elapsed >= budget_seconds else 'accepted_observation_limit'
                    invoke([frozen, 'stop-stages-job', job], output/'stop.log')
                    requested = time.monotonic()
                if requested is not None and time.monotonic()-requested >= grace_seconds:
                    child.terminate(); terminated = True; break
                time.sleep(.02)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                child.kill(); child.wait(timeout=5); terminated = True
        finally:
            if child.poll() is None:
                child.kill(); child.wait(timeout=5)
    name, run = snapshot(job)
    status = read(job/'status.json')
    # A stale status after termination is retained verbatim; it is not liveness.
    report = dict(schema='soil-fem-low-preload-nc-probe/1', passed=False,
        classification='未検証', p3b04_complete=False, p6_complete=False,
        solver_complete=False if run is None else run['complete'],
        stop_reason=reason, owned_child_terminated=terminated, exit_code=child.returncode,
        elapsed_seconds=time.monotonic()-start, job_state=status['state'], snapshot=name,
        resume_command=resume_command,
        solver_sha256=solver_sha, input_sha256=fingerprints,
        input_condition=dict(seed_p_kpa=100., seed_pc_kpa=100., preload_kpa=5., force_kn=[1., .5, 0.]),
        missing_acceptance_evidence=['独立塑性3D全場', '両実装の刻み感度', '全工程途中再開の連続経路一致'],
        observation=None)
    if run is not None:
        report['observation'] = inspect_run(*expected_inputs(), run)
        write(output/'observed-run.json', run)
        if not run['complete']:
            report['missing_acceptance_evidence'].insert(0, '全上部撤去の完了・全解放')
    else:
        report['missing_acceptance_evidence'].insert(0, '受理済み保存点')
    if sha(solver) != solver_sha or check_inputs(inputs) != fingerprints:
        raise ValueError('solver/input changed during observation')
    report['evidence_sha256'] = {p.relative_to(output).as_posix(): sha(p) for p in output.rglob('*')
                                 if p.is_file() and p != frozen}
    write(output/'probe.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--budget-seconds', type=float, default=60.)
    parser.add_argument('--grace-seconds', type=float, default=15.)
    parser.add_argument('--resume-from', type=Path)
    parser.add_argument('--stop-after-accepted', type=int)
    args = parser.parse_args()
    result = probe(args.solver.resolve(), args.output.resolve(), args.budget_seconds,
                   args.grace_seconds, args.resume_from.resolve() if args.resume_from else None,
                   args.stop_after_accepted)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    # Successful observation is deliberately distinct from numerical acceptance.
