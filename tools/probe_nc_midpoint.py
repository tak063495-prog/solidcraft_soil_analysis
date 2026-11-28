"""Bounded NC experiment with the EXISTING midpoint integration option.

No solver/source patch, tolerance change, pivot change or checkpoint resealing.
The numerical integration scheme changes, not the MCC parameters or loading.
This diagnostic cannot grant full-removal or independent-field acceptance.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time
from probe_low_preload_nc import NAMES, expected_inputs, inspect_run, invoke, read, sha, snapshot, write


def inputs(scheme):
    assert scheme in ('backward_euler', 'midpoint')
    values = expected_inputs(excavation_increment=1/128, minimum_increment=1/4194304,
                             maximum_newton_iterations=60, dense_lu=True, tangent_perturbation=1e-9)
    if scheme == 'midpoint':
        values[2]['settings']['mcc_integration'] = scheme
    return values


def resume_inputs(job, scheme, binary):
    """Check physical/numerical inputs before the CLI verifies every seal."""
    for name, expected in zip(NAMES, inputs(scheme)):
        path = job/'inputs'/('case.json' if name == 'stages' else name+'.json')
        if read(path) != expected:
            raise ValueError('5 kPa NC integration input mismatch: '+name)
    binding = read(job/'job.json')['binding']
    frozen = job/'solver/soil-fem-3d.exe'
    if binding['identity']['solver_binary_sha256'] != binary or sha(frozen) != binary:
        raise ValueError('resume requires the identical sealed solver binary')


def restart_observation(before, native):
    from probe_low_preload_nc import accepted_count, deterministic
    # Optional audit arrays are omitted by older serialized Run schemas.
    for key in ('attempts', 'frames', 'excavation_audits', 'transition_audits'):
        old = before.get(key, [])
        if deterministic(native.get(key, [])[:len(old)]) != deterministic(old):
            raise ValueError('resumed history prefix changed: '+key)
    if deterministic(native['geostatic']) != deterministic(before['geostatic']):
        raise ValueError('resumed geostatic history changed')
    return dict(accepted_before=accepted_count(before), accepted_after=accepted_count(native),
                new_accepted=accepted_count(native)-accepted_count(before),
                starting_progress=before['last_committed']['stage_progress'],
                preserved_attempt_and_frame_prefix=True, uninterrupted_full_path_compared=False)


def run(solver, output, scheme='midpoint', budget=180., resume_from=None):
    assert 0 < budget <= 180
    binary = sha(solver)
    before = None; resume_command = None
    if resume_from is not None:
        resume_inputs(resume_from, scheme, binary)
        _, before = snapshot(resume_from)
        if before is None or before['complete']:
            raise ValueError('resume requires an incomplete sealed snapshot')
        resume_command = ('prepare-stages-resume' if read(resume_from/'status.json')['state'] == 'checkpointed'
                          else 'prepare-stages-recovery')
    output.mkdir(parents=True, exist_ok=False); folder = output/'input'; folder.mkdir()
    values = inputs(scheme)
    for name, value in zip(NAMES, values):
        write(folder/(name+'.json'), value)
    fingerprints = {p.name: sha(p) for p in folder.iterdir()}
    job = output/'job'
    command = ([solver, resume_command, resume_from, job] if resume_from is not None else
               [solver, 'prepare-durable-stages-job', *[folder/(n+'.json') for n in NAMES], 1, job])
    invoke(command, output/'prepare.log')
    frozen = job/'solver/soil-fem-3d.exe'; assert sha(frozen) == binary
    if resume_from is not None:
        # A failed result may include rejected trials newer than the rotating
        # recovery checkpoint. Compare with the sealed state actually imported.
        before = read(job/'inputs/checkpoint.json')['run']
    started = time.monotonic(); requested = None; terminated = False
    with (output/'solve.log').open('x') as log:
        child = subprocess.Popen([str(frozen), 'run-stages-job', str(job)], stdout=log, stderr=subprocess.STDOUT)
        try:
            while child.poll() is None:
                if requested is None and time.monotonic()-started >= budget:
                    invoke([frozen, 'stop-stages-job', job], output/'stop.log'); requested = time.monotonic()
                if requested is not None and time.monotonic()-requested >= 20:
                    child.terminate(); terminated = True; break
                time.sleep(.05)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                child.kill(); child.wait(timeout=5); terminated = True
        finally:
            if child.poll() is None:
                child.kill(); child.wait(timeout=5); terminated = True
    name, native = snapshot(job)
    assert sha(solver) == binary and all(sha(folder/n) == h for n, h in fingerprints.items())
    # Persist process facts before postprocessing, so an audit error never
    # destroys the measured runtime or requires an expensive solver rerun.
    process = dict(budget_seconds=budget, stop_requested=requested is not None, owned_child_terminated=terminated,
                   elapsed_seconds=time.monotonic()-started, exit_code=child.returncode,
                   stop_reason='wall_time_budget' if requested is not None else 'solver_returned',
                   job_state=read(job/'status.json')['state'], resume_command=resume_command)
    write(output/'process-observation.json', process)
    report = dict(schema='soil-fem-nc-existing-midpoint-experiment/1', passed=False, production_adopted=False,
                  scheme=scheme, source_patched=False, solver_sha256=binary, input_sha256=fingerprints,
                  **process,
                  snapshot=name, complete=False if native is None else native['complete'],
                  input_condition=dict(pc_kpa=100., preload_kpa=5., force_kn=[1., .5, 0.]),
                  missing_acceptance=['full removal and post-release freeze', 'independent plastic 3D full path',
                                      'both implementations step sensitivity', 'restart through full removal'])
    if native is not None:
        write(output/'observed-run.json', native)
        report.update(failure_reason=native['failure_reason'], observation=inspect_run(*values, native),
                      last_attempt=native['attempts'][-1] if native['attempts'] else None)
    if before is not None:
        # Preserve and check all old attempts/frames; do not claim a direct
        # uninterrupted-path comparison from a single resumed observation.
        report['restart_observation'] = restart_observation(before, native)
        report['restart_observation']['imported_checkpoint_sha256'] = sha(job/'inputs/checkpoint.json')
    write(output/'experiment.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scheme', choices=('backward_euler', 'midpoint'), default='midpoint')
    parser.add_argument('--budget-seconds', type=float, default=180.)
    parser.add_argument('--resume-from', type=Path, help='finished job with the identical sealed binary and inputs')
    args = parser.parse_args()
    result = run(args.solver.resolve(), args.output.resolve(), args.scheme, args.budget_seconds,
                 args.resume_from.resolve() if args.resume_from else None)
    print(json.dumps({k: result.get(k) for k in ('passed', 'complete', 'failure_reason', 'elapsed_seconds')}, ensure_ascii=False, indent=2))
