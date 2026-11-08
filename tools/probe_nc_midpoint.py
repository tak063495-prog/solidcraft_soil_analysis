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


def run(solver, output, scheme='midpoint', budget=180.):
    assert 0 < budget <= 180
    output.mkdir(parents=True, exist_ok=False); folder = output/'input'; folder.mkdir()
    values = inputs(scheme)
    for name, value in zip(NAMES, values):
        write(folder/(name+'.json'), value)
    fingerprints = {p.name: sha(p) for p in folder.iterdir()}; binary = sha(solver)
    job = output/'job'
    invoke([solver, 'prepare-durable-stages-job', *[folder/(n+'.json') for n in NAMES], 1, job], output/'prepare.log')
    frozen = job/'solver/soil-fem-3d.exe'; assert sha(frozen) == binary
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
            child.wait(timeout=5)
        finally:
            if child.poll() is None:
                child.kill(); child.wait(timeout=5); terminated = True
    name, native = snapshot(job)
    assert sha(solver) == binary and all(sha(folder/n) == h for n, h in fingerprints.items())
    report = dict(schema='soil-fem-nc-existing-midpoint-experiment/1', passed=False, production_adopted=False,
                  scheme=scheme, source_patched=False, solver_sha256=binary, input_sha256=fingerprints,
                  budget_seconds=budget, stop_requested=requested is not None, owned_child_terminated=terminated,
                  elapsed_seconds=time.monotonic()-started, exit_code=child.returncode,
                  snapshot=name, complete=False if native is None else native['complete'],
                  input_condition=dict(pc_kpa=100., preload_kpa=5., force_kn=[1., .5, 0.]),
                  missing_acceptance=['full removal and post-release freeze', 'independent plastic 3D full path',
                                      'both implementations step sensitivity', 'restart through full removal'])
    if native is not None:
        write(output/'observed-run.json', native)
        report.update(failure_reason=native['failure_reason'], observation=inspect_run(*values, native),
                      last_attempt=native['attempts'][-1] if native['attempts'] else None)
    write(output/'experiment.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scheme', choices=('backward_euler', 'midpoint'), default='midpoint')
    parser.add_argument('--budget-seconds', type=float, default=180.)
    args = parser.parse_args()
    result = run(args.solver.resolve(), args.output.resolve(), args.scheme, args.budget_seconds)
    print(json.dumps({k: result.get(k) for k in ('passed', 'complete', 'failure_reason', 'elapsed_seconds')}, ensure_ascii=False, indent=2))
