"""Existing fixed-history multistart Newton diagnostic at new saved anchors.

Same force/linear/root-correction/strain gates and line search as branch_search.
Search bounds, root clustering and candidate continuations are numerical probes;
a root reached by a finite increment is not an admissible monotone removal jump.
"""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from branch_path import Adapter, parent_binding, digest, read, write, restore, follow, ROOT
from native_arc import tangent
from continue_path import append


def run(args):
    started = time.monotonic()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    m = Adapter(args.map, args.evidence, args.binary, args.reference, args.arc, out, args.material_steps)
    expected = parent_binding(args, m)
    seed = restore(args.seed, expected)
    if seed['contacts'] != m.contacts:
        raise ValueError('seed contact history differs')
    old = m.state(seed['states'])
    u = np.array(seed['displacement_m'])
    lam = seed['progress']
    t = np.array(seed['tangent'])
    B = m.m.rb if args.map == 'native' else m.m.B
    f0, K0, s0 = m.evaluate(np.zeros(len(m.L)), old, 'root-seed-zero')
    seed_rel = float(np.linalg.norm(m.R.T@m.m.load(lam)-f0)/max(np.linalg.norm(m.m.load(lam)), 1.))
    if seed_rel >= 1e-8 or m.serial(s0) != m.serial(old):
        raise ValueError('root seed equilibrium/freeze gate')
    eig, modes = np.linalg.eigh((K0+K0.T)/2)
    directions = [('saved_tangent', t[:-1]/np.linalg.norm(t[:-1])), ('symmetric_mode_0', modes[:, 0]), ('symmetric_mode_1', modes[:, 1])]
    guesses = [-4e-4, -2e-4, 0., 1e-4, 2e-4, 3e-4, 4e-4]
    targets = [lam, 1.]
    binding = dict(parent_binding=expected, seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
                   driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write(out/'controls.json', dict(binding=binding, targets=targets, initial_amplitudes_m=guesses,
        directions={k:v.tolist() for k,v in directions}, seed_symmetric_eigenvalues_kn_m=eig.tolist(),
        force_tolerance=1e-8, true_linear_tolerance=1e-10, root_correction_check_m=1e-11,
        maximum_total_strain_l1=.001, maximum_newton_iterations=20,
        line_search=[1., .5, .25, .125, .0625, .03125, .015625],
        cluster_distance_m=1e-8, budget_seconds=args.budget_seconds, production_adopted=False,
        seed_equilibrium_relative_residual=seed_rel, seed_zero_increment_state_exact=True,
        input_only=args.map == 'independent', stability_acceptance=False))
    trials, solutions, distinct = [], [], []
    exhausted = False
    for target in targets:
        for direction_name, direction in directions:
            for guess in guesses:
                if time.monotonic()-started >= args.budget_seconds:
                    exhausted = True
                    break
                du = guess*direction
                steps = []
                converged = False
                failure = None
                try:
                    for it in range(21):
                        if time.monotonic()-started >= args.budget_seconds:
                            raise ValueError('wall time budget; inconclusive')
                        f, K, new = m.evaluate(du, old, 'local-root')
                        rhs = m.R.T@m.m.load(target)-f
                        rel = float(np.linalg.norm(rhs)/max(np.linalg.norm(m.m.load(target)), 1.))
                        dx = np.linalg.solve(K, rhs)
                        lr = float(np.linalg.norm(K@dx-rhs)/max(np.linalg.norm(rhs), 1e-30))
                        append(out/'linear-systems.jsonl', dict(trial=len(trials), iteration=it, matrix=K.tolist(), rhs=rhs.tolist(), solution=dx.tolist(), true_relative_residual=lr))
                        steps.append(dict(iteration=it, relative_force_residual=rel, correction_norm_m=float(np.linalg.norm(dx)), true_relative_linear_residual=lr, reduced_increment_m=du.tolist()))
                        if lr > 1e-10:
                            raise ValueError('fixed-history root true linear gate; inconclusive')
                        if it >= 1 and rel < 1e-8 and np.linalg.norm(dx) < 1e-11:
                            row = dict(progress=target, initial_amplitude_m=guess, direction=direction_name,
                                reduced_increment_m=du.tolist(), displacement_m=(u+(m.R@du).reshape(u.shape)).tolist(),
                                states=m.serial(new), contacts=m.contacts, relative_force_residual=rel,
                                audit=m.audit(old, new, du, lam, target, K), tangent=tangent(K, m.L, t).tolist())
                            solutions.append(row)
                            if not any(np.linalg.norm(du-np.array(x['reduced_increment_m'])) < 1e-8 and target == x['progress'] for x in distinct):
                                distinct.append(row)
                            converged = True
                            break
                        if it == 20:
                            raise ValueError('bounded root Newton limit; inconclusive')
                        accepted = False
                        for alpha in [1., .5, .25, .125, .0625, .03125, .015625]:
                            candidate = du+alpha*dx
                            if np.abs(np.einsum('qai,i->qa', B, candidate)).sum(axis=1).max() > .001:
                                continue
                            ff, _, _ = m.evaluate(candidate, old, 'root-line-search')
                            if np.linalg.norm(m.R.T@m.m.load(target)-ff) < np.linalg.norm(rhs):
                                du = candidate
                                accepted = True
                                break
                        if not accepted:
                            raise ValueError('bounded root line search no decreasing admissible trial; inconclusive')
                except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
                    failure = str(exc)
                trials.append(dict(target_progress=target, direction=direction_name, initial_amplitude_m=guess,
                    converged=converged, stop_reason=failure, iterations=steps))
                write(out/'trials.json', trials)
                write(out/'solutions.json', solutions)
                write(out/'distinct-roots.json', distinct)
            if exhausted:
                break
        if exhausted:
            break
    if args.map == 'native':
        write(out/'material-processes.json', m.m.processes)
    report = dict(passed=False, history_converged=False, static_limit_confirmed=False, production_adopted=False,
        native_committed_increments=0, trials=len(trials), converged_trials=sum(x['converged'] for x in trials),
        distinct_local_root_clusters=len(distinct), nonzero_increment_clusters=sum(int(np.linalg.norm(x['reduced_increment_m']) >= 1e-8) for x in distinct),
        full_removal_root_found=any(x['progress'] == 1. for x in distinct),
        stop_reason='wall time budget; inconclusive' if exhausted else 'bounded multistart inventory exhausted; inconclusive',
        elapsed_seconds=time.monotonic()-started,
        limitations=['fixed saved material history and bounded strain only', 'symmetric eigenvectors are search directions, not physical stability criteria',
                     'not global/disconnected branch exclusion', 'incremental roots are not accepted monotone removal jumps'])
    write(out/'report.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    for name in ['evidence', 'binary', 'reference', 'arc']:
        p.add_argument('--'+name, type=Path)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--parent-step-size', type=float, default=1e-4)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    p.add_argument('--budget-seconds', type=float, default=600.)
    a = p.parse_args()
    if not 0 < a.budget_seconds <= 1800:
        p.error('bounded budget required')
    print(json.dumps(run(a), indent=2))
