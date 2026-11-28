"""Existing hyperplane arc control on the late bounded holding-load segment.

An archived accepted state is an explicitly identified diagnostic anchor, not
an edited/resealed checkpoint. Both material maps and the 1 m arc scaling are
the existing ones. This script never grants numerical acceptance.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE), str(ROOT/'tools'), str(ROOT/'interim/2026-11-10-nc-mean-stress-release'),
                str(ROOT/'interim/2026-11-09-nc-pseudo-arclength')]
from diagnose import checksum, read, rows, spectrum, write
from native_arc import tangent
from bounded_study import Model

spec = importlib.util.spec_from_file_location('bounded_linear_refinement', ROOT/'interim/2026-11-10-nc-bounded-refinement/continuation.py')
refinement = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refinement)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def component_bounds(model, s):
    if not .5 <= s <= 1.:
        raise ValueError('late arc load coordinate outside [.5,1]; no extrapolation')
    holding = model.load(s)-model.final
    allowance = 5e-14*np.maximum(1., np.abs(model.holding))
    if np.any(np.abs(holding) > np.abs(model.holding)+allowance) or np.any(holding*model.holding < -allowance):
        raise ValueError('original signed component holding-load bound exceeded')
    return holding


def run(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    deadline = start+args.budget_seconds
    folder = args.evidence / (args.map+'-complement-refined128')
    samples = list(rows(folder/'samples.jsonl.gz'))
    anchor = max(i for i, s in enumerate(samples) if s['coordinate'] <= args.anchor)
    if anchor < 1 or samples[anchor]['coordinate'] < .5:
        raise ValueError('accepted late-segment anchor unavailable')
    sample, previous = samples[anchor], samples[anchor-1]
    controls = SimpleNamespace(map=args.map, evidence=args.evidence/'prepared/origin', binary=args.binary,
                               reference=args.evidence/'prepared/reference400', material_steps=args.material_steps,
                               order='complement-first')
    m = Model(controls, out)
    archived = read(folder/'controls.json')
    for key, value in [('B', m.B), ('transform', m.R), ('weights', m.w), ('holding_kn', m.holding)]:
        if not np.array_equal(np.asarray(archived[key]), value):
            raise ValueError('anchor model differs: '+key)
    old = m.state(sample['states'])
    u = np.asarray(sample['displacement_m'])
    s = sample['coordinate']
    n = m.R.shape[1]
    reduced = np.linalg.lstsq(m.R, (u-np.asarray(previous['displacement_m'])).ravel(), rcond=None)[0]
    t = np.r_[reduced, s-previous['coordinate']]
    t /= np.linalg.norm(t)
    L = -2*m.R.T@m.vertical
    binding = dict(starting_commit='113d83b8a37b3b68ce46eb00b6664f042334eae0',
                   map=args.map, step_size_m=args.step_size, arc_lambda_scale_m=1., material_steps=args.material_steps,
                   archive_controls_sha256=checksum(folder/'controls.json'), archived_samples_sha256=checksum(folder/'samples.jsonl.gz'),
                   anchor_sample_sha256=digest(sample), source_sha256=checksum(Path(__file__)),
                   model_source=m.source, anchor_coordinate=s, force_tolerance=1e-8, true_linear_tolerance=1e-10,
                   arc_tolerance_m=1e-12, correction_check_m=1e-11, maximum_total_increment_strain_l1=.001,
                   maximum_newton_iterations=20, coordinate_bounds=[.5, 1.], production_adopted=False,
                   independent_input_only_history=args.map == 'independent',
                   diagnostic_anchor_is_not_a_resealed_checkpoint=True)
    write(out/'controls.json', binding)
    if args.resume:
        checkpoint = read(args.resume)
        if checkpoint['binding'] != binding or digest(checkpoint['payload']) != checkpoint['payload_sha256']:
            raise ValueError('arc checkpoint binding/payload mismatch')
        p = checkpoint['payload']
        old, u, s, t, first = m.state(p['states']), np.asarray(p['displacement_m']), p['coordinate'], np.asarray(p['tangent']), p['step']+1
        if p['contacts'] != m.contacts:
            raise ValueError('arc contact checkpoint mismatch')
    else:
        first = 0
    original_solve = np.linalg.solve
    linear = []
    def solve(a, b):
        if time.monotonic() >= deadline:
            raise ValueError('wall time budget; inconclusive')
        x, history = refinement.refined_solve(original_solve, a, b)
        lr = float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b), 1e-30))
        linear.append(dict(matrix=a.tolist(), rhs=b.tolist(), solution=x.tolist(), true_relative_residual=lr,
                           original_residual=history[0], refinement_candidates=history[1:]))
        if lr > 1e-10:
            raise ValueError('existing true linear residual gate; inconclusive')
        return x
    np.linalg.solve = solve
    events = []
    saved = []
    stop = None
    seed_s = s
    try:
        if args.map == 'independent':
            original_rate = m.material.rate
            def rate(*a, **kw):
                if time.monotonic() >= deadline:
                    raise ValueError('wall time budget inside material rate; inconclusive')
                return original_rate(*a, **kw)
            m.material.rate = rate
        component_bounds(m, s)
        ff, _, zero = m.evaluate(np.zeros(n), old, 'arc-anchor-zero')
        seed = dict(coordinate=s, relative_force_residual=float(np.linalg.norm(m.R.T@m.load(s)-ff)/max(np.linalg.norm(m.load(s)), 1.)),
                    zero_increment_entire_material_exact=m.serial(zero) == m.serial(old))
        write(out/'anchor-check.json', seed)
        if seed['relative_force_residual'] >= 1e-8 or not seed['zero_increment_entire_material_exact']:
            raise ValueError('arc anchor zero-increment equilibrium/state check failed')
        for index in range(first, args.count):
            if time.monotonic() >= deadline:
                raise ValueError('wall time budget; inconclusive')
            before = digest(m.serial(old))
            du = args.step_size*t[:-1]
            target = s+args.step_size*t[-1]
            step_events = []
            try:
                for it in range(21):
                    component_bounds(m, target)
                    f, K, new = m.evaluate(du, old, f'arc-{index}-{it}')
                    rhs = m.R.T@m.load(target)-f
                    rel = float(np.linalg.norm(rhs)/max(np.linalg.norm(m.load(target)), 1.))
                    gap = args.step_size-t[:-1]@du-t[-1]*(target-s)
                    A = np.zeros((n+1, n+1))
                    A[:n, :n], A[:n, n], A[n, :n], A[n, n] = K, -L, t[:-1]*1000, t[-1]*1000
                    dx = solve(A, np.r_[rhs, gap*1000])
                    step_events.append(dict(iteration=it, coordinate=float(target), relative_force_residual=rel,
                                            arc_residual_m=float(gap), correction_norm_m=float(np.linalg.norm(dx)),
                                            reduced_increment_m=du.tolist()))
                    if it >= 1 and rel < 1e-8 and abs(gap) < 1e-12 and np.linalg.norm(dx) < 1e-11:
                        next_t = tangent(K, L, t)
                        u = u+(m.R@du).reshape(u.shape)
                        holding = component_bounds(m, target)
                        payload = dict(step=index, coordinate=float(target), previous_coordinate=s,
                                       coordinate_increment=float(target-s), tangent=next_t.tolist(), states=m.serial(new),
                                       displacement_m=u.tolist(), contacts=m.contacts, reduced_increment_m=du.tolist(),
                                       relative_force_residual=rel, arc_residual_m=float(gap),
                                       removed_material_frozen=m.frozen(new), contact_history_frozen=True,
                                       remaining_holding_norm_kn=float(np.linalg.norm(holding)),
                                       original_component_bounds_verified=True)
                        saved.append(payload)
                        cp = dict(binding=binding, payload=payload, payload_sha256=digest(payload))
                        write(out/'restart.json', cp)
                        if index == args.count//2-1:
                            write(out/'middle.json', cp)
                        old, s, t = new, float(target), next_t
                        events.append(dict(step=index, accepted=True, events=step_events))
                        break
                    if it == 20:
                        raise ValueError('existing bounded arc Newton limit; inconclusive')
                    du += dx[:-1]
                    target += dx[-1]
                if s == 1.:
                    break
            except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
                events.append(dict(step=index, accepted=False, events=step_events, stop_reason=str(exc),
                                   rejected_committed_state_unchanged=before == digest(m.serial(old))))
                raise
    except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
        stop = str(exc)
    finally:
        np.linalg.solve = original_solve
        with gzip.open(out/'samples.jsonl.gz', 'wt') as stream:
            for item in saved:
                stream.write(json.dumps(item, allow_nan=False)+'\n')
        with gzip.open(out/'linear-systems.jsonl.gz', 'wt') as stream:
            for item in linear:
                stream.write(json.dumps(item, allow_nan=False)+'\n')
        write(out/'events.json', events)
    report = dict(passed=False, full_removal_reached=s == 1., full_release_acceptance=False,
                  static_limit_confirmed=False, history_converged=False, production_adopted=False,
                  independent_full_field_to_completion=False, restart_to_completion_verified=False,
                  accepted_increments=len(saved), final_coordinate=s,
                  maximum_coordinate=max([seed_s]+[p['coordinate'] for p in saved]),
                  minimum_coordinate=min([seed_s]+[p['coordinate'] for p in saved]),
                  negative_coordinate_increments=sum(p['coordinate_increment'] < 0 for p in saved),
                  stop_reason=stop or 'bounded arc point count; inconclusive', elapsed_seconds=time.monotonic()-start)
    write(out/'report.json', report)
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--binary', type=Path)
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    p.add_argument('--step-size', type=float, choices=[1e-4, 5e-5], required=True)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    p.add_argument('--anchor', type=float, default=.9744)
    p.add_argument('--count', type=int, default=128)
    p.add_argument('--budget-seconds', type=float, default=300)
    p.add_argument('--resume', type=Path)
    a = p.parse_args()
    if not 1 <= a.count <= 256 or not 0 < a.budget_seconds <= 600 or not .5 <= a.anchor < .9745:
        p.error('bounded existing-case diagnostic controls required')
    run(a)
