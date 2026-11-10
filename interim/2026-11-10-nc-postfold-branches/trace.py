"""Longer existing hyperplane arc traces and explicit transverse weak-mode probes.

All paths are diagnostic. No material, contact, geometry or acceptance gate changes.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from paths import HERE, ROOT, model
from branch_math import component_bounds, diagnose, probe_tangent
from diagnose import checksum, read, rows, write
from native_arc import tangent
import importlib.util

spec = importlib.util.spec_from_file_location('branch_refinement', ROOT/'interim/2026-11-10-nc-bounded-refinement/continuation.py')
refinement = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refinement)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def run(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    end = started+args.budget_seconds
    m = model(args, out)
    source = args.anchor_samples or args.evidence/(args.map+'-complement-refined128')/'samples.jsonl.gz'
    samples = list(rows(source))
    model_file = source.parent/'model.json'
    archived = read(model_file if model_file.exists() else source.parent/'controls.json')
    for key, value in [('B',m.B),('transform',m.R),('weights',m.w),('holding_kn',m.holding)]:
        if key not in archived or not np.array_equal(np.asarray(archived[key]),value):
            raise ValueError('anchor model differs: '+key)
    for names, value in [(('first_group_kn','horizontal_holding_kn'),m.horizontal),(('second_group_kn','vertical_holding_kn'),m.vertical)]:
        key=next((key for key in names if key in archived),None)
        if key is None or not np.array_equal(np.asarray(archived[key]),value):
            raise ValueError('anchor component order differs')
    if args.anchor_step is not None:
        anchor = next(i for i, p in enumerate(samples) if p['step'] == args.anchor_step)
    else:
        anchor = max(i for i, p in enumerate(samples) if p['coordinate'] <= args.anchor)
    if anchor < 1:
        raise ValueError('two accepted samples required to define anchor direction')
    p, previous = samples[anchor], samples[anchor-1]
    old, u, s = m.state(p['states']), np.asarray(p['displacement_m']), p['coordinate']
    n = m.R.shape[1]
    secant = np.linalg.lstsq(m.R, (u-np.asarray(previous['displacement_m'])).ravel(), rcond=None)[0]
    t = np.r_[secant, s-previous['coordinate']]
    t /= np.linalg.norm(t)
    if 'tangent' in p:
        t = np.asarray(p['tangent'])
    if args.probe:
        if 'branch_diagnostic' not in p:
            raise ValueError('weak mode must come from the actual saved accepted increment')
        t = probe_tangent(t, p['branch_diagnostic']['weak_right_mode'], args.probe)
    lower = .5 if s >= .5 else 0.
    upper = 1. if s >= .5 else .5
    L = -2*m.R.T@(m.vertical if s >= .5 else m.horizontal)
    if args.order != 'complement-first' and args.anchor_samples is None:
        raise ValueError('different order requires its own computed history')
    binding = dict(starting_commit='b235824ea55a0e70b8788e586ee8562549f8baf5',
                   source_sha256={f.name: checksum(f) for f in [Path(__file__), HERE/'paths.py', HERE/'branch_math.py']},
                   map=args.map, order=args.order, material_steps=args.material_steps, model_source=m.source,
                   anchor_samples_sha256=checksum(source), anchor_sample_sha256=digest(p), anchor_step=p['step'],
                   step_size_m=args.step_size, count=args.count, probe=args.probe,
                   force_tolerance=1e-8, true_linear_tolerance=1e-10, correction_check_m=1e-11,
                   arc_tolerance_m=1e-12, maximum_total_increment_strain_l1=.001,
                   maximum_newton_iterations=20, coordinate_bounds=[lower, upper], arc_coordinate_scale_m=1.,
                   anchor_is_not_edited_checkpoint=True, production_adopted=False,
                   numpy_version=np.__version__, longdouble_mantissa_bits=np.finfo(np.longdouble).nmant)
    write(out/'controls.json', binding)
    write(out/'model.json', dict(B=m.B.tolist(), transform=m.R.tolist(), weights=m.w.tolist(),
          final_load_kn=m.final.tolist(), holding_kn=m.holding.tolist(), first_group_kn=m.horizontal.tolist(),
          second_group_kn=m.vertical.tolist(), active_gp=m.pi.tolist(), initial_states=m.initial_old, contacts=m.contacts))
    def bounds(coordinate):
        if not lower <= coordinate <= upper:
            raise ValueError('current piecewise segment boundary; inconclusive')
        remainder = m.load(coordinate)-m.final
        component_bounds(m.holding, remainder)
        return remainder
    original = np.linalg.solve
    streams = {name: gzip.open(out/(name+'.jsonl.gz'), 'wt', compresslevel=3) for name in ['samples', 'linear-systems', 'events', 'tangents']}
    context = {}
    def append(name, row):
        streams[name].write(json.dumps(row, allow_nan=False)+'\n')
    def solve(a, b):
        if time.monotonic() >= end:
            raise ValueError('wall time budget; inconclusive')
        x, history = refinement.refined_solve(original, a, b)
        lr = float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b), 1e-30))
        append('linear-systems', dict(**context, matrix=a.tolist(), rhs=b.tolist(), solution=x.tolist(),
                                     true_relative_residual=lr, refinement_history=history))
        if lr > 1e-10:
            raise ValueError('existing true linear residual gate; inconclusive')
        return x
    count, negatives, turns = 0, 0, []
    maximum, minimum, cumulative, stop = s, s, 0., None
    ds = args.step_size
    if args.resume:
        cp = read(args.resume)
        if cp['binding'] != binding or digest(cp['payload']) != cp['payload_sha256']:
            raise ValueError('actual checkpoint binding/payload mismatch')
        p = cp['payload']
        old, u, s, t = m.state(p['states']), np.asarray(p['displacement_m']), p['coordinate'], np.asarray(p['tangent'])
        first, ds, cumulative = p['step']+1, p['next_step_size_m'], p['cumulative_arc_m']
        if p['contacts'] != m.contacts:
            raise ValueError('contact history checkpoint mismatch')
    else:
        first = 0
    maximum, minimum = s, s
    np.linalg.solve = solve
    try:
        bounds(s)
        f, _, zero = m.evaluate(np.zeros(n), old, 'anchor-zero')
        seed = dict(relative_force_residual=float(np.linalg.norm(m.R.T@m.load(s)-f)/max(np.linalg.norm(m.load(s)), 1.)),
                    zero_increment_entire_material_exact=m.serial(zero) == m.serial(old), coordinate=s)
        write(out/'anchor-check.json', seed)
        if seed['relative_force_residual'] >= 1e-8 or not seed['zero_increment_entire_material_exact']:
            raise ValueError('anchor equilibrium/state check failed')
        index = first
        attempt = 0
        while index < args.count:
            if time.monotonic() >= end:
                raise ValueError('wall time budget; inconclusive')
            attempt += 1
            before = digest(m.serial(old))
            du, target = ds*t[:-1], s+ds*t[-1]
            endpoint = bool(upper == 1. and target >= 1. and t[-1] > 0.)
            if endpoint:
                du = ((1.-s)/t[-1])*t[:-1]
                target = 1.
            events = []
            try:
                for it in range(21):
                    context.update(step=index, iteration=it, attempt=attempt, purpose='corrector')
                    bounds(target)
                    f, K, new = m.evaluate(du, old, f'arc-{index}-{it}')
                    rhs = m.R.T@m.load(target)-f
                    rel = float(np.linalg.norm(rhs)/max(np.linalg.norm(m.load(target)), 1.))
                    gap = 1.-target if endpoint else ds-t[:-1]@du-t[-1]*(target-s)
                    A = np.zeros((n+1, n+1))
                    A[:n, :n], A[:n, n] = K, -L
                    A[n, :n] = 0. if endpoint else t[:-1]*1000
                    A[n, n] = 1000. if endpoint else t[-1]*1000
                    dx = solve(A, np.r_[rhs, gap*1000])
                    events.append(dict(iteration=it, coordinate=float(target), relative_force_residual=rel,
                                       arc_residual_m=float(gap), correction_norm_m=float(np.linalg.norm(dx)), reduced_increment_m=du.tolist()))
                    if it >= 1 and rel < 1e-8 and abs(gap) < 1e-12 and np.linalg.norm(dx) < 1e-11:
                        context['purpose'] = 'tangent'
                        next_t = tangent(K, L, t)
                        diagnostic = diagnose(K, L, next_t, t)
                        new_u = u+(m.R@du).reshape(u.shape)
                        cumulative += float(np.linalg.norm(np.r_[du, target-s]))
                        payload = dict(step=index, coordinate=float(target), previous_coordinate=s, coordinate_increment=float(target-s),
                            tangent=next_t.tolist(), displacement_m=new_u.tolist(), states=m.serial(new), contacts=m.contacts,
                            reduced_increment_m=du.tolist(), relative_force_residual=rel, arc_residual_m=float(gap),
                            removed_material_frozen=m.frozen(new), contact_history_frozen=True,
                            remaining_holding_norm_kn=float(np.linalg.norm(bounds(target))), original_component_bounds_verified=True,
                            endpoint_load_control=endpoint, branch_diagnostic=diagnostic, cumulative_arc_m=cumulative, next_step_size_m=min(args.step_size, ds*2))
                        append('samples', payload)
                        append('tangents', dict(step=index, matrix=K.tolist(), load_direction=L.tolist(), diagnostic=diagnostic))
                        cp = dict(binding=binding, payload=payload, payload_sha256=digest(payload))
                        write(out/'restart.json', cp)
                        if diagnostic['turning_candidate']:
                            turns.append(dict(step=index, coordinate=float(target), **diagnostic))
                            write(out/f'turning-{index}.json', cp)
                        if index == args.count//2-1:
                            write(out/'middle.json', cp)
                        count += 1
                        negatives += target < s
                        maximum, minimum = max(maximum, target), min(minimum, target)
                        old, u, s, t, ds = new, new_u, float(target), next_t, payload['next_step_size_m']
                        append('events', dict(step=index, attempt=attempt, accepted=True, events=events))
                        index += 1
                        if index % 128 == 0:
                            print(json.dumps(dict(step=index, coordinate=s, max_coordinate=maximum, arc_m=cumulative)), flush=True)
                        break
                    if it == 20:
                        raise ValueError('existing arc Newton limit; inconclusive')
                    du += dx[:-1]
                    target += dx[-1]
            except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
                unchanged = before == digest(m.serial(old))
                append('events', dict(step=index, attempt=attempt, accepted=False, events=events, stop_reason=str(exc),
                                     rejected_committed_state_unchanged=unchanged, attempted_step_size_m=ds))
                if not unchanged:
                    raise ValueError('rejected committed state changed') from exc
                if 'wall time' in str(exc) or ds/2 < args.minimum_step:
                    raise
                ds /= 2
                continue
            if s == 1.:
                f, _, frozen = m.evaluate(np.zeros(n), old, 'full-release-actual-zero')
                write(out/'full-release-zero.json', dict(material_zero_increment_exact=m.serial(frozen) == m.serial(old),
                    contacts=m.contacts, remaining_holding_kn=bounds(s).tolist(),
                    relative_force_residual=float(np.linalg.norm(m.R.T@m.final-f)/max(np.linalg.norm(m.final), 1.))))
                stop = 'full release diagnostic reached; acceptance evidence still required'
                break
    except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
        stop = str(exc)
    finally:
        np.linalg.solve = original
        for stream in streams.values():
            stream.close()
    report = dict(passed=False, full_removal_reached=s == 1., full_release_acceptance=False, static_limit_confirmed=False,
        history_converged=False, independent_full_field_to_completion=False, restart_to_completion_verified=False,
        accepted_increments=count, maximum_coordinate=float(maximum), minimum_coordinate=float(minimum), final_coordinate=s,
        negative_coordinate_increments=int(negatives), turning_candidates=turns, bifurcation_confirmed=False,
        cumulative_arc_m=cumulative, remaining_holding_norm_kn=float(np.linalg.norm(m.load(s)-m.final)),
        stop_reason=stop or 'bounded arc point count; inconclusive', elapsed_seconds=time.monotonic()-started)
    write(out/'report.json', report)
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['evidence', 'output', 'binary', 'anchor-samples', 'resume']:
        p.add_argument('--'+name, type=Path, required=name in ['evidence', 'output'])
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    p.add_argument('--order', choices=['mean-first', 'complement-first', 'x-first', 'x-last', 'y-first', 'y-last', 'z-first', 'z-last'], default='complement-first')
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    p.add_argument('--anchor', type=float, default=.9744)
    p.add_argument('--anchor-step', type=int)
    p.add_argument('--probe', type=float, choices=[-1., -.25, 0., .25, 1.], default=0.)
    p.add_argument('--step-size', type=float, default=.0002)
    p.add_argument('--minimum-step', type=float, default=1e-9)
    p.add_argument('--count', type=int, default=2048)
    p.add_argument('--budget-seconds', type=float, default=600.)
    args = p.parse_args()
    if not 0 < args.minimum_step <= args.step_size <= .001 or not 1 <= args.count <= 8192 or not 0 < args.budget_seconds <= 1800:
        p.error('bounded numerical controls required')
    run(args)
