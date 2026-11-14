"""Extend the existing two 3D material-map diagnostics; never certify a job.

Same moving hyperplane, material maps, units and residual gates as the prior
diagnostic. Bounds stop the trace, not a physical failure criterion. JSONL keeps
every accepted full field and every actual solved bordered system auditable.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = HERE.parent / '2026-11-09-nc-pseudo-arclength'
sys.path.insert(0, str(PRIOR))
from native_arc import Native
from independent_map_arc import Model, read_arc_seed, seed_failure


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def append(path, value):
    with path.open('a') as stream:
        stream.write(json.dumps(value, allow_nan=False) + '\n')


def validate_controls(ds, count, budget):
    if not (np.isfinite(ds) and 0 < ds <= 1e-4 and type(count) is int
            and 0 < count <= 4096 and np.isfinite(budget) and 0 < budget <= 1800):
        raise ValueError('bounded diagnostic controls required')


def restore(path, binding):
    data = read(path)
    if data['binding'] != binding or data['payload_sha256'] != digest(data['payload']):
        raise ValueError('checkpoint binding or payload differs')
    return data['payload']


def checkpoint(path, binding, payload):
    write(path, dict(binding=binding, payload=payload, payload_sha256=digest(payload)))


class Adapter:
    def __init__(self, kind, evidence, binary, reference, arc, out, material_steps):
        self.kind = kind
        if kind == 'native':
            self.m = m = Native(evidence, binary, out / 'material')
            self.old = m.cp['material_states']
            self.u = np.array(m.cp['displacement_m'])
            self.lam = m.lam0
            t = np.r_[np.linalg.solve(m.k0, m.L), 1.]
            self.t = t / np.linalg.norm(t)
            self.contacts = m.cp['joints']
            self.source = dict(seed=m.cp, holding=m.holding.tolist(), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
        else:
            self.m = m = Model(reference, arc, material_steps)
            if m.system.live.any() or np.any(m.system.kjoint):
                raise ValueError('independent contact branch is not inactive')
            self.old = m.oldfull.copy()
            self.u = np.array(m.seed['displacement_m'])
            self.lam = m.seed['progress']
            self.t = np.array(m.seed['tangent'])
            self.contacts = m.seed['contacts']
            self.source = dict(seed=m.seed, unload=read(reference / 'unload.json'), inputs=[m.mesh, m.ground, m.case])
        self.R = m.R
        self.L = m.L

    def evaluate(self, du, old, label):
        if self.kind == 'native':
            return self.m.evaluate(du, old, label)
        return self.m.evaluate(du, np.array(old))

    def audit(self, old, new, du, lam0, lam, K):
        return self.m.audit(old, new, du, lam0, lam, K)

    def serial(self, old):
        return old if self.kind == 'native' else np.asarray(old).tolist()

    def state(self, old):
        return old if self.kind == 'native' else np.array(old)


def follow(model, out, ds, count, budget, binding, resume=None):
    started = time.monotonic()
    old, u, lam, t = model.old, model.u.copy(), model.lam, model.t.copy()
    index = 0
    arc_length = 0.
    if resume:
        saved = restore(resume, binding)
        old = model.state(saved['states'])
        u = np.array(saved['displacement_m'])
        lam = saved['progress']
        t = np.array(saved['tangent'])
        index = saved['step'] + 1
        arc_length = saved['arc_length_m']
        if saved['contacts'] != model.contacts:
            raise ValueError('checkpoint contact history differs')
    initial = dict(step=index-1, progress=float(lam), tangent=t.tolist(), states=model.serial(old),
                   displacement_m=u.tolist(), contacts=model.contacts, arc_length_m=arc_length)
    checkpoint(out / 'initial.json', binding, initial)
    samples = 0
    stop = None
    n = len(model.L)
    try:
        f, _, seed_state = model.evaluate(np.zeros(n), old, 'seed-zero')
        seed_residual = float(np.linalg.norm(model.R.T@model.m.load(lam)-f)/max(np.linalg.norm(model.m.load(lam)), 1.))
        seed_exact = model.serial(seed_state) == model.serial(old)
        write(out / 'seed-check.json', dict(relative_force_residual=seed_residual, zero_increment_state_exact=seed_exact))
        if seed_residual >= 1e-8 or not seed_exact:
            raise ValueError('seed equilibrium or zero increment freeze gate; inconclusive')
        for step in range(index, index+count):
            if time.monotonic() - started >= budget:
                raise ValueError('wall time budget; inconclusive')
            du = ds * t[:-1]
            next_lam = lam + ds*t[-1]
            endpoint = next_lam >= 1 and t[-1] > 0
            if endpoint:
                next_lam = 1.
            for it in range(21):
                if not (0 <= next_lam <= 1):
                    raise ValueError('trial outside diagnostic lambda [0,1]; inconclusive')
                f, K, new = model.evaluate(du, old, f'step-{step}:newton-{it}')
                target = model.m.load(next_lam)
                rhs = model.R.T @ target - f
                rel = float(np.linalg.norm(rhs) / max(np.linalg.norm(target), 1.))
                g = float(1-next_lam if endpoint else ds - t[:-1] @ du - t[-1]*(next_lam-lam))
                append(out / 'events.jsonl', dict(step=step, iteration=it, progress=float(next_lam),
                       relative_force_residual=rel, arc_residual_m=g, reduced_increment_m=du.tolist()))
                A = np.zeros((n+1, n+1))
                A[:n, :n], A[:n, n] = K, -model.L
                A[n, :] = t*1000
                if endpoint:
                    A[n, :] = 0.
                    A[n, n] = 1000.
                accepted = it >= 1 and rel < 1e-8 and abs(g) < 1e-12
                if it == 20 and not accepted:
                    raise ValueError('Newton limit; inconclusive')
                b = np.r_[np.zeros(n), 1000.] if accepted else np.r_[rhs, g*1000]
                x = np.linalg.solve(A, b)
                lr = float(np.linalg.norm(A@x-b)/max(np.linalg.norm(b), 1e-30))
                append(out / 'linear-systems.jsonl', dict(step=step, iteration=it,
                       kind='tangent' if accepted else 'corrector', matrix=A.tolist(), rhs=b.tolist(),
                       solution=x.tolist(), true_relative_residual=lr))
                if lr > 1e-10:
                    raise ValueError('true linear residual gate; inconclusive')
                if accepted:
                    next_t = x/np.linalg.norm(x)
                    if t@next_t < 0:
                        next_t = -next_t
                    u += (model.R@du).reshape(u.shape)
                    arc_length += ds
                    row = dict(step=step, arc_length_m=arc_length, progress=float(next_lam),
                               previous_progress=float(lam), progress_increment=float(next_lam-lam),
                               displacement_m=u.tolist(), states=model.serial(new), contacts=model.contacts,
                               tangent=next_t.tolist(), relative_force_residual=rel, arc_residual_m=g,
                               audit=model.audit(old, new, du, lam, next_lam, K))
                    append(out / 'samples.jsonl', row)
                    checkpoint(out / 'restart.json', binding, row)
                    old, lam, t = new, next_lam, next_t
                    samples += 1
                    break
                du += x[:-1]
                next_lam += x[-1]
            if endpoint:
                f, _, frozen = model.evaluate(np.zeros(n), old, 'full-release-zero')
                write(out / 'full-release.json', dict(progress=float(lam),
                      relative_force_residual=float(np.linalg.norm(model.R.T@model.m.load(lam)-f)/max(np.linalg.norm(model.m.load(lam)), 1.)),
                      zero_increment_material_exact=model.serial(frozen) == model.serial(old), contacts=model.contacts))
                stop = 'lambda=1 diagnostic endpoint; acceptance evidence still required'
                break
        else:
            stop = 'bounded step count; inconclusive'
    except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
        stop = str(exc)
    report = dict(passed=False, history_converged=False, static_limit_confirmed=False,
                  production_adopted=False, native_committed_increments=0, equilibrium_samples=samples,
                  stop_reason=stop, last_progress=float(lam), arc_length_m=arc_length,
                  full_removal_reached=bool(lam == 1.), full_release_freeze_verified=False,
                  restart_to_completion_verified=False, elapsed_seconds=time.monotonic()-started)
    write(out / 'report.json', report)
    if model.kind == 'native':
        write(out / 'material-processes.json', model.m.processes)
    return report


def run(args):
    validate_controls(args.step_size, args.count, args.budget_seconds)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    if args.map == 'independent':
        try:
            read_arc_seed(args.arc)
        except ValueError as exc:
            return seed_failure(out, exc)
    model = Adapter(args.map, args.evidence, args.binary, args.reference, args.arc, out, args.material_steps)
    sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in [Path(__file__), PRIOR/'native_arc.py', PRIOR/'independent_map_arc.py',
                         PRIOR/'independent_branches.py', ROOT/'tools/mcc_rate_fem.py', ROOT/'tools/mcc_rate_rk4.py']}
    binding = dict(map=args.map, source_sha256=sources, model_sha256=digest(model.source),
                   step_size_m=args.step_size, material_steps=args.material_steps)
    write(out / 'controls.json', dict(binding=binding, count=args.count, budget_seconds=args.budget_seconds,
          force_tolerance=1e-8, true_linear_tolerance=1e-10, arc_constraint_tolerance_m=1e-12,
          maximum_newton_iterations=20, minimum_correctors=1, arc_lambda_scale_m=1.,
          diagnostic_only=True, input_only=args.map == 'independent', production_adopted=False))
    return follow(model, out, args.step_size, args.count, args.budget_seconds, binding, args.resume)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    for name in ['evidence', 'binary', 'reference', 'arc', 'resume']:
        p.add_argument('--'+name, type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--step-size', type=float, required=True)
    p.add_argument('--count', type=int, required=True)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    p.add_argument('--budget-seconds', type=float, default=300.)
    result = run(p.parse_args())
    print(json.dumps(result, indent=2))
    sys.exit(1 if result.get('diagnostic_status') == 'failed' or result['equilibrium_samples'] == 0 else 0)
