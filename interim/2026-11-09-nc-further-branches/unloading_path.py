"""Bounded loading/unloading rays using the actual zero-increment material K.

Keep full saved histories; do not reverse plastic variables or symmetrize K for
solving. This is an additional branch predictor, not a production solver fix.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from branch_path import Adapter, parent_binding, restore, follow, write, digest, validate_controls, ROOT


def run(args):
    validate_controls(args.step_size, args.count, args.budget_seconds)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    m = Adapter(args.map, args.evidence, args.binary, args.reference, args.arc, out, args.material_steps)
    expected = parent_binding(args, m)
    seed = restore(args.seed, expected)
    if seed['contacts'] != m.contacts:
        raise ValueError('seed contact history differs')
    m.old = m.state(seed['states'])
    m.u = np.array(seed['displacement_m'])
    m.lam = seed['progress']
    f, K, frozen = m.evaluate(np.zeros(len(m.L)), m.old, 'alternative-seed-zero')
    rel = float(np.linalg.norm(m.R.T@m.m.load(m.lam)-f)/max(np.linalg.norm(m.m.load(m.lam)), 1.))
    if rel >= 1e-8 or m.serial(frozen) != m.serial(m.old):
        raise ValueError('alternative seed equilibrium/freeze gate')
    rate = np.linalg.solve(K, m.L)
    lr = float(np.linalg.norm(K@rate-m.L)/max(np.linalg.norm(m.L), 1e-30))
    write(out/'seed-tangent-system.json', dict(matrix=K.tolist(), rhs=m.L.tolist(), solution=rate.tolist(), true_relative_residual=lr))
    if lr > 1e-10:
        raise ValueError('alternative predictor linear gate')
    m.t = args.orientation*np.r_[rate, 1.]/np.linalg.norm(np.r_[rate, 1.])
    sources = {**expected['source_sha256']}
    for p in [Path(__file__), Path(__file__).with_name('branch_path.py')]:
        sources[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    binding = dict(map=args.map, source_sha256=sources, model_sha256=expected['model_sha256'],
                   step_size_m=args.step_size, material_steps=args.material_steps,
                   parent_binding=expected, branch_seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
                   predictor='actual zero-increment material tangent', orientation=args.orientation)
    write(out/'controls.json', dict(binding=binding, count=args.count, budget_seconds=args.budget_seconds,
        seed_progress=m.lam, seed_parent_arc_length_m=seed['arc_length_m'], force_tolerance=1e-8,
        true_linear_tolerance=1e-10, arc_constraint_tolerance_m=1e-12, maximum_newton_iterations=20,
        minimum_correctors=1, arc_lambda_scale_m=1., diagnostic_only=True, production_adopted=False,
        input_only=args.map == 'independent', predictor_is_stability_evidence=False,
        seed_equilibrium_relative_residual=rel, seed_zero_increment_state_exact=True))
    return follow(m, out, args.step_size, args.count, args.budget_seconds, binding, args.resume)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    for name in ['evidence', 'binary', 'reference', 'arc', 'resume']:
        p.add_argument('--'+name, type=Path)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--orientation', type=int, choices=[-1, 1], required=True)
    p.add_argument('--parent-step-size', type=float, default=1e-4)
    p.add_argument('--step-size', type=float, required=True)
    p.add_argument('--count', type=int, required=True)
    p.add_argument('--budget-seconds', type=float, default=900.)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    print(json.dumps(run(p.parse_args()), indent=2))
