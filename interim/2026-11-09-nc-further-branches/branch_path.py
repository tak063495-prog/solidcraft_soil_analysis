"""Diagnostic forward/reverse paths from verified old full-field checkpoints.

Keep the saved material history. Reversing a predictor is not reversing plastic
history. No claim of disconnected branch exclusion or physical stability.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BEYOND = HERE.parent/'2026-11-09-nc-beyond-fold'
PRIOR = HERE.parent/'2026-11-09-nc-pseudo-arclength'
sys.path.insert(0, str(BEYOND))
from continue_path import Adapter, digest, read, write, restore, follow, validate_controls


def parent_binding(args, model):
    paths = [BEYOND/'continue_path.py', PRIOR/'native_arc.py', PRIOR/'independent_map_arc.py',
             PRIOR/'independent_branches.py', ROOT/'tools/mcc_rate_fem.py', ROOT/'tools/mcc_rate_rk4.py']
    return dict(map=args.map, source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                model_sha256=digest(model.source), step_size_m=args.parent_step_size, material_steps=args.material_steps)


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
    m.t = args.orientation*np.array(seed['tangent'])
    binding = dict(map=args.map, source_sha256={**expected['source_sha256'], str(Path(__file__).relative_to(ROOT)): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                   model_sha256=expected['model_sha256'], step_size_m=args.step_size, material_steps=args.material_steps,
                   branch_seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
                   branch_seed_binding=expected, orientation=args.orientation)
    write(out/'controls.json', dict(binding=binding, count=args.count, budget_seconds=args.budget_seconds,
        seed_progress=m.lam, seed_parent_arc_length_m=seed['arc_length_m'], force_tolerance=1e-8,
        true_linear_tolerance=1e-10, arc_constraint_tolerance_m=1e-12, maximum_newton_iterations=20,
        minimum_correctors=1, arc_lambda_scale_m=1., diagnostic_only=True, production_adopted=False,
        input_only=args.map == 'independent', note='Saved material histories retained; predictor orientation is a numerical branch probe, not a stability certificate'))
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
    p.add_argument('--budget-seconds', type=float, default=600.)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    print(json.dumps(run(p.parse_args()), indent=2))
