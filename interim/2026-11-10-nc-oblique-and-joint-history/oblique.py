"""Diagnostic two-segment release using existing bounded Newton and material maps.

The waypoint changes only the numerical withdrawal schedule of the original H.
It is neither an added load nor a production construction-stage specification.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'interim/2026-11-10-nc-postfold-branches')]
from paths import model as previous_model
from branch_math import component_bounds
from trace import refinement


def fractions(s, waypoint):
    if not 0. <= s <= 1. or not all(0. < x < 1. for x in waypoint):
        raise ValueError('bounded withdrawal coordinates required')
    w = np.asarray(waypoint)
    result = 2*s*w if s <= .5 else w+(2*s-1)*(1-w)
    return tuple(float(x) for x in result)


def run(args):
    args.output = args.output.resolve()
    module = importlib.import_module('bounded_study')
    streams = refinement.Streams()
    module.append = streams.append
    module.fractions = lambda s, order: fractions(s, args.waypoint)
    original_write = module.write
    def write(path, value):
        if path.name == 'controls.json':
            value['grouping'] = 'existing clipped mean/complement groups, complement first group; oblique simultaneous withdrawal'
            value['withdrawal_waypoint_at_half'] = args.waypoint
            value['coordinate_is_not_original_proportional_lambda'] = True
            value['invocation_commit'] = args.invocation_commit
        original_write(path, value)
    module.write = write
    def build(controls, out):
        m = previous_model(controls, out)
        m.source.update(oblique_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        withdrawal_waypoint_at_half=args.waypoint,
                        invocation_commit=args.invocation_commit)
        def load(s):
            a, b = fractions(s, args.waypoint)
            result = m.final+(1-a)*m.horizontal+(1-b)*m.vertical
            component_bounds(m.holding, result-m.final)
            return result
        m.load = load
        write(out/'model.json', dict(B=m.B.tolist(), transform=m.R.tolist(), weights=m.w.tolist(),
              final_load_kn=m.final.tolist(), holding_kn=m.holding.tolist(),
              first_group_kn=m.horizontal.tolist(), second_group_kn=m.vertical.tolist(),
              active_gp=m.pi.tolist(), initial_states=m.initial_old, contacts=m.contacts,
              waypoint=args.waypoint, map=args.map))
        return m
    module.Model = build
    original = np.linalg.solve
    def solve(a, b):
        x, history = refinement.refined_solve(original, a, b)
        streams.append(args.output/'linear-refinement.jsonl', dict(
            original_relative_residual=history[0], refinement_candidates=history[1:],
            returned_relative_residual=float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b), 1e-30))))
        return x
    np.linalg.solve = solve
    try:
        module.run(args)
    finally:
        np.linalg.solve = original
        streams.close(args.output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('evidence', 'output', 'binary', 'reference', 'resume'):
        p.add_argument('--'+name, type=Path, required=name in ('evidence', 'output'))
    p.add_argument('--map', choices=('native', 'independent'), required=True)
    p.add_argument('--waypoint', type=float, nargs=2, required=True)
    p.add_argument('--invocation-commit', required=True)
    p.add_argument('--order', choices=('complement-first',), default='complement-first')
    p.add_argument('--material-steps', type=int, choices=(4, 8), default=4)
    p.add_argument('--step-size', type=float, choices=(1/128, 1/256), default=1/128)
    p.add_argument('--minimum-step', type=float, default=2**-22)
    p.add_argument('--count', type=int, default=2048)
    p.add_argument('--budget-seconds', type=float, default=180.)
    a = p.parse_args()
    if not 0 < a.budget_seconds <= 1800 or not 0 < a.count <= 4096 or not 0 < a.minimum_step <= a.step_size:
        p.error('bounded controls required')
    fractions(.5, a.waypoint)
    run(a)
