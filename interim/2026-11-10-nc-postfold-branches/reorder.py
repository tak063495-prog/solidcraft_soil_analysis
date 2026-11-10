"""Existing Newton/load control with a different bounded Cartesian component order."""
import argparse
import gzip
import importlib
import json
from pathlib import Path
import time
import numpy as np
from paths import HERE, ROOT, model
from branch_math import component_bounds
from trace import refinement
from diagnose import checksum


def run(args):
    module = importlib.import_module('bounded_study')
    streams = refinement.Streams()
    module.append = streams.append
    def build(controls, out):
        m = model(controls, out)
        m.source['reorder_driver_sha256'] = checksum(Path(__file__))
        m.source['path_driver_sha256'] = checksum(HERE/'paths.py')
        m.source['numpy_version'] = np.__version__
        m.source['longdouble_mantissa_bits'] = np.finfo(np.longdouble).nmant
        original_load = m.load
        def load(s):
            if not 0. <= s <= 1.:
                raise ValueError('release coordinate outside [0,1]')
            result = original_load(s)
            component_bounds(m.holding, result-m.final)
            return result
        m.load = load
        (out/'model.json').write_text(json.dumps(dict(B=m.B.tolist(), transform=m.R.tolist(), weights=m.w.tolist(),
            final_load_kn=m.final.tolist(), holding_kn=m.holding.tolist(), first_group_kn=m.horizontal.tolist(),
            second_group_kn=m.vertical.tolist(), active_gp=m.pi.tolist(), initial_states=m.initial_old, contacts=m.contacts),
            indent=2,allow_nan=False)+'\n',encoding='utf-8')
        return m
    module.Model = build
    original = np.linalg.solve
    def solve(a, b):
        x, history = refinement.refined_solve(original, a, b)
        streams.append(args.output/'linear-refinement.jsonl', dict(original_relative_residual=history[0],
            refinement_candidates=history[1:], returned_relative_residual=float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b), 1e-30))))
        return x
    np.linalg.solve = solve
    try:
        module.run(args)
    finally:
        np.linalg.solve = original
        streams.close(args.output)
    # The original run saves every accepted field and last 128 actual matrices.
    report = json.loads((args.output/'report.json').read_text())
    coordinates = [row['coordinate'] for row in (json.loads(line) for line in gzip.open(args.output/'samples.jsonl.gz', 'rt'))]
    report['maximum_coordinate'] = max(coordinates)
    report['component_order'] = args.order
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['evidence', 'output', 'binary', 'reference', 'resume']:
        p.add_argument('--'+name, type=Path, required=name in ['evidence', 'output'])
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    p.add_argument('--order', choices=['x-first', 'x-last', 'y-first', 'y-last', 'z-first', 'z-last'], required=True)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    p.add_argument('--step-size', type=float, choices=[1/128, 1/256], default=1/128)
    p.add_argument('--minimum-step', type=float, default=2**-22)
    p.add_argument('--count', type=int, default=2048)
    p.add_argument('--budget-seconds', type=float, default=300.)
    a = p.parse_args()
    if not 0 < a.budget_seconds <= 1800 or not 0 < a.count <= 4096 or not 0 < a.minimum_step <= a.step_size:
        p.error('bounded controls required')
    run(a)
