"""Numerical release ordering of the original holding vector, without extra loads."""
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'interim/2026-11-10-nc-tangent-branch-diagnosis'),
                str(ROOT/'interim/2026-11-10-nc-mean-stress-release'),
                str(ROOT/'interim/2026-11-09-nc-pseudo-arclength')]
from bounded_study import Model


def model(args, out):
    controls = SimpleNamespace(map=args.map, evidence=args.evidence/'prepared/origin', binary=args.binary,
                               reference=args.evidence/'prepared/reference400', material_steps=args.material_steps,
                               order=args.order if args.order in ('mean-first', 'complement-first') else 'complement-first')
    m = Model(controls, out)
    if args.order not in ('mean-first', 'complement-first'):
        axis, sequence = args.order.split('-')
        selected = np.arange(len(m.holding)) % 3 == 'xyz'.index(axis)
        if sequence == 'last':
            selected = ~selected
        m.horizontal = np.where(selected, m.holding, 0.)
        m.vertical = m.holding-m.horizontal
        m.order = args.order
    return m
