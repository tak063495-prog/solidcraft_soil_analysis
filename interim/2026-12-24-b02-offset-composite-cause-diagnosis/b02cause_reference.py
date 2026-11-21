"""Input-only whole-system reference quadrature/order sensitivity, in a child PID."""
import argparse
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from b02cause_binding import legacy


@contextmanager
def volume_control(order=5, reverse=False):
    module = legacy['reference']
    original = module.ref

    def volume(mesh, ei, independent=False, independent_order=5):
        q = original.volume(mesh, ei, independent, order if independent else independent_order)
        if independent and reverse:
            q = {k: (v[::-1].copy() if k != 'nodes' else v) for k, v in q.items()}
        return q

    module.ref = SimpleNamespace(**{**vars(original), 'volume': volume})
    try:
        yield
    finally:
        module.ref = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--order', type=int, choices=(5, 7), default=5)
    p.add_argument('--reverse', action='store_true')
    a = p.parse_args()
    with volume_control(a.order, a.reverse):
        legacy['reference'].solve(a.folder, a.output)
