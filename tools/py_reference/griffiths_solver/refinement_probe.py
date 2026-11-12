"""Small elastic h sensitivity on a unit block; NOT a Griffiths h series.

Fixed depth ny=1, x/z subdivisions 2/4/8, only 24/96/384 Tet10 elements.
This isolates geometric conditioning without a nonlinear constitutive path.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from .systems import assemble
from .solvers import compare, estimates
from plastic_spatial_convergence import inputs, lift


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists(): p.error('output already exists')
    rows = []
    for n in (2, 4, 8):
        mesh, _, case = inputs(n, depth_divisions=1)
        mesh, _ = lift(mesh, [[0., 0., 0.]]*len(mesh['nodes']))
        xyz = np.asarray(mesh['nodes'])
        mesh['node_sets'] = dict(all=list(range(len(xyz))),
                                sides=np.flatnonzero((xyz[:, 0] == 0) | (xyz[:, 0] == 1)).tolist(),
                                bottom=np.flatnonzero(xyz[:, 2] == 0).tolist())
        a, b, free, meta = assemble(mesh, case=case)
        rows.append(dict(xz_subdivisions=n, depth_divisions=1, relative_h=2/n,
                         metadata=meta, estimates=estimates(a), solvers=compare(a, b, free)))
    report = dict(schema='soil-fem-small-elastic-conditioning-probe/1',
                  domain='unit block; not Griffiths', griffiths_h_series=False,
                  p6_complete=False, nonlinear_solve_performed=False, levels=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    raise SystemExit(0 if all(r['accepted'] for level in rows for r in level['solvers']) else 1)


if __name__ == '__main__':
    main()
