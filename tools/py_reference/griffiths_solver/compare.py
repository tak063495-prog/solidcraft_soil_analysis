"""Compare frozen or reconstructed Griffiths tangent systems with SciPy only.

Run from the repository root: python -m tools.py_reference.griffiths_solver.compare
"""
import argparse
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from .systems import ROOT, read, sha, assemble, griffiths, save_system, load_system
from .solvers import METHODS, compare, estimates


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--system', type=Path, help='frozen matrix.npz/vectors.npz/metadata.json directory')
    p.add_argument('--mesh', type=Path)
    p.add_argument('--ground', type=Path)
    p.add_argument('--case', type=Path)
    p.add_argument('--gp-tangents', type=Path, help='NPZ key tangent: [element, GP, 6, 6], original algorithmic tensors')
    p.add_argument('--methods', nargs='+', choices=METHODS, default=METHODS)
    p.add_argument('--maximum-iterations', type=int, default=3000)
    p.add_argument('--restart', type=int, default=60)
    args = p.parse_args()
    if args.output.exists():
        p.error('output already exists')
    if args.system and any((args.mesh, args.ground, args.case, args.gp_tangents)):
        p.error('--system cannot be combined with reconstruction inputs')
    if any((args.mesh, args.ground, args.case)) and not all((args.mesh, args.ground, args.case)):
        p.error('--mesh, --ground and --case must be supplied together')
    if args.gp_tangents and not args.mesh:
        p.error('--gp-tangents requires explicit mesh/ground/case')
    if args.system:
        a, b, free, meta = load_system(args.system)
    else:
        mesh = read(args.mesh) if args.mesh else griffiths()
        case = read(args.case) if args.case else None
        material = dict(young_kpa=100000., poisson=.3)
        gamma = 20.
        if args.ground:
            ground = read(args.ground)
            init = ground['initialization']
            if (len(ground['materials']) != 1 or init.get('water_level_m') is not None
                    or init.get('surface_surcharge_kpa', 0.) != 0. or init['method'] != 'equilibrate'
                    or set(mesh['material_ids']) != {ground['materials'][0]['id']}):
                p.error('single-material dry gravity-only equilibrate slope required')
            if not args.gp_tangents:
                seed = init['geostatic']['seed']
                if seed['kind'] != 'body_states' or any(
                        np.any(state['effective_stress_compression_positive_kpa'])
                        or np.any(state['plastic_strain_compression_positive'])
                        for state in seed['states'].values()):
                    p.error('elastic seed reconstruction requires zero-stress virgin material')
            material = ground['materials'][0]['mohr_coulomb']
            gamma = ground['materials'][0]['unit_weight_kn_m3']
        tensors = None
        if args.gp_tangents:
            with np.load(args.gp_tangents, allow_pickle=False) as data:
                tensors = data['tangent']
        a, b, free, meta = assemble(mesh, material['young_kpa'], material['poisson'], gamma, case, tensors)
        meta['input_sha256'] = {str(v): sha(v) for v in (args.mesh, args.ground, args.case, args.gp_tangents) if v}
        if not args.mesh:
            meta['input_sha256'] = {'tests/data/griffiths_corner_geometry.json': sha(ROOT/'tests/data/griffiths_corner_geometry.json')}
        meta['rhs_kind'] = 'full_gravity_load; not a failed Newton correction RHS'
    args.output.mkdir(parents=True)
    if not args.system:
        save_system(args.output/'system', a, b, free, meta)
    rows = compare(a, b, free, args.methods, args.maximum_iterations, args.restart, meta.get('external_force_norm'))
    report = dict(schema='soil-fem-griffiths-linear-diagnostic/1',
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                  amg_status='not_run: pyamg absent from requirements-reference.lock.txt; SciPy-only comparison',
                  metadata=meta, estimates=estimates(a), solvers=rows,
                  p6_complete=False, xy2_failure_reproduced=False,
                  nonlinear_force_convergence_proven=False)
    (args.output/'comparison.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for row in rows:
        print(row['method'], row['status'], row['iterations'], row.get('true_linear_relative_residual'), flush=True)
    # Failures remain measured rows. Exit 1 if any requested method fails.
    raise SystemExit(0 if all(row['accepted'] for row in rows) else 1)


if __name__ == '__main__':
    main()
