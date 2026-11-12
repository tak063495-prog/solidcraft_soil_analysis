"""Prepare exactly one pinned xy1/xy2 case without running the solver."""
import argparse
import json
from pathlib import Path

import numpy as np
from .systems import read, ROOT
from prepare_slope_h_series import slope_mesh, GEOMETRY
from plastic_spatial_convergence import inputs
from griffiths_2d import load_reference, digest, GEOMETRY_DIGEST


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xy', type=int, choices=(1, 2), required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reference-2d', type=Path, help='explicit archived F=1.35 material; otherwise unreduced F=1')
    args = p.parse_args()
    if args.output.exists(): p.error('output already exists')
    geometry = read(GEOMETRY)
    if digest(geometry) != GEOMETRY_DIGEST: raise ValueError('pinned geometry changed')
    mesh, _ = slope_mesh(geometry, args.xy, 1)
    _, ground, case = inputs(2)
    zero = ground['initialization']['geostatic']['seed']['states']['Soil']
    ground['initialization']['geostatic']['seed']['states'] = {'Slope': zero}
    ground['bodies'] = [dict(body='Slope', material_id=1)]
    ground['materials'][0]['mohr_coulomb'].update(young_kpa=100000., poisson=.3, cohesion_kpa=10., friction_deg=20.)
    factor = 1.
    if args.reference_2d:
        reference = load_reference(args.reference_2d)
        ground['materials'][0]['mohr_coulomb'].update({k: v for k, v in reference['material'].items() if k != 'unit_weight_kn_m3'})
        ground['materials'][0]['unit_weight_kn_m3'] = reference['material']['unit_weight_kn_m3']
        zero['void_ratio'] = .8
        factor = reference['material_reduction_factor']
    ground['materials'][0]['name'] = f'Griffiths diagnostic; fixed reduction factor {factor}'
    case['boundary'] = dict(constraints=[dict(node_set=name, direction=np.eye(3)[axis].tolist(), value_m=0.)
        for name, axis in [('all', 1), ('fixed_x', 0), ('fixed_z', 2)]])
    case['settings']['initial_increment'] = 1/32
    case['stages'] = [dict(kind='boundary', name='retain_full_gravity_state', target_factor=0., initial_increment=1.)]
    args.output.mkdir(parents=True)
    for name, value in [('mesh', mesh), ('ground', ground), ('case', case)]:
        (args.output/f'{name}.json').write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(prepared=True, solved=False, elements=len(mesh['tetrahedra10']),
                         nodes=len(mesh['nodes']), material_reduction_factor=factor, p6_complete=False)))


if __name__ == '__main__':
    main()
