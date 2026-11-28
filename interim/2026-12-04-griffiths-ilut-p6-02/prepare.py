"""New inputs only. Preserves the fixed domain/material; no reference fitting."""
import argparse
import copy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from plastic_spatial_convergence import inputs, read, write
from prepare_slope_h_series import GEOMETRY, slope_mesh


def prepare(out, stride=5, ilut=False):
    geometry = read(GEOMETRY)
    if stride not in (1, 2, 5):
        raise ValueError('only fixed 20x10 grid strides 1/2/5')
    q = geometry['quads']
    geometry['quads'] = [[q[r*20+c][0], q[r*20+c+stride-1][1],
                          q[(r+stride-1)*20+c+stride-1][2],
                          q[(r+stride-1)*20+c][3]]
                         for r in range(0, 10, stride) for c in range(0, 20, stride)]
    mesh, _ = slope_mesh(geometry, 1, 1)
    _, ground, case = inputs(2)
    seed = ground['initialization']['geostatic']['seed']['states']['Soil']
    ground['initialization']['geostatic']['seed']['states'] = {'Slope': seed}
    ground['bodies'] = [dict(body='Slope', material_id=1)]
    ground['materials'][0]['name'] = 'Griffiths F=1, ILUT comparison only'
    ground['materials'][0]['mohr_coulomb'].update(
        young_kpa=100000., poisson=.3, cohesion_kpa=10., friction_deg=20.)
    case['boundary'] = dict(constraints=[dict(node_set=n, direction=d, value_m=0.)
        for n, d in [('all', [0.,1.,0.]), ('fixed_x',[1.,0.,0.]), ('fixed_z',[0.,0.,1.])]])
    case['settings']['initial_increment'] = 1/32
    case['stages'] = [dict(kind='boundary', name='retain_full_gravity_state',
                           target_factor=0., initial_increment=.25)]
    if ilut:
        case['settings']['linear_solver']['options']['ilut'] = dict(drop_tolerance=1e-4, fill_limit=200)
    out.mkdir(parents=True, exist_ok=False)
    for name, value in [('mesh',mesh), ('ground',ground), ('case',case)]:
        write(out/(name+'.json'), value)
    return mesh, ground, case


def add_ilut(source, target, drop=1e-4, fill=200):
    case = copy.deepcopy(read(source))
    options = case['settings']['linear_solver']
    if options['method'] != 'sparse_gmres_ilu0':
        raise ValueError('preserve the existing sparse GMRES method')
    options['options']['ilut'] = dict(drop_tolerance=drop, fill_limit=fill)
    write(target, case)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='mode', required=True)
    a = sub.add_parser('case'); a.add_argument('output', type=Path)
    a.add_argument('--stride', type=int, default=5); a.add_argument('--ilut', action='store_true')
    a = sub.add_parser('add-ilut'); a.add_argument('source', type=Path); a.add_argument('target', type=Path)
    a.add_argument('--drop', type=float, default=1e-4); a.add_argument('--fill', type=int, default=200)
    a = p.parse_args()
    if a.mode == 'case': prepare(a.output, a.stride, a.ilut)
    else: add_ilut(a.source, a.target, a.drop, a.fill)
