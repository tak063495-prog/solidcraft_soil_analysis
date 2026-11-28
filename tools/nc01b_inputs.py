"""Generate one explicitly selected Fourie first-lift subproblem, in m/kN/kPa."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import fem_reference as fem
from make_free_contact import make
from refine_tet10 import refine


ROOT = Path(__file__).resolve().parents[1]
PROFILES = ('k0', 'roller', 'pressure')
DEPTH, CUT, HALF_WIDTH = 13.5, 1., 5.5
GAMMA, WATER = 18.5, 9.80665
SIN_PHI = math.sin(math.radians(31))
K0 = 1 - SIN_PHI
M = 6 * SIN_PHI / (3 - SIN_PHI)
# Vertical OCR = 1: place the actual anisotropic stress on the MCC ellipse.
Q_OVER_P = 3 * (1 - K0) / (1 + 2 * K0)
PC_OVER_P = 1 + (Q_OVER_P / M)**2


def mesh_for(elements):
    if elements not in (12, 96):
        raise ValueError('NC-01b supports only 12 or 96 Tet10')
    mesh = json.loads((ROOT / 'examples/seepage-column-tet10-mesh.json').read_bytes())
    xyz = np.array(mesh['nodes'])
    z = xyz[:, 2].copy()
    xyz[:, 0] *= HALF_WIDTH
    xyz[:, 2] = np.where(z <= .5, 2*z*(DEPTH-CUT), DEPTH-CUT + (2*z-1)*CUT)
    mesh['nodes'] = xyz.tolist()
    for region in mesh['regions']:
        lower, upper = (0., DEPTH-CUT) if region['body'] == 'LowerClay' else (DEPTH-CUT, DEPTH)
        region['bounds_m'] = [[0., 0., lower], [HALF_WIDTH, 1., upper]]
        region['volume_m3'] = HALF_WIDTH * (upper-lower)
        region['source_sha256'] = hashlib.sha256(json.dumps(region['bounds_m']).encode()).hexdigest()
    for name, axis, value in [('xmin', 0, 0.), ('xmax', 0, HALF_WIDTH),
                               ('ymin', 1, 0.), ('ymax', 1, 1.)]:
        mesh['node_sets'][name] = np.flatnonzero(xyz[:, axis] == value).tolist()
    if elements == 96:
        mesh, _ = refine(mesh)
    return mesh


def supports(profile):
    axes = {'xmin': [0], 'ymin': [1], 'ymax': [1], 'bottom': [0, 1, 2]}
    if profile != 'pressure':
        axes['xmax'] = [0]
    if profile == 'k0':
        axes['bottom'] = [2]
    return [dict(node_set=name, direction=np.eye(3)[axis].tolist(), value_m=0.)
            for name, fixed in axes.items() for axis in fixed]


def inputs(profile, elements):
    if profile not in PROFILES:
        raise ValueError('unknown NC-01b physical profile')
    mesh = mesh_for(elements)
    _, ground, case = make(mcc=True, pc=100., preload_kpa=5.)
    case.pop('joint_system')  # Both lifts belong to one continuous clay deposit.
    case['monitor_sets'] = {'top': 'top'}
    case['monitors'] = {'floor': DEPTH-CUT}
    case['boundary']['constraints'] = supports(profile)
    init = ground['initialization']
    init.update(surface_surcharge_kpa=0., water_level_m=DEPTH, water_unit_weight_kn_m3=WATER)
    bodies = ground['bodies']
    for material in ground['materials']:
        material.update(name='Fourie uniform soft clay', unit_weight_kn_m3=GAMMA,
                        mcc=dict(critical_state_m=M, **{'lambda': .115}, kappa=.026,
                                 poisson=.25, reference_pressure_kpa=1.,
                                 specific_volume_ncl_at_reference=2.85))
    void = {'mode': 'from_mcc_reference'}
    if profile == 'k0':
        init.pop('geostatic')
        init['method'] = 'horizontal_k0'
        ground['layers'] = [dict(**body, k0_x=K0, k0_y=K0, void_ratio=void,
                                 consolidation=dict(definition='mcc_pc_over_mean_effective',
                                                    ratio=PC_OVER_P)) for body in bodies]
        del ground['bodies']
    else:
        geo = init['geostatic']
        effective_weight = GAMMA-WATER
        gradient = [[0., 0., -ratio*effective_weight] for ratio in (K0, K0, 1., 0., 0., 0.)]
        geo['seed'] = dict(kind='mcc_affine_fields', fields=[
            dict(body=b['body'], reference_position_m=[0., 0., DEPTH],
                 effective_stress_at_reference_compression_positive_kpa=[0.]*6,
                 effective_stress_gradient_kpa_per_m=gradient,
                 pc_over_mean_effective=PC_OVER_P, void_ratio=void) for b in bodies])
        pressure = dict(kind='pressure', face_set='xmax', pressure_kpa=0.,
                        gradient_kpa_per_m=[0., 0., -(K0*effective_weight + WATER)],
                        reference_point_m=[0., 0., DEPTH])
        geo['base_loads'] = dict(surface_loads=[pressure] if profile == 'pressure' else [])
    return mesh, ground, case


def initial_balance(mesh, ground, case):
    """Independently integrate the prescribed total stresses and physical loads."""
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    internal = np.zeros_like(mesh['nodes'], dtype=float)
    physical = np.zeros_like(internal)
    for ei in range(len(active)):
        gp = fem.volume(mesh, ei, independent=True)
        depth = DEPTH - gp['position'][:, 2]
        sv = GAMMA*depth
        sh = (K0*(GAMMA-WATER)+WATER)*depth
        stress = np.column_stack((sh, sh, sv, *[np.zeros_like(sv)]*3))
        internal[gp['nodes']] -= np.einsum('q,qai,qa->i', gp['weight'], gp['b'], stress).reshape(-1, 3)
        physical[gp['nodes'], 2] -= GAMMA * (gp['weight'] @ gp['shape'])
    loads = ground['initialization'].get('geostatic', {}).get('base_loads', {})
    physical += fem.conditions_load(mesh, loads, active)
    fixed = np.zeros_like(internal, dtype=bool)
    for bc in case['boundary']['constraints']:
        fixed[mesh['node_sets'][bc['node_set']], bc['direction'].index(1.)] = True
    residual = internal - physical
    reaction = np.where(fixed, residual, 0.).sum(axis=0)
    applied = physical.sum(axis=0)
    ratio = float(np.linalg.norm(residual[~fixed]) / np.linalg.norm(physical))
    tol = ground['initialization']['relative_equilibrium_tolerance']
    if ratio > tol or np.linalg.norm(reaction + applied) > tol*np.linalg.norm(physical):
        raise ValueError('prescribed physical seed is not in equilibrium')
    return dict(free_residual_norm_kn=float(np.linalg.norm(residual[~fixed])),
                relative_free_residual=ratio, applied_force_kn=applied.tolist(),
                support_reaction_kn=reaction.tolist(), imbalance_kn=(reaction+applied).tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=PROFILES, required=True)
    parser.add_argument('--elements', choices=(12, 96), type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mesh, ground, case = inputs(args.profile, args.elements)
    audit = initial_balance(mesh, ground, case)
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in [('mesh', mesh), ('ground', ground), ('stages', case), ('initial-balance', audit)]:
        (args.output / (name+'.json')).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n',
                                                 encoding='utf-8', newline='\n')


if __name__ == '__main__':
    main()
