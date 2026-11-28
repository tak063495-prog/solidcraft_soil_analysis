"""Frozen, explicit constant fresh stresses; no generated K0 or depth profile."""
import copy
import itertools
import math
import numpy as np
from a03shear_binding import study, reference, norm

HEIGHTS = (.5,)
# Voigt is xx, yy, zz, xy, yz, zx (compression positive).
STRESSES = {'zero': (None, 0.)}
for axis, index in (('xz', 5), ('yz', 4)):
    for percent, sign in itertools.product((1, 5), (-1, 1)):
        STRESSES[f'{axis}{percent}{"m" if sign < 0 else "p"}'] = (index, sign*percent*.05)


def plan():
    specs = [dict(name=f'h5-{s}-n{n}-gl5', height=.5, stress=s, n=n, volume=5,
                  minimum=1, reverse=False, element=False, face=False, engines=['native', 'independent'])
             for s, n in itertools.product(STRESSES, (4, 8, 16))]
    specs += [dict(name=f'h5-{s}-n8-gl7', height=.5, stress=s, n=8, volume=7,
                   minimum=1, reverse=False, element=False, face=False, engines=['independent']) for s in STRESSES]
    return dict(passed=False, a03_complete=False, diagnostic_only=True, planned_before_solves=True,
                formal_paths=63, formal_states=2304, common_coordinates='4 phases x i/4, i=1..4; no interpolation',
                shear_applied_to='fresh UpperClay only; LowerClay and inactive UpperClay diagonal constants', specs=specs)


def fixture(n, height, stress):
    if n not in (4, 8, 16) or height not in HEIGHTS or stress not in STRESSES:
        raise ValueError('outside fixed diagnostic matrix; do not solve')
    m, g, c = study.fixture(n, 2.)
    # All six components are explicit in the input; only this fresh entry differs.
    sigma = [2.5, 2.5, 5., 0., 0., 0.]
    index, tau = STRESSES[stress]
    if index is not None:
        sigma[index] = tau
    c['stages'][0]['states']['UpperClay']['effective_stress_compression_positive_kpa'] = sigma
    audit(m, g, c)
    return m, g, c


def audit(m, g, c):
    base = study.fixture(4, 2.)
    assert m == base[0] and g == base[1]
    if [stage.get('kind') for stage in c['stages']] != ['birth', 'joint_birth', 'pressure']:
        raise ValueError('stage sequence outside fixed elastic diagnostic scope; do not solve')
    assert c['settings']['kinematics'] == 'side_rollers' and c.get('boundary') is None
    assert c['joint_system'] == base[2]['joint_system'] and c['initially_inactive_bodies'] == ['UpperClay']
    assert c['stages'][0]['hydraulic_initialization'] == dict(kind='hydrostatic', level_m=2.)
    fresh = c['stages'][0]['states']['UpperClay']
    sigma = fresh['effective_stress_compression_positive_kpa']
    allowed = []
    for index, tau in STRESSES.values():
        v = [2.5, 2.5, 5., 0., 0., 0.]
        if index is not None:
            v[index] = tau
        allowed.append(v)
    if sigma not in allowed:
        raise ValueError('explicit fresh stress outside fixed matrix; do not solve')
    expected = dict(base[2]['stages'][0]['states']['UpperClay'])
    expected['effective_stress_compression_positive_kpa'] = sigma
    assert fresh == expected
    canonical = copy.deepcopy(c)
    canonical['stages'][0]['states']['UpperClay'] = base[2]['stages'][0]['states']['UpperClay']
    canonical['settings']['initial_increment'] = canonical['settings']['minimum_increment'] = .25
    for stage in canonical['stages']: stage['initial_increment'] = .25
    assert canonical == base[2]
    from mc_reference import normalized_yield
    initial_yield = {name: float(normalized_yield(g['materials'][i]['mohr_coulomb'], value['effective_stress_compression_positive_kpa']))
                     for i, (name, value) in enumerate(g['initialization']['geostatic']['seed']['states'].items())}
    initial_yield['fresh UpperClay'] = float(normalized_yield(g['materials'][1]['mohr_coulomb'], sigma))
    if max(initial_yield.values()) > 0.:
        raise ValueError('initial stress outside elastic MC scope; do not solve')
    return dict(allowed=True, nodes=54, tet10=12, connectivity_unchanged=True, height_m=.5,
                explicit_stresses=g['initialization']['geostatic']['seed']['states'], fresh=fresh,
                initial_normalized_mc_yield=initial_yield, fresh_mc=True, single_continuous_b=True,
                side_rollers=True, automatic_k0=False, depth_continuous_stress=False,
                trajectory_scope_checked_at_every_accepted_point=True)


def pre_solve_forces(inputs):
    """Assembly only: no inverse, solve, displacement or new force correction."""
    m = reference.Model(*inputs)
    stress = np.asarray(inputs[2]['stages'][0]['states']['UpperClay']['effective_stress_compression_positive_kpa'])
    diagonal = stress.copy(); diagonal[3:] = 0.
    total = np.zeros(m.nd); delta = np.zeros(m.nd); elements = []
    for e, (q, d, ids, ke) in enumerate(m.ops):
        if m.names[e] != 'UpperClay':
            continue
        force = -np.einsum('qai,a,q->i', q['b'], stress, q['weight'])
        shear = -np.einsum('qai,a,q->i', q['b'], stress-diagonal, q['weight'])
        total[ids] += force; delta[ids] += shear
        elements.append(dict(element=e, global_dof_ids=ids.tolist(), fresh_seed_force_kn=force.tolist(),
                             shear_only_force_kn=shear.tolist(), fresh_l2_kn=norm(force), shear_l2_kn=norm(shear)))
    t, _ = reference.constraint_space(m.mesh, m.case, np.ones(len(m.ops), bool), m.system, True)
    return dict(assembled_before_any_path_solve=True, global_solve_count=0,
                ordering='element input order / local DOF; global node xyz C-order; math.fsum norms',
                constant_stress_element_interior_divergence_kn_m3=[0., 0., 0.],
                note='Constant tensor has zero interior divergence; boundary/interface nodal force and mismatch are nonzero. These are seed contributions at zero strain, not the full placed BIRTH residual.',
                elements=elements, fresh_seed_force_kn=total.reshape(-1, 3).tolist(),
                shear_only_force_kn=delta.reshape(-1, 3).tolist(), fresh_l2_kn=norm(total), shear_l2_kn=norm(delta),
                free_fresh_l2_kn=norm(t.T@total), free_shear_l2_kn=norm(t.T@delta),
                fresh_resultant_kn=[math.fsum(map(float, total.reshape(-1, 3)[:, i])) for i in range(3)],
                treatment='unchanged production seed internal force / existing preparation holding-force release; no redistribution, correction, regularization')
