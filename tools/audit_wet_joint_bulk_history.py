"""Reassemble every accepted state of the existing small wet JOINT reference.

This extends contact-history replay to bulk compatibility, input seed response,
fixed water, prescribed constraints, reactions and free equilibrium. It does
not certify convergence, native integration, bulk plasticity or new JOINT laws.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import fem_reference as fem
from audit_wet_joint_accepted_history import audit as contact_audit, history
from hydrostatic_cut_volume_reference import force as pore_force
from joint_water_reference import Water
from mc_reference import normalized_yield
from verify_high_order import physical


def audit(folder, result, values):
    assert result['volume_rule'] == 'affine_cut_hydrostatic' and result['volume_order'] == 7, 'fixed affine-cut GL7 reference only'
    report = contact_audit(folder, result, values)
    mesh, ground, case = values
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    water = Water(*values); system = water.contact
    level = ground['initialization']['water_level_m']
    cfg = ground['initialization']['geostatic']
    materials = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
    regions = {r['id']: r['body'] for r in mesh['regions']}
    interface = water.hydrostatic(active, level)
    base, _ = physical(mesh, ground, active, ground['initialization']['surface_surcharge_kpa'],
                       case.get('pressure_on_original_surface_only', False), joint_system=system)
    base = (base+interface['load']).ravel()
    boundary = fem.conditions_load(mesh, case['boundary'], active).ravel()
    operators = []; seed_force = np.zeros_like(base); hydro = np.zeros_like(base); start = 0
    for ei, mid in enumerate(mesh['material_ids']):
        q = fem.volume(mesh, ei)
        ids = (3*q['nodes'][:, None]+np.arange(3)).ravel()
        mat = materials[mid]; d = fem.elastic_tensor(mat['young_kpa'], mat['poisson'])
        seed = cfg['seed']['states'][regions[mesh['region_ids'][ei]]]
        s0 = np.asarray(seed['effective_stress_compression_positive_kpa'])
        count = len(q['weight']); section = slice(start, start+count); start += count
        pf = pore_force(mesh, ei, level, water.gamma)['force'].ravel()
        hydro[ids] += pf
        seed_force[ids] -= np.einsum('qai,a,q->i', q['b'], s0, q['weight'])
        seed_force[ids] += pf
        operators.append((q, ids, mat, d, seed, s0, section))
    holding = seed_force-base
    xyz = np.asarray(mesh['nodes'])
    maximum = dict(force=0., constraint=0., stress=0., strain=0., void=0., pore=0., reaction=0.,
                   physical_pressure=0., water_load=0.)
    origin = 0.; stages = {'geostatic': (0., 0.)}
    for stage in case['stages']:
        stages[stage['name']] = (origin, stage['target_factor']); origin = stage['target_factor']
    rows = []
    for row, fields in history.accepted_states(folder):
        origin, target = stages[row['stage']]
        factor = origin+(target-origin)*row['progress']
        assert row['factor'] == factor, 'input boundary-factor schedule differs'
        u = fields['u'].ravel(); internal = hydro.copy()
        def compare(key, observed, expected, atol, rtol=0.):
            maximum[key] = max(maximum[key], float(np.max(np.abs(observed-expected))))
            np.testing.assert_allclose(observed, expected, atol=atol, rtol=rtol,
                                       err_msg=f'accepted bulk/water {key} at {row["stage"]}:{row["progress"]}')
        for q, ids, mat, d, seed, s0, section in operators:
            strain = np.einsum('qai,i->qa', q['b'], u[ids])
            stress = s0-strain@d.T
            compare('strain', fields['strain'][section], strain, 1e-10)
            compare('stress', fields['stress'][section], stress, 1e-7, 1e-9)
            compare('void', fields['void'][section], (1+seed['void_ratio'])*np.exp(strain[:, :3].sum(axis=1))-1, 1e-8)
            compare('pore', fields['pore'][section], water.gamma*np.maximum(level-q['position'][:, 2], 0), 1e-8)
            assert np.max(normalized_yield(mat, fields['stress'][section])) <= 1e-8, 'bulk plasticity outside reference'
            internal[ids] -= np.einsum('qai,qa,q->i', q['b'], fields['stress'][section], q['weight'])
        compare('physical_pressure', fields['physical_pressure'], interface['pressure'], 1e-8)
        compare('water_load', fields['water_load'], interface['load'], 1e-8)
        internal += np.einsum('qag,qa,q->g', system.c, fields['traction'], system.weight)
        load = base+(1-row['progress'])*holding if row['stage'] == 'geostatic' else base+factor*boundary
        reaction = internal-load
        compare('reaction', fields['reaction'], reaction.reshape(-1, 3), 1e-7, 1e-9)
        state = dict(active_elements=active, boundary_factor=factor, displacement_m=fields['u'])
        basis, _ = fem.staged_constraint_space(mesh, case, state)
        reduced = np.concatenate([b.T@reaction.reshape(-1, 3)[n] for n, b in basis.items()])
        residual = float(np.linalg.norm(reduced)/max(np.linalg.norm(load), 1.))
        assert residual < 1e-10, 'accepted reassembled free equilibrium'
        maximum['force'] = max(maximum['force'], residual)
        for condition in case['boundary']['constraints']:
            nodes = mesh['node_sets'][condition['node_set']]
            expected = factor*(condition['value_m']+(xyz[nodes]-condition.get('reference_point_m', [0., 0., 0.]))@np.asarray(condition.get('gradient', [0., 0., 0.])))
            actual = fields['u'][nodes]@np.asarray(condition['direction'])
            compare('constraint', actual, expected, 1e-12)
        rows.append(dict(stage=row['stage'], progress=row['progress'], factor=factor,
                         reassembled_relative_force_residual=residual))
    report.update(bulk_water_consistent=True, input_seed_holding_kn=holding.tolist(),
                  bulk_water_maximum_absolute_differences=maximum, reassembled_rows=rows,
                  independent_bulk_plastic_solution=False, native_comparison=False)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history', type=Path, required=True)
    p.add_argument('--result', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error('new output required')
    report = audit(args.history, json.loads(args.result.read_text(encoding='utf-8')), history.inputs())
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
