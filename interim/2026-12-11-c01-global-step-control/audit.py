"""Previous all-field audit, with adaptive coordinates and the original force gates."""
from study import operators, states, STAGES, np, probe, fem

def audit(folder, mesh, ground, case, rule, steps):
    op = operators(mesh, ground, case, rule); s = op.system
    count = 0; residual = 0.; pressure_delta = 0.; pore_delta = 0.; stress_delta = 0.; prev = None
    for row, f in states(folder):
        phase = STAGES.index(row['stage'])
        origin = 0. if phase <= 1 else case['stages'][phase-2]['target_factor']
        expected_factor = 0. if phase == 0 else origin+(case['stages'][phase-1]['target_factor']-origin)*row['progress']
        assert row['factor'] == expected_factor, 'input boundary-factor schedule differs'
        count += 1
        np.testing.assert_allclose(f['physical_pressure'], op.interface['pressure'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(f['water_load'], op.interface['load'], atol=1e-8, rtol=0)
        for condition in case['boundary']['constraints']:
            nodes = mesh['node_sets'][condition['node_set']]
            xyz = np.asarray(mesh['nodes'])[nodes]
            expected = expected_factor*(condition['value_m']+(xyz-condition.get('reference_point_m',[0.,0.,0.]))@np.asarray(condition.get('gradient',[0.,0.,0.])))
            np.testing.assert_allclose(f['u'][nodes]@np.asarray(condition['direction']),expected,atol=1e-12,rtol=0)
        u = f['u'].ravel(); q = f['traction']; gap = (f['jump']*s.normal).sum(axis=1); qn = (q*s.normal).sum(axis=1)
        np.testing.assert_allclose(f['jump'], np.einsum('qag,g->qa', s.c, u), atol=1e-10, rtol=0)
        np.testing.assert_allclose(qn, s.kn*np.minimum(gap, 0), atol=1e-7, rtol=1e-9)
        tau = q-qn[:, None]*s.normal
        assert np.all(np.linalg.norm(tau, axis=1) <= s.mu*np.maximum(-qn, 0)+1e-7)
        np.testing.assert_allclose(f['endpoint_work'], f['elastic_energy']+f['plastic_dissipation']+f['numerical_dissipation'], atol=1e-9, rtol=1e-8)
        assert not f['plastic_strain'].any(), 'elastic bulk branch only'
        assert np.min(f['accumulated']) >= 0 and np.all(f['accumulated']+1e-10 >= np.linalg.norm(f['plastic_slip'], axis=1))
        if prev is not None:
            assert np.all(f['accumulated']+1e-12 >= prev['accumulated'])
            pore_delta = max(pore_delta, float(np.max(abs(f['pore']-prev['pore']))))
            pressure_delta = max(pressure_delta, float(np.max(abs(f['physical_pressure']-prev['physical_pressure']))))
        force = op.k@u+op.seed+np.einsum('qag,qa,q->g', s.c, q, s.weight)
        load = op.base+(1-row['progress'])*op.holding if row['stage'] == 'geostatic' else op.base+row['factor']*op.boundary
        residual = max(residual, float(np.linalg.norm(op.reduce(force-load))/max(np.linalg.norm(load), 1.)))
        start = 0
        for ei, ((volume, ids), mid) in enumerate(zip(op.ops, mesh['material_ids'])):
            mat = next(m['mohr_coulomb'] for m in ground['materials'] if m['id'] == mid)
            d = fem.elastic_tensor(mat['young_kpa'], mat['poisson']); n = len(volume['weight'])
            # Both bodies have the same input isotropic effective seed; no history is injected.
            body = next(r['body'] for r in mesh['regions'] if r['id'] == mesh['region_ids'][ei])
            seed = np.asarray(ground['initialization']['geostatic']['seed']['states'][body]['effective_stress_compression_positive_kpa'])
            expected_stress = seed-f['strain'][start:start+n]@d.T
            np.testing.assert_allclose(f['stress'][start:start+n],expected_stress,atol=1e-7,rtol=1e-9)
            np.testing.assert_allclose(f['strain'][start:start+n],np.einsum('qai,i->qa',volume['b'],u[ids]),atol=1e-10,rtol=0)
            initial=ground['initialization']['geostatic']['seed']['states'][body]
            np.testing.assert_allclose(f['void'][start:start+n],(1+initial['void_ratio'])*np.exp(f['strain'][start:start+n,:3].sum(axis=1))-1,atol=1e-8,rtol=0)
            np.testing.assert_allclose(f['pore'][start:start+n],op.water.gamma*np.maximum(op.level-volume['position'][:,2],0),atol=1e-8,rtol=0)
            assert np.max(probe.normalized_yield(mat,f['stress'][start:start+n]))<=1e-8
            stress_delta = max(stress_delta, float(np.max(abs(f['stress'][start:start+n]-expected_stress))))
            start += n
        prev = f
    assert count >= 7*steps, 'incomplete path'
    native = folder.name == 'native'
    tolerance = 1e-8 if native else 1e-10
    assert residual < tolerance, ('existing force audit', residual, tolerance)
    return dict(accepted_states=count, maximum_reassembled_force_residual=residual, existing_force_audit_tolerance=tolerance,
                maximum_pressure_change_kpa=pressure_delta, maximum_pore_change_kpa=pore_delta,
                maximum_elastic_stress_compatibility_error_kpa=stress_delta,
                all_saved_states_audited=True, numerical_acceptance=False)
