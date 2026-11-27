"""Explicit unittest suite; deliberately outside existing test_*.py discovery."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from reference import Column, Darcy, cases
from numerical import integrate_forces, solve_head, head_errors, series, fitted_mesh
from compare import diagnose, synthetic, write_csv_bundle, read_export, read_json
from policies import trace, reservoir_case, mcc_invariants


class ClosedForms(unittest.TestCase):
    def test_hydrostatic_bottom_hand_value(self):
        c = Column(water=3.8)
        f = c.fields(0.)
        self.assertAlmostEqual(float(f['pressure']), 37.278)
        self.assertAlmostEqual(float(f['total_vertical']), 72.848)
        self.assertAlmostEqual(float(f['effective_vertical']), 35.57)

    def test_partial_new_layer_hand_values(self):
        c = Column()
        f = c.fields(1.2)
        self.assertAlmostEqual(float(f['pressure']), 8.829)
        self.assertAlmostEqual(float(f['total_vertical']), 38.3)
        self.assertAlmostEqual(float(f['effective_horizontal']), .55*(38.3-8.829))

    def test_top_pressure_and_surcharge(self):
        for c, flow in cases().values():
            p = c.pressure(c.top) if flow is None else flow.pressure(c.top)
            self.assertAlmostEqual(float(c.fields(c.top, p)['effective_vertical']), c.surcharge)

    def test_pressure_continuity_at_water_table(self):
        c = Column()
        self.assertEqual(c.pressure(c.water), 0.)
        self.assertAlmostEqual(float(c.total_vertical(c.water)), 20.3)

    def test_k0_is_explicit_not_empirical(self):
        c = Column(k0=.72)
        self.assertAlmostEqual(float(c.fields(1.)['effective_horizontal']), .72*float(c.fields(1.)['effective_vertical']))

    def test_two_layer_resistance_hand_value(self):
        d = Darcy()
        self.assertAlmostEqual(d.flux, .8/150000.)
        self.assertAlmostEqual(float(d.head(1.)), 3.866666666666667)

    def test_darcy_boundaries_and_interfaces(self):
        for c, d in cases().values():
            if d is None:
                continue
            self.assertAlmostEqual(float(d.head(0)), d.bottom_head)
            self.assertAlmostEqual(float(d.head(3)), d.top_head)
            for z in d.edges[1:-1]:
                self.assertLess(abs(float(d.head(z-1e-10)-d.head(z+1e-10))), 1e-9)

    def test_seepage_and_pressure_gradient(self):
        d = Darcy()
        for z in (.5, 2.):
            derivative = float((d.pressure(z+1e-5)-d.pressure(z-1e-5))/2e-5)
            self.assertAlmostEqual(derivative, d.pressure_gradient(z), places=8)
            self.assertAlmostEqual(float(-derivative-d.gamma_w), d.seepage(z), places=8)

    def test_no_negative_pressure_clipping(self):
        with self.assertRaises(ValueError):
            Darcy(bottom_head=0., top_head=0.).pressure(2.)

    def test_bad_layers_rejected(self):
        with self.assertRaises(ValueError):
            Darcy(conductivity=(0., 1e-5))


class IndependentNumerics(unittest.TestCase):
    def test_fitted_is_exact_not_claimed_second_order(self):
        for row in series(Darcy()):
            if row['fitted']:
                self.assertLess(row['head_l2'], 1e-11)
                self.assertIsNone(row['order'])

    def test_unfitted_first_order(self):
        rows = [r for r in series(Darcy()) if not r['fitted']]
        self.assertTrue(all(a['head_l2'] > b['head_l2'] for a, b in zip(rows[:-1], rows[1:])))
        self.assertTrue(all(.95 < r['order'] < 1.05 for r in rows[1:]))

    def test_three_layer_sparse_independent(self):
        d = cases()['darcy_three_layers'][1]
        z = fitted_mesh(d)
        result = solve_head(d, z)
        self.assertLess(head_errors(d, result)['head_l2'], 1e-11)
        self.assertLess(abs(result['water_defect']), 1e-15)

    def test_near_duplicate_geometry_rejected(self):
        with self.assertRaises(ValueError):
            solve_head(Darcy(), np.array([0., 1., 1.+1e-15, 3.]))

    def test_flux_conservation_unfitted_not_only_analytic(self):
        d = Darcy()
        result = solve_head(d, np.linspace(0, 3, 18))
        self.assertLess(abs(result['water_defect']), 1e-15)
        self.assertLess(np.ptp(result['flux']), 1e-15)
        self.assertLess(result['outward'][0], 0.)
        self.assertGreater(result['outward'][1], 0.)

    def test_all_node_force_identities(self):
        for c, flow in cases().values():
            metrics = integrate_forces(c, flow)['metrics']
            self.assertLess(metrics['pore_identity_max'], 1e-11)
            self.assertLess(metrics['nodal_equilibrium_max'], 1e-11)
            self.assertLess(max(abs(np.array(metrics['global_force']))), 1e-11)
            self.assertLess(abs(metrics['global_moment']), 1e-11)

    def test_continuous_interface_before_after(self):
        c, flow = cases()['partial_new_layer']
        result = integrate_forces(c, flow)
        m = result['metrics']
        self.assertEqual(m['continuous_pair_max'], 0.)
        np.testing.assert_allclose(m['before_old_exposed_top'], [[0., -4.4145], [0., -4.4145]], atol=1e-13)
        np.testing.assert_allclose(m['before_old_exposed_top'], m['after_old_interface'], atol=1e-13, rtol=0.)
        self.assertGreater(np.linalg.norm(result['forces']['pressure_interface']), 0.)

    def test_buoyancy_plus_seepage(self):
        c, d = cases()['darcy_upward']
        m = integrate_forces(c, d)['metrics']
        self.assertAlmostEqual(m['pressure_boundary_resultant'][1], 37.278)
        self.assertAlmostEqual(m['seepage_resultant'][1], 7.848)

    def test_numerical_pressure_force_balance(self):
        c, d = cases()['darcy_upward']
        solution = solve_head(d, np.linspace(0, 3, 18))
        m = integrate_forces(c, d, numerical_result=solution)['metrics']
        self.assertLess(m['pore_identity_max'], 1e-11)
        self.assertLess(m['nodal_equilibrium_max'], 1e-11)

    def test_split_waterline_quadrature_independent(self):
        c = Column(water=2.137)
        a = integrate_forces(c, quadrature=2)
        b = integrate_forces(c, quadrature=4)
        for key in a['forces']:
            np.testing.assert_allclose(a['forces'][key], b['forces'][key], atol=1e-12, rtol=0.)

    def test_balanced_corruption_must_fail_all_nodes(self):
        c, d = cases()['partial_new_layer']
        result = integrate_forces(c, d)
        force = result['forces']['pressure_exterior'].copy()
        force[0, 1] += 1.
        force[1, 1] -= 1.
        np.testing.assert_allclose(force.sum(axis=0), result['forces']['pressure_exterior'].sum(axis=0))
        self.assertEqual(np.max(abs(force-result['forces']['pressure_exterior'])), 1.)


class ExportDiagnostics(unittest.TestCase):
    def setUp(self):
        self.data, self.before = synthetic('partial_new_layer')

    def check(self, data=None):
        return diagnose(self.data if data is None else data, self.before, 'partial_new_layer')

    def test_roundtrip_never_accepted(self):
        report = self.check()
        self.assertTrue(report['valid_export'])
        self.assertIsNone(report['accepted'])
        self.assertFalse(report['a03_complete'])
        self.assertEqual(report['metrics']['pressure_kpa']['max_abs'], 0.)

    def test_csv_json_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_csv_bundle(self.data, tmp)
            self.assertEqual(self.check(read_export(tmp)), self.check())

    def test_undefined_pressure_rejected(self):
        del self.data['nodes'][0]['pressure_kpa']
        with self.assertRaises(KeyError):
            self.check()

    def test_filled_pressure_rejected(self):
        self.data['pressure_filled'] = True
        with self.assertRaises(ValueError):
            self.check()

    def test_missing_provenance_rejected(self):
        del self.data['pressure_filled']
        with self.assertRaises(KeyError):
            self.check()

    def test_origin_interpolation_rejected(self):
        self.data['nodes'][0]['pressure_origin'] = 'old_element_interpolation'
        with self.assertRaises(ValueError):
            self.check()

    def test_old_history_transfer_rejected(self):
        self.data['materials']['new'][0]['source_id'] = 'old-gp-0'
        with self.assertRaises(ValueError):
            self.check()

    def test_new_aged_history_rejected(self):
        self.data['materials']['new'][0]['plastic_strain'][0] = .001
        with self.assertRaises(ValueError):
            self.check()

    def test_old_history_rewrite_rejected(self):
        # deepcopy avoids changing the separately trusted before snapshot.
        self.data = copy.deepcopy(self.data)
        self.data['materials']['old'][0]['plastic_work_kpa'] = 0.
        with self.assertRaises(ValueError):
            self.check()

    def test_mcc_remains_spec_undefined(self):
        self.data['materials']['new'][0]['model'] = 'MCC'
        self.data['materials']['new'][0].update(effective_stress_kpa=[20., 20., 20., 0., 0., 0.],
                                               pc_kpa=50., void_ratio=.8, M=1.2)
        report = self.check()
        self.assertTrue(any('MCC' in x and '仕様未定義' in x for x in report['pending']))
        self.assertIsNone(report['accepted'])

    def test_timeout_cancelled_incomplete_rejected(self):
        for status in ('timeout', 'cancelled', 'running'):
            self.data['status'] = status
            with self.assertRaises(ValueError):
                self.check()

    def test_duplicate_node_rejected(self):
        self.data['nodes'].append(self.data['nodes'][0])
        with self.assertRaises(ValueError):
            self.check()

    def test_nan_rejected(self):
        self.data['nodes'][0]['pressure_kpa'] = float('nan')
        with self.assertRaises(ValueError):
            self.check()

    def test_invented_unsaturated_head_rejected(self):
        self.data['nodes'][-1]['head_m'] = 2.1
        with self.assertRaises(ValueError):
            self.check()

    def test_nodal_cancellation_not_hidden_by_global_sum(self):
        self.data['forces'][0]['pressure_exterior_z'] += 1.
        self.data['forces'][1]['pressure_exterior_z'] -= 1.
        report = self.check()
        self.assertEqual(report['metrics']['pressure_exterior_kn']['max_abs'], 1.)
        self.assertGreater(report['nodal_equilibrium_max_kn'], .99)
        self.assertIsNone(report['accepted'])

    def test_water_imbalance_reported(self):
        self.data['fluxes'][0]['outward_m3_s'] = .01
        report = self.check()
        self.assertEqual(report['water_defect_m3_s'], .01)
        self.assertIsNone(report['accepted'])

    def test_duplicate_json_key_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'x.json'
            path.write_text('{"x":1,"x":2}', encoding='utf-8')
            with self.assertRaises(ValueError):
                read_json(path)


class PolicyDiagnostics(unittest.TestCase):
    def test_soil_structure_internal_pressure_is_zero(self):
        row = trace('soil', 1.2, [0., 1.], head_a=4., soil_side='a')
        self.assertEqual(row['pore_b_kpa'], 0.)
        self.assertAlmostEqual(row['pressure_kpa'], 27.468)
        np.testing.assert_array_equal(np.array(row['traction_a'])+row['traction_b'], [0., 0.])

    def test_continuous_conflict_not_averaged(self):
        with self.assertRaises(ValueError):
            trace('continuous', 1., [0., 1.], head_a=3., head_b=4.)

    def test_prescribed_conflict_rejected(self):
        with self.assertRaises(ValueError):
            trace('prescribed', 1., [0., 1.], head_a=3., level=4.)

    def test_reservoir_exchange_included(self):
        result = reservoir_case()
        self.assertNotEqual(result['bottom_out']+result['top_out'], 0.)
        self.assertAlmostEqual(result['water_defect'], 0.)
        self.assertAlmostEqual(result['reservoir_out'], -5.555555555555556e-6)

    def test_mcc_invariant_only_no_birthing_rule(self):
        row = mcc_invariants(dict(effective_stress_kpa=[20., 20., 20., 0., 0., 0.],
                                 pc_kpa=50., void_ratio=.8, M=1.2))
        self.assertEqual(row['p_kpa'], 20.)
        self.assertEqual(row['q_kpa'], 0.)
        self.assertAlmostEqual(row['yield_function_kpa2'], -864., places=11)
        self.assertIsNone(row['accepted'])

    def test_mcc_zero_stress_is_not_auto_initialized(self):
        with self.assertRaises(ValueError):
            mcc_invariants(dict(effective_stress_kpa=[0.]*6, pc_kpa=50., void_ratio=.8, M=1.2))


if __name__ == '__main__':
    unittest.main(verbosity=2)
