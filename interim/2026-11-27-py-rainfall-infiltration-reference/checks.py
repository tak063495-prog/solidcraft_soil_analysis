"""Explicit opt-in tests: outside existing test_*.py discovery and CI counts."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.integrate import quad
from reference import (Gardner, analytical, step_u, step_q, gardner_fv,
                       vg, vg_column, mechanical_fv, mechanics)
from compare import compare, validate, read
from study import make_packet, water_level_exact
from pc_rainfall import solve as cartesian_rainfall


class Physics(unittest.TestCase):
    def test_initial_state_and_pulse_causality(self):
        u, q = analytical(np.array([.1, 1., 4.]), 0.)
        np.testing.assert_array_equal(u, [.2]*3)
        np.testing.assert_allclose(q, [2e-5]*3, rtol=1e-15)
        self.assertEqual(step_u(.1, -1., 2e-5), 0.)

    def test_analytical_darcy_relation(self):
        p, dx = Gardner(), 1e-4
        for x, t in ((.2, 120.), (.7, 450.)):
            u, q = analytical(np.array([x-dx, x, x+dx]), t)
            darcy = p.ks*(u[1]-(u[2]-u[0])/(2*dx*p.alpha))
            self.assertLess(abs(darcy-q[1]), 2e-12)

    def test_analytical_local_conservation(self):
        p, dx, dt = Gardner(), 1e-4, .01
        for x, t in ((.2, 120.), (.7, 450.)):
            before, _ = analytical(np.array([x]), t-dt)
            after, _ = analytical(np.array([x]), t+dt)
            _, q = analytical(np.array([x-dx, x+dx]), t)
            residual = p.delta*(after[0]-before[0])/(2*dt)+(q[1]-q[0])/(2*dx)
            self.assertLess(abs(residual), 1e-11)

    def test_analytical_rain_boundary(self):
        p, dx = Gardner(), 1e-5
        for t in (100., 450.):
            u, _ = analytical(np.array([0., dx, 2*dx]), t)
            derivative = (-3*u[0]+4*u[1]-u[2])/(2*dx)
            rain = 4e-5 if t < 300 else 2e-5
            self.assertLess(abs(p.ks*(u[0]-derivative/p.alpha)-rain), 2e-12)

    def test_analytical_continuous_mass(self):
        p = Gardner()
        integral = quad(lambda x: p.delta*(analytical(np.array([x]), 600.)[0][0]-.2),
                        0, 8, epsabs=2e-11)[0]
        self.assertLess(abs(integral-2e-5*300), 1e-10)

    def test_analytical_laplace_transform(self):
        p, x, s, jump = Gardner(), .2, .01, 2e-5
        numerical = quad(lambda t: np.exp(-s*t)*step_u(x, t, jump),
                         0, 5000, epsabs=1e-9)[0]
        root = np.sqrt(p.v*p.v+4*p.diffusion*s)
        exact = 2*jump/(p.delta*(p.v+root))*np.exp((p.v-root)*x/(2*p.diffusion))/s
        self.assertLess(abs(numerical-exact), 2e-8)

    def test_quadrature_tolerance_sensitivity(self):
        self.assertLess(abs(step_u(.1, 200, 2e-5, eps=2e-8)-step_u(.1, 200, 2e-5, eps=2e-12)), 1e-10)

    def test_fv_no_pulse_uniform_solution(self):
        r = gardner_fv(16, 10, end=100, jumps=())
        np.testing.assert_allclose(r['u'], .2, rtol=0, atol=5e-14)
        self.assertLess(r['max_mass_m'], 1e-12)

    def test_fv_mass_and_space_convergence(self):
        errors = []
        for n in (16, 32, 64):
            r = gardner_fv(n, .5, method='CN')
            u, _ = analytical(r['x'], 600.)
            errors.append(np.linalg.norm(r['psi']-np.log(u))/np.sqrt(n))
            self.assertLess(r['max_mass_m'], 1e-11)
        self.assertGreater(np.log2(errors[-2]/errors[-1]), 1.8)

    def test_fv_time_convergence(self):
        r = [gardner_fv(64, dt)['psi'] for dt in (20., 10., 5.)]
        self.assertGreater(np.log2(np.linalg.norm(r[0]-r[1])/np.linalg.norm(r[1]-r[2])), .85)

    def test_fv_stops_at_pulse_edge(self):
        r = gardner_fv(16, 17.)
        self.assertIn(300., r['history'][:, 0])
        self.assertAlmostEqual(r['history'][-1, 2], .018, places=13)

    def test_invalid_resolution_rejected(self):
        for n, dt in ((2, 1.), (16, -1.)):
            with self.assertRaises(ValueError):
                gardner_fv(n, dt)

    def test_gardner_saturation_branch_rejected(self):
        with self.assertRaises(ValueError):
            gardner_fv(16, 10, jumps=((0., 1.),))

    def test_vg_endpoints_and_storage(self):
        theta, k, storage = vg(np.array([-1., 0., 1.]))
        self.assertLess(theta[0], .45)
        np.testing.assert_allclose(theta[1:], .45)
        np.testing.assert_allclose(k[1:], 1e-6)
        self.assertAlmostEqual(storage[2], 1e-5)

    def test_vg_conservative_analogue_is_not_original_acceptance(self):
        r = vg_column(12, 100.)
        self.assertIsNone(r['accepted'])
        self.assertIsNone(r['original_3d_gate_value'])
        # Unit-area numerical conservation guard; not acceptance of the 3D case.
        self.assertLess(r['max_mass_m'], 1e-9)
        self.assertEqual([f['time_s'] for f in r['frames']], [3600., 7200., 10800.])
        self.assertGreater(r['frames'][-1]['rejected_rain_m'], 0.)

    def test_recharge_causes_water_level_rise(self):
        level = water_level_exact(600.)
        self.assertGreater(level, .2)
        self.assertLess(level, .2+.006/.2)

    def test_mechanical_equilibrium_and_exact_integral(self):
        for n in (4, 8, 16):
            r = mechanical_fv(n, .21)
            self.assertLess(r['max_displacement_error'], 1e-14)
            self.assertLess(r['max_stress_error'], 1e-10)
            self.assertLess(r['equilibrium_residual'], 1e-10)

    def test_hydraulic_stress_error_separation(self):
        a, b = mechanical_fv(8, .2), mechanical_fv(8, .201)
        np.testing.assert_allclose(b['effective_stress']-a['effective_stress'], -.00981, atol=1e-12)

    def test_mechanical_receptor_refuses_unsaturated_extension(self):
        with self.assertRaises(ValueError):
            mechanics(np.array([1.]), .2)

    def test_cartesian_strip_mass_symmetry_and_rain_geometry(self):
        r = cartesian_rainfall(dim=2, n=8, end=200.)
        h = np.asarray(r['frames'][-1]['head_m'])
        self.assertEqual(r['rain_patch_area'], .5)
        self.assertFalse(r['original_square_geometry'])
        np.testing.assert_allclose(h, h[:, ::-1], rtol=0, atol=1e-12)
        self.assertGreater(h[0, 3], h[0, 0])
        self.assertLess(r['maximum_mass_defect'], 1e-10)

    def test_cartesian_3d_requires_pc_opt_in(self):
        with self.assertRaisesRegex(ValueError, 'user-pc'):
            cartesian_rainfall(dim=3, n=4, end=100.)

    def test_cartesian_patch_requires_aligned_grid(self):
        with self.assertRaises(ValueError):
            cartesian_rainfall(dim=2, n=6, end=100.)


class Contract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = make_packet(8)

    def bad(self, mutate):
        candidate = copy.deepcopy(self.packet)
        mutate(candidate)
        with self.assertRaises((ValueError, KeyError, TypeError)):
            compare(self.packet, candidate)

    def test_zero_difference_is_undefined(self):
        r = compare(self.packet, self.packet)
        self.assertIsNone(r['accepted'])
        self.assertFalse(r['production_comparison'])
        self.assertTrue(all(v['maximum_absolute'] == 0 for v in r['metrics'].values()))

    def test_nonzero_difference_reports_separately(self):
        candidate = copy.deepcopy(self.packet)
        candidate['rows'][0]['value'] += .01
        r = compare(self.packet, candidate)
        self.assertAlmostEqual(r['metrics']['head']['maximum_absolute'], .01)
        self.assertEqual(r['metrics']['flux']['maximum_absolute'], 0)

    def test_units_rejected(self):
        self.bad(lambda p: p['metadata']['quantities']['head'].update(unit='mm'))

    def test_sign_rejected(self):
        self.bad(lambda p: p['metadata']['quantities']['flux'].update(sign='up'))

    def test_pressure_head_rejected(self):
        self.bad(lambda p: p['metadata']['quantities']['head'].update(meaning='pressure_head'))

    def test_effective_saturation_rejected(self):
        self.bad(lambda p: p['metadata']['quantities']['saturation'].update(meaning='effective_saturation'))

    def test_total_stress_rejected(self):
        self.bad(lambda p: p['metadata']['quantities']['effective_stress'].update(meaning='total'))

    def test_integrated_flow_rejected(self):
        self.bad(lambda p: p['metadata']['quantities']['flux'].update(meaning='integrated'))

    def test_missing_id_time_rejected(self):
        self.bad(lambda p: p['rows'].pop())

    def test_duplicate_id_time_rejected(self):
        self.bad(lambda p: p['rows'].append(p['rows'][0].copy()))

    def test_unknown_id_rejected(self):
        self.bad(lambda p: p['rows'][0].update(id='other'))

    def test_time_shift_rejected(self):
        self.bad(lambda p: p['rows'][0].update(time_s=.000001))

    def test_progress_clock_rejected(self):
        self.bad(lambda p: p['metadata'].update(clock='progress'))

    def test_left_switch_boundary_trace_rejected(self):
        self.bad(lambda p: p['metadata'].update(event_trace='left_limit'))

    def test_interpolation_rejected(self):
        self.bad(lambda p: p['metadata'].update(interpolation='linear'))

    def test_missing_coordinate_rejected(self):
        self.bad(lambda p: p['rows'][0].pop('coordinate_m'))

    def test_coordinate_shift_rejected(self):
        self.bad(lambda p: p['rows'][0].update(coordinate_m=123.))

    def test_nonfinite_value_rejected(self):
        self.bad(lambda p: p['rows'][0].update(value=float('nan')))

    def test_incomplete_rejected(self):
        self.bad(lambda p: p['metadata'].update(complete=False))

    def test_parameter_change_rejected(self):
        self.bad(lambda p: p['metadata']['parameters'].update(ks=.001))

    def test_empty_quantity_rejected(self):
        self.bad(lambda p: p['metadata']['locations'].update(head=[]))

    def test_json_and_csv_same_contract(self):
        import csv
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'p.csv'
            with path.open('w', newline='', encoding='utf-8') as f:
                w = csv.DictWriter(f, fieldnames=list(self.packet['rows'][0]))
                w.writeheader()
                w.writerows(self.packet['rows'])
            path.with_suffix('.metadata.json').write_text(json.dumps(self.packet['metadata']), encoding='utf-8')
            loaded = read(path)
            self.assertEqual(validate(loaded), validate(self.packet))


if __name__ == '__main__':
    unittest.main(verbosity=2)
