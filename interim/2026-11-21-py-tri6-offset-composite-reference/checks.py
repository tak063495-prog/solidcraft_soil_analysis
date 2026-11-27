"""Explicit independent checks. Filename deliberately outside test_*.py collection."""
import argparse
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest
import numpy as np
import reference as r
import compare as c
from run_reference import packet


class Mechanics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xy, cls.cells = r.mesh(1, 1)
        cls.p = cls.xy[cls.cells[0]]

    def test_direct_thickness_energy_matches_transformed_matrix_both_signs(self):
        q = np.random.default_rng(91).normal(size=36)*1.e-4
        for e in [-.1, 0., .1]:
            K, _ = r.element(self.p, 2.e5, .3, .2, e)
            self.assertAlmostEqual(r.direct_energy(self.p, q, 2.e5, .3, .2, e), .5*q@K@q, delta=1e-15)

    def test_rigid_six_and_no_extra_zero_modes(self):
        for e in [-.1, 0., .1]:
            result = r.rigid_diagnostics(e)
            self.assertEqual(result['nullity'], 6); self.assertEqual(result['extra_zero_modes'], 0)
            self.assertLess(result['rigid_relative_residual'], 1e-13)

    def test_assembled_free_patch_has_only_six_rigid_modes(self):
        xy, cells = r.mesh(2, 2)
        K, _ = r.assemble(xy, cells, 2.e5, .3, .2, .1)
        eig = np.linalg.eigvalsh(K.toarray())
        self.assertEqual(np.sum(abs(eig) < max(abs(eig))*1e-10), 6)
        self.assertGreater(eig[6], 0.)

    def test_axial_membrane_to_bending_closed_form_both_signs(self):
        for e in [-.1, 0., .1]:
            report, *_ = r.strip_run(2, e=e)
            self.assertLess(report['errors']['displacement_max_m'], 1e-12)
            self.assertLess(report['errors']['rotation_max_rad'], 1e-12)

    def test_simple_supported_axial_and_bending_closed_form(self):
        for e in [-.1, .1]:
            report, *_ = r.strip_run(2, e=e, N=3., M=.2, support='simple')
            self.assertLess(report['errors']['displacement_max_m'], 1e-12)
            self.assertLess(report['errors']['rotation_max_rad'], 1e-12)

    def test_pure_bending_and_membrane_are_energy_conjugate(self):
        for N, M in [(0., .2), (3., 0.), (3., .2)]:
            report, *_ = r.strip_run(2, N=N, M=M)
            self.assertLess(report['diagnostics']['energy_work_absolute_kn_m'], 1e-13)

    def test_section_resultants_have_eccentric_moment_transfer(self):
        for e in [-.1, .1]:
            _, xy, cells, d, _ = r.strip_run(2, e=e)
            for cell in cells:
                v = r.section(xy, cell, d, 2.e5, 0., .2, e)
                np.testing.assert_allclose(v[:3], [3., 0., 0.], atol=1e-10)
                np.testing.assert_allclose(v[3:6], [e*3., 0., 0.], atol=1e-10)
                np.testing.assert_allclose(v[6:9], [0., 0., 0.], atol=1e-10)

    def test_strips_balance_all_free_and_constrained_nodes(self):
        report, *_ = r.strip_run(4)
        self.assertLess(report['free_residual_relative'], 1e-10)
        self.assertLess(report['diagnostics']['all_node_force_residual_max_kn'], 1e-10)
        np.testing.assert_allclose(report['diagnostics']['force_sum_kn'], 0., atol=1e-10)
        np.testing.assert_allclose(report['diagnostics']['moment_sum_kn_m'], 0., atol=1e-10)

    def test_winkler_analytic_amplitude_satisfies_two_field_equilibrium(self):
        x, y = .37, .29
        w, bx, by = r.winkler_exact(x, y)
        a, b = np.pi/2, np.pi; k2 = a*a+b*b
        W = w/(np.sin(a*x)*np.sin(b*y))
        R = -bx/(a*np.cos(a*x)*np.sin(b*y))
        _, Dm, S = r.constitutive(2.e5, .3, .2)
        self.assertAlmostEqual(S*(W-R), Dm[0, 0]*k2*R, delta=1e-15)
        self.assertAlmostEqual(S*k2*(W-R)+100*W, 1., delta=1e-13)
        self.assertAlmostEqual(by, -b*R*np.sin(a*x)*np.cos(b*y), delta=1e-15)

    def test_winkler_h_errors_are_measured_separately_from_quadrature(self):
        values = [r.winkler_run(n)[0] for n in [1, 2, 4]]
        self.assertTrue(all(b['w_relative_L2'] < a['w_relative_L2'] for a, b in zip(values, values[1:])))
        for report in values:
            self.assertLess(report['free_residual_relative'], 1e-10)
            np.testing.assert_allclose(report['diagnostics']['force_sum_kn'], 0., atol=1e-10)

    def test_single_linear_truss_tip_solution_and_load_balance(self):
        a = r.strip_run(2, N=0., anchor=50.)[0]
        b = r.strip_run(4, N=0., anchor=50.)[0]
        compliance = 2**3/(3*(2.e5*.2**3/12))+2/(5/6*1.e5*.2)
        self.assertAlmostEqual(a['errors']['w_tip_exact_m'], compliance/(1+50*compliance))
        self.assertLess(b['errors']['w_tip_average_absolute_m'], a['errors']['w_tip_average_absolute_m'])
        np.testing.assert_allclose(b['diagnostics']['force_sum_kn'], 0., atol=1e-10)

    def test_hydrostatic_closed_polynomials_and_cut_wet_polygon(self):
        p = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
        for h in [-.1, 0., .4, 1., 1.3]:
            np.testing.assert_allclose(r.pressure_nodes(p, h)[2::6], r.pressure_closed(h), atol=2e-14)

    def test_water_normal_force_has_no_offset_couple_or_uniform_lumping(self):
        p = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
        f0 = r.pressure_nodes(p, .4, e=0.)
        for e in [-.1, .1]:
            np.testing.assert_array_equal(r.pressure_nodes(p, .4, e=e), f0)
        self.assertGreater(np.ptp(f0[2::6]), .01)

    def test_water_virtual_work_exact_quadratic_field(self):
        p = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
        # delta w=x^2: exact integral over x=0..1-y, y=0..h.
        h = .4
        exact = -9.81/3*(h*(h-1.5*h*h+h**3-h**4/4)-(h*h/2-h**3+.75*h**4-h**5/5))
        self.assertAlmostEqual((p[:, 0]**2)@r.pressure_nodes(p, h)[2::6], exact, delta=1e-14)

    def test_out_of_range_offset_and_nonfinite_rejected(self):
        for e in [.100000001, -.100000001, float('nan')]:
            with self.assertRaises(ValueError): r.element(self.p, 2.e5, .3, .2, e)

    def test_curved_edge_and_inverted_jacobian_rejected(self):
        p = self.p.copy(); p[3, 1] += .01
        with self.assertRaises(ValueError): r.geometry(p)
        with self.assertRaises(ValueError): r.geometry(self.p[[0, 2, 1, 5, 4, 3]])

    def test_mcc_existing_invariants_only_no_coupled_acceptance(self):
        v = r.mcc_invariants(np.eye(3)*100, 150., 1.2)
        self.assertEqual(v['q_kpa'], 0.); self.assertAlmostEqual(v['F_kpa2'], -7200.)
        self.assertIsNone(v['accepted']); self.assertEqual(v['coupled_offset_soil'], '仕様未定義')
        with self.assertRaises(ValueError): r.mcc_invariants(-np.eye(3), 150., 1.2)

    def test_installation_reference_uses_relative_displacement_before_offset(self):
        K, _ = r.element(self.p, 2.e5, .3, .2, .1)
        d0 = np.random.default_rng(2).normal(size=36)*.001
        np.testing.assert_array_equal(r.installed_force(K, d0, d0), np.zeros(36))
        delta = np.random.default_rng(3).normal(size=36)*.0001
        np.testing.assert_allclose(r.installed_force(K, d0+delta, d0), K@delta, atol=1e-13)

    def test_death_saved_force_release_keeps_input_immutable(self):
        old, new, removed = np.array([10., -5., 6.]), np.array([2., -1., 3.]), np.array([3., 1., -2.])
        frozen = removed.copy()
        for x in [0., .25, 1.]:
            load, holding = r.death_load(old, new, removed, x)
            np.testing.assert_array_equal(holding, (1-x)*np.array([5., -5., 5.]))
            np.testing.assert_array_equal(load, new+holding)
        np.testing.assert_array_equal(removed, frozen)
        with self.assertRaises(ValueError): r.death_load(old, new, removed, 1.001)


class Contracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = packet(r.strip_run(2))

    def bad(self, change):
        other = copy.deepcopy(self.base); change(other)
        with self.assertRaises(ValueError): c.compare(self.base, other)

    def test_identical_packet_still_has_undefined_acceptance(self):
        v = c.compare(self.base, self.base)
        self.assertIsNone(v['accepted']); self.assertEqual(v['tolerance'], '未定義')
        self.assertFalse(v['production_export'])

    def test_csv_json_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'csv'; c.write_csv(self.base, p)
            self.assertEqual(c.read(p), self.base)

    def test_wrong_units_refused(self):
        self.bad(lambda x: x['metadata']['units'].update(length='mm'))

    def test_wrong_signed_offset_refused(self):
        self.bad(lambda x: x['metadata'].update(e_m=-.06))

    def test_outside_offset_refused(self):
        self.bad(lambda x: x['metadata'].update(e_m=.1000000001))

    def test_basis_normal_and_reference_mismatch_refused(self):
        for update in [dict(normal=[0., 0., -1.]), dict(coordinate_system='global-assumed'),
                       dict(displacement_surface='neutral'), dict(section_moment_surface='attachment-only'),
                       dict(install_reference='unknown'), dict(offset_definition='opposite')]:
            self.bad(lambda x: x['metadata'].update(update))

    def test_silent_legacy_fallback_refused_even_zero_offset(self):
        for formulation in ['cst_dkt', 'mindlin_tri6', None]:
            self.bad(lambda x: x['metadata'].update(actual_formulation=formulation, e_m=0.))

    def test_unfinished_or_cancelled_refused(self):
        self.bad(lambda x: x['metadata'].update(completed=False))
        self.bad(lambda x: x['metadata'].update(cancelled=True))

    def test_missing_metadata_and_duplicate_missing_points_refused(self):
        self.bad(lambda x: x['metadata'].pop('actual_formulation'))
        self.bad(lambda x: x['nodes'].append(copy.deepcopy(x['nodes'][0])))
        self.bad(lambda x: x['nodes'].pop())

    def test_sample_coordinates_cannot_be_silently_interpolated(self):
        self.bad(lambda x: x['nodes'][0].update(x_m=.0001))
        self.bad(lambda x: x['sections'][0].update(r=.25))

    def test_relabelled_opposite_moment_transfer_refused(self):
        self.bad(lambda x: x['sections'][0].update(Msxx_kn=x['sections'][0]['Mpxx_kn']+.06*x['sections'][0]['Nxx_kn_m']))

    def test_nonfinite_and_boolean_values_refused(self):
        self.bad(lambda x: x['nodes'][0].update(ux_m=float('nan')))
        self.bad(lambda x: x['nodes'][0].update(ux_m=True))

    def test_numeric_mismatch_recorded_without_invented_threshold(self):
        other = copy.deepcopy(self.base); other['nodes'][0]['ux_m'] += .001
        result = c.compare(self.base, other)
        self.assertAlmostEqual(result['tables']['nodes']['fields']['ux_m']['max_absolute'], .001)
        self.assertIsNone(result['accepted'])

    def test_duplicate_json_keys_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'bad.json'; p.write_text('{"nodes":[],"nodes":[]}', encoding='utf-8')
            with self.assertRaises(ValueError): c.read(p)


def main():
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    suite = unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    def flatten(s):
        for t in s:
            if isinstance(t, unittest.TestSuite): yield from flatten(t)
            else: yield t.id()
    ids = list(flatten(suite)); started = time.monotonic(); log = io.StringIO()
    result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    text = log.getvalue(); print(text)
    status = int(not result.wasSuccessful() or bool(result.skipped) or bool(result.expectedFailures) or bool(result.unexpectedSuccesses))
    report = dict(test_ids=ids, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
                  skipped=len(result.skipped), expected_failures=len(result.expectedFailures), unexpected_successes=len(result.unexpectedSuccesses),
                  exit_code=status, elapsed_seconds=time.monotonic()-started)
    (args.output/'checks.log').write_text(text, encoding='utf-8')
    (args.output/'validation.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return status


if __name__ == '__main__':
    sys.exit(main())
