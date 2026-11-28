"""Light B-03 contracts; sealed matrix solves are never repeated in CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

from diagnostic_source import diagnostic_bytes

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-26-b03-offset-kink-cause-diagnosis'
_path = sys.path[:]
try:
    sys.path.insert(0, str(HERE))
    from b03cause_evidence import verify
    from b03cause_math import cross, document_text, dot, frame, norm, observed_order, shared_guard
    from b03cause_report import prior_assessment
finally:
    sys.path[:] = _path


class B03OffsetKinkCauseDiagnosis(unittest.TestCase):
    def measurements(self):
        return json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))

    def test_fixed_order_arithmetic_keeps_cancelled_component(self):
        self.assertEqual(dot([1e16, 1., -1e16], [1., 1., 1.]), 1.)
        self.assertEqual(norm([[3., 4.], [0., 12.]]), 13.)
        with self.assertRaises(ValueError):
            dot([math.nan], [1.])
        with self.assertRaises(ValueError):
            dot([1.], [1., 2.])

    def test_signed_arm_is_invariant_only_for_simultaneous_reversal(self):
        points = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]
        n = frame(points); reversed_n = frame([points[i] for i in (0, 2, 1)])
        self.assertEqual((.05*n).tolist(), (-.05*reversed_n).tolist())
        self.assertEqual(norm(.05*n-.05*reversed_n), .1)
        self.assertEqual(norm(.05*n-(-.05*n)), .1)

    def test_work_conjugacy_requires_force_and_arm_couple(self):
        import numpy as np
        arm = np.array([0., 0., .05]); force = np.array([2., -3., 1.])
        u = np.array([.1, .2, -.1]); theta = np.array([.03, -.02, .01])
        neutral = u+cross(arm)@theta
        moment = -np.cross(arm, force)
        self.assertAlmostEqual(dot(force, neutral), dot(force, u)+dot(moment, theta), places=15)
        self.assertNotEqual(dot(force, neutral), dot(force, u))

    def test_actual_range_guard_does_not_clip_one_ulp_outside(self):
        rows = self.measurements()['production_contract_observations']['model_guards']
        data = {r['id']:r['production_model_guard']['accepted'] for r in rows}
        self.assertTrue(data['range-0.05'])
        self.assertTrue(data['range--0.05'])
        self.assertFalse(data['range-'+str(math.nextafter(.05, math.inf))])
        self.assertFalse(data['range-'+str(-math.nextafter(.05, math.inf))])
        self.assertFalse(data['range-0.051'])

    def test_actual_order_dependence_is_recorded_without_changing_guard(self):
        rows = self.measurements()['production_contract_observations']['model_guards']
        data = {r['id']:r['production_model_guard']['accepted'] for r in rows}
        self.assertTrue(data['three-order-012'])
        self.assertTrue(data['three-order-210'])
        self.assertFalse(data['three-order-021'])
        self.assertTrue(data['thickness-order-False'])
        self.assertFalse(data['thickness-order-True'])
        tau = 64*sys.float_info.epsilon*.1
        self.assertTrue(shared_guard([[0, 0, 0], [0, 0, .75*tau], [0, 0, 1.5*tau]], [.1]*3)['accepted'])
        self.assertFalse(shared_guard([[0, 0, 0], [0, 0, 1.5*tau], [0, 0, .75*tau]], [.1]*3)['accepted'])

    def test_internal_rigid_modes_do_not_override_production_refusal(self):
        a = self.measurements(); rows = {r['id']:r for r in a['formal_patch_observations']['cases']}
        self.assertEqual(rows['fold90-1/base']['inherited_algebraic_zero_count'], 6)
        self.assertEqual(rows['split90/base']['inherited_algebraic_zero_count'], 12)
        guards = {r['id']:r['production_model_guard']['accepted'] for r in a['production_contract_observations']['model_guards']}
        self.assertFalse(guards['fold-90-positive'])
        self.assertIsNone(rows['fold90-1/base']['original_general_joint_acceptance'])
        self.assertFalse(a['passed'])
        self.assertEqual(a['formal_patch_observations']['global_equilibrium_solves'], 46)

    def test_original_controls_and_headless_preparation_refusals_preserved(self):
        a = self.measurements()['production_contract_observations']
        self.assertTrue(a['all_cli_contracts_confirmed'])
        self.assertTrue(a['predeclared_model_expectations_confirmed'])
        self.assertTrue(a['preparation_and_legacy']['refusal_contract_passed'])
        self.assertTrue(a['preparation_and_legacy']['legacy_controls_completed'])
        self.assertFalse(a['gui_gpu_real_cad_native_display_executed'])
        self.assertEqual(sum(r['expected_diagnostic'] is None for r in a['cli_controls']), 6)

    def test_two_spatial_levels_do_not_create_order_or_reference_acceptance(self):
        a = prior_assessment()
        self.assertEqual(a['original_b03_spatial_levels'], 2)
        self.assertIsNone(a['spatial_observed_order'])
        self.assertIsNone(observed_order([.2, .1]))
        self.assertFalse(a['continuum_reference_certified'])
        self.assertFalse(a['literature_values_adopted_as_acceptance'])
        self.assertFalse(a['b02_old_failure_and_negative_orders_replaced'])

    def test_prior_record_and_production_source_bytes_remain_bound(self):
        binding = json.loads((HERE/'prior-record-binding.json').read_text())
        for name, expected in binding.items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(), expected, name)
        binding = json.loads((HERE/'source-path-diagnosis.json').read_text())
        # Bind the diagnosed code to the base in the B-03 report; current sources evolve.
        commit = '3747630ab615b4be2117ca1aaa13e798cd436374'
        for name, value in binding.items():
            self.assertEqual(hashlib.sha256(diagnostic_bytes(ROOT, commit, name)).hexdigest(), value['sha256'], name)

    def test_document_suffix_normalizes_real_crlf_only_for_text(self):
        suffix = 'B-03は未達、passed=false。\n'
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'README.ja.md'
            p.write_bytes(('記録\r\n'+suffix.replace('\n','\r\n')).encode())
            self.assertTrue(document_text(p).endswith(suffix))
        self.assertTrue(document_text(HERE/'README.ja.md').endswith(suffix))

    def test_sealed_evidence_false_status_and_sparse_checkout_contract(self):
        a = json.loads((HERE/'acceptance.json').read_text())
        self.assertFalse(a['passed'])
        self.assertFalse(a['b03_complete'])
        self.assertFalse(a['production_scope_expanded'])
        self.assertIsNone(a['general_joint_acceptance'])
        result = verify(HERE)
        self.assertTrue(result['all_sha256_verified'])
        self.assertEqual(result['numerical_reruns'], 0)
        workflow = (ROOT/'.github/workflows/fem-ci.yml').read_text()
        self.assertIn("matrix.group == 'python-c01-cause'", workflow)
        self.assertNotIn('!/interim/2026-12-26-b03-offset-kink-cause-diagnosis/', workflow)


if __name__ == '__main__':
    unittest.main()
