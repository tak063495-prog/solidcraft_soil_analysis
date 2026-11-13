"""Explicit unittest entry point; intentionally outside test_*.py collection."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import compare as cmp
import reference as ref


class Projection(unittest.TestCase):
    def setUp(self):
        self.xy, self.w = ref.quadrature()
        self.target = ref.exact_nodal([2, 1, -10])

    def test_partial_area_exact(self):
        self.assertAlmostEqual(sum(self.w), float(ref.integral({(0, 0): F(1)})), places=14)

    def test_polynomial_integrals_independent_of_GP(self):
        for i in range(5):
            for j in range(5-i):
                self.assertAlmostEqual(float(sum(self.w*self.xy[:, 0]**i*self.xy[:, 1]**j)),
                                       float(ref.integral({(i, j): F(1)})), places=14)

    def test_closed_form_projection(self):
        p = ref.project(self.xy, self.w, self.target)
        self.assertEqual(p["outcome"], "representable_dry")
        self.assertEqual(p["rank"], 5)
        np.testing.assert_allclose(p["traction"], ref.traction(self.xy, [2, 1, -10]), rtol=0, atol=5e-13)

    def test_minimum_weighted_norm_stationarity(self):
        p = ref.project(self.xy, self.w, self.target)
        c = ref.assembly(self.xy, self.w)/np.sqrt(self.w)
        _, _, vt = np.linalg.svd(c, full_matrices=True)
        null = vt[5:].T
        z = np.sqrt(self.w)[:, None]*p["traction"]
        np.testing.assert_allclose(null.T@z, 0, atol=5e-14)
        self.assertGreater(np.linalg.norm(z+null[:, :1]), np.linalg.norm(z))

    def test_all_node_force_moment_virtual_work(self):
        t = ref.project(self.xy, self.w, self.target)["traction"]
        b = ref.balances(self.xy, self.w, t, self.target)
        self.assertLess(b["nodal_max_kN"], ref.FORCE_ATOL)
        self.assertLess(max(abs(x) for x in b["side_a_moment_defect_kNm"]), 1e-13)
        self.assertLess(abs(b["virtual_work_defect_kNm"]), 1e-16)

    def test_random_permutation_and_repeat_exact(self):
        ids = np.random.default_rng(1947).permutation(len(self.w))
        a = ref.project(self.xy, self.w, self.target)["traction"]
        b = ref.project(self.xy[ids], self.w[ids], self.target, ids)["traction"]
        self.assertTrue(np.array_equal(a[ids], b))
        self.assertTrue(np.array_equal(a, ref.project(self.xy, self.w, self.target)["traction"]))

    def test_tension_is_not_clipped_into_admissibility(self):
        p = ref.project(self.xy, self.w, ref.exact_nodal([2, 1, 3]))
        self.assertTrue(p["reproducible"])
        self.assertFalse(p["admissible"])
        self.assertEqual(p["outcome"], "rejected")
        self.assertGreater(np.max(p["traction"][:, 2]), 0)

    def test_outside_friction_is_rejected(self):
        self.assertEqual(ref.project(self.xy, self.w, ref.exact_nodal([8, 6, -10]))["outcome"], "rejected")

    def test_shared_tip_force_unrepresentable(self):
        self.target[0, 2] = 1
        p = ref.project(self.xy, self.w, self.target)
        self.assertFalse(p["reproducible"])
        self.assertEqual(p["outcome"], "rejected")
        self.assertEqual(p["residual"][0, 2], -1)

    def test_same_side_unbalanced_force_unrepresentable(self):
        self.target[1, 0] += .1
        self.assertEqual(ref.project(self.xy, self.w, self.target)["outcome"], "rejected")

    def test_bad_weight_refused(self):
        self.w[0] = 0
        with self.assertRaises(ValueError):
            ref.project(self.xy, self.w, self.target)

    def test_nan_refused(self):
        self.target[1, 0] = float("nan")
        with self.assertRaises(ValueError):
            ref.project(self.xy, self.w, self.target)

    def test_duplicate_GP_refused(self):
        with self.assertRaises(ValueError):
            ref.project(self.xy, self.w, self.target, np.zeros(len(self.w), int))

    def test_clip_diagnostic_preserves_failure(self):
        r = ref.projection_cases()
        for name in ["tension", "outside"]:
            c = r[name]["corrections"]["clip"]
            self.assertGreater(c["nodal_max_kN"], .5)
            self.assertGreater(abs(c["virtual_work_defect_kNm"]), 1e-4)
            self.assertGreater(abs(c["fixed_jump_work_change_kNm"]), 1e-4)
            self.assertIsNone(r[name]["accepted"])

    def test_smoothing_cannot_hide_behind_global_resultant(self):
        c = ref.projection_cases()["inside"]["corrections"]["smooth_then_clip"]
        self.assertLess(max(abs(x) for x in c["action_reaction_kN"]), 1e-13)
        self.assertGreater(c["nodal_max_kN"], .02)
        self.assertGreater(max(abs(x) for x in c["side_a_moment_defect_kNm"]), .01)


class NumericalAndBirth(unittest.TestCase):
    def test_separated_series_and_refinement(self):
        r = ref.numerical_series()
        self.assertEqual([x["surface_elements"] for x in r["gp_series"]], [1]*4)
        self.assertEqual([x["gp_per_element"] for x in r["element_series"]], [8]*4)
        for name in ["gp_series", "element_series"]:
            self.assertGreater(r[name][-1]["observed_order"], 1.9)  # test of manufactured smooth problem only
            self.assertLess(r[name][-1]["observed_order"], 2.1)
        self.assertIsNone(r["accepted"])

    def test_independent_penalty_weak_equilibrium_and_work(self):
        for r in ref.numerical_series()["element_series"]:
            self.assertGreater(r["minimum_pressure_kpa"], 0)
            self.assertLess(r["weak_force_residual_kN"], 1e-13)
            self.assertLess(r["work_balance_defect_kNm"], 1e-16)

    def test_joint_only_death_keeps_all_owners_and_balances(self):
        r = ref.birth_ledger()["joint_only"]
        self.assertEqual(r["owner_active"], [True, True])
        self.assertEqual(r["death_residual_kN"], [0.]*4)
        self.assertEqual(r["rebirth_B_preparation_residual_kN"], [0.]*4)
        self.assertEqual((r["frozen_generation"], r["new_generation"]), (1, 2))
        self.assertEqual(r["new_slip_m"], 0)
        self.assertEqual(r["new_work_kNm"], 0)

    def test_body_birth_changes_full_node_work_and_is_undefined(self):
        r = ref.birth_ledger()["crossing_body_birth"]
        self.assertEqual(r["joint_only_retained_load_reuse_defect_kN"], [0., 0., -9., 11.])
        self.assertAlmostEqual(r["virtual_work_defect_kNm"], -.015)
        self.assertEqual(r["generation_policy"], "仕様未定義")
        self.assertIsNone(r["accepted"])


class Comparison(unittest.TestCase):
    def setUp(self):
        self.data = cmp.fixture()

    def reject(self):
        r = cmp.compare(self.data)
        self.assertFalse(r["contract_valid"])
        self.assertEqual(r["exit_code"], 1)
        self.assertIsNone(r["accepted"])

    def test_clean_is_not_acceptance(self):
        r = cmp.compare(self.data)
        self.assertTrue(r["contract_valid"], r["errors"])
        self.assertEqual(r["exit_code"], 2)
        self.assertIsNone(r["accepted"])
        self.assertEqual(r["source_kind"], "synthetic_reference")

    def test_json_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            for suffix in [".json", ".csv"]:
                path = Path(tmp)/("sample"+suffix)
                cmp.save(self.data, path)
                self.assertTrue(cmp.compare(cmp.load(path))["contract_valid"])

    def test_row_order_irrelevant(self):
        self.data["rows"].reverse()
        self.assertTrue(cmp.compare(self.data)["contract_valid"])

    def test_traction_clip_detected_even_when_node_force_forged(self):
        for r in self.data["rows"]:
            if r["kind"] == "contact":
                r["tn"] = max(r["tn"], -5.)
        self.reject()

    def test_smoothing_detected(self):
        self.data["rows"][0]["tx"] *= .9
        self.reject()

    def test_generation_missing_refused(self):
        del self.data["rows"][0]["generation"]
        with self.assertRaises(ValueError):
            cmp.compare(self.data)

    def test_generation_regression_and_skipping_refused(self):
        for value in [0, 1, 3]:
            with self.subTest(value=value):
                data = copy.deepcopy(self.data)
                next(r for r in data["rows"] if r["frame"] == 2)["generation"] = value
                self.assertFalse(cmp.compare(data)["contract_valid"])

    def test_nonreproducible_projection_detected(self):
        next(r for r in self.data["rows"] if r["run"] == "reversed")["tn"] += 1e-8
        self.reject()

    def test_out_of_scope_silent_acceptance_refused(self):
        for scope in cmp.REFUSALS:
            with self.subTest(scope=scope):
                data = copy.deepcopy(self.data)
                next(r for r in data["rows"] if r["status"] == scope+":rejected")["status"] = "accepted"
                self.assertFalse(cmp.compare(data)["contract_valid"])

    def test_missing_or_duplicate_physical_ID_refused(self):
        missing = copy.deepcopy(self.data)
        missing["rows"].pop()
        self.assertFalse(cmp.compare(missing)["contract_valid"])
        self.data["rows"].append(self.data["rows"][0].copy())
        with self.assertRaises(ValueError):
            cmp.compare(self.data)

    def test_slip_work_reference_and_old_frozen_force_tamper_refused(self):
        for field, frame in [("slip_x", 2), ("work", 2), ("reference_z", 2), ("tn", 1)]:
            with self.subTest(field=field):
                data = copy.deepcopy(self.data)
                next(r for r in data["rows"] if r["frame"] == frame)[field] += .001
                self.assertFalse(cmp.compare(data)["contract_valid"])
        data = copy.deepcopy(self.data)
        next(r for r in data["rows"] if r["frame"] == 2)["reference_z"] = 0.
        self.assertFalse(cmp.compare(data)["contract_valid"])  # keep the nonzero 1e-12 m reference

    def test_global_force_cancellation_does_not_replace_all_node_check(self):
        for r in self.data["rows"]:
            if r["kind"] == "node" and r["frame"] == 0:
                r["fz"] = 0
        self.reject()
        data = cmp.fixture()
        # A simultaneous internal/external edit preserves residual but changes the reference.
        for kind in ["ledger_internal", "ledger_external"]:
            next(r for r in data["rows"] if r["kind"] == kind)["fx"] += 1.
        self.assertFalse(cmp.compare(data)["contract_valid"])

    def test_nan_inf_bool_and_fractional_generation_refused(self):
        for field, value in [("tx", float("nan")), ("tn", float("inf")), ("active", True), ("generation", 1.5)]:
            with self.subTest(field=field):
                data = copy.deepcopy(self.data)
                data["rows"][0][field] = value
                with self.assertRaises(ValueError):
                    cmp.compare(data)

    def test_units_water_completion_input_and_acceptance_refused(self):
        for value in [None, [], "invalid"]:
            data = copy.deepcopy(self.data)
            data["metadata"] = value
            with self.assertRaises(ValueError):
                cmp.compare(data)
        for key, value in [("units", "N-mm"), ("water", True), ("complete", False),
                           ("input_sha256", "wrong"), ("accepted", True)]:
            with self.subTest(key=key):
                data = copy.deepcopy(self.data)
                data["metadata"][key] = value
                self.assertFalse(cmp.compare(data)["contract_valid"])

    def test_duplicate_JSON_key_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"bad.json"
            path.write_text('{"metadata":{},"metadata":{},"rows":[]}', encoding="utf-8")
            with self.assertRaises(ValueError):
                cmp.load(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
