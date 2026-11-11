"""Explicit unittest entry; filename deliberately excluded from test_*.py discovery."""
import ast
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.special import erf

import compare
from diagnose import study
from finite_volume import backward_euler, faces, residual
from reference import (ALPHA, CASES, DELTA_THETA, KS, RA, RB, THETA_R,
                       fields, primary, richards_modes, richards_u, storage)


class AnalyticalChecks(unittest.TestCase):
    def test_terzaghi_early_time_independent_erf_limit(self):
        # Image from the far impermeable wall is exponentially small here.
        h, _ = fields("terzaghi", .2, .02)
        self.assertAlmostEqual(float(h), float(erf(.2/(2*np.sqrt(.02)))), places=12)

    def test_terzaghi_drain_and_impermeable_end(self):
        for t in (.02, .1, .2):
            h, q = fields("terzaghi", np.array([0., 1.]), t)
            self.assertEqual(h[0], 0.)
            self.assertLess(abs(q[1]), 1e-14)
            self.assertTrue(0 < h[1] < 1)

    def test_moving_level_boundaries(self):
        for t in (0., .1, .2):
            h, _ = fields("moving_level", np.array([0., 1.]), t)
            np.testing.assert_allclose(h, [1+.1*np.sin(2*np.pi*t), 1.], atol=1e-14)

    def test_richards_eigenvalues_and_projection(self):
        k, b = richards_modes()
        self.assertLess(np.max(abs(k*np.cos(k)+.5*np.sin(k))), 2e-11)
        # Srivastava-Yeh (10),(11), dimensionless L=1: residue denominator.
        residue_form = 4*(RA-RB)*np.exp(.5)*np.sin(k)/(1+.5+2*k*k)
        np.testing.assert_allclose(b, residue_form, rtol=2e-12, atol=1e-14)
        for kk, bb in zip(k[:6], b[:6]):
            projection = quad(lambda z: 2*(RA-RB)*np.sinh(z/2)*np.sin(kk*z), 0, 1,
                              epsabs=1e-13)[0]
            norm = quad(lambda z: np.sin(kk*z)**2, 0, 1, epsabs=1e-13)[0]
            self.assertAlmostEqual(bb, projection/norm, places=13)

    def test_richards_water_table_surface_flux_and_wetting(self):
        x = np.linspace(0, 1, 31)
        before = richards_u(x, .5)
        after = richards_u(x, 5.)
        self.assertTrue(np.all(after >= before-1e-14))
        self.assertTrue(np.all((after > 0) & (after <= 1+1e-14)))
        self.assertAlmostEqual(float(fields("richards", 0., 5.)[0]), 0., places=14)
        self.assertAlmostEqual(float(fields("richards", 1., 5.)[1]), -KS*RB, places=14)
        np.testing.assert_allclose(richards_u(x, 0), RA+(1-RA)*np.exp(-x), atol=1e-14)

    def test_series_truncation_control(self):
        x = np.linspace(0, 1, 33)
        for name in ("terzaghi", "richards"):
            for a, b in zip(fields(name, x, CASES[name].start, 64),
                            fields(name, x, CASES[name].start, 128)):
                self.assertLess(np.max(abs(a-b)), 1e-13)

    def test_continuum_pde_without_manufactured_source(self):
        dx, dt = 1e-5, 1e-5
        for name, c in CASES.items():
            t = (c.start+c.end)/2
            x = np.array([.2, .5, .8])
            splus = storage(name, primary(name, x, t+dt))[0]
            sminus = storage(name, primary(name, x, t-dt))[0]
            qplus = fields(name, x+dx, t)[1]
            qminus = fields(name, x-dx, t)[1]
            np.testing.assert_allclose((splus-sminus)/(2*dt)+(qplus-qminus)/(2*dx),
                                       0., atol=2e-8, err_msg=name)

    def test_continuum_integrated_mass_independent_quadrature(self):
        for name, c in CASES.items():
            delta = quad(lambda x: float(storage(name, primary(name, x, c.end))[0]
                                          -storage(name, primary(name, x, c.start))[0]), 0, 1,
                         epsabs=1e-12)[0]
            net = quad(lambda t: float(fields(name, 0., t)[1]-fields(name, 1., t)[1]),
                       c.start, c.end, epsabs=1e-12)[0]
            self.assertLess(abs(delta-net), 1e-11, name)

    def test_invalid_domains_are_rejected(self):
        for x, t in ((-1, 1), (1.1, 1), (.5, -1), (np.nan, 1), (.5, np.inf)):
            with self.assertRaises(ValueError):
                fields("richards", x, t)
        with self.assertRaises(ValueError):
            fields("terzaghi", 0, 0)


class NumericalChecks(unittest.TestCase):
    def test_nonlinear_jacobian_against_central_difference(self):
        h = np.array([-.1, -.2, -.3, -.4, -.5])
        old = storage("richards", h-.01)[0]
        _, band, _ = residual("richards", h, old, 1., .1)
        jac = np.diag(band[1])+np.diag(band[0, 1:], 1)+np.diag(band[2, :-1], -1)
        numerical = np.empty((5, 5))
        for j in range(5):
            d = np.zeros(5); d[j] = 1e-6
            numerical[:, j] = (residual("richards", h+d, old, 1., .1)[0]
                                -residual("richards", h-d, old, 1., .1)[0])/2e-6
        np.testing.assert_allclose(jac, numerical, atol=1e-10, rtol=1e-8)

    def test_backward_euler_mass_is_closed_for_each_case(self):
        for name, c in CASES.items():
            _, r = backward_euler(name, 20, c.start, c.end, 20)
            self.assertTrue(r["completed"])
            self.assertLess(abs(r["mass_defect"]), 1e-10)
            self.assertLess(r["mass_relative"], 1e-8)

    def test_richards_whole_flux_step_mass_including_initial_boundary_jump(self):
        n = 20
        initial = np.log(richards_u((np.arange(n)+.5)/n, 0))/ALPHA
        h, r = backward_euler("richards", n, 0., 5., 80, initial)
        self.assertLess(r["mass_relative"], 1e-8)
        self.assertAlmostEqual(r["inflow"], KS*RB*5, places=12)
        self.assertTrue(np.all(h <= 0))

    def test_separated_refinement_orders_and_temporal_contamination(self):
        for name in CASES:
            report = study(name)
            self.assertTrue(all(report["reference_checks"].values()), name)

    def test_hydrostatic_richards_faces_have_zero_flux(self):
        x = (np.arange(20)+.5)/20
        q = faces("richards", -x, 1.)[0]
        # Imposed top influx remains nonzero; all interior and bottom faces vanish.
        np.testing.assert_allclose(q[:-1], 0., atol=1e-15)

    def test_no_application_imports_or_process_execution(self):
        folder = Path(__file__).resolve().parent
        allowed = {"argparse", "csv", "dataclasses", "functools", "json", "pathlib", "platform",
                   "time", "numpy", "scipy", "reference", "finite_volume", "diagnose"}
        for name in ("reference.py", "finite_volume.py", "diagnose.py", "compare.py"):
            tree = ast.parse((folder/name).read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    self.assertTrue(all(n.name.split('.')[0] in allowed for n in node.names), name)
                elif isinstance(node, ast.ImportFrom):
                    self.assertIn(node.module.split('.')[0], allowed, name)


class AdapterChecks(unittest.TestCase):
    def fixture(self, name="terzaghi"):
        with tempfile.TemporaryDirectory() as temp:
            compare.demo(name, Path(temp))
            return compare.read_input(Path(temp)/f"{name}.json")

    def test_exact_json_and_csv_demo_are_identical(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            for name in CASES:
                compare.demo(name, folder)
                a = compare.compare(compare.read_input(folder/f"{name}.json"))
                b = compare.compare(compare.read_input(folder/f"{name}.csv", folder/f"{name}.metadata.json"))
                self.assertEqual(a, b)
                self.assertEqual(a["pde"]["head_m"]["rms"], 0.)
                self.assertIsNone(a["accepted"])
                self.assertEqual(a["A-04"], "未達")

    def test_perturbed_pde_output_is_measured_without_becoming_accepted(self):
        data = self.fixture()
        data["samples"][0]["head"] += .1
        result = compare.compare(data)
        self.assertGreater(result["pde"]["head_m"]["rms"], 0)
        self.assertIsNone(result["pde"]["accepted"])

    def test_units_sign_constitutive_benchmark_and_time_semantics_are_strict(self):
        for key in ("time_unit", "head_kind", "flux_sign", "flux_kind", "benchmark", "coordinate",
                    "complete", "run_status"):
            data = self.fixture()
            data[key] = "wrong"
            with self.assertRaises(ValueError):
                compare.compare(data)

    def test_incomplete_cancelled_or_boolean_aliased_status_cannot_be_read_as_complete(self):
        for value in (False, 1, "true", None):
            data = self.fixture(); data["complete"] = value
            with self.assertRaises(ValueError):
                compare.compare(data)
        data = self.fixture(); data["run_status"] = "cancelled"
        with self.assertRaises(ValueError):
            compare.compare(data)

    def test_nan_inf_duplicate_and_out_of_window_samples_are_rejected(self):
        for key, val in (("head", np.nan), ("flux", np.inf), ("time", -1), ("x", 2), ("head", True)):
            data = self.fixture(); data["samples"][0][key] = val
            with self.assertRaises(ValueError):
                compare.compare(data)
        data = self.fixture(); data["samples"].append(copy.deepcopy(data["samples"][0]))
        with self.assertRaises(ValueError):
            compare.compare(data)

    def test_single_time_and_missing_flux_are_rejected(self):
        data = self.fixture()
        data["samples"] = data["samples"][:5]
        with self.assertRaises(ValueError):
            compare.compare(data)
        data = self.fixture(); del data["samples"][0]["flux"]
        with self.assertRaises(KeyError):
            compare.compare(data)

    def test_mass_missing_is_unevaluated_and_leak_is_reported(self):
        data = self.fixture()
        self.assertIn("未評価", compare.compare(data)["mass_balance"]["status"])
        data["balance"] = [{"start": .02, "end": .2, "inflow_m": 0., "outflow_m": .3,
                            "storage_start_m": 1., "storage_end_m": .8}]
        self.assertAlmostEqual(compare.compare(data)["mass_balance"]["intervals"][0]["defect_m"], .1)
        data["balance"][0]["outflow_m"] = -1
        with self.assertRaises(ValueError):
            compare.compare(data)

    def test_mechanical_path_and_pde_propagation_are_separate(self):
        data = self.fixture()
        row = data["mechanics"][0]
        row["settlement_m"] += .002
        result = compare.compare(data)["mechanics"]["rows"][0]
        self.assertAlmostEqual(result["mechanical_path_error_m"], .002)
        self.assertAlmostEqual(result["pde_propagation_error_m"], 0.)
        self.assertAlmostEqual(result["total_error_m"], sum(result[k] for k in
            ("mechanical_path_error_m", "pde_propagation_error_m", "sampling_quadrature_error_m")))
        row["head"] = [v+.1 for v in row["head"]]
        changed = compare.compare(data)["mechanics"]["rows"][0]
        self.assertNotEqual(changed["pde_propagation_error_m"], 0.)

    def test_wet_joint_mechanics_cannot_use_linear_settlement_reference(self):
        data = self.fixture("richards")
        data["mechanics"] = self.fixture()["mechanics"]
        with self.assertRaises(ValueError):
            compare.compare(data)


if __name__ == "__main__":
    unittest.main()
