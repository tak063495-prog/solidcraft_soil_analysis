"""CLI contract for the h-refinement report. Numerical gates live in the Rust test."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "target" / "release" / ("soil-fem-3d.exe" if os.name == "nt" else "soil-fem-3d")


class SpatialConvergenceCli(unittest.TestCase):
    def test_ci_profile_keeps_p6_open_and_separates_uniform_srm_from_general_convergence(self):
        self.assertTrue(BINARY.is_file(), f"release binary missing: {BINARY}")
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "report.json"
            completed = subprocess.run(
                [str(BINARY), "spatial-convergence", "ci", str(output)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(output.read_text(encoding="utf-8"))
            refused = subprocess.run(
                [str(BINARY), "spatial-convergence", "ci", str(output)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(refused.returncode, 0)
        self.assertEqual(report["schema"], "soil-fem-spatial-convergence/1")
        self.assertFalse(report["p6_complete"])
        self.assertFalse(report["general_3d_spatial_convergence_demonstrated"])
        self.assertTrue(report["smooth_linear_elasticity_orders_demonstrated"])
        self.assertTrue(report["uniform_srm_factor_within_search_spacing"])
        self.assertAlmostEqual(report["analytical_uniform_srm_factor"], 4 / 3)
        self.assertLess(report["manufactured_balance_error"], 1e-8)
        self.assertIn("5%", " ".join(report["limitations"]))
        gravity = next(series for series in report["elasticity"] if series["id"] == "gravity-tet4")
        sine = next(series for series in report["elasticity"] if series["id"] == "sine-tet10")
        self.assertGreater(gravity["displacement_orders"][-1], 1.6)
        self.assertLess(gravity["displacement_orders"][-1], 2.4)
        self.assertGreater(sine["displacement_orders"][-1], 2.5)
        self.assertLess(abs(gravity["monitor"]["estimate_minus_reference"]), 1e-6)
        midpoints = [level["midpoint"] for level in report["srm"]]
        self.assertTrue(all(level["analytical_factor_inside_bracket"] for level in report["srm"]))
        self.assertTrue(all(abs(value - midpoints[0]) <= 1e-12 for value in midpoints))
        self.assertLess(abs(midpoints[0] - 4 / 3), report["srm_search_spacing"])
