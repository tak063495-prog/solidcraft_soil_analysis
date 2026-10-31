"""CLI contract for the plane-strain comparison with the original 2D column."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "target" / "release" / ("soil-fem-3d.exe" if os.name == "nt" else "soil-fem-3d")


class PlaneStrainOriginal2d(unittest.TestCase):
    def test_ci_profile_keeps_the_slope_failure_and_does_not_close_p6(self):
        self.assertTrue(BINARY.is_file(), f"release binary missing: {BINARY}")
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "report.json"
            completed = subprocess.run(
                [str(BINARY), "plane-strain-2d", "ci", str(output)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(output.read_text(encoding="utf-8"))
            refused = subprocess.run(
                [str(BINARY), "plane-strain-2d", "ci", str(output)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(refused.returncode, 0)
            series = {
                "schema": "soil-fem-plane-strain-2d-series/1",
                "quantity": "archived slope is not an in-plane h series",
                "reference": -1.0,
                "nested_h_refinement": False,
                "values": [-0.8441800129308275, -0.8411430807048378],
            }
            spec = Path(folder) / "series.json"
            judged = Path(folder) / "judged.json"
            spec.write_text(json.dumps(series), encoding="utf-8")
            decision = subprocess.run(
                [str(BINARY), "plane-strain-2d", "judge", str(spec), str(judged)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(decision.returncode, 0, decision.stderr)
            judged_report = json.loads(judged.read_text(encoding="utf-8"))
        self.assertEqual(report["schema"], "soil-fem-plane-strain-2d/1")
        self.assertFalse(report["p6_complete"])
        self.assertFalse(report["griffiths_slope_five_percent_closed"])
        self.assertEqual(report["displacement_sample_relative_l2_limit"], 0.05)
        self.assertEqual(report["column_finest_mesh_five_percent"], "fail")
        self.assertEqual(report["column_continuum_five_percent_after_separation"], "pass")
        self.assertTrue(all(row["five_percent_verdict"] == "fail" for row in report["archived_slope"]))
        self.assertTrue(all(not row["in_plane_h_series"] for row in report["archived_slope"]))
        tet4 = next(item for item in report["series"] if item["element"] == "tet4")
        tet10 = next(item for item in report["series"] if item["element"] == "tet10")
        self.assertEqual(tet4["finest_mesh_five_percent"], "fail")
        self.assertEqual(tet4["continuum_five_percent_after_separation"], "pass")
        self.assertEqual(tet10["finest_mesh_five_percent"], "pass")
        self.assertGreater(
            tet4["levels"][-1]["stages"][0]["displacement_sample_relative_l2"],
            0.05,
        )
        self.assertLess(
            tet4["field_relative_l2"][0]["continuum_relative_error"],
            0.05,
        )
        self.assertFalse(judged_report["p6_complete"])
        self.assertFalse(judged_report["nested_h_refinement"])
        self.assertEqual(judged_report["judgment"]["continuum_verdict"], "indeterminate")
        self.assertTrue(all(verdict == "fail" for verdict in judged_report["judgment"]["discrete_verdicts"]))
