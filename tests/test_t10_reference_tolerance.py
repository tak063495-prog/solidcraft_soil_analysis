"""T10 compares only the approved probes; an unmet gate is a reportable result."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "interim/2026-11-23-impl-t10-reference-tolerance"
sys.path.insert(0, str(HERE))
import t10_compare as study
import t10_history as history


class ReferenceWidths(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binaries = ROOT / "target/release"
        cls.suffix = ".exe" if os.name == "nt" else ""
        examples = ("t10_layered_cantilever", "joint_birth_points")
        if any(
            not (cls.binaries / "examples" / (name + cls.suffix)).exists()
            for name in examples
        ):
            subprocess.run(
                [
                    "cargo",
                    "build",
                    "--locked",
                    "--release",
                    "--example",
                    examples[0],
                    "--example",
                    examples[1],
                ],
                cwd=ROOT,
                check=True,
            )

    def test_signed_relative_comparison_and_boundary(self):
        self.assertTrue(study.relative(-97.0, -100.0, 0.03)["passed"])
        self.assertFalse(study.relative(-96.999, -100.0, 0.03)["passed"])
        self.assertFalse(study.relative(100.0, -100.0, 0.03)["passed"])
        for v in (float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                study.relative(v, -100.0, 0.03)

    def test_only_cylinder_is_gated_and_missing_rows_survive(self):
        rows = study.shells()
        self.assertEqual(len(rows), 18)
        for r in rows:
            if r["case"] == "pinched_cylinder" and r["measured"] is not None:
                self.assertEqual(
                    (r["reference"], r["relative_width"]), (-1.8248e-5, 0.03)
                )
                self.assertEqual(r["passed"], abs(r["signed_relative"]) <= 0.03)
            else:
                self.assertIsNone(r["passed"])
        fine_hemisphere = [
            r for r in rows if r["case"] == "hemisphere_hole" and r["n"] == 8
        ]
        self.assertEqual(len(fine_hemisphere), 2)
        self.assertTrue(
            all(r["measured"] is None and r["diagnostic"] for r in fine_hemisphere)
        )

    def test_cantilever_probes_are_layer_resolved(self):
        r = study.cantilever(
            self.binaries / "examples" / ("t10_layered_cantilever" + self.suffix)
        )
        self.assertEqual(r["triangles"], 200)
        self.assertEqual(r["tip_dz_m"]["relative_width"], 0.005)
        self.assertEqual({v["offset_m"] for v in r["layers"]}, {-0.2, 0.2})
        for layer in r["layers"]:
            for key, ref in (("qx_n_m", -500.0), ("mxx_nm_m", 5000.0)):
                gate = layer[key]
                self.assertEqual(
                    (gate["reference"], gate["relative_width"]), (ref, 0.001)
                )
                self.assertEqual(
                    gate["passed"], abs(gate["measured"] - ref) <= abs(ref) * 0.001
                )

    def test_joint_laws_against_analytic_patch_and_multipliers(self):
        with tempfile.TemporaryDirectory() as tmp:
            laws = study.joints(
                self.binaries / ("soil-fem-3d" + self.suffix), Path(tmp)
            )
        for name in ("rigid", "hinge", "spring"):
            self.assertTrue(laws[name]["passed"])
            self.assertEqual(
                {c["group"] for c in laws[name]["multipliers"]["cuts"].values()},
                {"wing-0", "wing-1"},
            )

    def test_work_from_prescribed_slip_and_effective_pressure(self):
        self.assertTrue(all(r["passed"] for r in study.edf_force_points()))
        patch = study.contact_patch(
            self.binaries / "examples" / ("joint_birth_points" + self.suffix)
        )
        self.assertTrue(patch["open_passed"])
        self.assertEqual(
            [r["state"] for r in patch["points"]], ["stick", "slip", "slip"]
        )
        for r in patch["points"]:
            self.assertTrue(r["force"]["passed"])
            self.assertTrue(r["work_passed"])
            self.assertEqual(r["pressure_kpa"], 1000.0)
            self.assertEqual(r["reference_work_kpa_m"], 300.0 * r["reference_slip_m"])

    def test_full_history_preserves_boundary_labels(self):
        r = history.wet_history()
        self.assertEqual(
            (len(r["stages"]), r["accepted"], r["contact_gp"]), (13, 416, 6)
        )
        self.assertEqual(sum(r["counts"].values()), 2496)
        self.assertEqual(r["counts"]["2/1"], 133)
        self.assertIsNone(r["passed"])
        self.assertIsNone(r["relative_width"])
        bad = [{"row": {"stage": "a", "progress": i / 32}} for i in range(1, 33)]
        history.accepted(bad, ["a"])
        bad[3] = bad[4]
        with self.assertRaises(ValueError):
            history.accepted(bad, ["a"])

    def test_placement_and_both_physical_time_budgets(self):
        r = study.water()
        self.assertTrue(all(p["passed"] for p in r["placements"]))
        self.assertIsNone(r["steady_domain_balance"])
        self.assertTrue(r["excavation"]["passed"])
        self.assertEqual(r["excavation"]["hand_m3"], 0.5)
        for sample in ("fixed_soil", "lifecycle"):
            for name, (step, cumulative, floor) in study.LIMITS.items():
                profile = r[sample][name]
                self.assertTrue(profile["passed"])
                for row in profile["global_rows"]:
                    self.assertEqual(row["step_limit"], max(floor, step * row["D"]))
                    self.assertTrue(
                        max(abs(row["E"]), row["U"])
                        <= max(floor, cumulative * row["sum_D"])
                    )
