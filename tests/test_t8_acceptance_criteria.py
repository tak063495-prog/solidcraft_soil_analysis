import io
import json
from pathlib import Path
import sys
import unittest
import zipfile

import numpy as np

STUDY = Path(__file__).resolve().parents[1] / "interim/2026-11-22-impl-t8-acceptance-criteria"
sys.path.insert(0, str(STUDY))
from branches import paired
from refinement import richardson
from two_mesh import metrics
from observations import verify


class AcceptanceCriteria(unittest.TestCase):
    def estimate(self, values, widths=(1., .5, .25), **changes):
        conditions = dict(similar=True, complete=True, reference_converged=True,
                          other_axes_fixed=True, smooth=True)
        return richardson(values, widths, **(conditions | changes))

    def test_manufactured_second_order_scalar_two_successive_triples(self):
        for widths in ((1., .5, .25), (.5, .25, .125)):
            with self.subTest(widths=widths):
                result = self.estimate([2 + h*h for h in widths], widths)
                self.assertEqual(result["holds"], [])
                self.assertEqual(result["observed_order"], 2)
                self.assertEqual(result["extrapolated"], 2)
                self.assertAlmostEqual(result["gci_fine"], 1.25 * widths[-1]**2 / (2+widths[-1]**2))

    def test_hold_conditions_keep_extrapolation_empty(self):
        cases = [((1, 1.5, 1.25), {}, "reverse"), ((1, 2, 4), {}, "non-positive"),
                 ((1, 1, 1), {}, "roundoff"), ((1, 2, None), {}, "missing"),
                 ((1, 2), {}, "three levels"), ((1, 2, 2.5), {"similar":False}, "non-similar"),
                 ((1, 2, 2.5), {"widths":(1., .4, .25)}, "unequal"),
                 ((1, 2, 2.5), {"complete":False}, "incomplete"),
                 ((1, 2, 2.5), {"reference_converged":False}, "unresolved"),
                 ((1, 2, 2.5), {"other_axes_fixed":False}, "mixed"),
                 ((1, 2, 2.5), {"smooth":False}, "non-smooth")]
        for q, options, reason in cases:
            with self.subTest(reason=reason):
                result = self.estimate(q, **options)
                self.assertIn(reason, " ".join(result["holds"]))
                self.assertIsNone(result["extrapolated"])
                self.assertIsNone(result["gci_fine"])

    def test_zero_fine_observable_retains_dimensional_correction(self):
        result = self.estimate([3., 1., 0.])
        self.assertAlmostEqual(result["dimensional_correction"], 1)
        self.assertIsNone(result["gci_fine"])

    def test_branch_correspondence_precedes_field_difference(self):
        fields = dict(stage="open", generations=[1,1], event=["stick","slip"])
        point = dict(actual=.02, bracket=[.01,.03], matched_side=True, fields=fields)
        self.assertIsNone(paired(point, point, .02))
        for altered in (dict(matched_side=False), dict(bracket=[.03,.01]), dict(actual=.025),
                        dict(fields=fields|{"stage":"birth"}), dict(fields=fields|{"generations":[1,2]}),
                        dict(fields=fields|{"event":["stick","stick"]})):
            with self.subTest(altered=altered):
                self.assertIsNotNone(paired(point, point | altered, .02))

    def test_integral_cancellation_preserves_point_error_and_location(self):
        x = np.array([[[1.,0.],[-1.,0.]]])
        r = metrics(x, np.zeros_like(x), np.array([.5,.5]), [("open",1)],
                    np.array([[0.,0.,0.],[1.,0.,0.]]))
        self.assertEqual(r["surface_integral"]["absolute_l2"], 0)
        self.assertEqual(r["area_weighted"]["absolute_l2"], 1)
        self.assertEqual(r["point"]["maximum_absolute"], 1)
        self.assertEqual(r["maximum_at"]["accepted"], ("open",1))


    def test_all_original_rows_and_verified_limited_scopes(self):
        self.assertGreater(verify(STUDY), 0)
        result = json.loads((STUDY/"acceptance.json").read_bytes())
        self.assertEqual(set(result["items"]), {f"A-0{i}" for i in range(1,5)} |
                         {f"B-0{i}" for i in range(1,4)} | {f"C-0{i}" for i in range(1,5)} | {"NC-01"})
        self.assertTrue(all(r["original"] == "未達" for r in result["items"].values()))
        self.assertFalse(result["p6_complete"])

    def test_finite_references_and_strict_category_inventory(self):
        with zipfile.ZipFile(STUDY/"observations.zip") as z:
            c1 = json.loads(z.read("c01-results.json"))
            c2 = json.loads(z.read("c02-extended.json"))
            c3 = json.loads(z.read("c03/results.json"))
            c4 = json.loads(z.read("c04-results.json"))
        self.assertLess(c1["directions"]["independent"]["plastic_dissipation"]["observed_order"], 0)
        self.assertIsNone(c1["gci_fine"])
        self.assertEqual(c2["common_points"], 160)
        self.assertEqual(len(c2["time_mesh_exchange"]), 12)
        self.assertIsNone(c2["observed_order"])
        sequence = [c3["sensitivity"][f"m0c{i}-native--m0c{i+1}-native"]["all_fields"]["plastic_slip"]["absolute_l2"]
                    for i in range(3)]
        self.assertGreater(sequence[1],sequence[0])
        self.assertFalse(c3["reference_assessed_as_exact"])
        self.assertEqual(c4["preserved_mismatch_counts"], [36,63,133])
        fine = next(r for r in c4["measured"] if (r["mechanical_steps"],r["water_steps"]) == (32,32))
        self.assertEqual(fine["strict_mismatches"], 133)
        self.assertEqual(len(fine["boundary"]), 133)
        self.assertEqual(sum(fine["confusion"].values()), 2496)
        self.assertEqual(sum(r["zero_saved_plastic_increment_both"] for r in fine["boundary"]), 133)
        self.assertFalse(c4["boundary_reclassification"])

    def test_admissible_projection_and_moving_water_plane(self):
        with zipfile.ZipFile(STUDY/"observations.zip") as z:
            projection = json.loads(z.read("a02/results.json"))
            placements = json.loads(z.read("t8-water/placement.json"))
            water = json.loads(z.read("t8-water/affine.json"))
        self.assertTrue(projection["projection/good"]["accepted"])
        for name in ("tension", "outside-circle", "balanced-unreproducible-tip"):
            self.assertFalse(projection["projection/"+name]["accepted"])
        self.assertEqual([r["level_m"] for r in placements], [.625,.75,.875])
        for row in placements:
            self.assertIsNone(row["whole_domain_budget"])
            for key, value in row["expected"].items():
                self.assertLess(abs(value-row["measured"][key]), 1e-12)
        self.assertEqual(len(water), 6)
        for row in water:
            self.assertTrue(all(v["passed"] for v in row["water"]["profiles"].values()))
            for point in row["traces"]:
                self.assertLess(abs(point["a_flow_error_m3_s"]), 1e-14)
                self.assertLess(abs(point["b_flow_error_m3_s"]), 1e-14)
                self.assertLess(point["max_nodal_force_error_kn"], 1e-10)
                self.assertEqual(point["expected_a_outflow_m3_s"] == 0, row["slope"] == 0)

    def test_shell_references_and_reassembly_failures_remain_visible(self):
        with zipfile.ZipFile(STUDY/"observations.zip") as z:
            standard = json.loads(z.read("b01/standard/measurements.json"))
            joints = json.loads(z.read("b03/controls.json"))
        fine = [r for r in standard if r["n"] == 8]
        self.assertEqual(len(fine), 6)
        self.assertEqual(sum("displacement" not in r for r in fine), 2)
        self.assertTrue(all(not r["passed"] for r in standard))
        for r in joints:
            self.assertEqual(r["modes"]["free"]["nullity"], {"rigid":6,"hinge":7,"spring":6}[r["law"]])
        table = json.loads((STUDY/"comparisons.json").read_bytes())
        self.assertEqual(sum(r["native"]["max_relative"]>1e-10 for r in table["thin_original"].values()), 4)
        self.assertIsNone(table["nc01"]["physical_collapse_load"])
        self.assertEqual(table["nc01"]["physical_classification"], "D")

    def test_status_suffix_under_lf_and_crlf(self):
        ending = "A-01～C-04、NC-01、P6 の元項目は未達。\n"
        for name in ("README.ja.md","REPRODUCE.ja.md"):
            data = (STUDY/name).read_bytes().replace(b"\r\n", b"\n")
            for newline in (b"\n", b"\r\n"):
                with io.TextIOWrapper(io.BytesIO(data.replace(b"\n",newline)), encoding="utf-8") as f:
                    self.assertTrue(f.read().endswith(ending), name)


if __name__ == "__main__":
    unittest.main()
