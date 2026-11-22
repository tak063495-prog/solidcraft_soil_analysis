import json
from pathlib import Path
import sys
import unittest

import numpy as np

STUDY = Path(__file__).resolve().parents[1] / "interim/2026-11-22-impl-t8-acceptance-criteria"
sys.path.insert(0, str(STUDY))
from branches import paired
from refinement import richardson
from two_mesh import metrics


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


if __name__ == "__main__":
    unittest.main()
