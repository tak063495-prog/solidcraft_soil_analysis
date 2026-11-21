"""Keep sampled historical evidence distinct from new affine cut assembly."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from diagnose_cut_volume_history import inputs
from diagnose_wet_joint_cut_volume import components
import fem_reference as fem
from hydrostatic_volume_audit import correction, enabled


class HydrostaticVolumeAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, cls.ground, _ = inputs()
        cls.q = [fem.volume(cls.mesh, i) for i in range(12)]

    def test_independent_cut_minus_sampled_correction_keeps_original_gp_layout(self):
        q_before = [g['position'].copy() for g in self.q]
        sampled = components(self.mesh, self.ground, 'gj', 3)
        exact = components(self.mesh, self.ground, 'cut', 7)
        actual = correction(self.mesh, [True]*12, self.q, .75, 10.)
        np.testing.assert_allclose(actual.ravel(), exact['pore_force']-sampled['pore_force'], atol=1e-14, rtol=1e-12)
        self.assertGreater(abs(actual).max(), .003)
        self.assertEqual(sum(len(q) for q in q_before), 324)
        for a, b in zip(q_before, self.q):
            np.testing.assert_array_equal(a, b['position'])

    def test_removed_and_structure_cells_receive_no_correction(self):
        full = correction(self.mesh, [True]*12, self.q, .75, 10.)
        mask = [i < 6 for i in range(12)]
        # Lower cells are fully wet; the partial upper cells are excluded.
        self.assertFalse(correction(self.mesh, mask, self.q, .75, 10.).any())
        self.assertFalse(correction(self.mesh, [True]*12, self.q, .75, 10., mask).any())
        self.assertTrue(full.any())

    def test_dry_fully_wet_or_unspecified_water_retains_declared_rule(self):
        for level in (None, -1., 0., 1., 2.):
            self.assertFalse(correction(self.mesh, [True]*12, self.q, level, 10.).any())

    def test_curved_tet10_is_not_flattened_by_audit(self):
        import copy
        mesh = copy.deepcopy(self.mesh)
        node = mesh['tetrahedra10'][6][4]
        mesh['nodes'][node][2] += .001
        active = [i == 6 for i in range(12)]
        self.assertFalse(correction(mesh, active, self.q, .75, 10.).any())

    def test_old_missing_marker_stays_sampled_and_unknown_marker_is_rejected(self):
        self.assertFalse(enabled({}))
        self.assertFalse(enabled({'hydrostatic_volume_integration': None}))
        self.assertTrue(enabled({'hydrostatic_volume_integration': 'affine_cut_partial'}))
        with self.assertRaisesRegex(AssertionError, 'unknown'):
            enabled({'hydrostatic_volume_integration': 'invented'})


if __name__ == '__main__':
    unittest.main()
