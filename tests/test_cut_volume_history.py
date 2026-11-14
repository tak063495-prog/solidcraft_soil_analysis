"""Full-path/field corruption must not masquerade as cut-volume convergence."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import diagnose_cut_volume_history as study
from wet_joint_cycle_reference import cut_volume_reference, differences


class CutVolumeHistory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = study.inputs()
        cls.coarse = cut_volume_reference(*cls.inputs, 4, 8, 7)
        cls.fine = cut_volume_reference(*cls.inputs, 8, 8, 7)

    def test_existing_full_path_and_all_fields_are_checked(self):
        rows = study.validate_path(self.coarse, *self.inputs, 4, 8)
        self.assertEqual(len(rows), 7)
        self.assertEqual(rows[2]['open_points'], 32)
        measured = study.summary(differences(self.coarse, self.fine))
        self.assertEqual(set(measured), set(study.FIELDS))
        self.assertGreater(measured['plastic_slip']['maximum_absolute']['value'], 0.)

    def test_changed_schedule_missing_fields_and_nonfinite_history_are_rejected(self):
        for kind in ('schedule', 'missing', 'nonfinite'):
            bad = copy.deepcopy(self.coarse)
            if kind == 'schedule':
                bad['path'][0]['progress'] += .01
            elif kind == 'missing':
                bad['frames'].pop()
            else:
                bad['frames'][2]['plastic_slip'][0, 0] = np.nan
            with self.assertRaises(AssertionError):
                study.validate_path(bad, *self.inputs, 4, 8)

    def test_corrupt_contact_energy_history_and_changed_physical_scope_are_rejected(self):
        bad = copy.deepcopy(self.coarse); bad['frames'][1]['endpoint_work'][0] += .01
        with self.assertRaises(AssertionError):
            study.validate_path(bad, *self.inputs, 4, 8)
        mesh, ground, case = copy.deepcopy(self.inputs); ground['initialization']['water_level_m'] = 1.
        with self.assertRaisesRegex(AssertionError, 'existing wet JOINT inputs only'):
            study.validate_path(self.coarse, mesh, ground, case, 4, 8)

    def test_old_256_path_is_checked_and_read_without_new_solve(self):
        value, provenance = study.recorded_256()
        self.assertFalse(provenance['rerun'])
        self.assertEqual(value['steps'], 256)
        self.assertEqual(len(provenance['entries']), 4)
        self.assertEqual(value['volume_rule'], 'affine_cut_hydrostatic')


if __name__ == '__main__':
    unittest.main()
