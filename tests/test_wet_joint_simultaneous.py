"""Keep pressure alignment and simultaneous path measurements fail-closed."""
import copy
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import diagnose_wet_joint_simultaneous as study
from diagnose_cut_volume_history import validate_path
from wet_joint_cycle_reference import cut_volume_reference


class SimultaneousWetJoint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = study.inputs()
        cls.names = ['separated/g256-m5e-06/'+n+'.json' for n in ('mesh', 'ground', 'stages', 'stage-results')]
        cls.old, cls.provenance = study.recorded(study.NATIVE, cls.names)
        cls.native_run = cls.old[cls.names[-1]]
        cls.alignment = study.force_alignment(*cls.inputs, cls.native_run)
        cls.coarse = cut_volume_reference(*cls.inputs, 4, 8, 7)
        cls.fine = cut_volume_reference(*cls.inputs, 8, 16, 7)

    def test_pressure_trace_can_match_while_volume_equilibrium_is_rejected(self):
        self.assertFalse(self.alignment['exact_volume_equilibrium_accepted'])
        self.assertFalse(self.alignment['native_result_modified'])
        self.assertFalse(self.alignment['independent_plastic_solution'])
        self.assertEqual(len(self.alignment['frames']), 7)
        for row in self.alignment['frames']:
            self.assertLess(row['sampled_relative_force_residual'], 1.06e-8)
            self.assertLess(row['physical_pressure']['maximum_absolute'], 1e-10)
            self.assertLess(row['interface_load']['maximum_absolute'], 1e-9)
        self.assertGreater(self.alignment['volume_force_difference']['maximum_absolute'], 1e-6)

    def test_changed_pressure_load_history_and_missing_frame_are_rejected(self):
        for kind in ('pressure', 'load', 'missing', 'history'):
            bad = copy.deepcopy(self.native_run)
            if kind == 'pressure':
                bad['frames'][0]['state']['joint_water']['physical_pressure_kpa'][0] += .1
            elif kind == 'load':
                bad['frames'][0]['state']['joint_water']['nodal_load_kn'][0] += .1
            elif kind == 'missing':
                bad['frames'].pop()
            else:
                bad['frames'][0]['state']['joints']['points'][0]['endpoint_work_kpa_m'] += .1
            with self.assertRaises(AssertionError):
                study.force_alignment(*self.inputs, bad)

    def test_every_axis_keeps_all_fields_and_does_not_invent_history_gate(self):
        paths = {key: self.coarse if key.startswith('g512') else self.fine
                 for pair in study.PAIRS for key in pair}
        rows = study.compare_paths(paths)
        self.assertEqual(len(rows), 5)
        for value in rows.values():
            self.assertEqual(len(value['fields']), 7)
            self.assertEqual(len(value['maxima']), 18)
        self.assertGreater(rows['g512-l8__g1024-l32']['maxima']['plastic_slip']['maximum_absolute']['value'], 0.)

    def test_incomplete_or_nonfinite_new_path_does_not_enter_comparisons(self):
        for kind in ('path', 'residual', 'field'):
            bad = copy.deepcopy(self.fine)
            if kind == 'path':
                bad['path'].pop()
            elif kind == 'residual':
                bad['path'][0]['relative_force_residual'] = 1e-8
            else:
                bad['frames'][0]['accumulated'][0] = np.nan
            with self.assertRaises(AssertionError):
                validate_path(bad, *self.inputs, 8, 16)

    def test_changed_archive_is_rejected_without_solving(self):
        self.assertFalse(self.provenance['rerun'])
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            (folder/'archive-sha256.json').write_text('{"archive_sha256":"wrong"}')
            (folder/'evidence.zip').write_bytes(b'changed')
            with patch.object(study, 'cut_volume_reference', side_effect=RuntimeError('must not solve')):
                with self.assertRaisesRegex(AssertionError, 'archive changed'):
                    study.recorded(folder, ['anything'])


if __name__ == '__main__':
    unittest.main()
