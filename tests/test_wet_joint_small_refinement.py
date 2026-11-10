"""Keep refined small JOINT evidence and software success distinct from acceptance."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-10-wet-joint-small-refinement'


def module(name):
    spec = importlib.util.spec_from_file_location('small_joint_'+name, HERE/(name+'.py'))
    value = importlib.util.module_from_spec(spec); spec.loader.exec_module(value)
    return value


study = module('study')
replay = module('replay')


class SmallJointRefinement(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'evidence'
        cls.manifest = replay.extract(cls.folder)
        cls.measured = study.measure(cls.folder/'paths')

    def test_every_refined_accepted_field_passes_existing_history_and_force_audits(self):
        for steps in (16, 32, 64, 128):
            value = self.measured['paths'][f'g{steps}-l8']
            self.assertEqual(value['accepted_states'], 7*steps)
            self.assertTrue(value['bulk_water_consistent'])
            self.assertTrue(value['history_consistent'])
            self.assertLess(value['bulk_water_maximum_absolute_differences']['force'], 1e-10)
        self.assertEqual(sum(p['accepted_states'] for p in self.measured['paths'].values()), 1680)
        for flag in ('passed', 'full_history_convergence_accepted', 'independent_bulk_plastic_solution',
                     'completion_restart_verified', 'p3b_p4_complete', 'structure_water_joint_complete', 'p6_complete'):
            self.assertIs(self.measured[flag], False)

    def test_prior_16_32_measurement_is_preserved_and_no_unsampled_point_is_certified(self):
        saved = study.read(self.folder/'paths/measurements.json')
        prior = study.read(ROOT/'interim/2026-11-10-nc-oblique-and-joint-history/joint-measurements.json')
        self.assertEqual(saved['adjacent_comparisons']['16__32'], prior['comparison'])
        for pair, points in (('16__32', 112), ('32__64', 224), ('64__128', 448)):
            result = self.measured['adjacent_comparisons'][pair]
            self.assertEqual(result['common_accepted_states'], points)
            self.assertFalse(result['full_history_acceptance'])
            self.assertFalse(result['finer_unmatched_states_compared'])
        # The saved finer difference is not monotone: software success must
        # never promote these fields to a convergence certificate.
        field = lambda pair: saved['adjacent_comparisons'][pair]['fields']['plastic_slip']['pooled_relative_l2']
        self.assertGreater(field('64__128'), field('32__64'))

    def test_forged_completion_flag_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            forged = copy.deepcopy(self.measured); forged['passed'] = True
            study.write(folder/'measurements.json', forged)
            with patch.object(study, 'measure', return_value=self.measured):
                with self.assertRaisesRegex(ValueError, 'acceptance flags differ'):
                    study.replay(folder)

    def test_missing_interior_accepted_chunk_is_rejected(self):
        folder = self.folder/'paths/g128-l8/accepted'
        manifest = study.read(folder/'manifest.json')
        path = folder/manifest['chunks'][1]['file']; raw = path.read_bytes()
        try:
            path.unlink()
            with self.assertRaises(FileNotFoundError):
                study.measure(self.folder/'paths')
        finally:
            path.write_bytes(raw)


if __name__ == '__main__':
    unittest.main()
