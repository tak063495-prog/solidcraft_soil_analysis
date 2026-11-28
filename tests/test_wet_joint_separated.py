"""Independent input paths and rejected/changed histories cannot certify convergence."""
import copy
import inspect
import os
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import diagnose_wet_joint_convergence as study
import probe_nc_midpoint as nc
from wet_joint_cycle_reference import discrete_volume_reference, native_fields, compare


class SeparatedWetHistory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'study'
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.report = study.run(solver, cls.folder, steps=4)
        cls.native = study.read(cls.folder/'g16-m5e-06/stage-results.json')

    def test_separate_axes_and_all_frames_are_measurements(self):
        self.assertFalse(self.report['passed'])
        self.assertTrue(self.report['diagnostic_completed'])
        for key in ('native_global', 'native_local', 'independent_global', 'independent_local'):
            self.assertEqual(len(self.report[key]), 2)
            for rows in self.report[key]:
                self.assertEqual([r['frame'] for r in rows], [f['name'] for f in self.native['frames']])
                for row in rows:
                    for field in ('plastic_slip', 'accumulated', 'plastic_dissipation', 'numerical_dissipation', 'endpoint_work'):
                        self.assertTrue(np.isfinite(list(row[field].values())).all())
        self.assertEqual(len(self.report['volume_discretization']), 7)
        self.assertEqual(len(self.report['independent_volume_order']), 7)
        self.assertTrue(any(r['u']['maximum_absolute'] > 0 for r in self.report['volume_discretization']))

    def test_input_only_discrete_reference_keeps_physical_pressure_and_all_points(self):
        self.assertEqual(list(inspect.signature(discrete_volume_reference).parameters),
                         ['mesh', 'ground', 'case', 'steps', 'contact_steps'])
        a = study.read(self.folder/'g16-l8-discrete.json'); b = study.read(self.folder/'g16-l8-gl7.json')
        self.assertFalse(a['reference_uses_native_results'] or b['reference_uses_native_results'])
        self.assertEqual(a['volume_rule'], 'input_gauss_jacobi')
        for x, y in zip(a['frames'], b['frames']):
            self.assertEqual(np.shape(x['stress']), (324, 6))
            self.assertEqual(np.shape(x['traction']), (32, 3))
            self.assertEqual(x['pore'], y['pore'])
            self.assertEqual(x['physical_pressure'], [2.5]*32)
            self.assertEqual(x['physical_pressure'], y['physical_pressure'])

    def test_changed_history_is_visible_even_when_existing_u_traction_gates_pass(self):
        expected = study.read(self.folder/'g16-l8-discrete.json')
        bad = copy.deepcopy(self.native)
        bad['frames'][2]['state']['joints']['points'][0]['accumulated_plastic_slip_m'] += .1
        # Demonstrate why the old two gates alone cannot certify histories.
        rows = compare(bad, expected)
        clean = study.all_fields(native_fields(self.native), expected)
        self.assertGreater(rows[2]['accumulated']['maximum_absolute'], clean[2]['accumulated']['maximum_absolute']+.05)
        self.assertFalse(self.report['passed'])
        bad['frames'][2]['state']['joints']['points'][0]['plastic_slip_m'][0] = float('nan')
        with self.assertRaises(AssertionError):
            study.all_fields(native_fields(bad), expected)

    def test_incomplete_or_changed_global_path_cannot_enter_sensitivity(self):
        bad = copy.deepcopy(self.native); bad['complete'] = False
        with self.assertRaisesRegex(AssertionError, 'incomplete'):
            study.schedule(bad, 16)
        bad = copy.deepcopy(self.native); bad['attempts'][0]['progress_to'] += .001
        with self.assertRaisesRegex(AssertionError, 'global schedule changed'):
            study.schedule(bad, 16)
        bad = copy.deepcopy(self.native); bad['attempts'][0]['linear_solves'][0]['true_relative_residual'] = 1.
        with self.assertRaises(AssertionError):
            study.schedule(bad, 16)

    def test_refined_every_phase_restart_and_rejected_same_path(self):
        self.assertEqual(len(self.report['restart_checks']), 7)
        for row in self.report['restart_checks']:
            self.assertTrue(row['whole_run_exact_except_elapsed_seconds'])
            self.assertTrue(0 < row['progress'] < 1)
        self.assertTrue(self.report['rejected_history_exact_same_accepted_path'])
        self.assertEqual(sum(r['rejected'] for r in self.report['rejection_schedule']), 3)
        self.assertTrue(self.report['actual_open_slip']['all_rejections_preserve_committed_state'])

    def test_nc_experiment_changes_only_existing_integration_control(self):
        backward = nc.inputs('backward_euler'); midpoint = nc.inputs('midpoint')
        self.assertEqual(midpoint[2]['settings'].pop('mcc_integration'), 'midpoint')
        self.assertEqual(backward, midpoint)
        original = study.read(ROOT/'examples/stages-joint-water-cycle.json')
        saved = copy.deepcopy(original); refined = study.configuration(original, 64, 5e-6, True)
        self.assertEqual(original, saved)
        self.assertEqual(refined['joint_system']['joints'], original['joint_system']['joints'])
        self.assertEqual(refined['boundary'], original['boundary'])


if __name__ == '__main__':
    unittest.main()
