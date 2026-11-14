"""Actual opening/sliding evidence must not become a full-history certificate."""
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
import fem_reference as fem
import verify_wet_joint_cycle as driver
from wet_joint_cycle_reference import reference, compare, native_fields


class WetJointCycle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'cycle'
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.report = driver.run(solver, cls.folder, (8, 16))
        cls.mesh = driver.read(cls.folder/'original/mesh.json')
        cls.ground = fem.normalize_ground(driver.read(cls.folder/'original/ground.json'))
        cls.case = fem.normalize_stage_case(driver.read(cls.folder/'original/stages.json'))
        cls.native = driver.read(cls.folder/'original/stage-results.json')
        cls.expected = driver.read(cls.folder/'reference-16.json')

    def test_existing_geometry_input_and_actual_open_slip(self):
        self.assertEqual(self.mesh, driver.read(ROOT/'examples/joint-water-column-tet10-mesh.json'))
        self.assertEqual(driver.read(self.folder/'original/stages.json'), driver.read(ROOT/'examples/stages-joint-water-cycle.json'))
        self.assertEqual(len(self.mesh['nodes']), 54)
        self.assertEqual(len(self.mesh['tetrahedra']), 12)
        rows = self.report['native_routes']['original']['frames']
        self.assertEqual(next(x for x in rows if x['frame'] == 'shear')['status_counts']['slip'], 32)
        self.assertEqual(next(x for x in rows if x['frame'] == 'open')['status_counts']['open'], 32)

    def test_reference_signature_and_independent_schedule(self):
        self.assertEqual(list(inspect.signature(reference).parameters), ['mesh', 'ground', 'case', 'steps', 'contact_steps', 'order'])
        self.assertFalse(self.expected['reference_uses_native_results'])
        self.assertEqual(len(self.expected['path']), 7*16)
        for name in ['geostatic']+[s['name'] for s in self.case['stages']]:
            self.assertEqual([p['progress'] for p in self.expected['path'] if p['stage'] == name],
                             [i/16 for i in range(1, 17)])

    def test_all_fields_all_frames_and_independent_work_identity(self):
        self.assertEqual(len(self.expected['frames']), 7)
        for frame in self.expected['frames']:
            self.assertEqual(np.shape(frame['u']), (54, 3))
            self.assertEqual(np.shape(frame['stress']), (324, 6))
            self.assertEqual(np.shape(frame['traction']), (32, 3))
            self.assertEqual(np.shape(frame['plastic_slip']), (32, 3))
            for key in ('u', 'stress', 'strain', 'pore', 'void', 'jump', 'traction', 'reaction',
                        'plastic_slip', 'accumulated', 'plastic_dissipation', 'numerical_dissipation', 'endpoint_work'):
                self.assertTrue(np.isfinite(frame[key]).all(), key)
            np.testing.assert_allclose(frame['endpoint_work'], np.asarray(frame['elastic_energy'])+
                                       frame['plastic_dissipation']+np.asarray(frame['numerical_dissipation']),
                                       rtol=1e-8, atol=1e-9)
        self.assertGreater(max(p['slipping_points'] for p in self.expected['path']), 0)
        self.assertEqual(max(p['open_points'] for p in self.expected['path']), 32)

    def test_full_history_and_quadrature_differences_are_not_certified(self):
        self.assertFalse(self.report['passed'])
        self.assertTrue(self.report['independent_reference_completed'])
        self.assertFalse(self.report['p3b_p4_rows_complete'] or self.report['p6_complete'] or self.report['structure_water_joint_complete'])
        self.assertFalse(self.report['nc']['executed'] or self.report['nc']['passed'])
        self.assertEqual(len(self.report['independent_step_sensitivity']), 7)
        self.assertEqual(len(self.report['native_step_sensitivity']), 7)
        self.assertEqual(len(self.report['independent_volume_order_sensitivity']), 7)
        self.assertTrue(any(r['accumulated']['maximum_absolute'] > 0 for r in self.report['independent_step_sensitivity']))
        self.assertTrue(any(r['u']['maximum_absolute'] > 0 for r in self.report['independent_volume_order_sensitivity']))

    def test_every_phase_interior_restart_matches_whole_run(self):
        self.assertEqual(len(self.report['restart_checks']), 7)
        for row in self.report['restart_checks']:
            self.assertGreater(row['progress'], 0); self.assertLess(row['progress'], 1)
            self.assertTrue(row['whole_run_exact_except_elapsed_seconds'])
        self.assertTrue(self.report['native_routes']['original']['all_rejections_preserve_committed_state'])

    def test_missing_frames_nonfinite_and_changed_displacement_fail(self):
        missing = copy.deepcopy(self.native); missing['frames'].pop()
        with self.assertRaisesRegex(AssertionError, 'missing native frames'):
            native_fields(missing)
        bad = copy.deepcopy(self.native)
        bad['frames'][2]['state']['displacement_m'][10][0] += 1
        with self.assertRaisesRegex(AssertionError, 'displacement gate'):
            compare(bad, self.expected)
        bad = copy.deepcopy(self.native)
        bad['frames'][2]['state']['joints']['points'][0]['plastic_slip_m'][0] = float('nan')
        with self.assertRaises(AssertionError):
            compare(bad, self.expected)

    def test_unsupported_water_domain_and_bulk_plasticity_refused(self):
        bad = copy.deepcopy(self.case); bad['stages'][0]['kind'] = 'excavate'
        with self.assertRaisesRegex(AssertionError, 'boundary cycle'):
            reference(self.mesh, self.ground, bad, 2, 1)
        bad = copy.deepcopy(self.case); bad['joint_system']['joints'][0]['water']['kind'] = 'soil'
        with self.assertRaisesRegex(AssertionError, 'CONTINUOUS'):
            reference(self.mesh, self.ground, bad, 2, 1)
        bad = copy.deepcopy(self.ground)
        next(iter(bad['initialization']['geostatic']['seed']['states'].values()))['equivalent_plastic_strain'] = 1e-14
        with self.assertRaises(AssertionError):
            reference(self.mesh, bad, self.case, 2, 1)
        bad = copy.deepcopy(self.native); bad['frames'][0]['state']['material_states'][0]['equivalent_plastic_strain'] = 1e-14
        with self.assertRaisesRegex(AssertionError, 'bulk plasticity'):
            native_fields(bad)


if __name__ == '__main__':
    unittest.main()
