"""Explicit PRESCRIBED scope, independent 3D fields and refusal gates."""
import copy
import os
from pathlib import Path
import sys
import tempfile
import unittest

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
import verify_wet_structures as wet


class PrescribedWetStructures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.output = Path(cls.temp.name)/'prescribed'
        cls.report = wet.run(solver, cls.output, 'prescribed')
        cls.mesh = wet.read(cls.output/'mesh.json')
        cls.ground = fem.normalize_ground(wet.read(cls.output/'ground.json'))
        cls.case = fem.normalize_stage_case(wet.read(cls.output/'stages.json'))
        cls.native = wet.read(cls.output/'stage-results.json')
        cls.expected = wet.reference(cls.mesh, cls.ground, cls.case, water_modes=('prescribed',))

    def test_full_fields_quadrature_and_all_stage_restarts(self):
        self.assertTrue(self.report['passed'])
        self.assertEqual(self.report['water_modes'], ['prescribed'])
        self.assertEqual(len(self.report['frames']), 8)
        self.assertEqual(len(self.report['checkpoints']), 8)
        self.assertLess(self.report['independent_gl5_to_gl7_relative_difference'], 1e-7)
        self.assertFalse(self.report['p6_complete'])
        self.assertFalse(self.report['p3b_p4_rows_complete'])

    def test_physical_pressure_both_side_loads_and_member_reference_tampering_fail(self):
        for target in ('physical_pressure', 'side_load', 'placement', 'source'):
            with self.subTest(target=target):
                changed = copy.deepcopy(self.native)
                state = changed['frames'][4]['state']
                if target == 'physical_pressure':
                    state['joint_water']['physical_pressure_kpa'][0] += 1.
                elif target == 'side_load':
                    state['joint_water']['nodal_load_kn'][0] += 1.
                elif target == 'placement':
                    state['beams']['plate_reference_global_kinematics'][0][0] += .001
                else:
                    state['joint_water']['sources'][0] = 'continuous'
                with self.assertRaises(AssertionError):
                    wet.compare(self.mesh, self.ground, self.case, changed, self.expected)

    def test_unscoped_fixed_head_and_soil_members_are_refused(self):
        with self.assertRaisesRegex(AssertionError, 'explicit reference scope'):
            wet.reference(self.mesh, self.ground, self.case)
        fixed = copy.deepcopy(self.case)
        fixed['joint_system']['joints'][0]['water']['level_m'] = 1.25
        with self.assertRaisesRegex(AssertionError, 'fixed prescribed head'):
            wet.reference(self.mesh, self.ground, fixed, water_modes=('prescribed',))
        soil = copy.deepcopy(self.case)
        source = wet.read(ROOT/'examples/stages-joint-water-soil.json')['joint_system']
        soil['joint_system']['joints'][0]['water'] = source['joints'][0]['water']
        soil['joint_system']['structure_bodies'] = source['structure_bodies']
        with self.assertRaises(AssertionError):
            wet.reference(self.mesh, self.ground, soil, water_modes=('soil',))


if __name__ == '__main__':
    unittest.main()
