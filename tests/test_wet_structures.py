"""Real small CLI + independent 3D reference and adversarial evidence checks."""
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


class WetStructures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.output = Path(cls.temp.name)/'wet'
        cls.report = wet.run(solver, cls.output)
        cls.mesh = wet.read(cls.output/'mesh.json')
        cls.ground = fem.normalize_ground(wet.read(cls.output/'ground.json'))
        cls.case = fem.normalize_stage_case(wet.read(cls.output/'stages.json'))
        cls.native = wet.read(cls.output/'stage-results.json')
        cls.expected = wet.reference(cls.mesh, cls.ground, cls.case)

    def test_small_integrated_3d_and_eight_separate_process_restarts(self):
        self.assertTrue(self.report['passed'])
        self.assertEqual(len(self.report['frames']), 8)
        self.assertEqual(len(self.report['checkpoints']), 8)
        self.assertEqual({c['stage'] for c in self.report['checkpoints']},
                         {'geostatic'} | {s['name'] for s in self.case['stages']})
        self.assertFalse(self.report['p6_complete'])
        self.assertFalse(self.report['p3b_p4_rows_complete'])

    def test_wrong_physical_water_and_contact_traction_are_rejected(self):
        for target in ('water', 'traction', 'missing_flow'):
            native = copy.deepcopy(self.native)
            state = native['frames'][4]['state']
            if target == 'water':
                state['joint_water']['physical_pressure_kpa'][0] += 1.
            elif target == 'traction':
                state['joints']['points'][0]['traction_kpa'][2] += 1.
            else:
                native['flow_audits'] = []
            with self.assertRaises(AssertionError):
                wet.compare(self.mesh, self.ground, self.case, native, self.expected)

    def test_wrong_structural_placement_and_soil_field_are_rejected(self):
        for target in ('placement', 'soil'):
            native = copy.deepcopy(self.native)
            state = native['frames'][3]['state']
            if target == 'placement':
                state['beams']['plate_reference_global_kinematics'][0][0] += .001
            else:
                state['displacement_m'][8][2] += .001
            with self.assertRaises(AssertionError):
                wet.compare(self.mesh, self.ground, self.case, native, self.expected)

    def test_unverified_birth_richards_partial_wetting_and_sliding_cannot_pass(self):
        for kind in ('birth', 'richards'):
            case = copy.deepcopy(self.case); case['stages'][0]['kind'] = kind
            with self.assertRaisesRegex(AssertionError, 'unsupported stage'):
                wet.reference(self.mesh, self.ground, case)
        ground = copy.deepcopy(self.ground)
        ground['initialization']['water_level_m'] = .75
        with self.assertRaisesRegex(AssertionError, 'partial volume wetting'):
            wet.reference(self.mesh, ground, self.case)
        case = copy.deepcopy(self.case)
        case['joint_system']['joints'][0]['material']['friction_coefficient'] = 0.
        with self.assertRaisesRegex(AssertionError, 'sliding reference refused'):
            wet.reference(self.mesh, self.ground, case)

    def test_soil_and_prescribed_with_members_are_not_continuous_acceptance(self):
        for name in ('soil', 'prescribed'):
            with self.subTest(water=name):
                source = wet.read(ROOT/'examples'/('stages-joint-water-'+name+'.json'))
                case = copy.deepcopy(self.case)
                case['joint_system']['joints'][0]['water'] = source['joint_system']['joints'][0]['water']
                with self.assertRaises(AssertionError):
                    wet.reference(self.mesh, self.ground, case)

    def test_plastic_seed_cannot_be_certified_by_elastic_member_reference(self):
        ground = copy.deepcopy(self.ground)
        seed = ground['initialization']['geostatic']['seed']['states']
        next(iter(seed.values()))['equivalent_plastic_strain'] = .001
        with self.assertRaises(AssertionError):
            wet.reference(self.mesh, ground, self.case)


if __name__ == '__main__':
    unittest.main()
