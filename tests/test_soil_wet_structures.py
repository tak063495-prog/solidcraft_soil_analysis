"""SOIL continuum role, inherited flow boundary and independent member fields."""
import copy
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
import verify_wet_structures as wet


class SoilWetStructures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.output = Path(cls.temp.name)/'soil'
        cls.report = wet.run(solver, cls.output, 'soil')
        cls.mesh = wet.read(cls.output/'mesh.json')
        cls.ground = fem.normalize_ground(wet.read(cls.output/'ground.json'))
        cls.case = fem.normalize_stage_case(wet.read(cls.output/'stages.json'))
        cls.native = wet.read(cls.output/'stage-results.json')
        cls.expected = wet.reference(cls.mesh, cls.ground, cls.case, water_modes=('soil',))

    def test_inherited_roles_flow_and_limited_full_field_gate(self):
        source = wet.read(ROOT/'examples/stages-joint-water-soil.json')
        self.assertEqual(self.case['joint_system']['structure_bodies'], source['joint_system']['structure_bodies'])
        flow = next(s for s in wet.read(self.output/'stages.json')['stages'] if s['kind'] == 'steady_seepage')
        self.assertEqual(flow['flow'], source['stages'][0]['flow'])
        self.assertTrue(self.report['passed'])
        self.assertEqual(len(self.report['frames']), 8)
        self.assertEqual(len(self.report['checkpoints']), 8)
        self.assertFalse(self.report['p6_complete'])
        self.assertFalse(self.report['p3b_p4_rows_complete'])

    def test_structure_pore_water_force_and_placement_corruption_fail(self):
        water = wet.Water(self.mesh, self.ground, self.case)
        structure_point = next(i for i, p in enumerate(self.native['integration_points'])
                               if not water.soil[p['element_index']])
        for target in ('pore', 'pressure', 'side_load', 'placement', 'soil', 'head'):
            with self.subTest(target=target):
                native = copy.deepcopy(self.native)
                state = native['frames'][4]['state']
                if target == 'pore':
                    state['pore_pressure_kpa'][structure_point] = 1.
                elif target == 'pressure':
                    state['joint_water']['physical_pressure_kpa'][0] += 1.
                elif target == 'side_load':
                    state['joint_water']['face_nodal_force_a_kn'][0][0][2] += 1.
                elif target == 'placement':
                    state['beams']['plate_reference_global_kinematics'][0][0] += .001
                elif target == 'soil':
                    state['displacement_m'][8][2] += .001
                else:
                    native['flow_audits'][0]['joint_flow']['original_field']['defined_nodes'][0] = False
                with self.assertRaises(AssertionError):
                    wet.compare(self.mesh, self.ground, self.case, native, self.expected)

    def test_missing_role_wrong_soil_side_and_structure_head_boundary_fail(self):
        for target in ('role', 'side', 'top_head'):
            with self.subTest(target=target):
                case = copy.deepcopy(self.case)
                if target == 'role':
                    case['joint_system']['structure_bodies'] = []
                elif target == 'side':
                    case['joint_system']['joints'][0]['water']['soil_side'] = 'b'
                else:
                    flow = next(s for s in case['stages'] if s['kind'] == 'steady_seepage')
                    flow['flow']['heads'].append(dict(face_set='top', head_m=3.))
                with self.assertRaises(AssertionError):
                    wet.reference(self.mesh, self.ground, case, water_modes=('soil',))

    def test_plastic_and_unscoped_soil_are_not_accepted(self):
        with self.assertRaisesRegex(AssertionError, 'explicit reference scope'):
            wet.reference(self.mesh, self.ground, self.case)
        ground = copy.deepcopy(self.ground)
        next(iter(ground['initialization']['geostatic']['seed']['states'].values()))['equivalent_plastic_strain'] = .001
        with self.assertRaises(AssertionError):
            wet.reference(self.mesh, ground, self.case, water_modes=('soil',))


if __name__ == '__main__':
    unittest.main()
