"""Nonzero SOIL flow uses the existing geometry, Darcy source and member laws."""
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
from verify_soil_source_column import audit


class SoilSourceWetStructures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.temp.cleanup)
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.output = Path(cls.temp.name)/'soil-source'
        cls.report = wet.run(solver, cls.output, 'soil-source')
        cls.mesh = wet.read(cls.output/'mesh.json')
        cls.ground = fem.normalize_ground(wet.read(cls.output/'ground.json'))
        cls.case = fem.normalize_stage_case(wet.read(cls.output/'stages.json'))
        cls.native = wet.read(cls.output/'stage-results.json')
        cls.expected = wet.reference(cls.mesh, cls.ground, cls.case, water_modes=('soil',))

    def test_same_small_geometry_roles_and_source_only_flow_difference(self):
        self.assertEqual(self.mesh, wet.read(ROOT/'examples/joint-water-column-tet10-mesh.json'))
        source = wet.read(ROOT/'examples/stages-joint-water-soil.json')
        self.assertEqual(self.case['joint_system']['structure_bodies'], source['joint_system']['structure_bodies'])
        raw = wet.read(self.output/'stages.json')
        flow = copy.deepcopy(next(s for s in raw['stages'] if s['kind'] == 'steady_seepage')['flow'])
        for material in flow['materials']:
            self.assertEqual(material.pop('volume_source_per_s'), 2*material['hydraulic_conductivity_m_s'][2][2])
        self.assertEqual(flow, source['stages'][0]['flow'])
        self.assertEqual(len(raw['beam_system']['plates']), 2)
        self.assertEqual(len(raw['beam_system']['members']), 1)

    def test_nonzero_analytic_flow_full_3d_fields_and_restart_scope(self):
        self.assertTrue(self.report['passed'] and self.report['soil_volume_source'])
        row = self.report['analytic_soil_source']['frames'][0]
        self.assertEqual(row['head_at_seam_m'], 2.25)
        self.assertEqual(row['physical_joint_pressure_kpa'], 17.5)
        self.assertGreater(row['source_and_bottom_outflow_m3_s'], 0)
        flow = self.native['flow_audits'][0]
        self.assertTrue(any(abs(p['darcy_flux_m_s'][2]) > 1e-8 for p in flow['solution']['points']))
        self.assertEqual(len(self.report['frames']), 8)
        self.assertEqual(len(self.report['checkpoints']), 8)
        self.assertTrue(self.report['clean_rejected_path_exact'])
        self.assertFalse(self.report['p6_complete'] or self.report['p3b_p4_rows_complete'])

    def test_flux_head_pressure_and_balance_corruption_are_rejected(self):
        for target in ('flux', 'head', 'pressure', 'source_balance', 'outflow'):
            with self.subTest(target=target):
                native = copy.deepcopy(self.native)
                flow = native['flow_audits'][0]
                if target == 'flux': flow['solution']['points'][0]['darcy_flux_m_s'][2] += 1e-5
                elif target == 'head': flow['joint_flow']['original_field']['head_m'][0] += .01
                elif target == 'pressure': native['frames'][flow['stage_index']+1]['state']['joint_water']['physical_pressure_kpa'][0] += 1.
                elif target == 'source_balance': flow['solution']['balances'][0]['source_m3_s'] = 0.
                else: flow['solution']['balances'][0]['dirichlet_reaction_outflow_m3_s'] = 0.
                with self.assertRaises(AssertionError): audit(self.mesh, self.ground, self.case, native)

    def test_zero_source_structure_pressure_and_side_force_not_nonzero_acceptance(self):
        case = copy.deepcopy(self.case)
        next(s for s in case['stages'] if s['kind'] == 'steady_seepage')['flow']['materials'][0]['volume_source_per_s'] = 0.
        with self.assertRaises(AssertionError): audit(self.mesh, self.ground, case, self.native)
        water = wet.Water(self.mesh, self.ground, self.case)
        index = self.native['flow_audits'][0]['stage_index']+1
        point = next(i for i, p in enumerate(self.native['integration_points']) if not water.soil[p['element_index']])
        for target in ('pore', 'side_force', 'placement'):
            with self.subTest(target=target):
                native = copy.deepcopy(self.native); state = native['frames'][index]['state']
                if target == 'pore': state['pore_pressure_kpa'][point] = 1.
                elif target == 'side_force': state['joint_water']['face_nodal_force_a_kn'][0][0][2] += 1.
                else: state['beams']['plate_reference_global_kinematics'][0][0] += .001
                with self.assertRaises(AssertionError): wet.compare(self.mesh, self.ground, self.case, native, self.expected)

    def test_plastic_unscoped_and_structural_head_are_refused(self):
        with self.assertRaises(AssertionError): wet.reference(self.mesh, self.ground, self.case)
        ground = copy.deepcopy(self.ground)
        next(iter(ground['initialization']['geostatic']['seed']['states'].values()))['equivalent_plastic_strain'] = .001
        with self.assertRaises(AssertionError): wet.reference(self.mesh, ground, self.case, water_modes=('soil',))
        case = copy.deepcopy(self.case)
        next(s for s in case['stages'] if s['kind'] == 'steady_seepage')['flow']['heads'].append(dict(face_set='top', head_m=3.))
        with self.assertRaises(AssertionError): wet.reference(self.mesh, self.ground, case, water_modes=('soil',))


if __name__ == '__main__':
    unittest.main()
