"""Existing fixed-reservoir policy with small members, plus refusal gates."""
import copy
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
import fem_reference as fem
import verify_wet_structures as wet


class FixedPrescribedWetStructures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.output = Path(cls.temp.name)/'fixed'
        cls.report = wet.run(cls.solver, cls.output, 'prescribed-fixed')
        cls.mesh = wet.read(cls.output/'mesh.json')
        cls.ground = fem.normalize_ground(wet.read(cls.output/'ground.json'))
        cls.case = fem.normalize_stage_case(wet.read(cls.output/'stages.json'))
        cls.native = wet.read(cls.output/'stage-results.json')
        cls.expected = wet.reference(cls.mesh, cls.ground, cls.case,
                                     water_modes=('prescribed',), fixed_prescribed=True)

    def test_existing_inputs_and_limited_full_fields_restarts(self):
        self.assertEqual(self.mesh, wet.read(ROOT/'examples/joint-water-column-tet10-mesh.json'))
        self.assertEqual(wet.read(self.output/'ground.json'), wet.read(ROOT/'examples/ground-wet-plate-beam.json'))
        source = wet.read(ROOT/'examples/stages-wet-plate-beam.json')
        actual = wet.read(self.output/'stages.json')
        source['joint_system']['joints'][0]['water'] = dict(kind='prescribed', level_m=1.)
        source['stages'] = [s for s in source['stages'] if s['name'] != 'rise']
        self.assertEqual(actual, source)
        self.assertTrue(self.report['passed'])
        self.assertTrue(self.report['fixed_prescribed'])
        self.assertEqual(len(self.report['frames']), 7)
        self.assertEqual(len(self.report['checkpoints']), 7)
        self.assertFalse(self.report['p6_complete'])
        self.assertFalse(self.report['p3b_p4_rows_complete'])

    def test_fixed_interface_head_is_not_the_exterior_or_average(self):
        water = wet.Water(self.mesh, self.ground, self.case)
        pressures = np.array([water.gamma*(1.-p['position_m'][2]) for p in water.contact.points])
        for frame in self.native['frames']:
            np.testing.assert_allclose(frame['state']['joint_water']['physical_pressure_kpa'], pressures, rtol=0, atol=1e-9)
        flow = self.native['flow_audits'][0]
        self.assertEqual(next(s for s in self.case['stages'] if s['kind'] == 'steady_seepage')['exterior_water_level_m'], .75)
        field = flow['joint_flow']['original_field']['head_m']
        for point in water.faces:
            for side in ('a', 'b'):
                np.testing.assert_allclose(np.asarray(field)[point['nodes_'+side]], 1., rtol=0, atol=1e-9)
        row = next(r for r in self.report['frames'] if r['name'] == 'flow')
        self.assertLess(row['independent_darcy']['maximum_flux_error_m_s'], 1e-12)
        self.assertTrue(any(abs(v) > 1e-8 for p in flow['solution']['points'] for v in p['darcy_flux_m_s']))

    def test_flow_flux_side_pressure_and_placement_tampering_fail(self):
        for target in ('flux', 'head', 'pressure', 'side', 'placement'):
            with self.subTest(target=target):
                native = copy.deepcopy(self.native)
                frame = next(f for f in native['frames'] if f['name'] == 'flow')
                if target == 'flux':
                    native['flow_audits'][0]['solution']['points'][0]['darcy_flux_m_s'][2] += 1e-5
                elif target == 'head':
                    native['flow_audits'][0]['joint_flow']['original_field']['head_m'][0] += .01
                elif target == 'pressure':
                    frame['state']['joint_water']['physical_pressure_kpa'][0] += 1.
                elif target == 'side':
                    frame['state']['joint_water']['face_nodal_force_a_kn'][0][0][2] += 1.
                else:
                    frame['state']['beams']['plate_reference_global_kinematics'][0][0] += .001
                with self.assertRaises(AssertionError):
                    wet.compare(self.mesh, self.ground, self.case, native, self.expected)

    def test_unscoped_mismatched_reservoir_and_plastic_seed_refused(self):
        with self.assertRaisesRegex(AssertionError, 'fixed prescribed head'):
            wet.reference(self.mesh, self.ground, self.case, water_modes=('prescribed',))
        for target in ('initial', 'hydrostatic', 'plastic'):
            with self.subTest(target=target):
                ground = copy.deepcopy(self.ground); case = copy.deepcopy(self.case)
                if target == 'initial':
                    case['joint_system']['joints'][0]['water']['level_m'] = 1.25
                elif target == 'hydrostatic':
                    next(s for s in case['stages'] if s['kind'] == 'hydrostatic_water')['target_level_m'] = 1.25
                else:
                    next(iter(ground['initialization']['geostatic']['seed']['states'].values()))['equivalent_plastic_strain'] = .001
                with self.assertRaises(AssertionError):
                    wet.reference(self.mesh, ground, case, water_modes=('prescribed',), fixed_prescribed=True)

    def test_native_conflicting_hydrostatic_level_refused_before_output(self):
        case = wet.read(self.output/'stages.json')
        reset = next(s for s in case['stages'] if s['kind'] == 'hydrostatic_water')
        reset['target_level_m'] = 1.25
        path = self.output/'conflicting-stages.json'
        path.write_text(wet.json.dumps(case), encoding='utf-8')
        result = self.output/'conflicting-result.json'
        wet.invoke([str(self.solver), 'solve-stages', str(self.output/'mesh.json'),
                    str(self.output/'ground.json'), str(path), str(result)], self.output/'conflict.log', 1)
        self.assertFalse(result.exists())
        self.assertIn('prescribed', (self.output/'conflict.log').read_text().lower())


if __name__ == '__main__':
    unittest.main()
