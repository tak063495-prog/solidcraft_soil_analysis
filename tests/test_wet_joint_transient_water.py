"""C-04 observation: fixed mechanical path, conservative PDE and refusals."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-22-wet-joint-transient-water'
spec = importlib.util.spec_from_file_location('c04_study_tests', HERE/'study.py')
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
spec = importlib.util.spec_from_file_location('c04_measure_tests', HERE/'measure.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class TransientWater(unittest.TestCase):
    def test_only_water_stage_steps_change_in_three_level_matrix(self):
        cases = [s.mechanical_case(32, n) for n in (8, 16, 32)]
        canonical = []
        for case, n in zip(cases, (8, 16, 32)):
            phases = s.schedule(case)
            self.assertEqual([p['n'] for p in phases if p['kind'] != 'hydrostatic_water'], [32]*7)
            self.assertEqual([p['n'] for p in phases if p['kind'] == 'hydrostatic_water'], [n]*6)
            for stage in case['stages']:
                if stage['kind'] == 'hydrostatic_water': stage['initial_increment'] = None
            canonical.append(case)
        self.assertEqual(canonical[0], canonical[1]); self.assertEqual(canonical[1], canonical[2])

    def test_hydrostatic_exterior_buoyancy_and_interface_pair(self):
        mesh = s.read(ROOT/'examples/joint-water-column-mesh.json')
        ground = s.fem.normalize_ground(s.read(ROOT/'examples/ground-joint-water-mc.json'))
        ground['initialization']['geostatic']['hydrostatic_face_sets'] = ['exterior']
        model = s.mechanics.Model(mesh, ground, s.mechanical_case(32, 8))
        for level in (.25, .75, 1.25):
            model.set_level(level)
            # Unit cube, wet height clipped to [0,1]: Archimedes' resultant.
            np.testing.assert_allclose(model.exterior_water.sum(axis=0), [0, 0, 10*min(level, 1)], rtol=0, atol=1e-12)
            np.testing.assert_allclose(model.interface['load'].sum(axis=0), 0, rtol=0, atol=1e-12)
            np.testing.assert_allclose(model.interface['pressure'], 10*max(level-.5, 0), rtol=0, atol=1e-12)

    def test_cut_water_force_is_reassembled_at_intermediate_level(self):
        mesh = s.read(ROOT/'examples/joint-water-column-mesh.json')
        ground = s.fem.normalize_ground(s.read(ROOT/'examples/ground-joint-water-mc.json'))
        model = s.mechanics.Model(mesh, ground, s.mechanical_case(32, 8))
        model.set_level(.25); low = model.pore_force.copy()
        model.set_level(1.25); high = model.pore_force.copy()
        model.set_level(.75); actual = model.pore_force.copy()
        self.assertGreater(np.linalg.norm(actual-(low+high)/2), .1)
        model.set_level(.25); np.testing.assert_array_equal(model.pore_force, low)

    def test_hybrid_saturated_affine_head_patch_has_exact_darcy_flux(self):
        mesh = s.read(ROOT/'examples/seepage-column-tet4-mesh.json')
        case = s.hydraulic_case(4)
        gradient = np.array([.04, -.02, .03]); affine = dict(head_m=2., gradient=gradient.tolist())
        case['initial_head'] = affine
        topology = {}
        for cell in mesh['tetrahedra']:
            for i in range(4):
                key = tuple(sorted(np.delete(cell, i))); topology[key] = topology.get(key, 0)+1
        mesh['face_sets']['all-exterior'] = [list(k) for k, n in topology.items() if n == 1]
        period = dict(name='patch', end_time_s=.25, boundaries=[dict(face_set='all-exterior', condition=dict(kind='head', value=affine))])
        case['analysis']['periods'] = [period]
        model = s.hydraulics.Model(mesh, case); result = model.solve()['rows'][0]
        vertices = np.asarray(mesh['nodes'])[mesh['tetrahedra']]
        np.testing.assert_allclose(result['cell_head_m'], 2.+vertices.mean(axis=1)@gradient, rtol=0, atol=1e-13)
        for e, cell in enumerate(vertices):
            for i in range(4):
                face = np.delete(cell, i, axis=0); area = np.cross(face[1]-face[0], face[2]-face[0])/2
                if area@(cell[i]-face[0]) > 0: area = -area
                self.assertAlmostEqual(float(result['outflow_m3_s'][e, i]), float(-.02*gradient@area), places=14)

    def test_closed_saturated_source_stores_exact_water_without_boundary_flux(self):
        mesh = s.read(ROOT/'examples/seepage-column-tet4-mesh.json'); case = s.hydraulic_case(4)
        case['initial_head'] = dict(head_m=2.)
        for mat in case['materials']: mat['volume_source_per_s'] = .001
        case['analysis']['periods'] = [dict(name='closed', end_time_s=.25, boundaries=[])]
        result = s.hydraulics.Model(mesh, case).solve()['rows'][0]
        np.testing.assert_allclose(result['cell_head_m'], 2.+.001/.02*.25, rtol=0, atol=1e-13)
        self.assertLess(abs(result['boundary_outflow_m3_s']), 1e-14)
        self.assertLess(abs(result['global_mass_defect_m3']), 1e-14)
        self.assertLess(result['maximum_scaled_residual'], .01)

    def test_vg_reference_refuses_unimplemented_atmospheric_policy(self):
        mesh = s.read(ROOT/'examples/seepage-column-tet4-mesh.json'); case = s.hydraulic_case(4)
        case['analysis']['periods'][0]['boundaries'][0]['condition'] = dict(kind='atmospheric_no_ponding', rainfall_vertical_m_s=1e-6)
        with self.assertRaisesRegex(ValueError, 'head/outward-flux'): s.hydraulics.Model(mesh, case).solve()

    def test_missing_interior_state_or_changed_full_field_is_exposed(self):
        with tempfile.TemporaryDirectory() as directory:
            a = Path(directory)/'a'; b = Path(directory)/'b'; a.mkdir(); b.mkdir()
            def rows(value, drop=False):
                values = [dict(row=dict(stage='water', progress=p, factor=0., time_s=p), fields=dict(u=[[value, 0., 0.]], plastic_slip=[[0., 0., 0.]])) for p in (.5, 1.)]
                if drop: values.pop(0)
                return ''.join(json.dumps(v)+'\n' for v in values)
            (a/'fields.jsonl').write_text(rows(0.)); (b/'fields.jsonl').write_text(rows(1e-6))
            compared = m.pair(a, b)
            self.assertFalse(compared['numerical_acceptance'])
            self.assertEqual(compared['fields']['u']['maximum_absolute'], 1e-6)
            (b/'fields.jsonl').write_text(rows(0., True))
            with self.assertRaisesRegex(AssertionError, 'missing common'): m.pair(a, b)

    def test_independent_mechanical_scope_does_not_silently_accept_richards(self):
        case = s.mechanical_case(32, 8); case['stages'][1]['kind'] = 'richards'
        with self.assertRaisesRegex(ValueError, 'Richards JOINT unsupported'):
            s.mechanics.Model(s.read(ROOT/'examples/joint-water-column-mesh.json'), s.read(ROOT/'examples/ground-joint-water-mc.json'), case)

    def test_real_one_way_cli_refuses_joint_before_producing_result(self):
        case = s.read(ROOT/'examples/one-way-mc.json'); wet = s.mechanical_case(32, 8)
        case['flow'] = s.hydraulic_case(4)
        case['exterior_periods'] = [dict(name=p['name'], reservoirs=[]) for p in case['flow']['analysis']['periods']]
        case['mechanics'].update(stages=[], boundary=wet['boundary'], joint_system=wet['joint_system'], pressure_on_original_surface_only=True)
        case['mechanics']['settings']['kinematics'] = 'explicit'
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if sys.platform == 'win32' else 'soil-fem-3d')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory); s.write(path/'case.json', case)
            result = subprocess.run([str(solver), 'solve-one-way', str(ROOT/'examples/joint-water-column-mesh.json'),
                                     str(ROOT/'examples/ground-one-way-mc.json'), str(path/'case.json'), str(path/'result.json')], capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('ONE_WAY requires dry-reference weights', result.stderr)
            self.assertFalse((path/'result.json').exists())


if __name__ == '__main__': unittest.main()
