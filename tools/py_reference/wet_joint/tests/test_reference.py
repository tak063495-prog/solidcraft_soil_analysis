import copy
import importlib
import json
from pathlib import Path
import tempfile
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tools.py_reference.wet_joint.geometry import (basis, tet_geometry, hydro_volume, triangle_rule,
                                                  pressure_face, refine, wet_simplices)
from tools.py_reference.wet_joint.solver import Model, run, coulomb, elastic, linear_correction
from tools.py_reference.wet_joint.compare import (pair, metric, observed, load_series, spatial_series,
                                                 save, load, diagnose, matching_nodes, validate)
from tools.py_reference.wet_joint.native import configuration, input_subset


def inputs():
    return [json.loads((ROOT/'examples'/name).read_text(encoding='utf-8')) for name in
            ('joint-water-column-tet10-mesh.json', 'ground-joint-water-mc.json', 'stages-joint-water-cycle.json')]


def unit_tet():
    from tools.py_reference.wet_joint.geometry import EDGES
    vertices = np.vstack((np.zeros(3), np.eye(3)))
    return np.vstack((vertices, [(vertices[i]+vertices[j])/2 for i, j in EDGES]))


class GeometryTests(unittest.TestCase):
    def test_partition_and_linear_completeness(self):
        points = unit_tet(); lam = np.array([[.1, .2, .3, .4], [.25]*4])
        n, grad = basis(lam)
        np.testing.assert_allclose(n.sum(axis=1), 1, atol=1e-14)
        np.testing.assert_allclose(grad.sum(axis=1), 0, atol=1e-14)
        np.testing.assert_allclose(n@points, lam@points[:4], atol=1e-14)

    def test_affine_patch_strain_and_energy(self):
        points = unit_tet(); _, b, _, volume = tet_geometry(points)
        affine = np.array([[.01, .02, .03], [.04, .05, .06], [.07, .08, .09]])
        u = (points@affine.T+np.array([2, -3, 4])).ravel()
        expected = np.array([.01, .05, .09, .06, .14, .10])
        strain = np.einsum('qai,i->qa', b, u)
        np.testing.assert_allclose(strain, np.tile(expected, (4, 1)), atol=1e-13)
        d = elastic(1000, .25)
        energy = .5*np.einsum('qa,ab,qb->', strain, d, strain)*volume/4
        self.assertAlmostEqual(energy, .5*expected@d@expected/6, places=12)

    def test_rigid_rotation_has_zero_strain(self):
        p = unit_tet(); _, b, _, _ = tet_geometry(p)
        u = (np.cross(np.array([.2, -.3, .1]), p)+[4, 2, 7]).ravel()
        np.testing.assert_allclose(np.einsum('qai,i->qa', b, u), 0, atol=1e-13)

    def test_cut_hydrostatic_virtual_work_exact(self):
        p = unit_tet(); force = hydro_volume(p, .5, 10).reshape(-1, 3)
        # For u=x, div(u)=3. Integral p over z<1/2 is 10*17/384.
        integral_p = 10*17/384
        self.assertAlmostEqual(np.sum(force*p), -3*integral_p, places=13)
        np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-13)

    def test_dry_and_fully_wet_pressure(self):
        p = unit_tet()
        np.testing.assert_array_equal(hydro_volume(p, -.1, 10), np.zeros(30))
        f = hydro_volume(p, 2, 10).reshape(-1, 3)
        self.assertAlmostEqual(np.sum(f*p), -3*10*(2-.25)/6, places=12)

    def test_clip_volume_and_product_triangle_moments(self):
        p = unit_tet()[:4]; pieces = wet_simplices(p, .5)
        volume = sum(abs(np.linalg.det((s[1:]-s[0]).T))/6 for s in pieces)
        self.assertAlmostEqual(volume, 7/48, places=14)
        l, w = triangle_rule()
        self.assertAlmostEqual(w.sum(), .5, places=14)
        np.testing.assert_allclose(w@l, [1/6]*3, atol=1e-14)

    def test_constant_pressure_resultant_and_moment(self):
        xyz = np.array([[0,0,0], [1,0,0], [0,1,0], [.5,0,0], [.5,.5,0], [0,.5,0]])
        f = pressure_face(xyz, list(range(6)), 12)
        np.testing.assert_allclose(f.sum(axis=0), [0, 0, -6], atol=1e-14)
        np.testing.assert_allclose(np.cross(xyz, f).sum(axis=0), [-2, 2, 0], atol=1e-14)

    def test_refinement_preserves_volume_constraints_and_seam_identity(self):
        mesh, ground, case = inputs(); fine = refine(mesh)
        self.assertEqual((len(fine['nodes']), len(fine['tetrahedra'])), (250, 96))
        total = sum(tet_geometry(np.array(fine['nodes'])[cell])[3] for cell in fine['tetrahedra10'])
        self.assertAlmostEqual(total, 1., places=12)
        m = Model(fine, ground, case)
        self.assertEqual(len(m.weight), 128)
        self.assertAlmostEqual(m.weight.sum(), 1., places=12)
        np.testing.assert_allclose(np.array(fine['nodes'])[fine['node_sets']['top'], 2], 1.)
        np.testing.assert_allclose(np.array(fine['nodes'])[fine['node_sets']['bottom'], 2], 0.)

    def test_curved_cell_is_rejected(self):
        p = unit_tet(); p[4, 2] += .01
        with self.assertRaisesRegex(ValueError, 'affine'):
            tet_geometry(p)


class ArchiveTests(unittest.TestCase):
    def make_archive(self, folder, name='field.txt'):
        import hashlib
        import zipfile
        from tools.py_reference.wet_joint.compare import digest
        payload = b'preserved evidence\n'; part = folder/'evidence.part001'
        with zipfile.ZipFile(part, 'w') as z:
            z.writestr(name, payload)
        index = dict(parts=[dict(file=part.name, bytes=part.stat().st_size, sha256=digest(part))],
                     archive_sha256=digest(part), entries=[dict(file=name, bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())])
        path = folder/'index.json'; path.write_text(json.dumps(index))
        return path, part

    def test_verified_archive_and_corrupt_part(self):
        from tools.py_reference.wet_joint.archive import extract
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); index, part = self.make_archive(root)
            report = extract(index, root/'good')
            self.assertTrue(report['software_integrity_passed'])
            self.assertFalse(report['numerical_acceptance'])
            self.assertEqual((root/'good/data/field.txt').read_bytes(), b'preserved evidence\n')
            with part.open('ab') as f:
                f.write(b'tamper')
            with self.assertRaisesRegex(ValueError, 'SHA'):
                extract(index, root/'bad')
            self.assertTrue((root/'bad/failure.json').exists())

    def test_parent_path_in_archive_is_rejected(self):
        from tools.py_reference.wet_joint.archive import extract
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); index, _ = self.make_archive(root, '../escape.txt')
            with self.assertRaisesRegex(ValueError, 'unsafe'):
                extract(index, root/'bad')
            self.assertFalse((root/'escape.txt').exists())


class ContactTests(unittest.TestCase):
    def update(self, start, traction, end, count=1):
        return coulomb(np.array([start], dtype=float), np.array([traction], dtype=float), np.array([end], dtype=float),
                       np.array([[0.,0.,1.]]), np.array([1000.]), np.array([100.]), np.array([.3]), count)

    def test_closed_stick_then_slip_closed_form(self):
        q, dp, acc = self.update([0,0,-.1], [0,0,-100], [.1,0,-.1])
        np.testing.assert_allclose(q[0,0], [10,0,-100], atol=1e-13)
        np.testing.assert_allclose(dp, 0, atol=1e-14)
        q, dp, acc = self.update([0,0,-.1], [0,0,-100], [.5,0,-.1])
        np.testing.assert_allclose(q[0,0], [30,0,-100], atol=1e-13)
        np.testing.assert_allclose(dp[0,0], [.2,0,0], atol=1e-14)
        self.assertAlmostEqual(acc[0,0], .2, places=14)

    def test_opening_and_reclosure_closed_fraction(self):
        q, dp, acc = self.update([0,0,.1], [0,0,0], [.1,0,-.1])
        np.testing.assert_allclose(q[0,0], [5,0,-100], atol=1e-13)
        q, dp, acc = self.update([0,0,-.1], [0,0,-100], [.1,0,.1])
        np.testing.assert_allclose(q, 0, atol=1e-14)
        self.assertAlmostEqual(acc[0,0], .05, places=14)

    def test_open_path_and_trial_state_are_immutable(self):
        start = np.array([[.2,0,.1]]); traction = np.zeros((1,3)); before = start.copy()
        q, dp, acc = coulomb(start, traction, np.array([[.8,.3,.2]]), np.array([[0.,0.,1.]]), np.array([1000.]), np.array([100.]), np.array([.3]), 8)
        np.testing.assert_array_equal(start, before); np.testing.assert_array_equal(traction, np.zeros((1,3)))
        np.testing.assert_array_equal(q, np.zeros((1,1,3))); np.testing.assert_array_equal(dp, q)

    def test_accumulated_length_is_not_net_vector(self):
        q, dp, acc = self.update([0,0,-.1], [0,0,-100], [1.,0,-.1], 8)
        q2, dp2, acc2 = self.update([1.,0,-.1], q[0,0], [-1.,0,-.1], 8)
        self.assertGreater(acc[0,0]+acc2[0,0], np.linalg.norm(dp[0,0]+dp2[0,0]))


class SolverAndComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, cls.ground, cls.case = inputs()
        cls.model = Model(cls.mesh, cls.ground, cls.case)
        cls.values = [run(cls.model, n) for n in (4,8,16)]

    def test_full_uniform_path_force_and_water_invariance(self):
        value = self.values[1]; audit = diagnose(value, self.model)
        self.assertEqual(audit['accepted_states'], 56)
        self.assertLess(audit['maximum_reassembled_force_residual'], 1e-10)
        self.assertEqual(set(audit['invariant_water_maximum_changes'].values()), {0.})
        self.assertFalse(value['metadata']['reference_uses_native_results'])
        self.assertGreater(np.max(value['fields']['gap']), 0)
        self.assertGreater(value['fields']['accumulated'][-1].max(), 0)

    def test_roundtrip_sha_and_missing_history_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'result'; save(path, self.values[0]); restored = load(path)
            self.assertEqual(pair(self.values[0], restored)['fields']['plastic_slip']['absolute_l2'], 0)
            with (path/'fields.npz').open('ab') as f:
                f.write(b'changed')
            with self.assertRaisesRegex(ValueError, 'SHA'):
                load(path)
        changed = copy.deepcopy(self.values[0]); changed['rows'].pop()
        with self.assertRaises(ValueError):
            validate(changed)

    def test_load_order_uses_same_coordinates(self):
        result = load_series(self.values)
        self.assertEqual(result['observed_orders'][0]['common_coordinates'], 28)
        self.assertFalse(result['numerical_acceptance'])
        x = np.array([4., 0]); y = np.array([2., 0])
        self.assertAlmostEqual(observed(x, y)['order'], 1)
        self.assertAlmostEqual(observed(y, x)['order'], -1)
        self.assertIsNone(observed(x*0, y)['order'])

    def test_zero_denominator_and_nonfinite_are_explicit(self):
        self.assertTrue(metric(np.ones(2), np.zeros(2))['relative_l2_undefined'])
        self.assertIsNone(metric(np.ones(2), np.zeros(2))['relative_l2'])
        with self.assertRaises(ValueError):
            metric(np.array([np.nan]), np.array([0.]))

    def test_duplicate_seam_coordinates_do_not_merge_bodies(self):
        a = self.values[0]['metadata']; ia, ib = matching_nodes(a, a)
        np.testing.assert_array_equal(ia, ib)
        self.assertEqual(len(ia), len(self.mesh['nodes']))

    def test_space_two_levels_has_no_order_and_rejects_mixed_time(self):
        fine = run(Model(refine(self.mesh), self.ground, self.case), 4)
        triangles = [np.array(self.mesh['nodes'])[f] for f in self.mesh['face_sets']['Seam:a']]
        value = spatial_series([self.values[0], fine], triangles)
        self.assertEqual(value['observed_orders'], [])
        self.assertEqual(value['pairs'][0]['common_nodes'], 54)
        self.assertFalse(value['pointwise_spatial_convergence_certified'])
        with self.assertRaises(ValueError):
            spatial_series([self.values[1], fine], triangles)

    def test_spatial_three_level_known_second_order(self):
        values = [copy.deepcopy(self.values[0]) for _ in range(3)]
        for i, value in enumerate(values):
            value['metadata']['elements'] = 12*8**i
            value['fields']['plastic_slip'][:] = 4.**(-i)
        triangles = [np.array(self.mesh['nodes'])[f] for f in self.mesh['face_sets']['Seam:a']]
        result = spatial_series(values, triangles)
        self.assertAlmostEqual(result['observed_orders'][0]['parent_area_means']['plastic_slip']['order'], 2.)
        self.assertFalse(result['pointwise_spatial_convergence_certified'])

    def test_true_linear_residual_does_not_hide_small_rhs(self):
        from unittest.mock import patch
        matrix = np.diag([2., 3.]); rhs = np.array([2e-8, 3e-8])
        np.testing.assert_allclose(linear_correction(matrix, rhs), [1e-8, 1e-8], rtol=1e-15)
        inaccurate = np.array([1e-8*(1+1e-8), 1e-8])
        with patch('tools.py_reference.wet_joint.solver.solve', return_value=inaccurate):
            with self.assertRaisesRegex(RuntimeError, 'true linear residual'):
                linear_correction(matrix, rhs)

    def test_python_only_study_does_not_require_native_step_defaults(self):
        from types import SimpleNamespace
        from tools.py_reference.wet_joint.study import execute
        with tempfile.TemporaryDirectory() as tmp:
            report = execute(SimpleNamespace(output=Path(tmp)/'new', steps=[2,4], levels=[0], user_pc=False,
                                            native_solver=None, native_steps=[16,32,64], native_timeout_seconds=900))
            self.assertFalse(report['passed'])
            self.assertEqual(report['load_series']['h0']['pairs'][0]['steps'], [2,4])

    def test_existing_output_is_never_changed_by_failure_recording(self):
        from unittest.mock import patch
        from tools.py_reference.wet_joint.study import main
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); (path/'prior.txt').write_text('preserved')
            with patch('sys.argv', ['study', '--output', str(path)]):
                with self.assertRaises(FileExistsError):
                    main()
            self.assertEqual({p.name for p in path.iterdir()}, {'prior.txt'})
            self.assertEqual((path/'prior.txt').read_text(), 'preserved')

    def test_force_tolerance_may_only_be_tightened(self):
        with self.assertRaisesRegex(ValueError, 'tightened'):
            run(self.model, 4, tolerance=1e-8)
        tight = run(self.model, 4, tolerance=1e-11)
        self.assertLess(max(r['relative_force_residual'] for r in tight['rows']), 1e-11)

    def test_unsupported_water_and_plastic_seed_fail(self):
        case = copy.deepcopy(self.case); case['joint_system']['joints'][0]['water']['kind'] = 'soil'
        with self.assertRaisesRegex(ValueError, 'CONTINUOUS'):
            Model(self.mesh, self.ground, case)
        ground = copy.deepcopy(self.ground)
        ground['initialization']['geostatic']['seed']['states']['LowerClay']['equivalent_plastic_strain'] = .01
        with self.assertRaisesRegex(ValueError, 'seed'):
            Model(self.mesh, ground, self.case)

    def test_bulk_yield_is_rejected_before_acceptance(self):
        ground = copy.deepcopy(self.ground)
        for material in ground['materials']:
            material['mohr_coulomb']['cohesion_kpa'] = 0.
            material['mohr_coulomb']['friction_deg'] = 0.
        model = Model(self.mesh, ground, self.case)
        with self.assertRaisesRegex(ValueError, 'bulk yielded'):
            model.check_elastic(np.arange(self.model.xyz.size)*1e-3)

    def test_native_configuration_does_not_change_tolerances_or_inputs(self):
        before = copy.deepcopy(self.case); cfg = configuration(self.case, 32)
        self.assertEqual(self.case, before)
        for key in ('relative_force_tolerance', 'relative_linear_tolerance'):
            self.assertEqual(cfg['settings'][key], before['settings'][key])
        self.assertEqual(cfg['settings']['initial_increment'], 1/32)
        self.assertTrue(all(s['initial_increment'] == 1/32 for s in cfg['stages']))

    def test_comparison_rejects_different_factors_or_contact_controls(self):
        wrong = copy.deepcopy(self.values[1]); wrong['rows'][1]['factor'] = 99.
        with self.assertRaisesRegex(ValueError, 'factors'):
            pair(self.values[0], wrong)
        wrong = copy.deepcopy(self.values[1]); wrong['metadata']['local'] = 16
        with self.assertRaisesRegex(ValueError, 'contact integration'):
            load_series([self.values[0], wrong])

    def test_historical_record_is_read_without_rewriting(self):
        from tools.py_reference.wet_joint.compare import historical_comparison, digest
        path = ROOT/'interim/2026-11-15-wet-joint-convergence/measurements.json'
        before = digest(path)
        current = dict(pairs=[dict(steps=[16, 32], fields=dict(plastic_slip=dict(relative_l2=.7)))])
        result = historical_comparison(current, path)
        self.assertFalse(result['historical_passed'])
        self.assertFalse(result['numerical_acceptance'])
        self.assertEqual(result['comparisons'][0]['new_relative_l2'], .7)
        self.assertAlmostEqual(result['comparisons'][0]['recorded_relative_l2'], .6762689785372746)
        self.assertEqual(digest(path), before)

    def test_native_physical_input_mutation_is_rejected(self):
        expected = configuration(self.case, 32); actual = copy.deepcopy(expected)
        actual['boundary']['constraints'][3]['value_m'] += .001
        with self.assertRaisesRegex(ValueError, 'native input value differs'):
            input_subset(expected, actual)
        actual = copy.deepcopy(expected); actual['extra_serializer_default'] = []
        input_subset(expected, actual)

    def test_solver_does_not_import_legacy_implementations(self):
        import ast
        for name in ('geometry', 'solver'):
            module = importlib.import_module('tools.py_reference.wet_joint.'+name)
            tree = ast.parse(Path(module.__file__).read_text())
            imports = [n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
            imports += [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
            self.assertTrue(all(m.startswith(('numpy', 'scipy.', 'geometry', 'copy', 'itertools')) for m in imports), imports)


if __name__ == '__main__':
    unittest.main()
