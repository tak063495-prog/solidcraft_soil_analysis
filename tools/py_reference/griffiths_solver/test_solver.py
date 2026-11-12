"""Assembly/solver gates: actual small Rust solve, frozen systems, failure paths."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import splu

from tools.py_reference.griffiths_solver import systems as s
from tools.py_reference.griffiths_solver import solvers as l
from plastic_spatial_convergence import inputs, lift


def small():
    mesh, ground, case = inputs(2)
    mesh, _ = lift(mesh, [[0., 0., 0.]]*len(mesh['nodes']))
    xyz = np.asarray(mesh['nodes'])
    mesh['node_sets'] = dict(all=list(range(len(xyz))),
                            sides=np.flatnonzero((xyz[:, 0] == 0) | (xyz[:, 0] == 1)).tolist(),
                            bottom=np.flatnonzero(xyz[:, 2] == 0).tolist())
    material = ground['materials'][0]['mohr_coulomb']
    material.update(young_kpa=100000., poisson=.3, cohesion_kpa=1e6, friction_deg=0.)
    case['stages'] = [dict(kind='boundary', name='hold', target_factor=0., initial_increment=1.)]
    return mesh, ground, case


class Assembly(unittest.TestCase):
    def test_tet10_partition_volume_and_affine_patch_energy(self):
        mesh, _, case = small()
        geometry = s.fem.Geometry(mesh)
        xyz = np.asarray(mesh['nodes'])
        gradient = np.array([[.01, .02, -.01], [.03, -.02, .01], [-.01, .02, .03]])
        u = xyz @ gradient.T + [.1, -.2, .3]
        expected = np.r_[np.diag(gradient), gradient[0, 1]+gradient[1, 0],
                         gradient[1, 2]+gradient[2, 1], gradient[2, 0]+gradient[0, 2]]
        integrated = 0.
        d = s.fem.elastic_tensor(100000., .3)
        for i, cell in enumerate(s.fem.cells(mesh)):
            q = geometry.volume(i)
            strain = np.einsum('qai,i->qa', q['b'], u[cell].ravel())
            np.testing.assert_allclose(strain, np.broadcast_to(expected, strain.shape), atol=1e-14)
            np.testing.assert_allclose(q['shape'].sum(axis=1), 1., atol=1e-14)
            integrated += np.einsum('qa,ab,qb,q->', strain, d, strain, q['weight'])
        self.assertAlmostEqual(integrated, expected @ d @ expected, places=10)
        a, b, _, meta = s.assemble(mesh, case=case)
        self.assertAlmostEqual(meta['volume_m3'], 1., places=13)
        self.assertLess(np.linalg.norm((a-a.T).data)/np.linalg.norm(a.data), 1e-14)
        self.assertTrue(np.all(np.linalg.eigvalsh(a.toarray()) > 0))
        self.assertGreater(np.linalg.norm(b), 0.)

    def test_actual_rust_elastic_geostatic_displacements_match(self):
        mesh, ground, case = small()
        binary = s.ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        self.assertTrue(binary.is_file(), f'release build required: {binary}')
        a, b, free, _ = s.assemble(mesh, case=case)
        x = splu(a.tocsc()).solve(b)
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            for name, value in [('mesh', mesh), ('ground', ground), ('case', case)]:
                (folder/f'{name}.json').write_text(json.dumps(value), encoding='utf-8')
            done = subprocess.run([str(binary), 'solve-stages', *(str(folder/f'{name}.json') for name in ('mesh', 'ground', 'case', 'result'))],
                                  capture_output=True, text=True, timeout=60)
            self.assertEqual(done.returncode, 0, done.stderr)
            result = s.read(folder/'result.json')
        self.assertTrue(result['complete'])
        self.assertTrue(result['geostatic']['complete'])
        actual = np.asarray(result['geostatic']['last_committed']['displacement_m']).ravel()[free]
        self.assertLess(np.linalg.norm(x-actual)/np.linalg.norm(x), 1e-8)
        self.assertLess(np.linalg.norm(a@actual-b)/np.linalg.norm(b), 1e-8)

    def test_full_xy1_geometry_and_gravity_inventory(self):
        mesh = s.griffiths()
        self.assertEqual(len(mesh['tetrahedra10']), 1200)
        a, b, free, meta = s.assemble(mesh)
        self.assertAlmostEqual(meta['volume_m3'], 220., places=9)
        self.assertEqual(meta['integration_points'], 32400)
        self.assertTrue(np.all(free % 3 != 1))
        self.assertEqual(a.shape, (len(free), len(free)))
        self.assertFalse(meta['plastic_failure_tangent_reproduced'])
        self.assertTrue(np.isfinite(b).all())

    def test_imported_nonsymmetric_gp_tensor_is_not_symmetrized(self):
        mesh, _, case = small()
        count = len(mesh['tetrahedra10'])
        tensor = np.broadcast_to(s.fem.elastic_tensor(100000., .3), (count, 27, 6, 6)).copy()
        tensor[:, :, 0, 2] *= .5
        a, _, _, meta = s.assemble(mesh, case=case, tangents=tensor)
        self.assertGreater(np.linalg.norm((a-a.T).data), 1.)
        self.assertEqual(meta['tangent_kind'], 'imported_algorithmic_GP_snapshot')
        with self.assertRaises(ValueError):
            s.assemble(mesh, case=case, tangents=tensor[:, :26])

    def test_noncartesian_nonzero_and_coupled_supports_refused(self):
        mesh, _, original = small()
        for change in ('direction', 'value', 'coupling', 'anchor'):
            case = copy.deepcopy(original)
            if change == 'direction': case['boundary']['constraints'][0]['direction'] = [1., 1., 0.]
            if change == 'value': case['boundary']['constraints'][0]['value_m'] = .1
            if change == 'coupling': case['ties'] = ['unsupported']
            if change == 'anchor': case['anchor_system'] = {'unsupported': True}
            with self.subTest(change=change), self.assertRaises(ValueError):
                s.assemble(mesh, case=case)

    def test_pc_preparation_creates_only_requested_level_without_solving(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'xy1'
            run = subprocess.run([sys.executable, '-m', 'tools.py_reference.griffiths_solver.prepare_pc',
                                  '--xy', '1', '--output', str(output)], cwd=s.ROOT,
                                 capture_output=True, text=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stderr)
            report = json.loads(run.stdout)
            self.assertFalse(report['solved'])
            self.assertFalse(report['p6_complete'])
            self.assertEqual(report['elements'], 1200)
            self.assertEqual(set(p.name for p in output.iterdir()), {'mesh.json', 'ground.json', 'case.json'})
            case = s.read(output/'case.json')
            self.assertEqual(case['settings']['relative_force_tolerance'], 1e-8)
            self.assertEqual(case['settings']['relative_linear_tolerance'], 1e-10)

    def test_xy2_and_invalid_material_refused_before_assembly(self):
        mesh, _, case = small()
        with self.assertRaises(ValueError): s.griffiths(2)
        mesh['tetrahedra10'] *= 51
        with self.assertRaises(ValueError): s.assemble(mesh, case=case)
        mesh, _, case = small()
        with self.assertRaises(ValueError): s.assemble(mesh, poisson=.5, case=case)


class LinearSolvers(unittest.TestCase):
    def setUp(self):
        self.a = csr_matrix([[4., 1., 0., 0.], [-.5, 5., 2., 0.], [0., 1., 4., 1.], [1., 0., -.5, 3.]])
        self.b = np.array([1., -2., 3., .5])
        self.free = np.array([0, 2, 3, 5])

    def test_all_methods_solve_a_nonsymmetric_system_with_original_gates(self):
        rows = l.compare(self.a, self.b, self.free)
        self.assertEqual(len(rows), len(l.METHODS))
        for row in rows:
            self.assertTrue(row['accepted'], row)
            self.assertLessEqual(row['true_linear_relative_residual'], 1e-10)
            self.assertLessEqual(row['force_normalized_linear_residual'], 1e-8)
            self.assertLessEqual(row['iterations'], 3000)

    def test_ilu0_no_fill_factor_has_expected_incomplete_product(self):
        apply, nnz = l.ilu0(self.a)
        self.assertLessEqual(nnz, self.a.nnz+self.a.shape[0])
        x = apply(self.b)
        self.assertTrue(np.isfinite(x).all())
        self.assertGreater(np.linalg.norm(self.a@x-self.b), .01)
        triangular = csr_matrix([[2., 1., 0.], [0., 3., 1.], [0., 0., 4.]])
        exact, _ = l.ilu0(triangular)
        np.testing.assert_allclose(triangular@exact(np.ones(3)), np.ones(3), atol=1e-14)

    def test_structural_zero_fill_slots_survive_row_scaling(self):
        # A full element clique includes zero-valued entries. Dropping them
        # before no-fill factorization changes the preconditioner itself.
        a = csr_matrix((self.a.toarray().ravel(), np.tile(np.arange(4), 4), np.arange(0, 17, 4)), shape=(4, 4))
        self.assertEqual(a.nnz, 16)
        row = l.compare(a, self.b, self.free, methods=['ilu0'])[0]
        self.assertTrue(row['accepted'], row)
        self.assertEqual(row['iterations'], 1)

    def test_zero_pivot_is_reported_without_shift_or_fallback(self):
        a = csr_matrix([[0., 1.], [1., 1.]])
        row = l.compare(a, np.ones(2), np.array([0, 2]), methods=['ilu0'])[0]
        self.assertEqual(row['status'], 'error')
        self.assertFalse(row['accepted'])
        self.assertIn('pivot', row['error'])

    def test_iteration_limit_counts_inner_iterations_and_is_not_success(self):
        row = l.compare(self.a, self.b, self.free, methods=['none'], maximum_iterations=3, restart=2)[0]
        self.assertEqual(row['iterations'], 3)
        self.assertEqual(row['solver_info'], 3)
        self.assertEqual(row['status'], 'iteration_limit')
        self.assertFalse(row['accepted'])

    def test_scaled_system_is_checked_on_original_matrix(self):
        a = diags([1e-5, 1., 1e3, 1e6]) @ self.a
        b = a @ np.array([1., 2., -1., 3.])
        row = l.compare(a, b, self.free, methods=['ilut'])[0]
        self.assertTrue(row['accepted'], row)
        self.assertLessEqual(row['true_linear_relative_residual'], 1e-10)

    def test_force_gate_keeps_max_F_one_and_is_distinct_from_linear_gate(self):
        r = l.original_residuals(csr_matrix([[1.]]), np.array([1e-12]), np.array([0.]), 0.)
        self.assertEqual(r['force_normalized_linear_residual'], 1e-12)
        self.assertEqual(r['true_linear_relative_residual'], 1.)
        self.assertEqual(l.FORCE_LIMIT, 1e-8)
        self.assertEqual(l.LINEAR_LIMIT, 1e-10)

    def test_zero_rhs_needs_no_factorization(self):
        for row in l.compare(self.a, np.zeros(4), self.free):
            self.assertTrue(row['accepted'], row)
            self.assertEqual(row['iterations'], 0)
            self.assertEqual(row['true_linear_relative_residual'], 0.)

    def test_condition_estimate_and_spd_extremes_on_known_spectrum(self):
        report = l.estimates(diags([1., 2., 4., 8., 16.]))
        self.assertAlmostEqual(report['condition_1_estimate'], 16., places=10)
        self.assertAlmostEqual(report['sampled_eigenvalue_magnitude_ratio'], 16., places=10)
        self.assertEqual(report['relative_asymmetry'], 0.)

    def test_singular_direct_and_estimate_remain_failures(self):
        a = csr_matrix([[1., 1.], [1., 1.]])
        row = l.compare(a, np.ones(2), np.array([0, 2]), methods=['direct'])[0]
        self.assertFalse(row['accepted'])
        self.assertEqual(row['status'], 'error')
        self.assertIn('estimate_error', l.estimates(a))

    def test_frozen_roundtrip_and_tampering_refusal(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)/'system'
            s.save_system(folder, self.a, self.b, self.free, {'tangent_kind': 'test'})
            a, b, free, meta = s.load_system(folder)
            np.testing.assert_array_equal(a.toarray(), self.a.toarray())
            np.testing.assert_array_equal(b, self.b)
            np.testing.assert_array_equal(free, self.free)
            self.assertEqual(meta['tangent_kind'], 'test')
            with (folder/'vectors.npz').open('ab') as stream: stream.write(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed'): s.load_system(folder)

    def test_cli_records_failure_returns_nonzero_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            s.save_system(folder/'frozen', self.a, self.b, self.free, {'tangent_kind': 'test'})
            command = [sys.executable, '-m', 'tools.py_reference.griffiths_solver.compare',
                       '--system', str(folder/'frozen'), '--output', str(folder/'result'),
                       '--methods', 'none', '--maximum-iterations', '1']
            run = subprocess.run(command, cwd=s.ROOT, capture_output=True, text=True, timeout=30)
            self.assertEqual(run.returncode, 1, run.stderr)
            report = s.read(folder/'result/comparison.json')
            self.assertFalse(report['p6_complete'])
            self.assertFalse(report['solvers'][0]['accepted'])
            again = subprocess.run(command, cwd=s.ROOT, capture_output=True, text=True, timeout=30)
            self.assertEqual(again.returncode, 2)

    def test_invalid_budget_nonfinite_system_and_unknown_method_rejected(self):
        with self.assertRaises(ValueError): l.compare(self.a, self.b, self.free, maximum_iterations=3001)
        with self.assertRaises(ValueError): l.compare(self.a, self.b*np.nan, self.free)
        with self.assertRaises(ValueError): l.compare(self.a, self.b, self.free, methods=['CG'])


if __name__ == '__main__':
    unittest.main()
