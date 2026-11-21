"""Input fidelity, independent coupled rate balance, and four real CLI restarts."""
import os
from pathlib import Path
import sys
import tempfile
import unittest

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
import probe_low_preload_nc as low
from low_preload_nc_reference import ClosedStickSystem
from mcc_rate_fem import seed
from check_low_preload_nc_restarts import check


class LowNCFollowup(unittest.TestCase):
    def test_refinement_changes_only_declared_existing_controls(self):
        mesh, ground, case = low.expected_inputs()
        refined = low.expected_inputs(.03125, 12, 1/4096, True)
        self.assertEqual(mesh, refined[0]); self.assertEqual(ground, refined[1])
        actual = refined[2]
        actual['stages'][2]['initial_increment'] = .5
        actual['settings']['maximum_newton_iterations'] = 60
        actual['settings']['minimum_increment'] = case['settings']['minimum_increment']
        actual['settings']['linear_solver'] = case['settings']['linear_solver']
        actual['settings']['maximum_free_dofs'] = case['settings']['maximum_free_dofs']
        self.assertEqual(actual, case)
        for value in (0., float('nan'), .75, 1e-10):
            with self.assertRaises(ValueError): low.expected_inputs(excavation_increment=value)
        for value in (0, 61, 1.5):
            with self.assertRaises(ValueError): low.expected_inputs(maximum_newton_iterations=value)
        for value in (0., float('inf'), .5):
            with self.assertRaises(ValueError): low.expected_inputs(minimum_increment=value)

    def test_refined_input_still_refuses_100_kpa_and_unselected_controls(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            for name, value in zip(low.NAMES, low.expected_inputs(.03125, 12)):
                low.write(folder/(name+'.json'), value)
            low.check_inputs(folder, excavation_increment=.03125, maximum_newton_iterations=12)
            with self.assertRaises(ValueError): low.check_inputs(folder)
            case = low.read(folder/'stages.json')
            case['boundary']['surface_loads'][0]['traction_kpa'][2] = -100.
            (folder/'stages.json').unlink(); low.write(folder/'stages.json', case)
            with self.assertRaises(ValueError):
                low.check_inputs(folder, excavation_increment=.03125, maximum_newton_iterations=12)

    def test_independent_full_3d_elastoplastic_rate_and_joint_balance(self):
        mesh, ground, case = low.expected_inputs()
        model = ClosedStickSystem(mesh, ground, [True]*12, case)
        state = seed(mesh, ground).reshape(-1, 14)
        force = np.random.default_rng(143).normal(size=len(mesh['nodes'])*3)
        for mask in (np.zeros(len(state), dtype=bool), np.ones(len(state), dtype=bool)):
            du, dy, _, _ = model.rate(state, force, mask)
            residual = model.reduce(model.total_force(dy, du)-force)
            self.assertLess(np.linalg.norm(residual)/np.linalg.norm(force), 1e-10)
            self.assertGreater(np.linalg.norm(du.reshape(-1, 3), axis=0).min(), 0.)
            if mask.all(): self.assertGreater(np.linalg.norm(dy[:, 8:]), 0.)

    def test_reference_refuses_opening_and_sliding_instead_of_extending_law(self):
        mesh, ground, case = low.expected_inputs()
        model = ClosedStickSystem(mesh, ground, [True]*12, case)
        nodes = np.unique([n for p in model.joint.points for n in p['nodes_b']])
        u = np.zeros_like(mesh['nodes'], dtype=float)
        u[nodes, 2] = .001
        with self.assertRaisesRegex(ValueError, 'closed-contact'): model.contacts(u.ravel())
        u[nodes, 2] = -.00001; u[nodes, 0] = .001
        with self.assertRaisesRegex(ValueError, 'sticking'): model.contacts(u.ravel())

    def test_real_separate_process_restarts_inside_every_phase(self):
        solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        with tempfile.TemporaryDirectory() as temp:
            report = check(solver, Path(temp)/'restarts')
            self.assertFalse(report['passed'])
            self.assertTrue(report['partial_restart_comparisons_exact'])
            self.assertEqual(len(report['comparisons']), 4)


if __name__ == '__main__':
    unittest.main()
