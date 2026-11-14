"""Actual 5 kPa NC CLI checkpoints and refusal to certify missing evidence."""
import copy
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import probe_low_preload_nc as low


class LowPreloadNC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)
        cls.solver = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
        cls.first = low.probe(cls.solver, cls.folder/'first', 60., 15., stop_after_accepted=9)
        cls.second = low.probe(cls.solver, cls.folder/'second', 60., 15.,
                               resume_from=cls.folder/'first/job', stop_after_accepted=1)
        cls.a = low.read(cls.folder/'first/observed-run.json')
        cls.b = low.read(cls.folder/'second/observed-run.json')
        count = low.accepted_count(cls.b)
        inputs = [cls.folder/'first/input'/(n+'.json') for n in low.NAMES]
        completed = subprocess.run([str(cls.solver), 'checkpoint-stages', *map(str, inputs),
                                    str(count), str(cls.folder/'direct.json')],
                                   capture_output=True, timeout=90, check=False)
        if completed.returncode:
            raise AssertionError(completed.stderr.decode('utf-8', errors='replace'))
        cls.direct = low.read(cls.folder/'direct.json')['run']

    def test_actual_nc_plastic_checkpoint_and_separate_process_prefix_match(self):
        self.assertGreater(low.accepted_count(self.b), low.accepted_count(self.a))
        self.assertTrue(low.compare_prefix(self.b, self.direct)['exact'])
        self.assertGreater(self.second['observation']['plastic_points'], 0)
        self.assertTrue(self.second['observation']['all_rejections_preserve_committed_state'])
        self.assertTrue(self.second['observation']['stored_stress_force_checks'])
        self.assertFalse(self.first['owned_child_terminated'])
        self.assertFalse(self.second['owned_child_terminated'])

    def test_partial_solver_and_restart_are_not_full_excavation_acceptance(self):
        for report in (self.first, self.second):
            self.assertFalse(report['passed'])
            self.assertFalse(report['solver_complete'])
            self.assertEqual(report['classification'], '未検証')
            self.assertFalse(report['p3b04_complete'])
            self.assertFalse(report['p6_complete'])
            self.assertFalse(report['observation']['independent_plastic_3d_displacement'])
            self.assertFalse(report['observation']['fully_released_and_frozen'])
            self.assertIn('独立塑性3D全場', report['missing_acceptance_evidence'])

    def test_100_kpa_preload_cannot_replace_5_kpa_original(self):
        folder = self.folder/'wrong'; folder.mkdir()
        mesh, ground, case = low.make(mcc=True, pc=100., preload_kpa=100., force=(1., .5, 0.))
        for name, value in zip(low.NAMES, (mesh, ground, case)):
            low.write(folder/(name+'.json'), value)
        with self.assertRaisesRegex(ValueError, 'input mismatch: stages'):
            low.check_inputs(folder)
        job = self.folder/'wrong-job'; (job/'inputs').mkdir(parents=True)
        for name, value in zip(('mesh', 'ground', 'case'), (mesh, ground, case)):
            low.write(job/'inputs'/(name+'.json'), value)
        with self.assertRaisesRegex(ValueError, 'input mismatch: stages'):
            low.probe(self.solver, self.folder/'not-started', resume_from=job)

    def test_corrupt_material_contact_and_displacement_cannot_pass_prefix(self):
        for kind in ('material', 'contact', 'displacement', 'attempt'):
            changed = copy.deepcopy(self.b)
            state = changed['last_committed']
            if kind == 'material':
                state['material_states'][0]['preconsolidation_pressure_kpa'] += .01
            elif kind == 'contact':
                state['joints']['points'][0]['traction_kpa'][0] += .01
            elif kind == 'displacement':
                state['displacement_m'][0][0] += .001
            else:
                changed['attempts'][0]['committed_state_unchanged_on_rejection'] = False
            with self.assertRaisesRegex(ValueError, 'prefix differs'):
                low.compare_prefix(changed, self.direct)

    def test_nonfinite_and_unbounded_runtime_budgets_are_rejected(self):
        for budget in (float('nan'), float('inf'), -1., 0., 601.):
            with self.assertRaisesRegex(ValueError, 'budget must'):
                low.probe(self.solver, self.folder/'not-created', budget)

    def test_actual_wall_time_stop_never_becomes_collapse_or_acceptance(self):
        report = low.probe(self.solver, self.folder/'timeout', .001, .001)
        self.assertEqual(report['stop_reason'], 'wall_time_budget')
        self.assertFalse(report['passed'])
        self.assertFalse(report['solver_complete'])
        self.assertEqual(report['classification'], '未検証')
        self.assertTrue(report['owned_child_terminated'])


if __name__ == '__main__':
    unittest.main()
