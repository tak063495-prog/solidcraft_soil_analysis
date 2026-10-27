"""The public runner must initialize independently and refuse another backend."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from prepare_mcc_structures import inputs
import run_projected_reference as runner
from assess_projected_mcc_job import reference_bindings


class ReferenceRunner(unittest.TestCase):
    def fixture(self, folder):
        mesh, ground, case = inputs(Path(__file__).resolve().parents[1])
        case.pop('beam_system')
        case['stages'] = [dict(kind='pressure', name='small_load', target_pressure_kpa=101.)]
        paths = [folder/(name+'.json') for name in ('mesh', 'ground', 'case', 'identity')]
        for path, value in zip(paths[:3], (mesh, ground, case), strict=True):
            path.write_text(json.dumps(value), encoding='utf-8')
        paths[3].write_text(json.dumps(dict(schema='soil-fem-mesh-identity/1',
            mesh_sha256=hashlib.sha256(paths[0].read_bytes()).hexdigest())), encoding='utf-8')
        return paths

    def test_both_integrators_create_bound_initialization_and_construction(self):
        with tempfile.TemporaryDirectory() as name, contextlib.redirect_stdout(io.StringIO()):
            folder = Path(name);paths = self.fixture(folder)
            for method in ('midpoint', 'rk4'):
                output = folder/method
                result = runner.run(*paths, output, steps=8, integration=method)
                self.assertTrue(result['complete'])
                initial = runner.read(output/'initial.json')
                cycle, _ = reference_bindings(output/'cycle.json')
                self.assertEqual(initial['load_fraction'], 1.)
                self.assertEqual(initial['reference_steps'], 8)
                self.assertEqual(len(cycle['frames']), 2)
                self.assertEqual(cycle['frames'][0]['material_state'], initial['material_state'])
                self.assertIn(str(Path(runner.__file__).resolve()), cycle['input_and_source_sha256'])
                with self.assertRaises(FileExistsError):
                    runner.run(*paths, output, steps=8, integration=method)

    def test_saved_initialization_rejects_another_integrator(self):
        with tempfile.TemporaryDirectory() as name, contextlib.redirect_stdout(io.StringIO()):
            folder = Path(name);paths = self.fixture(folder)
            runner.run(*paths, folder/'seed', steps=8, phase='initial')
            original = folder/'seed/initial.json';before = original.read_bytes()
            with self.assertRaisesRegex(ValueError, 'backend differs'):
                runner.run(*paths, folder/'wrong', steps=8, phase='cycle', integration='rk4', initial=original)
            self.assertEqual(original.read_bytes(), before)
            self.assertFalse((folder/'wrong/cycle.json').exists())
            runner.run(*paths, folder/'valid', steps=8, phase='cycle', initial=original)
            self.assertTrue((folder/'valid/cycle.json').is_file())

    def test_load_holds_preserve_all_histories_and_reaudit_original_gates(self):
        with tempfile.TemporaryDirectory() as name, contextlib.redirect_stdout(io.StringIO()):
            folder = Path(name);paths = self.fixture(folder)
            case = runner.read(paths[2])
            case['stages'] = [dict(kind='pressure', name='initial_hold', target_pressure_kpa=100.),
                case['stages'][0], dict(kind='pressure', name='loaded_hold', target_pressure_kpa=101.),
                dict(kind='boundary', name='boundary_hold', target_factor=0.)]
            paths[2].write_text(json.dumps(case), encoding='utf-8')
            for method in ('midpoint', 'rk4'):
                output = folder/method
                runner.run(*paths, output, steps=8, integration=method)
                frames = runner.read(output/'cycle.json')['frames']
                for previous, current in ((0, 1), (2, 3), (3, 4)):
                    for key in ('displacement_m', 'material_state', 'active_elements', 'pressure_dofs',
                                'volume_dofs', 'compatibility_offsets'):
                        self.assertEqual(frames[previous][key], frames[current][key], key)
                    self.assertEqual(frames[current]['diagnostics']['accepted_steps'], 0)
                    self.assertLess(frames[current]['relative_work_conjugate_force_residual'], 1e-8)
                    self.assertLess(frames[current]['relative_frozen_volume_constraint_residual'], 1e-8)

    def test_rk4_separate_cli_processes_share_the_same_bound_backend(self):
        with tempfile.TemporaryDirectory() as name:
            folder = Path(name);paths = self.fixture(folder)
            for phase in ('initial', 'cycle'):
                args = [sys.executable, str(Path(runner.__file__).resolve()), *map(str, paths),
                        str(folder/phase), '--phase', phase, '--steps', '8', '--integration', 'rk4']
                if phase == 'cycle':
                    args += ['--initial', str(folder/'initial/initial.json')]
                process = subprocess.run(args, capture_output=True, text=True, encoding='utf-8')
                self.assertEqual(process.returncode, 0, process.stdout+process.stderr)
            initial = runner.read(folder/'initial/initial.json')
            cycle, _ = reference_bindings(folder/'cycle/cycle.json')
            self.assertEqual(initial['reference_backend'], cycle['reference_backend'])
            self.assertEqual(cycle['reference_backend'], '__main__.Rk4MatrixFreeProjectedSystem')


if __name__ == '__main__':
    unittest.main()
