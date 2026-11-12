"""Original-system residual and retained failures of the input-only arc."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'interim/2026-11-09-nc-pseudo-arclength'
EVIDENCE = ROOT/'interim/2026-12-05-nc-independent-arc-residual'
sys.path[:0] = [str(ROOT/'tools'), str(PRIOR)]
import independent_arc as arc
from independent_map_arc import read_arc_seed


class IndependentArcResidual(unittest.TestCase):
    def fixture(self, directory):
        archive = EVIDENCE/'reference400-fixture.zip'
        manifest = json.loads((EVIDENCE/'reference400-fixture-sha256.json').read_text())
        self.assertEqual(hashlib.sha256(archive.read_bytes()).hexdigest(), manifest[archive.name])
        ref = directory/'reference'
        with zipfile.ZipFile(archive) as z:
            z.extractall(ref)
        return ref

    def matrix(self):
        data = json.loads((EVIDENCE/'bordered-regression.json').read_text())
        return np.array(data['matrix']), np.array(data['rhs'])

    def test_actual_bordered_matrix_improves_without_changing_norm_or_inputs(self):
        a, b = self.matrix(); old_a, old_b = a.copy(), b.copy()
        x, diagnostics = arc.solve_bordered(a, b)
        measured = float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b), 1e-30))
        self.assertEqual(measured, diagnostics['true_relative_residual'])
        self.assertLessEqual(measured, 1e-10)
        self.assertLessEqual(measured, diagnostics['unscaled_direct_relative_residual'])
        self.assertEqual(len(diagnostics['equilibrated_and_refined_relative_residuals']), 5)
        np.testing.assert_array_equal(a, old_a); np.testing.assert_array_equal(b, old_b)

    def test_decimal_corrections_when_longdouble_has_no_extra_precision(self):
        a, b = self.matrix()
        with patch.object(arc.np, 'longdouble', np.float64):
            x, diagnostics = arc.solve_bordered(a, b)
        self.assertEqual(diagnostics['correction_precision'], 'decimal50')
        self.assertLessEqual(float(np.linalg.norm(a@x-b)/np.linalg.norm(b)), 1e-10)
        # The correction residual sees cancellation that float64 A@x loses.
        with patch.object(arc.np, 'longdouble', np.float64):
            r = arc.correction_residual(np.array([[1e16, 1., -1e16]]), np.ones(3), np.zeros(1))
        self.assertEqual(r.tolist(), [-1.])

    def test_above_tolerance_is_saved_and_fails_with_zero_samples(self):
        def bad_solve(a, b):
            x = np.zeros(len(b))
            return x, dict(true_relative_residual=float(np.linalg.norm(a@x-b)/np.linalg.norm(b)))
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); ref = self.fixture(root); out = root/'failed'
            with patch.object(arc, 'solve_bordered', bad_solve):
                report = arc.run(ref, out, 2e-5, 1)
            self.assertFalse(report['passed']); self.assertEqual(report['equilibrium_samples'], 0)
            self.assertIn('existing 1e-10', report['stop_reason'])
            self.assertEqual(json.loads((out/'samples.json').read_text()), [])
            systems = json.loads((out/'linear-systems.json').read_text())
            self.assertEqual(len(systems), 1)
            a, b, x = (np.array(systems[0][key]) for key in ('matrix', 'rhs', 'solution'))
            self.assertEqual(float(np.linalg.norm(a@x-b)/np.linalg.norm(b)), 1.)
            self.assertTrue((out/'failure.log').is_file()); self.assertTrue((out/'segments.json').is_file())
            with self.assertRaisesRegex(ValueError, 'equilibrium_samples=0'):
                read_arc_seed(out)

    def test_later_failure_keeps_accepted_sample_and_exact_solved_vectors(self):
        integrate = arc.advance_rk4; calls = 0
        def fail_second(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise ValueError('intentional second-step failure')
            return integrate(*args, **kwargs)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); ref = self.fixture(root); out = root/'partial'
            with patch.object(arc, 'advance_rk4', fail_second):
                report = arc.run(ref, out, 2e-5, 2)
            self.assertEqual(report['equilibrium_samples'], 1); self.assertFalse(report['passed'])
            self.assertEqual(report['diagnostic_status'], 'failed')
            self.assertEqual(len(json.loads((out/'samples.json').read_text())), 1)
            self.assertTrue((out/'restart.json').is_file()); self.assertTrue((out/'failure.json').is_file())
            with self.assertRaisesRegex(ValueError, 'failed trace; retain partial samples'):
                read_arc_seed(out)
            for row in json.loads((out/'linear-systems.json').read_text()):
                a, b, x = (np.array(row[k]) for k in ('matrix', 'rhs', 'solution'))
                self.assertEqual(float(np.linalg.norm(a@x-b)/np.linalg.norm(b)), row['true_relative_residual'])
                self.assertLessEqual(row['true_relative_residual'], 1e-10)

    def test_missing_and_empty_arc_cli_fail_clearly_and_preserve_sources(self):
        scripts = [ROOT/'tools/nc_fold_path.py', PRIOR/'independent_map_arc.py',
                   ROOT/'interim/2026-11-09-nc-beyond-fold/continue_path.py']
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); source = root/'arc'; source.mkdir()
            (source/'old.log').write_text('retained failure\n')
            for mode in ('missing', 'empty'):
                if mode == 'empty':
                    (source/'samples.json').write_text('[]\n')
                    (source/'report.json').write_text('{"passed":false,"equilibrium_samples":0}\n')
                before = {p.name: p.read_bytes() for p in source.iterdir()}
                for index, script in enumerate(scripts):
                    with self.subTest(mode=mode, script=script.name):
                        out = root/f'{mode}-{index}'
                        args = [sys.executable, str(script), '--reference', str(root/'unused'), '--arc', str(source),
                                '--output', str(out), '--step-size', '0.00002', '--count', '1']
                        if script.name != 'independent_map_arc.py':
                            args += ['--map', 'independent']
                        child = subprocess.run(args, capture_output=True, text=True, cwd=ROOT)
                        self.assertEqual(child.returncode, 1, child.stderr)
                        self.assertNotIn('FileNotFoundError', child.stderr)
                        report = json.loads((out/'report.json').read_text())
                        self.assertFalse(report['passed']); self.assertFalse(report['nc01_complete'])
                        self.assertEqual(report['equilibrium_samples'], 0)
                        self.assertEqual(report['failure_stage'], 'independent_arc_seed')
                        self.assertIn('missing' if mode == 'missing' else 'empty', report['stop_reason'])
                        self.assertTrue((out/'failure.log').is_file())
                self.assertEqual(before, {p.name: p.read_bytes() for p in source.iterdir()})

    def test_empty_trace_cannot_be_reported_as_usable(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); ref = self.fixture(root); out = root/'zero'
            report = arc.run(ref, out, 2e-5, 0)
            self.assertEqual(report['diagnostic_status'], 'failed')
            self.assertFalse(report['passed']); self.assertEqual(report['equilibrium_samples'], 0)
            self.assertIn('no independent arc samples', report['stop_reason'])


if __name__ == '__main__':
    unittest.main()
