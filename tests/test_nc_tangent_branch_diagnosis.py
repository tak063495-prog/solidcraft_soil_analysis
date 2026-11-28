"""Refuse corrupted inputs and out-of-bound NC diagnostic paths."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-10-nc-tangent-branch-diagnosis'
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location('nc_tangent_diagnose_test', HERE/'diagnose.py')
diagnose = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnose)
spec = importlib.util.spec_from_file_location('nc_bounded_arc_test', HERE/'arc.py')
arc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(arc)


class TangentBranchDiagnostics(unittest.TestCase):
    def test_actual_nonsymmetric_spectrum_and_both_singular_equations(self):
        matrix = np.array([[2., 7.], [0., .01]])
        result = diagnose.spectrum(matrix)
        actual = sorted(x['real'] for x in result['eigenvalues_actual_matrix'])
        np.testing.assert_allclose(actual, [.01, 2.])
        self.assertLess(result['symmetric_part_eigenvalues'][0], 0.)
        self.assertGreater(result['nonsymmetry_relative'], 1.)
        self.assertLess(result['left_singular_equation_residual'], 1e-14)
        self.assertLess(result['right_singular_equation_residual'], 1e-14)
        self.assertFalse(result['static_limit_certified'])

    def test_nonfinite_or_nonsquare_tangents_cannot_be_spectra(self):
        for matrix in ([[float('nan')]], [[1., 2.]], [[float('inf')]]):
            with self.assertRaisesRegex(ValueError, 'finite square'):
                diagnose.spectrum(matrix)

    def test_late_arc_refuses_extrapolation_and_extra_holding(self):
        model = SimpleNamespace(final=np.array([2., -3., 0.]), holding=np.array([1., -2., 0.]))
        model.load = lambda s: model.final+(2-2*s)*model.holding
        np.testing.assert_array_equal(arc.component_bounds(model, 1.), np.zeros(3))
        np.testing.assert_array_equal(arc.component_bounds(model, .5), model.holding)
        for s in (.499, 1.000001):
            with self.assertRaisesRegex(ValueError, 'outside'):
                arc.component_bounds(model, s)
        model.load = lambda s: model.final+2*model.holding
        with self.assertRaisesRegex(ValueError, 'bound exceeded'):
            arc.component_bounds(model, .9)

    def test_changed_retained_evidence_is_rejected_before_diagnosis(self):
        with tempfile.TemporaryDirectory() as folder:
            evidence = Path(folder)/'evidence'
            evidence.mkdir()
            (evidence/'raw.json').write_bytes(b'changed')
            expected = dict(sha256=hashlib.sha256(b'original').hexdigest(), size=8)
            (evidence/'extraction.json').write_text(json.dumps(dict(retained_entries={'raw.json': expected})))
            with self.assertRaisesRegex(ValueError, 'diagnostic input changed'):
                diagnose.run(evidence, Path(folder)/'result')


if __name__ == '__main__':
    unittest.main()
