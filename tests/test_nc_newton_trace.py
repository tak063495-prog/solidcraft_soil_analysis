"""Trace association and nonsymmetric matrix orientation are evidence-critical."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from measure_nc_newton_trace import groups, tangent_report


class NCNewtonTrace(unittest.TestCase):
    def test_scoped_attempts_keep_rejections_and_separate_restarts(self):
        events = []
        for stage, progress, outcome in ((2, 1., 'accepted'), (3, .953125, 'rejected'),
                                         (3, .951296329498291, 'rejected')):
            events.append(dict(kind='progress', event=dict(kind='attempt_started', stage_index=stage,
                                                           progress=progress)))
            events += [dict(kind='iterate', iteration=0, projected_load=projected) for projected in (True, False)]
            events.append(dict(kind='progress', event=dict(kind=outcome, reason='unconverged')))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'trace.jsonl'
            path.write_text('\n'.join(json.dumps(e) for e in events), encoding='utf-8')
            measured = list(groups(path))
        self.assertEqual(len(measured), 1)
        self.assertEqual(measured[0]['outcome']['kind'], 'rejected')
        self.assertEqual([e['projected_load'] for e in measured[0]['events']], [True, False])

    def test_actual_nonsymmetric_matrix_columns_are_not_symmetrized(self):
        event = dict(matrix_columns=[[4., 1.], [2., 3.]], reduced_residual=[1., 2.],
                     iteration=7, projected_load=False)
        report = tangent_report(event)
        self.assertEqual(report['dimension'], 2)
        self.assertAlmostEqual(report['numpy_correction'][0], -.1)
        self.assertAlmostEqual(report['numpy_correction'][1], .7)
        self.assertLess(report['numpy_linear_relative_residual'], 1e-10)
        self.assertFalse(report['projected_load'])
        self.assertEqual(len(report['numpy_refinement_residuals']), 5)
        self.assertIn('no native acceptance', report['interpretation'])


if __name__ == '__main__':
    unittest.main()
