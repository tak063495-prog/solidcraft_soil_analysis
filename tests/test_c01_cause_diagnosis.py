"""Analytical checks for the isolated diagnostic, not new C-01 acceptance gates."""
import importlib.util
from pathlib import Path
import sys
import unittest
import tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-13-c01-cause-diagnosis'
sys.path.insert(0, str(HERE))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


diag = load('cause_diagnostic_tests',HERE/'c01_cause_diagnostics.py')
experiment = load('cause_experiment_tests',HERE/'c01_cause_experiment.py')
compare = load('cause_compare_tests',HERE/'c01_cause_compare.py')


class CauseDiagnosticTests(unittest.TestCase):
    def test_constant_pressure_return_and_work_balance(self):
        # kn=100, penetration=.1 -> p=10; kt=100, mu=.5 -> tau limit=5.
        start,q0,end = [0.,0.,-.1],[0.,0.,-10.],[.2,0.,-.1]
        q,dp,d,w,nd = diag.local_replay(start,q0,end,[0.,0.,1.],100.,100.,.5,1)
        np.testing.assert_allclose(q,[5.,0.,-10.],rtol=0,atol=1e-14)
        np.testing.assert_allclose(dp,[.15,0.,0.],rtol=0,atol=1e-14)
        self.assertAlmostEqual(d,.75)
        self.assertAlmostEqual(w,diag.energy(q,[0.,0.,1.],100.,100.)-.5+d+nd)
        self.assertGreaterEqual(nd,0.)

    def test_opening_endpoint_is_not_internal_dissipation_quadrature(self):
        # First half stays closed and yields; final endpoint is open with zero q.
        start,q0,end = [0.,0.,-.1],[0.,0.,-10.],[.2,0.,.1]
        q,dp,d,_,_ = diag.local_replay(start,q0,end,[0.,0.,1.],100.,100.,.5,4)
        self.assertEqual(q,[0.,0.,0.])
        self.assertGreater(d,0.)
        self.assertEqual(diag.dot(q,dp),0.)

    def test_global_vector_basis_rotation_invariance(self):
        rotation = np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
        start,q0,end,n = [0.,0.,-.1],[0.,0.,-10.],[.2,.05,-.12],[0.,0.,1.]
        a = diag.local_replay(start,q0,end,n,100.,100.,.5,4)
        b = diag.local_replay(*(rotation@x for x in (start,q0,end,n)),100.,100.,.5,4)
        np.testing.assert_allclose(rotation@a[0],b[0],rtol=0,atol=1e-14)
        np.testing.assert_allclose(rotation@a[1],b[1],rtol=0,atol=1e-14)
        np.testing.assert_allclose(a[2:],b[2:],rtol=0,atol=1e-14)

    def test_actual_coordinates_required_without_interpolation(self):
        data = ({('shear',.5):0},{'u':np.array([[1.]])})
        with self.assertRaisesRegex(ValueError,'no interpolation'):
            compare.select(data,[('shear',.25)])
        np.testing.assert_array_equal(compare.select(data,[('shear',.5)])['u'],[[1.]])

    def test_fixed_order_weighted_norm_has_analytic_value(self):
        a,b = np.array([[3.,4.],[0.,0.]]),np.array([[0.,0.],[0.,5.]])
        m = diag.metric(a,b,np.array([1.,2.]))
        self.assertEqual(m['difference_squared'],91.)
        self.assertEqual(m['reference_squared'],50.)
        self.assertAlmostEqual(m['pooled_relative_l2'],np.sqrt(91/50))
        self.assertEqual(m,diag.metric(a.copy(order='F'),b.copy(order='F'),np.array([1.,2.])))

    def test_contact_rules_preserve_area_and_quadratic_patch(self):
        old,_ = diag.dependencies()
        mesh = old.read(ROOT/'examples/joint-water-column-tet10-mesh.json')
        case = old.fem.normalize_stage_case(old.read(ROOT/'examples/stages-joint-water-cycle.json'))
        base = experiment.contact_system(mesh,case['joint_system'],4)
        for order in (4,5,6):
            system = experiment.contact_system(mesh,case['joint_system'],order)
            self.assertEqual(len(system.points),2*order*order)
            self.assertAlmostEqual(sum(system.weight),sum(base.weight),places=12)
            # Same physical quadratic field on the two duplicated sides gives zero jump.
            xyz = np.asarray(mesh['nodes'])
            u = np.column_stack((xyz[:,0]**2+xyz[:,1],xyz[:,1]*xyz[:,2],xyz[:,2]**2))
            np.testing.assert_allclose(np.einsum('qag,g->qa',system.c,u.ravel()),0.,rtol=0,atol=1e-14)
            for axis in range(3):
                moment = sum(p['weight_m2']*p['position_m'][axis]**2 for p in system.points)
                expected = sum(p['weight_m2']*p['position_m'][axis]**2 for p in base.points)
                self.assertAlmostEqual(moment,expected,places=12)

    def test_sealed_diagnostic_integrity_keeps_original_decision(self):
        evidence = load('cause_evidence_tests',HERE/'c01_cause_evidence.py')
        with tempfile.TemporaryDirectory() as tmp:
            result = evidence.replay(Path(tmp)/'evidence')
            self.assertTrue(result['measurements_exact'])
            self.assertFalse(result['numerical_solve'])
        report = diag.read(HERE/'acceptance.json')
        self.assertIs(report['passed'],False)
        self.assertEqual(report['c01_decision'],'未達')
        self.assertIs(report['original_metrics_exact'],True)

    def test_document_endings_accept_real_crlf(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text = (HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('C-01は未達、passed=false。\n'))
            crlf = text.replace('\n','\r\n').encode('utf-8')
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp)/name
                path.write_bytes(crlf)
                self.assertEqual(path.read_text(encoding='utf-8'),text)
                self.assertTrue(path.read_text(encoding='utf-8').endswith('C-01は未達、passed=false。\n'))


if __name__ == '__main__':
    unittest.main()
