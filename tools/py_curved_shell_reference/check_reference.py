"""Explicit standalone tests. Filename avoids existing test_*.py collectors."""
import contextlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

import reference as ref
import compare as cmp


CASES=ref.catalog()['benchmarks']


def flat_patch():
    x=np.array([[0.,0.,0.],[3.,0.,0.],[3.,2.,0.],[0.,2.,0.]])
    d=np.tile([0.,0.,1.],(4,1))
    b=np.tile(np.array([[1.,0.],[0.,1.],[0.,0.]]),(4,1,1))
    return x,d,b


class ElementChecks(unittest.TestCase):
    def test_shape_partition_and_derivatives(self):
        for r,s in ((0,0),(-.73,.42),(1,-1)):
            N,dN=ref.shape(r,s)
            self.assertAlmostEqual(sum(N),1.)
            np.testing.assert_allclose(dN.sum(axis=0),0.,atol=1e-15)

    def test_exact_flat_membrane_energy(self):
        x,d,b=flat_patch();E,nu,t=2.1e7,.25,.07
        K,_=ref.element(x,d,b,E,nu,t)
        q=np.zeros((4,5));eps=np.array([.002,-.0003,.0007])
        q[:,0]=eps[0]*x[:,0]+eps[2]*x[:,1]/2
        q[:,1]=eps[1]*x[:,1]+eps[2]*x[:,0]/2
        C=E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
        exact=.5*6*t*(eps@C@eps)
        self.assertAlmostEqual(float(q.ravel()@K@q.ravel()/2)/exact,1.,places=12)

    def test_flat_cylindrical_bending_energy(self):
        x,d,b=flat_patch();E,nu,t=2.1e7,.25,.07;kappa=.003
        K,_=ref.element(x,d,b,E,nu,t)
        q=np.zeros((4,5));q[:,2]=kappa*x[:,0]**2/2;q[:,3]=-kappa*x[:,0]
        exact=.5*6*E*t**3/(12*(1-nu**2))*kappa**2
        self.assertAlmostEqual(float(q.ravel()@K@q.ravel()/2)/exact,1.,places=9)

    def test_flat_rigid_modes_and_no_hourglass(self):
        x,d,b=flat_patch();K,_=ref.element(x,d,b,1.,.3,.1)
        modes=ref.rigid_modes(x,d,b)
        np.testing.assert_allclose(K@modes,0.,atol=2e-16)
        eig=np.linalg.eigvalsh(K)
        self.assertEqual(sum(abs(eig)<1e-12),6)
        self.assertGreater(eig[6],1e-8)

    def test_curved_rigid_modes_in_cylinder(self):
        self.assert_rigid(CASES[1])

    def test_curved_rigid_modes_in_sphere(self):
        self.assert_rigid(CASES[2])

    def assert_rigid(self,case):
        K,_,_,x,d,b=ref.assemble(case,3)
        R=ref.rigid_modes(x,d,b)
        self.assertLess(np.linalg.norm(K@R)/(np.linalg.norm(K.data)*np.linalg.norm(R)),1e-14)

    def test_roof_free_and_supported_spectrum(self):
        self.assert_audit(CASES[0])

    def test_cylinder_free_and_supported_spectrum(self):
        self.assert_audit(CASES[1])

    def test_sphere_free_and_supported_spectrum(self):
        self.assert_audit(CASES[2])

    def assert_audit(self,case):
        for n in (2,4):
            audit=ref.eigen_audit(case,n)
            self.assertTrue(audit['passed'])
            self.assertEqual(audit['detected_zero_modes'],6)
            self.assertEqual(audit['extra_zero_modes'],0)
            self.assertGreater(min(audit['supported_eigenvalues']),1e-9)

    def test_roof_total_surface_load(self):
        case=CASES[0];_,f,_,_,_,_=ref.assemble(case,8)
        analytic=90*25*25*math.radians(40)
        # Bilinear cylindrical geometry has the known chord/arc area ratio.
        ratio=math.sin(math.radians(40)/16)/(math.radians(40)/16)
        self.assertAlmostEqual(-sum(f[2::5])/analytic,ratio,places=12)

    def test_point_load_symmetry_factors(self):
        _,f,_,_,_,_=ref.assemble(CASES[1],2)
        self.assertEqual(f[2],-.25)
        self.assertEqual(np.count_nonzero(f),1)
        _,f,_,_,_,_=ref.assemble(CASES[2],2)
        self.assertEqual(f[0],1.)
        self.assertEqual(f[11],-1.)
        self.assertEqual(np.count_nonzero(f),2)

    def test_equilibrium_work_and_finite_displacements(self):
        for case in CASES:
            row,_,_=ref.solve(case,4)
            self.assertLess(row['equation_residual'],1e-7)
            self.assertAlmostEqual(2*row['strain_energy']/row['external_work'],1.,places=8)
            self.assertGreater(row['displacement']*case['reference']['signed_value'],0.)

    def test_catalog_hole_and_fixed_targets(self):
        self.assertEqual([c['reference']['signed_value'] for c in CASES],[-.3024,-1.8248e-5,.094])
        self.assertEqual(CASES[2]['geometry']['hole_polar_angle_deg'],18)
        for case in CASES:
            for source in case['input_source_ids']+case['reference']['source_ids']:
                self.assertIn(source,ref.catalog()['sources'])
        self.assertIsNone(ref.catalog()['acceptance']['relative_displacement_tolerance'])

    def test_mesh_rejects_invalid_levels(self):
        for n in (0,-1,1.5):
            with self.assertRaises(ValueError):
                ref.mesh(CASES[0],n)

    def test_refuse_large_mesh_without_pc_flag(self):
        with tempfile.TemporaryDirectory() as folder:
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as err:
                ref.main(['--levels','32','64','--out',str(Path(folder)/'out')])
            self.assertEqual(err.exception.code,2)

    def test_previous_evidence_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'results.json';path.write_text('old evidence')
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as err:
                ref.main(['--levels','4','8','16','--out',folder])
            self.assertEqual(err.exception.code,2)
            self.assertEqual(path.read_text(),'old evidence')

    def test_reference_has_no_production_import_or_process(self):
        # A physical solve still runs when process creation is forbidden.
        with mock.patch.object(subprocess,'Popen',side_effect=AssertionError('process forbidden')):
            row,_,_=ref.solve(CASES[2],2)
        self.assertTrue(math.isfinite(row['displacement']))
        self.assertNotIn('mcc_rate_fem',sys.modules)


class ConvergenceChecks(unittest.TestCase):
    def test_orders_for_known_second_order_sequence(self):
        rows=[dict(n=n,displacement=1-1/n**2,relative_error=1/n**2) for n in (4,8,16)]
        ref.convergence(rows)
        self.assertAlmostEqual(rows[-1]['reference_error_order'],2.)
        self.assertAlmostEqual(rows[-1]['difference_order'],2.)
        self.assertAlmostEqual(rows[-1]['richardson_limit'],1.)

    def test_negative_order_is_retained(self):
        rows=[dict(n=n,displacement=v,relative_error=e) for n,v,e in ((4,.09,.1),(8,.095,.01),(16,.097,.03))]
        ref.convergence(rows)
        self.assertLess(rows[-1]['reference_error_order'],0.)

    def test_oscillatory_or_nondecreasing_differences_no_extrapolation(self):
        for values in ([1.,2.,1.5],[1.,2.,4.],[1.,1.,1.]):
            rows=[dict(n=n,displacement=v,relative_error=.1) for n,v in zip((4,8,16),values)]
            ref.convergence(rows)
            self.assertIsNone(rows[-1]['difference_order'])
            self.assertIsNone(rows[-1]['richardson_limit'])


class ComparisonChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)/'nodes.json'
        self.row=dict(node_id=71,x=10,y=0,z=0,ux=.094,uy=0,uz=0)

    def read(self,rows):
        self.path.write_text(json.dumps(rows),encoding='utf-8')
        return cmp.read_nodes(self.path,'hemisphere_hole')

    def test_json_and_undefined_is_never_pass(self):
        result=cmp.compare(CASES[2],self.read([self.row]),ref.catalog()['acceptance'])
        self.assertEqual(result['relative_difference'],0.)
        self.assertEqual(result['status'],'未定義')
        self.assertIsNone(result['passed'])
        self.assertFalse(result['B_01_passed'])

    def test_csv_roundtrip_physical_probe_without_node_order(self):
        p=Path(self.tmp.name)/'nodes.csv'
        ref.write_nodes(p,np.array([[0.,0.,0.],[10.,0.,0.]]),np.array([[0.,0.,0.],[.093,0.,0.]]))
        result=cmp.compare(CASES[2],cmp.read_nodes(p,'hemisphere_hole'),ref.catalog()['acceptance'])
        self.assertAlmostEqual(result['signed_difference'],-.001)

    def test_wrong_sign_is_not_hidden_by_absolute_displacement(self):
        self.row['ux']=-.094
        result=cmp.compare(CASES[2],self.read([self.row]),ref.catalog()['acceptance'])
        self.assertEqual(result['relative_difference'],2.)

    def test_reject_nonfinite_missing_duplicate_and_boolean(self):
        for key,value in [('ux',float('nan')),('z',float('inf')),('ux',True)]:
            with self.assertRaises(ValueError):self.read([{**self.row,key:value}])
        with self.assertRaises(ValueError):self.read([{'x':0}])
        with self.assertRaises(ValueError):self.read([self.row,self.row])
        with self.assertRaises(ValueError):self.read([self.row,{**self.row,'x':11}])

    def test_reject_missing_or_ambiguous_probe(self):
        for xs in ([10.1],[10,10+1e-8]):
            data=self.read([{**self.row,'node_id':i,'x':x} for i,x in enumerate(xs)])
            with self.assertRaises(ValueError):cmp.compare(CASES[2],data,ref.catalog()['acceptance'])

    def test_reject_wrong_benchmark_metadata(self):
        with self.assertRaises(ValueError):self.read({'benchmark_id':'scordelis_lo','nodes':[self.row]})

    def test_cli_undefined_and_invalid_exit_codes(self):
        self.read([self.row])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(cmp.main(['--benchmark','hemisphere_hole','--input',str(self.path)]),2)
            self.assertEqual(cmp.main(['--benchmark','hemisphere_hole','--input',str(self.path)+'missing']),3)


if __name__=='__main__':
    unittest.main(verbosity=2)
