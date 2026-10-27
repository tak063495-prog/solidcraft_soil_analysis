"""Audit geometry must preserve values and isolate its reusable snapshot."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import fem_reference as ref
from test_mcc_factor_reference import specimen

class ReferenceGeometry(unittest.TestCase):
    def test_rotated_curved_quadratic_geometry_is_bit_exact_to_standalone_rules(self):
        mesh,_,_=specimen()
        mesh['nodes'][4]=(np.asarray(mesh['nodes'][4])+[.012,.003,-.008]).tolist()
        for order in (3,4,5):
            mesh['quadrature_order']=order;geometry=ref.Geometry(mesh)
            for independent in (False,True):
                expected=ref.volume(mesh,0,independent,5);actual=geometry.volume(0,independent,5)
                for key in expected:np.testing.assert_array_equal(actual[key],expected[key],err_msg=key)

    def test_tet4_affine_patch_and_volume(self):
        mesh=dict(nodes=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],tetrahedra=[[0,1,2,3]])
        q=ref.Geometry(mesh).volume(0)
        gradient=np.array([[.1,.2,.3],[.4,.5,.6],[.7,.8,.9]])
        u=np.asarray(mesh['nodes'])@gradient.T+[2.,-3.,4.]
        np.testing.assert_allclose(q['b']@u.reshape(-1),[[.1,.5,.9,.6,1.4,1.]],rtol=0,atol=1e-14)
        self.assertAlmostEqual(q['weight'].sum(),1/6)

    def test_snapshot_and_returned_arrays_cannot_corrupt_later_elements(self):
        mesh,_,_=specimen();geometry=ref.Geometry(mesh);expected=ref.volume(mesh,0)
        first=geometry.volume(0)
        for key in first:first[key].fill(0)
        mesh['nodes'][1][0]+=.05;mesh['quadrature_order']=5
        again=geometry.volume(0)
        for key in expected:np.testing.assert_array_equal(again[key],expected[key])
        updated=ref.Geometry(mesh).volume(0);standalone=ref.volume(mesh,0)
        for key in standalone:np.testing.assert_array_equal(updated[key],standalone[key])
        self.assertNotEqual(len(updated['natural']),len(expected['natural']))

    def test_common_quadrature_is_generated_once_per_context_and_rule(self):
        mesh,_,_=specimen();geometry=ref.Geometry(mesh)
        with patch.object(ref,'volume_rule',wraps=ref.volume_rule) as rules,patch.object(ref,'basis',wraps=ref.basis) as bases:
            for _ in range(3):geometry.volume(0)
            self.assertEqual((rules.call_count,bases.call_count),(1,1))
            for _ in range(3):geometry.volume(0,True,5)
            self.assertEqual((rules.call_count,bases.call_count),(2,2))

if __name__=='__main__':unittest.main()
