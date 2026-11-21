"""Check the continuous manufactured load independently of the discrete stiffness."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / 'interim/2026-11-21-impl-t5-curved-tri6-shell'


def load(name):
    spec = importlib.util.spec_from_file_location(name, SOURCE / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CurvedManufactured(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import unittest.mock
        with unittest.mock.patch.dict(sys.modules, standard=load('standard')):
            cls.m = load('manufactured')

    def test_planar_compatible_transverse_shear(self):
        p = np.array([[.22,.36],[.43,.13],[.19,.61]])
        u,du,theta = self.m.exact(p,0.,0.)
        eps,stress = self.m.continuum(p,0.,0.,0.)
        np.testing.assert_allclose(2*eps[:,0,2],du[:,2,0]+theta[:,1],rtol=1e-12)
        np.testing.assert_allclose(2*eps[:,1,2],du[:,2,1]-theta[:,0],rtol=1e-12)
        np.testing.assert_allclose(eps[:,0,0],du[:,0,0],rtol=1e-12)
        np.testing.assert_allclose(stress[:,2,2],0.,atol=1e-14)
        self.assertEqual(u.shape,(3,3))

    def test_body_and_boundary_balance_for_constant_translations(self):
        a,b=.18,-.12
        case,_=self.m.mesh_case(2,a,b,False)
        load=np.array(case['steps'][1]['nodal_load']).reshape(-1,6)
        # All natural edges were applied except r=0. Its outward conormal is -e_r.
        p=np.column_stack((np.zeros_like(self.m.RULE),self.m.RULE))
        omitted=-(self.m.WEIGHT[:,None]*self.m.fluxes(p,a,b)[0][:,:3,0]).sum(axis=0)
        np.testing.assert_allclose(load[:,:3].sum(axis=0)+omitted,0.,atol=2e-9)
