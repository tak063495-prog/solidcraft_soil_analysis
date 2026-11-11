"""Constitutive and assembly properties of the C-03 observation reference."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
SPEC=importlib.util.spec_from_file_location('c03_plastic_reference',ROOT/'interim/2026-11-20-wet-joint-plastic-bulk-integration/reference.py')
sys.path.insert(0,str(Path(SPEC.origin).parent))
ref=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(ref)
import mc_reference as mc
import fem_reference as fem


def inputs():
    read=lambda name:json.loads((ROOT/'examples'/name).read_text(encoding='utf-8'))
    return read('joint-water-column-mesh.json'),read('ground-joint-water-mc.json'),read('stages-joint-water-cycle.json')

def material():return dict(young_kpa=100000.,poisson=.3,cohesion_kpa=10.,friction_deg=30.,dilation_deg=5.,tension_cutoff_kpa=None,flow='non_associated')
def state():return dict(effective_stress_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],plastic_strain_compression_positive=[0.]*6,void_ratio=.7,equivalent_plastic_strain=0.,plastic_work_kpa=0.,last_return='elastic')
def voigt(t,engineering=False):return np.array([t[0,0],t[1,1],t[2,2],t[0,1]*(2 if engineering else 1),t[1,2]*(2 if engineering else 1),t[2,0]*(2 if engineering else 1)])


class PlasticReference(unittest.TestCase):
    def test_elastic_limit_matches_tensor(self):
        p=material();old=state();de=np.array([1e-6,-2e-6,3e-6,4e-6,0.,0.]);new=ref.integrate(p,old,de,3)
        np.testing.assert_allclose(new['effective_stress_compression_positive_kpa'],np.asarray(old['effective_stress_compression_positive_kpa'])+fem.elastic_tensor(p['young_kpa'],p['poisson'])@de,atol=1e-10,rtol=1e-12)
        self.assertEqual(new['equivalent_plastic_strain'],0.)

    def test_plastic_nonassociated_return_is_admissible_and_dissipative(self):
        p=material();new=ref.integrate(p,state(),[0.,0.,0.,.006,.002,0.],8)
        self.assertGreater(new['equivalent_plastic_strain'],0.);self.assertGreater(new['plastic_work_kpa'],0.)
        self.assertLess(float(mc.normalized_yield(p,new['effective_stress_compression_positive_kpa'])),1e-8)
        self.assertGreater(np.linalg.norm(new['plastic_strain_compression_positive']),0.)

    def test_elastic_unload_preserves_irreversible_history(self):
        p=material();old=ref.integrate(p,state(),[0.,0.,0.,.006,.002,0.],8);new=ref.integrate(p,old,[0.,0.,0.,-1e-5,-3e-6,0.],1)
        for key in ('equivalent_plastic_strain','plastic_work_kpa','plastic_strain_compression_positive'):np.testing.assert_allclose(new[key],old[key],atol=1e-14,rtol=0)

    def test_trial_does_not_mutate_committed_history(self):
        old=state();saved=copy.deepcopy(old);ref.integrate(material(),old,[0.,0.,0.,.006,.002,0.],8);self.assertEqual(old,saved)

    def test_plastic_tensor_is_rotation_covariant(self):
        p=material();angle=.37;r=np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])
        e=np.array([[.001,.003,.001],[.003,-.0004,0.],[.001,0.,.0001]]);a=ref.integrate(p,state(),voigt(e,True),8);b=ref.integrate(p,state(),voigt(r@e@r.T,True),8)
        np.testing.assert_allclose(mc.tensor(b['effective_stress_compression_positive_kpa']),r@mc.tensor(a['effective_stress_compression_positive_kpa'])@r.T,atol=1e-8,rtol=1e-9)
        pa=np.asarray(a['plastic_strain_compression_positive']);pb=np.asarray(b['plastic_strain_compression_positive']);pa[3:]/=2;pb[3:]/=2
        np.testing.assert_allclose(mc.tensor(pb),r@mc.tensor(pa)@r.T,atol=1e-11,rtol=1e-9)
        self.assertAlmostEqual(a['equivalent_plastic_strain'],b['equivalent_plastic_strain'],places=11)
        self.assertAlmostEqual(a['plastic_work_kpa'],b['plastic_work_kpa'],places=9)

    def test_affine_hydrostatic_internal_force_preserves_rigid_work(self):
        mesh,ground,case=inputs();m=ref.Model(mesh,ground,case);force=m.pore_force.reshape(-1,3)
        np.testing.assert_allclose(force.sum(axis=0),0,atol=1e-12,rtol=0)
        np.testing.assert_allclose(np.cross(m.xyz,force).sum(axis=0),0,atol=1e-12,rtol=0)
        self.assertGreater(np.linalg.norm(force),0.)
        np.testing.assert_allclose(m.strains(np.tile([.1,.2,.3],len(m.xyz))),0,atol=1e-15,rtol=0)

    def test_rankine_apex_keeps_tensile_surface_and_return_category(self):
        p=material();p['tension_cutoff_kpa']=2.;new=ref.integrate(p,state(),[-.005,-.005,-.005,0.,0.,0.],4)
        np.testing.assert_allclose(new['effective_stress_compression_positive_kpa'][:3],[-2.]*3,atol=1e-9,rtol=0)
        self.assertEqual(new['last_return'],'tension_apex');self.assertGreater(new['plastic_work_kpa'],0.)

    def test_material_subdivision_is_exercised_without_changing_yield(self):
        p=material();de=[.001,-.0004,.0001,.006,.002,0.]
        a=ref.integrate(p,state(),de,1);b=ref.integrate(p,state(),de,16)
        self.assertGreater(abs(a['plastic_work_kpa']-b['plastic_work_kpa']),1e-8)
        for result in (a,b):self.assertLess(float(mc.normalized_yield(p,result['effective_stress_compression_positive_kpa'])),1e-8)

    def test_mcc_is_explicitly_outside_this_reference(self):
        mesh,ground,case=inputs();ground['materials'][0]['mcc']={}
        with self.assertRaisesRegex(ValueError,'MC reference only'):ref.Model(mesh,ground,case)

if __name__=='__main__':unittest.main()
