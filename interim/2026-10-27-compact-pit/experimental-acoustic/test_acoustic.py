from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent));from mcc_acoustic_diagnostic import branch,acoustic_tensor,directions,sample,scan_batch,diagnose
root=next(parent for parent in Path(__file__).resolve().parents if (parent/'tools/mcc_rate_fem.py').is_file());sys.path.insert(0,str(root/'tools'));from mcc_rate_fem import terms
M=dict(critical_state_m=1.4,**{'lambda':.12},kappa=.025,poisson=.2)

def tensor(v):return np.array([[v[0],v[3],v[5]],[v[3],v[1],v[4]],[v[5],v[4],v[2]]])
def voigt(t):return np.array([t[0,0],t[1,1],t[2,2],t[0,1],t[1,2],t[0,2]])
def strain_matrix(n):
 x,y,z=n
 return np.array([[x,0,0],[0,y,0],[0,0,z],[y,x,0],[0,z,y],[z,0,x]],float)
class Acoustic(unittest.TestCase):
 def test_elastic_eigenvalues_for_rotated_directions(self):
  d=branch([120,90,60,8,4,3],190,.9,M)
  for n in directions(3,12):
   values=np.linalg.eigvalsh(acoustic_tensor(d,n,plastic=False));np.testing.assert_allclose(values,[d['shear_kpa'],d['shear_kpa'],d['bulk_kpa']+4*d['shear_kpa']/3],rtol=1e-12)
 def test_yield_gradient_matches_independent_tensor_finite_difference(self):
  stress=np.array([120.,90.,60.,8.,4.,3.]);d=branch(stress,190,.9,M)
  def phi(s):
   t=tensor(s);mean=np.trace(t)/3;dev=t-mean*np.eye(3)
   return 1.5*np.trace(dev@dev)+M['critical_state_m']**2*mean*(mean-190)
  gradient=[]
  for i in range(6):
   delta=np.zeros(6);delta[i]=1e-3;gradient.append((phi(stress+delta)-phi(stress-delta))/(2e-3))
  np.testing.assert_allclose(d['normal'],gradient,rtol=1e-9,atol=1e-8)
 def test_plastic_acoustic_matches_independent_engineering_voigt_rate(self):
  for pc in [150.,270.,400.]:
   stress=np.array([120.,90.,60.,8.,4.,3.]);d=branch(stress,pc,.9,M);state=np.r_[stress,np.log(pc),np.log(1.9),np.zeros(6)][None,:];parameters=np.array([[1.4**2,.12,.025,.2]])
   _,ce,n,cn,den=terms(state,parameters);cep=ce[0]-np.outer(cn[0],cn[0])/den[0]
   for direction in directions(3,12):
    b=strain_matrix(direction);np.testing.assert_allclose(acoustic_tensor(d,direction,plastic=True),b.T@cep@b,rtol=1e-11,atol=1e-8)
 def test_plastic_tensor_is_covariant_under_arbitrary_rotation(self):
  rng=np.random.default_rng(712);stress=np.array([120.,90.,60.,8.,4.,3.]);n=np.array([1.,2.,3.]);n/=np.linalg.norm(n);d=branch(stress,270,.9,M);q=acoustic_tensor(d,n,plastic=True)
  for _ in range(12):
   r,_=np.linalg.qr(rng.normal(size=(3,3)));rotated=branch(voigt(r@tensor(stress)@r.T),270,.9,M)
   np.testing.assert_allclose(acoustic_tensor(rotated,r@n,plastic=True),r@q@r.T,rtol=1e-11,atol=1e-8)
 def test_softening_witness_is_negative_but_not_global_stability_claim(self):
  p=100.;pc=400.;q=1.4*np.sqrt(p*(pc-p));d=branch([p+q/3,p+q/3,p-2*q/3,0,0,0],pc,.9,M);self.assertLess(d['consistency_H_kpa3'],0)
  result=sample(d,directions());self.assertLess(result['eigenvalue_kpa'],0);self.assertTrue(result['actual_loading_branch_not_determined']);self.assertTrue(result['is_not_global_or_algorithmic_tangent']);self.assertAlmostEqual(result['eigenvalue_kpa'],result['rayleigh_kpa'],places=8)
 def test_softening_by_itself_does_not_imply_a_negative_sampled_acoustic_mode(self):
  p=100.;pc=400.;q=1.4*np.sqrt(p*(pc-p));d=branch([p+2*q/3,p-q/3,p-q/3,0,0,0],pc,.9,M)
  self.assertLess(d['consistency_H_kpa3'],0);self.assertGreater(sample(d,directions())['eigenvalue_kpa'],0)
 def test_invalid_pressure_parameters_direction_and_denominator_are_rejected(self):
  for bad in [dict(M,kappa=.2),dict(M,poisson=.5)]:
   with self.assertRaises(ValueError):branch([100]*3+[0]*3,200,.9,bad)
  with self.assertRaises(ValueError):branch([0]*6,200,.9,M)
  d=branch([100]*3+[0]*3,200,.9,M)
  with self.assertRaises(ValueError):acoustic_tensor(d,[1,1,1],plastic=True)
  d['denominator_kpa3']=0
  with self.assertRaises(ValueError):acoustic_tensor(d,[1,0,0],plastic=True)
class Report(unittest.TestCase):
 def inputs(self):
  import copy
  p=100.;pc=400.;q=1.4*np.sqrt(p*(pc-p));states=[dict(effective_stress_compression_positive_kpa=[p+q/3,p+q/3,p-2*q/3,0.,0.,0.],preconsolidation_pressure_kpa=pc,void_ratio=.9),dict(effective_stress_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],preconsolidation_pressure_kpa=200.,void_ratio=.9)]
  mesh=dict(tetrahedra=[[0,1,2,3],[0,1,2,4]],material_ids=[1,1]);ground=dict(materials=[dict(id=1,mcc=M)])
  result=dict(complete=True,integration_points=[dict(element_index=0,position_m=[0.,0.,0.],weight_m3=2.),dict(element_index=1,position_m=[1.,0.,0.],weight_m3=3.)],frames=[dict(name='excavated',state=dict(active_elements=[True,False],material_states=states))])
  return mesh,ground,result
 def test_batch_matches_independent_single_point_scan(self):
  data=[branch([120.,90.,60.,8.,4.,3.],pc,.9,M) for pc in [150.,270.,400.]];normals=directions(3,12);batch=scan_batch(data,normals)
  for d,b in zip(data,batch):self.assertAlmostEqual(b['normalized_by_shear'],sample(d,normals)['normalized_by_shear'],places=11)
 def test_report_integrates_active_volume_preserves_input_and_reports_candidate_witness(self):
  import copy
  inputs=self.inputs();saved=copy.deepcopy(inputs);report=diagnose(*inputs,normals=directions(6,24));self.assertEqual(inputs,saved);__import__('json').dumps(report,allow_nan=False);row=report['frames'][0]
  self.assertEqual(row['active_volume_m3'],2.);self.assertEqual(row['sampled_negative_volume_m3'],2.);self.assertEqual(row['sampled_negative_acoustic_candidates'],1);self.assertLess(row['witness']['rayleigh_kpa'],0);self.assertTrue(report['actual_continuing_loading_branch_not_determined'])
 def test_bad_inventory_and_non_mcc_material_are_rejected(self):
  inputs=self.inputs();inputs[2]['complete']=False
  with self.assertRaises(ValueError):diagnose(*inputs)
  inputs=self.inputs();inputs[2]['integration_points'][0]['weight_m3']=-1.
  with self.assertRaises(ValueError):diagnose(*inputs)
  inputs=self.inputs();inputs[1]['materials'][0]['mcc']=None
  with self.assertRaises(ValueError):diagnose(*inputs)
if __name__=='__main__':unittest.main()
