from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import unittest
import numpy as np
from mcc_acoustic_diagnostic import branch,acoustic_tensor,directions
from mcc_acoustic_exact import criterion

M=dict(critical_state_m=1.4,lambda_=.18,kappa=.04,poisson=.3)
M['lambda']=M.pop('lambda_')
def data(p=100.,pc=400.,extension=True):
 q=1.4*np.sqrt(p*(pc-p))
 stress=[p+q/3,p+q/3,p-2*q/3,0.,0.,0.] if extension else [p+2*q/3,p-q/3,p-q/3,0.,0.,0.]
 return branch(stress,pc,.9,M)
def tensor(cn):return np.array([[cn[0],cn[3],cn[5]],[cn[3],cn[1],cn[4]],[cn[5],cn[4],cn[2]]])
class Exact(unittest.TestCase):
 def test_extension_has_independently_negative_rayleigh_witness(self):
  d=data();r=criterion(d);self.assertEqual(r['decision'],'negative_continuum_candidate')
  n,z=np.asarray(r['normal']),np.asarray(r['polarization']);self.assertLess(z@acoustic_tensor(d,n,plastic=True)@z,0)
 def test_softening_compression_is_positive_control(self):
  d=data(extension=False);self.assertLess(d['consistency_H_kpa3'],0)
  self.assertEqual(criterion(d)['decision'],'positive_continuum_candidate_all_directions')
 def test_hydrostatic_degenerate_axes_match_closed_form(self):
  d=branch([100.,100.,100.,0.,0.,0.],100.,.9,M);r=criterion(d);a=tensor(d['cn'])[0,0]
  expected=a*a/((d['bulk_kpa']+4*d['shear_kpa']/3)*d['denominator_kpa3'])
  self.assertAlmostEqual(r['maximum_rank_one_ratio'],expected,places=13)
 def test_dense_angular_search_cannot_exceed_analytic_maximum(self):
  rng=np.random.default_rng(385);normals=rng.normal(size=(12000,3));normals/=np.linalg.norm(normals,axis=1)[:,None]
  for extension in [True,False]:
   for pc in [250.,300.,400.,500.]:
    d=data(pc=pc,extension=extension);a=tensor(d['cn']);v=normals@a.T
    beta=(d['bulk_kpa']+d['shear_kpa']/3)/(d['bulk_kpa']+4*d['shear_kpa']/3)
    sampled=(np.sum(v*v,axis=1)-beta*np.sum(v*normals,axis=1)**2)/(d['shear_kpa']*d['denominator_kpa3'])
    exact=criterion(d)['maximum_rank_one_ratio']
    self.assertLessEqual(float(sampled.max()),exact+1e-12)
    self.assertLess(exact-float(sampled.max()),.002)
 def test_rotation_covariance_of_criterion(self):
  d=data();rng=np.random.default_rng(93);q,_=np.linalg.qr(rng.normal(size=(3,3)));a=q@tensor(d['cn'])@q.T
  changed=dict(d,cn=np.array([a[0,0],a[1,1],a[2,2],a[0,1],a[1,2],a[2,0]]))
  self.assertAlmostEqual(criterion(d)['maximum_rank_one_ratio'],criterion(changed)['maximum_rank_one_ratio'],places=12)
 def test_ratio_matches_full_positive_definiteness_for_many_directions(self):
  d=data();a=tensor(d['cn'])
  for n in directions(5,20):
   qe=acoustic_tensor(d,n,plastic=False);ratio=float((a@n)@np.linalg.solve(qe,a@n)/d['denominator_kpa3'])
   eigen=np.linalg.eigvalsh(acoustic_tensor(d,n,plastic=True))[0]
   if abs(ratio-1)>1e-9:self.assertEqual(eigen>0,ratio<1)
 def test_boundary_classification_is_indeterminate(self):
  d=data();r=criterion(d);d['denominator_kpa3']*=r['maximum_rank_one_ratio']
  self.assertEqual(criterion(d)['decision'],'near_boundary_indeterminate')
 def test_extreme_compressibility_and_stress_scaling(self):
  d=data();base=criterion(d)['maximum_rank_one_ratio']
  for scale in [1e-4,1e4]:
   changed=dict(d,bulk_kpa=d['bulk_kpa']*scale,shear_kpa=d['shear_kpa']*scale,denominator_kpa3=d['denominator_kpa3']*scale**3,cn=d['cn']*scale**2)
   self.assertAlmostEqual(criterion(changed)['maximum_rank_one_ratio'],base,places=12)
  for nu in [-.95,.4999]:
   m=dict(M,poisson=nu);criterion(branch([300.,100.,-100.,0.,0.,0.],400.,.9,m))
 def test_bad_inputs_are_rejected(self):
  for changes in [dict(denominator_kpa3=0),dict(shear_kpa=-1),dict(cn=[1,2,3]),dict(bulk_kpa=float('inf'))]:
   with self.assertRaises(ValueError):criterion(dict(data(),**changes))
if __name__=='__main__':unittest.main()
