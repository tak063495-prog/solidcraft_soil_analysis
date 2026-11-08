"""Compare production tetrahedral wedge integration with independent face-fan moments."""
import argparse, hashlib, json
from pathlib import Path
import subprocess, sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as fem
from hydrostatic_cut_volume_reference import force
from diagnose_cut_volume_history import inputs

def check(probe,out):
 out.mkdir(parents=True,exist_ok=False);requests=[];expected=[];labels=[];checks=[]
 mesh,_,_=inputs();rng=np.random.default_rng(20261108)
 specimens=[(np.asarray(mesh['nodes'])[nodes],.75,f'existing-{i}') for i,nodes in enumerate(mesh['tetrahedra10'])]
 for i in range(16):
  p=rng.normal(size=(4,3));jac=(p[1:]-p[0]).T
  if np.linalg.det(jac)<0:p[[1,2]]=p[[2,1]]
  if abs(np.linalg.det(jac))<.05:continue
  for count in (4,10):
   xyz=p if count==4 else np.vstack((p,[(p[a]+p[b])/2 for a,b in fem.EDGES]))
   z=sorted(p[:,2]);levels=[z[0]-.1,*[(a+b)/2 for a,b in zip(z,z[1:])],z[-1]+.1,*z]
   for j,h in enumerate(levels):specimens.append((xyz,h,f'random-{i}-tet{count}-{j}'))
 for xyz,h,label in specimens:
  m=dict(nodes=xyz.tolist(),tetrahedra=[[0,1,2,3]])
  if len(xyz)==10:m['tetrahedra10']=[list(range(10))]
  exact=force(m,0,h,10)['force'];gl=force(m,0,h,10,3)['force']
  np.testing.assert_allclose(exact,gl,atol=1e-11,rtol=1e-11)
  requests.append(dict(points=xyz.tolist(),level=float(h),gamma=10.));expected.append(exact);labels.append(label)
 # Curvature must be refused, even though the old material rule remains available.
 curved=np.asarray(mesh['nodes'])[mesh['tetrahedra10'][6]].copy();curved[4,2]+=1e-7
 requests.append(dict(points=curved.tolist(),level=.75,gamma=10.));expected.append(None);labels.append('curved-not-affine')
 (out/'requests.json').write_text(json.dumps(requests))
 subprocess.run([str(probe),str(out/'requests.json'),str(out/'production-forces.json')],check=True)
 actual=json.loads((out/'production-forces.json').read_text())
 for label,a,b in zip(labels,actual,expected,strict=True):
  if b is None:assert a is None;checks.append(dict(case=label,curved_shortcut_refused=True));continue
  np.testing.assert_allclose(a,b,atol=1e-11,rtol=1e-11)
  checks.append(dict(case=label,maximum_absolute_kn=float(np.max(abs(np.asarray(a)-b)))))
 report=dict(software_integration_check_passed=True,full_history_acceptance=False,requests=len(requests),checks=checks,
             maximum_absolute_kn=max(c.get('maximum_absolute_kn',0) for c in checks),probe_sha256=hashlib.sha256(probe.read_bytes()).hexdigest(),
             reference_source_sha256=hashlib.sha256((ROOT/'tools/hydrostatic_cut_volume_reference.py').read_bytes()).hexdigest())
 (out/'integration-check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ('requests','maximum_absolute_kn','software_integration_check_passed')}))
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();check(a.probe.resolve(),a.output.resolve())
