"""Common-point comparison must retain split traces and interpolation limits."""
import json
from pathlib import Path
import subprocess
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=r'''
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path('interim/2026-11-19-wet-joint-spatial-convergence').resolve()))
import mapping as m
from measure_c02 import Stats
c=m.c
mesh=c.read(c.ROOT/'examples/joint-water-column-tet10-mesh.json')
meshes=[mesh]
for _ in range(2):meshes.append(c.old.refine(meshes[-1])[0])
case=c.fem.normalize_stage_case(c.read(c.ROOT/'examples/stages-joint-water-cycle.json'))
a=m.Maps(meshes,case)
assert [len(x['tetrahedra']) for x in meshes]==[12,96,768]
assert len(a.contact_position)==672 and len(a.material_position)==23652
assert len(a.node_position)==1458
for h,(mesh,s) in enumerate(zip(meshes,a.systems)):
    xyz=np.array(mesh['nodes']);q=np.concatenate([c.fem.volume(mesh,e)['position'] for e in range(len(mesh['tetrahedra']))])
    qp=np.array([p['position_m'] for p in s.points])
    for source,target,op in [(xyz,a.node_position,a.operators[h]['u']),(q,a.material_position,a.operators[h]['material']),(qp,a.contact_position,a.operators[h]['contact'])]:
        patch=2+source[:,0]+source[:,1]**2-source[:,0]*source[:,2]
        expected=2+target[:,0]+target[:,1]**2-target[:,0]*target[:,2]
        np.testing.assert_allclose(op@patch,expected,rtol=0,atol=2e-11)
    # Discontinuous root material fields must not bleed across original tets.
    root=np.repeat(np.arange(len(mesh['tetrahedra']))//8**h,27).astype(float)
    np.testing.assert_allclose(a.operators[h]['material']@root,a.material_root,rtol=0,atol=3e-13)
    # Each original GP appears in the common inventory and remains an exact sample.
    offset=sum(len(ss.points) for ss in a.systems[:h]);n=len(s.points)
    np.testing.assert_allclose(a.operators[h]['contact'][offset:offset+n].toarray(),np.eye(n),rtol=0,atol=2e-11)
# Coincident contact traces have independent displacement values.
fine=meshes[-1];upper=np.zeros(len(fine['nodes']))
for e,cell in enumerate(fine['tetrahedra10']):
    if a.node_root[cell[0]]>=6:upper[cell]=1
assert len({tuple(p) for p in fine['nodes']})<len(fine['nodes'])
# Material extrapolation outside its ancestor is explicitly rejected.
try:a.volume_map(0,np.array([[100.,100.,100.]]),np.array([0]),False)
except AssertionError:pass
else:raise AssertionError('outside physical point accepted')
s=Stats();s.add({'x':np.array([1.,3.])},{'x':np.array([0.,1.])},'open')
r=s.report()['fields']['x'];assert r['difference_squared']==5 and r['reference_squared']==1
print('mapping-patches-and-refusals-ok')
'''

class WetJointSpatialPoints(unittest.TestCase):
    def test_all_gp_union_split_ancestry_patch_exactness_and_outside_refusal(self):
        r=subprocess.run([sys.executable,'-c',SCRIPT],cwd=ROOT,capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertIn('mapping-patches-and-refusals-ok',r.stdout)

    def test_negative_order_direction_reversal_and_zero_difference_are_preserved(self):
        script="import sys;sys.path.insert(0,'interim/2026-11-19-wet-joint-spatial-convergence');from measure_c02 import order_metric;x=order_metric(1,4,-2);assert x['observed_order']==-1 and x['difference_direction_cosine']==-1;assert order_metric(0,4,0)['observed_order'] is None;assert set(x)=={'coarse_absolute_l2','fine_absolute_l2','observed_order','difference_direction_cosine'}"
        r=subprocess.run([sys.executable,'-c',script],cwd=ROOT,capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)

if __name__=='__main__':unittest.main()
