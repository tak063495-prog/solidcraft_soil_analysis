"""Independent shared-tip projection gates; no saved JSON is used as a solve."""
from pathlib import Path
import subprocess
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-11-25-joint-birth-projection-a02'


class BirthProjectionA02(unittest.TestCase):
    def check(self,body):
        prefix=(f"import sys,copy,json,tempfile\nfrom pathlib import Path\nsys.path.insert(0,{str(HERE)!r})\n"
                "from study import fixture,Model,exact\nfrom projection import recover\nimport numpy as np\n"
                "m,g,c=fixture('transfer_force',2);s=Model(m,g,c).system\n"
                "def force(q):return np.einsum('qag,a,q->g',s.c,q,s.weight)\n")
        result=subprocess.run([sys.executable,'-c',prefix+body],cwd=ROOT,text=True,capture_output=True,timeout=90)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_all_nodes_reproduce_without_tip_forces(self):
        self.check("""
f=force([12.,5.,-100.]);q,gate=recover(s,f)
assert max(gate.values())<1e-10
np.testing.assert_allclose(np.einsum('qag,qa,q->g',s.c,q,s.weight),f,rtol=0,atol=1e-10)
for n in c['joint_system']['joints'][0]['crack_tip']['edges'][0]:assert np.all(s.c[:,:,3*n:3*n+3]==0)
""")

    def test_balanced_shared_tip_residual_is_rejected(self):
        self.check("""
f=force([12.,5.,-100.]);before=f.copy();tip=c['joint_system']['joints'][0]['crack_tip']['edges'][0]
d=np.array(m['nodes'][tip[1]])-m['nodes'][tip[0]]
f[3*tip[0]:3*tip[0]+3]+=d;f[3*tip[1]:3*tip[1]+3]-=d
delta=(f-before).reshape(-1,3)
assert np.max(abs(delta.sum(axis=0)))<1e-10 and np.max(abs(np.cross(m['nodes'],delta).sum(axis=0)))<1e-10
try:recover(s,f)
except ValueError as e:assert 'free nodal' in str(e)
else:raise AssertionError('balanced unreproducible force admitted')
""")

    def test_tensile_and_outside_circle_are_rejected_without_clipping(self):
        self.check("""
for q in ([0.,0.,100.],[40.,0.,-100.]):
 f=force(q);old=f.copy()
 try:recover(s,f)
 except ValueError as e:assert 'cone' in str(e)
 else:raise AssertionError('inadmissible projection admitted')
 np.testing.assert_array_equal(f,old)
""")

    def test_reference_keeps_scope_rejections(self):
        self.check("""
for mode in ('wet','inactive','body','multiple'):
 mm,gg,cc=copy.deepcopy((m,g,c))
 if mode=='wet':gg['initialization']['water_level_m']=1.5
 elif mode=='inactive':cc['initially_inactive_bodies']=['UpperClay']
 elif mode=='body':cc['stages'][0]['kind']='birth'
 else:cc['joint_system']['joints']*=2
 try:Model(mm,gg,cc)
 except ValueError:pass
 else:raise AssertionError(mode)
""")

    def test_comparison_rejects_missing_state_and_preserves_signed_zero(self):
        self.check("""
from fields import compare
assert exact({'elapsed_seconds':1,'u':[-0.]})==exact({'elapsed_seconds':2,'u':[-0.]})
assert exact({'u':[-0.]})!=exact({'u':[0.]})
with tempfile.TemporaryDirectory() as temp:
 a=Path(temp)/'a';b=Path(temp)/'b';v=json.dumps(dict(row=dict(stage='birth',progress=.5),fields=dict(u=[1.,2.,3.])))+'\\n'
 a.write_text(v);b.write_text(v+v)
 try:compare(a,b)
 except AssertionError as e:assert 'missing accepted state' in str(e)
 else:raise AssertionError('state omitted')
""")


if __name__=='__main__':unittest.main()
