"""A-01 observation/reference properties, isolated from historical imports."""
from pathlib import Path
import subprocess
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'interim/2026-11-23-joint-combo-extension'


class ComboExtension(unittest.TestCase):
    def check_code(self,body):
        prefix=("import sys,copy,json,tempfile\nfrom pathlib import Path\nimport numpy as np\n"
                f"sys.path.insert(0,{str(STUDY)!r})\n"
                "from fixtures import *\nfrom reference import Model,run\nfrom fields import compare,Sink\n")
        result=subprocess.run([sys.executable,'-c',prefix+body],cwd=ROOT,text=True,capture_output=True,timeout=90)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_two_nonintersecting_interfaces_have_no_shared_physical_nodes(self):
        self.check_code("""
m,g,c=separated();model=Model(m,g,c);s=model.system
a=set(n for p in s.points if p['joint_id']=='Seam' for n in p['nodes_a']+p['nodes_b'])
b=set(n for p in s.points if p['joint_id']=='Second' for n in p['nodes_a']+p['nodes_b'])
assert a and b and not a&b
assert len(m['tetrahedra'])==24 and len(s.points)==12
np.testing.assert_allclose(s.c.sum(axis=2),0,atol=1e-15,rtol=0)
""")

    def test_shared_edge_uses_the_same_side_nodes_and_distinct_opposite_dofs(self):
        self.check_code("""
m,g,c=shared();model=Model(m,g,c);s=model.system
assert len(set(m['node_sets']['First:a'])&set(m['node_sets']['Second:a']))==2
assert len(set(m['node_sets']['First:b'])&set(m['node_sets']['Second:b']))==2
assert not set(m['node_sets']['First:a'])&set(m['node_sets']['First:b'])
assert len(s.points)==6
""")

    def test_both_policies_form_two_generations_and_frozen_archives(self):
        self.check_code("""
for policy in ('transfer_force','stress_free'):
 rows=[];model=Model(*shared(policy,steps=2));r=run(model,2,lambda a,b:rows.append((a,b)))
 assert r['complete'] and len(rows)==24
 assert np.array_equal(rows[-1][1]['generation'],[2,2])
 assert np.array_equal(rows[-1][1]['archived_generation_count'],[1,1])
 selected=[f for row,f in rows if row['stage_index']>=8]
 for f in selected[1:]:
  for k in selected[0]:
   if k.startswith('frozen_'):np.testing.assert_array_equal(f[k],selected[0][k])
 for row,f in rows:
  if row['stage']=='rebirth' and policy=='stress_free':assert np.all(f['initial_traction']==0)
""")

    def test_frame_plate_and_tendon_are_coupled_and_nonzero(self):
        self.check_code("""
for kind in ('beam','plate','anchor'):
 rows=[];r=run(Model(*structure(kind,steps=2)),2,lambda a,b:rows.append(b))
 assert r['complete']
 key={'beam':'beam_Beam_local_end_actions','plate':'plate_Plate_membrane_energy_kn_m','anchor':'anchor_0_segment_forces_kn'}[kind]
 assert max(np.max(np.abs(f[key])) for f in rows)>0
 if kind=='plate':
  for f in rows:
   np.testing.assert_array_equal(f['plate_Plate_transverse_shear_strain'],np.zeros((3,2)))
   np.testing.assert_array_equal(f['plate_Plate_transverse_shear_force_kn_m'],np.zeros((3,2)))
""")

    def test_reference_refuses_undeclared_wet_or_tip_scope(self):
        self.check_code("""
for change in ('wet','tip'):
 m,g,c=specimen()
 if change=='wet':g['initialization']['water_level_m']=1.5
 else:c['joint_system']['joints'][0]['crack_tip']={}
 try:Model(m,g,c)
 except ValueError:pass
 else:raise AssertionError('unsupported reference scope admitted')
""")

    def test_comparator_detects_single_nodal_field_corruption(self):
        self.check_code(r"""
with tempfile.TemporaryDirectory() as tmp:
 a=Path(tmp)/'a';b=Path(tmp)/'b'
 v=dict(row=dict(stage='open',progress=.5),fields=dict(u=[[0.,0.,0.],[1.,2.,3.]]))
 a.write_text(json.dumps(v)+'\n');v['fields']['u'][1][2]+=1e-4;b.write_text(json.dumps(v)+'\n')
 r=compare(a,b);assert r['all_fields']['u']['maximum_absolute']>9e-5
 assert r['numerical_acceptance'] is False
""")

    def test_comparator_rejects_missing_accepted_state(self):
        self.check_code(r"""
with tempfile.TemporaryDirectory() as tmp:
 a=Path(tmp)/'a';b=Path(tmp)/'b';v=json.dumps(dict(row=dict(stage='open',progress=.5),fields=dict(u=[1.,2.,3.])))+'\n'
 a.write_text(v);b.write_text(v+v)
 try:compare(a,b)
 except AssertionError as e:assert 'missing accepted state' in str(e)
 else:raise AssertionError('missing state was ignored')
""")

    def test_restart_equality_preserves_signed_zero_and_all_nonclock_fields(self):
        self.check_code("""
from study import exact
a=dict(elapsed_seconds=1.,state=dict(generation=2,u=[-0.,1.]))
b=copy.deepcopy(a);b['elapsed_seconds']=5.;assert exact(a)==exact(b)
b['state']['u'][0]=0.;assert exact(a)!=exact(b)
b=copy.deepcopy(a);b['state']['generation']=3;assert exact(a)!=exact(b)
""")


if __name__=='__main__':unittest.main()
