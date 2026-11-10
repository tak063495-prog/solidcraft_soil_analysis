"""New BIRTH gates use a separate reference, with unchanged numerical limits."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from verify_joint_birth import fixture,shear_fixture,run_case
from joint_birth_reference import point_path,audit_elastic,elastic_reference,audit_hydrostatic_reassembly
BINARY=ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')

class JointBirth(unittest.TestCase):
    def test_independent_dry_modes_and_known_inherited_shear(self):
        with tempfile.TemporaryDirectory() as temp:
            for mode in ('transfer_force','stress_free','shear'):
                mesh,ground,case=shear_fixture() if mode=='shear' else fixture(mode,False)
                run,_,_=run_case(BINARY,Path(temp)/mode,mesh,ground,case)
                from joint_birth_io import export
                import meshio
                folder=Path(temp)/mode;export(folder/'mesh.json',folder/'run.json',folder/'exports')
                for fi,frame in enumerate(run['frames']):
                    data=meshio.read(folder/'exports'/f'joints-{fi:03}.vtu').point_data
                    life=frame['state']['joints']['lifecycle'];points=frame['state']['joints']['points']
                    self.assertTrue(np.array_equal(data['joint_generation'],np.full(len(points),life[0]['generation'],dtype=np.uint64)))
                    self.assertTrue(np.array_equal(data['birth_initial_traction_kpa'],np.asarray([(p.get('birth') or {}).get('initial_traction_kpa',[0.,0.,0.]) for p in points])))
                    self.assertEqual(json.loads((folder/'exports'/f'joint-generations-{fi:03}.json').read_text()),life)
                self.assertTrue(audit_elastic(mesh,ground,run)['passed'])
                prepared=run['transition_audits'][0]['prepared_state']['joints']['points']
                for p in prepared:
                    self.assertEqual(p['plastic_slip_m'],[0.,0.,0.]);self.assertEqual(p['plastic_dissipation_kpa_m'],0.)
                    self.assertEqual(p['numerical_dissipation_kpa_m'],0.);self.assertEqual(p['endpoint_work_kpa_m'],0.)
                    if mode=='shear':self.assertLess(np.linalg.norm(np.asarray(p['traction_kpa'])-[12.,0.,-100.]),1e-7)
    def test_independent_water_roles_reassemble_each_rhs_once(self):
        with tempfile.TemporaryDirectory() as temp:
            for policy in ('hydrostatic','flow','soil','prescribed'):
                mesh,ground,case=fixture(water=policy)
                run,_,_=run_case(BINARY,Path(temp)/policy,mesh,ground,case)
                self.assertTrue(audit_elastic(mesh,ground,run)['passed'])
                from joint_birth_io import export
                import meshio
                folder=Path(temp)/policy;export(folder/'mesh.json',folder/'run.json',folder/'exports')
                final=meshio.read(folder/'exports'/f'integration-points-{len(run["frames"])-1:03}.vtu').point_data
                self.assertIn('material_stress_reference',final)
                if policy=='soil':self.assertTrue(np.any(final['material_stress_reference']==1))
    def test_reference_refuses_unknown_birth_combinations_and_nonlinear_scope(self):
        mesh,ground,case=fixture()
        for stage in ('excavate','joint_death','birth'):
            changed=copy.deepcopy(case);changed['stages'][0]['kind']=stage
            with self.assertRaises(AssertionError):elastic_reference(mesh,ground,changed)
            with self.assertRaises(AssertionError):audit_hydrostatic_reassembly(mesh,ground,dict(complete=True,schema='soil-fem-stage-result/15',case=changed,frames=[]))
        changed=copy.deepcopy(case);changed['joint_system']['joints'][0]['crack_tip']={'edges':[[0,1,2]]}
        with self.assertRaises(AssertionError):elastic_reference(mesh,ground,changed)
    def test_scalar_contact_reference_enforces_initial_cone_and_energy(self):
        data=dict(material=dict(normal_stiffness_kpa_m=1e6,tangential_stiffness_kpa_m=2e5,friction_coefficient=.3),normal=[0,0,1],reference=[.02,-.01,.03],traction=[12.,5.,-100.],maximum_increment_m=1e-5,increments=[[.003,.001,-.0001],[.003,.001,.001],[.004,.002,-.0001]])
        for traction in ([0,0,1],[31,0,-100]):
            invalid={**data,'traction':traction}
            with self.assertRaises(AssertionError):point_path(invalid)
        for p in point_path(data):
            q=p['traction_kpa'];energy=q[2]**2/2e6+np.linalg.norm(q[:2])**2/4e5
            self.assertLess(abs(p['initial_energy']+p['endpoint_work_kpa_m']-energy-p['plastic_dissipation_kpa_m']-p['numerical_dissipation_kpa_m']),1e-10)
    def test_new_reference_does_not_expand_existing_wet_structure_allowlist(self):
        from verify_wet_structures import reference
        mesh,ground,case=fixture(water='hydrostatic')
        with self.assertRaisesRegex(AssertionError,'unsupported stage'):reference(mesh,ground,case)

if __name__=='__main__':unittest.main()
