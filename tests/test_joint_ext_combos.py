"""Extension references stay bounded and require actual completed CLI runs."""
import copy
import os
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from verify_joint_ext_combos import tip_fixture,wet_fixture,solve,invoke,dump
from joint_ext_reference import audit_tip_elastic,audit_wet_birth_elastic,audit_reassembly,tip_elastic,constraint_space,audit_wet_birth_holding,audit_clean_path
from joint_reference import System
BINARY=Path(os.environ.get('JOINT_EXT_SOLVER',str(ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d'))))

class JointExtensions(unittest.TestCase):
    def test_tip_modes_independent_elastic_all_nodes_and_gp(self):
        with tempfile.TemporaryDirectory() as temp:
            for policy in ('transfer_force','stress_free'):
                mesh,ground,case=tip_fixture(policy)
                run,_=solve(BINARY,Path(temp)/policy,mesh,ground,case)
                self.assertTrue(audit_tip_elastic(mesh,ground,run)['passed'])
                self.assertTrue(audit_reassembly(mesh,ground,run)['passed'])
                self.assertTrue(audit_clean_path(run,copy.deepcopy(run))['passed'])
                changed=copy.deepcopy(run)
                changed['transition_audits'][0]['initial_holding_load_kn'][0]+=1.
                with self.assertRaises(AssertionError):audit_clean_path(run,changed)
                initial=run['transition_audits'][0]['prepared_state']['joints']['points']
                if policy=='transfer_force':
                    self.assertGreater(max(abs(p['traction_kpa'][0]) for p in initial),1e-3)
                    self.assertGreater(np.ptp([p['traction_kpa'][2] for p in initial]),1.)
                else:self.assertTrue(all(p['traction_kpa']==[0.,0.,0.] for p in initial))

    def test_wet_birth_input_only_placement_and_elastic_reload(self):
        mesh,ground,case=wet_fixture();case['stages']=case['stages'][:3];case['verification_reject_once_stages']=[]
        case['settings']['kinematics']='side_rollers';case.pop('boundary',None)
        with tempfile.TemporaryDirectory() as temp:
            run,_=solve(BINARY,Path(temp)/'wet',mesh,ground,case)
            self.assertTrue(audit_wet_birth_elastic(mesh,ground,run)['passed'])
            self.assertTrue(audit_reassembly(mesh,ground,run)['passed'])
            self.assertTrue(audit_wet_birth_holding(mesh,ground,run)['passed'])

    def test_shared_tip_stays_a_single_free_physical_node_in_global_tie_nullspace(self):
        mesh,_,case=tip_fixture();case['boundary']={}
        system=System(mesh,case['joint_system'])
        t,_=constraint_space(mesh,case,[True]*len(mesh['tetrahedra']),system,True)
        for n in case['joint_system']['joints'][0]['crack_tip']['edges'][0]:
            for axis in range(3):
                u=np.zeros(t.shape[0]);u[3*n+axis]=1.
                self.assertLess(np.max(np.abs(t@(t.T@u)-u)),1e-10)

    def test_new_references_refuse_unsupported_and_incomplete_paths(self):
        mesh,ground,case=tip_fixture()
        for kind in ('birth','steady_seepage','richards'):
            bad=copy.deepcopy(case);bad['stages'][0]['kind']=kind
            with self.assertRaises(AssertionError):tip_elastic(mesh,ground,bad)
        for audit in (audit_tip_elastic,audit_reassembly,audit_wet_birth_elastic):
            with self.assertRaises(AssertionError):audit(mesh,ground,dict(complete=False,case=case,frames=[]))
        from joint_birth_reference import elastic_reference
        from verify_wet_structures import reference
        with self.assertRaises(AssertionError):elastic_reference(mesh,ground,case)
        with self.assertRaises(AssertionError):reference(mesh,ground,case)

    def test_missing_hydraulic_input_and_unsupported_flow_initialization_are_rejected(self):
        mesh,ground,case=wet_fixture()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for n,d in [('mesh',mesh),('ground',ground)]:dump(root/f'{n}.json',d)
            for i,kind in enumerate((None,'steady_seepage','richards')):
                bad=copy.deepcopy(case)
                if kind is None:bad['stages'][0].pop('hydraulic_initialization')
                else:bad['stages'][0]['hydraulic_initialization']['kind']=kind
                dump(root/'case.json',bad)
                invoke(BINARY,['solve-stages',root/'mesh.json',root/'ground.json',root/'case.json',root/f'refused-{i}.json'],root/f'refused-{i}.log',False)
                self.assertFalse((root/f'refused-{i}.json').exists())

if __name__=='__main__':unittest.main()
