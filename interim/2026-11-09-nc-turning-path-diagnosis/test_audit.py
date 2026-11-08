"""Corruption guards on the actual diagnostic evidence, not acceptance tests."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import unittest
import numpy as np
import audit

EVIDENCE=None
OUTPUT=None
class Guards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        e=EVIDENCE
        records=audit.read(e/'field-records.json')
        seed=next(r['event'] for r in records if r['label']=='terminal-rejected:0')
        cls.cp=seed['committed'];cls.q,cls.offsets=audit.prepare(audit.read(e/'prior/fresh-head-midpoint/input/mesh.json'))
        nd=len(seed['reduced_residual']);cls.R=np.zeros((len(seed['displacement_m'])*3,nd))
        for g,row in enumerate(seed['dof_rows']):
            for i,c in row:cls.R[g,i]=c
        run=audit.read(e/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')['run']
        cls.holding=np.array(run['excavation_audits'][0]['initial_holding_load_kn']).ravel()
        cls.load=np.array(cls.cp['applied_nodal_load_kn']).ravel()
        cls.folder=e/'path-diagnostic-coarse';cls.sample=audit.read(cls.folder/'samples.json')[0]
        cls.control=audit.read(cls.folder/'controls.json')
    def check(self,s):
        return audit.check_sample(s,self.cp,self.cp,self.R,self.q,self.offsets,self.load,self.holding,
            np.array(self.control['plane']),self.control['step_size_m'])
    def test_saved_sample_reassembles_under_original_force_gate(self):
        r,indices=self.check(self.sample)
        self.assertLess(r,1e-8);self.assertEqual(len(indices),162)
    def test_changed_active_stress_is_rejected(self):
        s=copy.deepcopy(self.sample);s['states'][0]['effective_stress_compression_positive_kpa'][0]+=1.
        with self.assertRaisesRegex(ValueError,'force gate'):self.check(s)
    def test_changed_removed_history_is_rejected(self):
        s=copy.deepcopy(self.sample);s['states'][162]['preconsolidation_pressure_kpa']+=1e-4
        with self.assertRaisesRegex(ValueError,'removed material'):self.check(s)
    def test_changed_displacement_control_is_rejected(self):
        s=copy.deepcopy(self.sample);d=self.R@np.array(self.control['plane'])*1e-5
        s['displacement_m']=(np.array(s['displacement_m'])+d.reshape(-1,3)).tolist()
        with self.assertRaisesRegex(ValueError,'displacement plane'):self.check(s)
    def test_changed_progress_bookkeeping_is_rejected(self):
        s=copy.deepcopy(self.sample);s['progress_increment']+=1e-5
        with self.assertRaisesRegex(ValueError,'bookkeeping'):self.check(s)
    def test_corrupted_linear_correction_is_rejected(self):
        m=audit.read(next(self.folder.glob('linear-*.json')))
        audit.check_linear(m);m['solution']=(np.array(m['solution'])*1.01).tolist()
        with self.assertRaisesRegex(ValueError,'linear gate'):audit.check_linear(m)
    def test_wrong_material_scheme_is_explicitly_refused(self):
        out=OUTPUT/'wrong-scheme';out.mkdir(parents=True,exist_ok=False)
        req=audit.read(EVIDENCE/'material-requests.json')[:1]
        req[0]['parameters']['integration']='backward_euler'
        audit.write(out/'request.json',req)
        binary=EVIDENCE/('target/release/examples/material_diagnostic'+('.exe' if os.name=='nt' else ''))
        with (out/'process.log').open('x') as log:
            p=subprocess.run([str(binary),str(out/'request.json'),str(out/'result.json')],stdout=log,stderr=subprocess.STDOUT,check=False)
        audit.write(out/'process.json',dict(exit_code=p.returncode,independent_reference=False,passed=False))
        self.assertEqual(p.returncode,0)
        point=audit.read(out/'result.json')['points'][0]
        self.assertIn('requires midpoint',point['error']);self.assertNotIn('updated_state',point)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--evidence',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args,rest=parser.parse_known_args();EVIDENCE=args.evidence.resolve();OUTPUT=args.output.resolve();OUTPUT.mkdir(parents=True,exist_ok=False)
    unittest.main(argv=[__file__,*rest])
