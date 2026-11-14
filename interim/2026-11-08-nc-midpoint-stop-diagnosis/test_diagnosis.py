"""Verify that the isolated diagnostic rejects changed seals and out-of-window work."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

class DiagnosticGuardTests(unittest.TestCase):
    def rejected(self,mesh,checkpoint,target):
        with tempfile.TemporaryDirectory(prefix='nc-diagnostic-guard-') as d:
            d=Path(d);(d/'mesh.json').write_text(json.dumps(mesh));(d/'checkpoint.json').write_text(json.dumps(checkpoint))
            r=subprocess.run([str(PROBE),str(d/'mesh.json'),str(d/'checkpoint.json'),str(target),str(d/'output')],capture_output=True,text=True,timeout=20)
            self.assertNotEqual(r.returncode,0)
            self.assertIn('diagnostic',r.stderr)
            self.assertFalse((d/'output/diagnostic.json').exists())
            self.assertEqual((d/'output/trace.jsonl').read_text(),'')
    def values(self):
        return (json.loads((PRIOR/'nc-resumed/input/mesh.json').read_text()),
                json.loads((PRIOR/'nc-resumed/job/recovery-checkpoint-0.json').read_text()))
    def test_changed_mesh_rejected(self):
        m,c=self.values();m['nodes'][0][0]+=.001;self.rejected(m,c,.9514002799987793)
    def test_changed_history_rejected(self):
        m,c=self.values();c['run']['last_committed']['material_states'][0]['preconsolidation_pressure_kpa']+=.001
        self.rejected(m,c,.9514002799987793)
    def test_changed_producer_identity_rejected(self):
        m,c=self.values();c['identity']['solver_binary_sha256']='0'*64;self.rejected(m,c,.9514002799987793)
    def test_outside_window_and_no_progress_rejected(self):
        m,c=self.values()
        for target in (1.,.9514,.9514000415802002):
            with self.subTest(target=target):self.rejected(m,c,target)
    def test_recorded_success_and_failure_remain_distinct(self):
        _,c=self.values()
        a=json.loads((TRACES/'last-accepted/diagnostic.json').read_text())
        b=json.loads((TRACES/'terminal-rejected/diagnostic.json').read_text())
        self.assertEqual(a['outcome']['state'],c['run']['last_committed'])
        self.assertTrue(a['outcome']['converged']);self.assertFalse(b['outcome']['converged'])
        self.assertEqual(b['outcome']['error'],'numerical solve failed: true linear residual 1.4875587579753504e-10')
        for row in (a,b):
            self.assertFalse(row['passed']);self.assertFalse(row['checkpoint_resume']);self.assertFalse(row['new_continuation_created'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',type=Path,required=True);p.add_argument('--prior',type=Path,required=True);p.add_argument('--traces',type=Path,required=True);a,rest=p.parse_known_args()
    PROBE=a.probe.resolve();PRIOR=a.prior.resolve();TRACES=a.traces.resolve()
    unittest.main(argv=[__file__]+rest)
