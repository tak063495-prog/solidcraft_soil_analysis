"""Fail-closed tests using the actually captured minimum-trial fields."""
import argparse
import copy
from pathlib import Path
import sys
import unittest
import audit_evidence as audit

ROOT=None

class AuditGuards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        prior=ROOT/'experiment/prior/evidence/prior/nc-resumed'
        cls.mesh,cls.ground,cls.case=[audit.read(prior/'input'/(n+'.json')) for n in audit.low.NAMES]
        cls.saved_run=audit.read(prior/'job/recovery-checkpoint-0.json')['run']
        cls.trace=audit.events(ROOT/'experiment/retain-trial/trace.jsonl')
        cls.first=next(e for e in cls.trace if e['kind']=='iterate')

    def measure(self,trace):
        return audit.force_history(self.mesh,self.ground,self.case,self.saved_run,trace,True)

    def test_missing_trace_and_missing_linear_data_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'missing Newton'):self.measure([])
        with self.assertRaisesRegex(ValueError,'missing independent linear'):audit.linear_audit([],[])

    def test_contact_activity_is_not_hidden_by_stale_stick_label(self):
        e=copy.deepcopy(self.first);e['joint_report']['points'][0]['active']=True
        with self.assertRaisesRegex(ValueError,'contact history or activity'):self.measure([e])

    def test_removed_material_freeze_is_checked(self):
        e=copy.deepcopy(self.first)
        index=next(i for i,a in enumerate(self.saved_run['last_committed']['active_elements']) if not a)*27
        e['states'][index]['preconsolidation_pressure_kpa']+=1.
        with self.assertRaisesRegex(ValueError,'removed material'):self.measure([e])

    def test_reassembly_detects_edited_active_stress(self):
        e=copy.deepcopy(self.first)
        index=next(i for i,a in enumerate(self.saved_run['last_committed']['active_elements']) if a)*27
        e['states'][index]['effective_stress_compression_positive_kpa']=[.999*x for x in e['states'][index]['effective_stress_compression_positive_kpa']]
        with self.assertRaisesRegex(ValueError,'force reassembly'):self.measure([e])

    def test_affine_release_is_not_reconstructed_from_trial_stress(self):
        e=copy.deepcopy(self.first);e['target'][0]+=1e-5
        with self.assertRaisesRegex(ValueError,'affine release'):self.measure([e])

    def test_repartition_must_retain_current_evaluated_iterate(self):
        trace=copy.deepcopy(self.trace)
        position=next(i for i,e in enumerate(trace) if e['kind']=='repartition')
        next_iterate=next(e for e in trace[position+1:] if e['kind']=='iterate')
        next_iterate['displacement_m'][0][0]+=1e-5
        with self.assertRaisesRegex(ValueError,'retain evaluated current'):self.measure(trace)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--evidence',type=Path,required=True)
    args,remaining=parser.parse_known_args();ROOT=args.evidence.resolve()
    unittest.main(argv=[sys.argv[0],*remaining])
