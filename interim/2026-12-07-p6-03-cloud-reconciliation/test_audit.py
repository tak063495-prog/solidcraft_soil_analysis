"""Mutation checks for original failure preservation and P6-01 consistency."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location('p603_audit_module',HERE/'audit_status.py')
audit=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(audit)

class P603AuditTests(unittest.TestCase):
    def setUp(self):self.s=audit.read(HERE/'status.json')
    def codes(self):return audit.audit(self.s,verify_files=False)['errors']
    def test_status_and_previous_ledger_agree(self):self.assertEqual(self.codes(),[])
    def test_overall_completion_rejected(self):
        self.s['p6_03_complete']=True;self.assertIn('COMPLETION:p6_03_complete',self.codes())
    def test_estimate_cannot_accept_tet4(self):
        self.s['original_mc_column']['estimate_accepts_tet4']=True;self.assertIn('MC_SEPARATION',self.codes())
    def test_old_rain_failure_cannot_be_deleted(self):
        self.s['old_rainfall_five'].pop();self.assertIn('RAIN_COUNT',self.codes())
    def test_rain_threshold_cannot_be_relaxed(self):
        self.s['old_rainfall_five'][0]['upper_bound']=.005;self.assertIn('RAIN_ORIGINAL:1',self.codes())
    def test_boolean_failure_cannot_be_replaced_by_zero(self):
        self.s['old_rainfall_five'][0]['passed']=0
        self.assertIn('RAIN_ORIGINAL:1',self.codes())
    def test_rain_strict_operator_preserved(self):
        self.s['old_rainfall_five'][0]['operator']='<=';self.assertIn('RAIN_ORIGINAL:1',self.codes())
    def test_rain_false_cannot_be_replaced(self):
        self.s['old_rainfall_five'][0]['passed']=True;self.assertIn('RAIN_ORIGINAL:1',self.codes())
    def test_missing_ledger_row_rejected(self):
        self.s['rows'][0]['ledger'].append('missing');self.assertIn('P601_ROW:missing',self.codes())
    def test_indeterminate_scope_not_broadened(self):
        self.s['rows'][1]['state']='限定受入';self.assertIn('ROW_STATE:extrusion-stress',self.codes())
    def test_pc_requirement_not_removed(self):
        self.s['rows'][1]['pc']='';self.assertIn('PC:extrusion-stress',self.codes())
    def test_duplicate_target_rejected(self):
        self.s['rows'][1]=copy.deepcopy(self.s['rows'][0]);self.assertIn('COVERAGE',self.codes())
    def test_historical_unexecuted_preserved(self):
        self.s['historical_python']['p601_unexecuted']=0;self.assertIn('HISTORICAL_COUNTS',self.codes())
    def test_ci_cannot_be_invented(self):
        self.s['windows_ci']['success_claimed']=True;self.assertIn('CI_NOT_CONFIRMED',self.codes())
    def test_pile_difference_is_not_exact_error(self):
        self.s['pile_metrics']['metric_is_exact_solution_error']=True
        self.assertIn('PILE_METRICS',self.codes())
    def test_legacy_monitor_cannot_be_replaced(self):
        self.s['reference_legacy_rows']=[]
        self.assertIn('REFERENCE_LEGACY_ROWS',self.codes())
    def test_exact_files_and_original_audit(self):self.assertTrue(audit.audit()['integrity_ok'],audit.audit()['errors'])
    def test_table_matches_json(self):self.assertEqual((HERE/'status.ja.md').read_text(),audit.table(self.s))
    def test_later_appendix_mutation_rejected_by_previous_audit(self):
        old=audit.load_old()
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);dest=root/'interim/2026-12-07-p6-03-cloud-reconciliation';dest.mkdir(parents=True)
            (dest/'document-appendices.json').write_text('{}\n')
            ledger={'base_commit':old.BASE,'sources':{},'rows':[],'evidence':{}}
            codes={r['code'] for r in old.audit(root,ledger,verify_git=False)['findings']}
            self.assertIn('LATER_APPENDIX_MANIFEST',codes)

if __name__=='__main__':unittest.main(verbosity=2)
