"""Rerun affected original IDs after environment/fixture recovery, without edits."""
import argparse
import json
from pathlib import Path
import sys
import time
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from run_ci_python_tests import discover

class Recorded(unittest.TextTestResult):
    def __init__(self,*a,**kw):super().__init__(*a,**kw);self.outcomes={}
    def addSuccess(self,test):super().addSuccess(test);self.outcomes[test.id()]='passed'
    def addFailure(self,test,err):super().addFailure(test,err);self.outcomes[test.id()]='failed'
    def addError(self,test,err):super().addError(test,err);self.outcomes[test.id()]='error'
    def addSkip(self,test,reason):super().addSkip(test,reason);self.outcomes[test.id()]='skipped'

def run(output):
    output.mkdir(parents=True,exist_ok=False);all_tests,partitions=discover(ROOT)
    selected=sorted(k for k in all_tests if k.startswith(('p603_reconciliation_cases.P603AuditTests.','p6_status_reconciliation_cases.FixedRepositoryTests.')) or k.endswith(('.test_powershell_gate_when_available','.test_real_rerun_saves_specific_failure_and_preserves_pc3')))
    assert len(selected)==23,selected
    started=time.monotonic()
    with (output/'python-tests.log').open('x',encoding='utf-8') as log:
        result=unittest.TextTestRunner(stream=log,verbosity=2,resultclass=Recorded).run(unittest.TestSuite(all_tests[k] for k in selected))
    report=dict(executed_test_ids=selected,outcomes=result.outcomes,tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),elapsed_seconds=time.monotonic()-started,
                discovered_total=len(all_tests),original_tests_unchanged=True,reason='restore exact inherited JSON fixture; install real PowerShell to execute previously skipped bodies',exit_code=int(not result.wasSuccessful() or bool(result.skipped)))
    (output/'retry.json').write_text(json.dumps(report,indent=2)+'\n');return report['exit_code']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();sys.exit(run(a.output.resolve()))
