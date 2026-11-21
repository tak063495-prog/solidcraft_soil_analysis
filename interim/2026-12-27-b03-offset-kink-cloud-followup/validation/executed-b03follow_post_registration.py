import json
from pathlib import Path
import sys
import time
import unittest

root=Path(__file__).resolve().parents[1]
output=root/'interim/2026-12-27-b03-offset-kink-cloud-followup/validation'
patterns=['test_b03_offset_kink_cloud_followup.py','test_b03_offset_kink_cause_diagnosis.py',
          'test_p6_status_reconciliation.py','test_ci_python_partition.py','test_ci_rust_partition.py']
suite=unittest.TestSuite(unittest.TestLoader().discover(str(root/'tests'),pattern=name) for name in patterns)
def flatten(suite):
    for t in suite:
        if isinstance(t,unittest.TestSuite):yield from flatten(t)
        else:yield t.id()
ids=list(flatten(suite));started=time.monotonic()
with (output/'post-registration-tests.log').open('w',encoding='utf-8') as stream:
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
report=dict(discovered=len(ids),test_ids=ids,executed=result.testsRun,
            passed=result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped),
            failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),
            expected_failures=len(result.expectedFailures),unexpected_successes=len(result.unexpectedSuccesses),
            unexecuted=len(ids)-result.testsRun,elapsed_seconds=time.monotonic()-started,
            included_in_full_unique_counts_not_added_again=True,
            full_unique_rust_passed=467,full_unique_python_passed=859)
(output/'post-registration.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in report.items() if k!='test_ids'}),flush=True)
if not result.wasSuccessful() or result.skipped or result.expectedFailures or result.unexpectedSuccesses or report['unexecuted']:
    sys.exit(1)
