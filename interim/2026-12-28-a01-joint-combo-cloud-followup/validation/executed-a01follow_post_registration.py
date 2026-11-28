import json,time,unittest
from pathlib import Path
root=Path.cwd();out=root/'interim/2026-12-28-a01-joint-combo-cloud-followup/validation'
patterns=['test_a01_joint_combo_cloud_followup.py','test_a01_joint_combo_cause_diagnosis.py','test_p6_status_reconciliation.py','test_ci_python_partition.py','test_ci_rust_partition.py']
suite=unittest.TestSuite(unittest.TestLoader().discover(str(root/'tests'),pattern=p) for p in patterns)
def flatten(s):
 for t in s:
  if isinstance(t,unittest.TestSuite):yield from flatten(t)
  else:yield t.id()
ids=list(flatten(suite));start=time.monotonic()
with (out/'post-registration-tests.log').open('x',encoding='utf-8') as log:result=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
report=dict(discovered=len(ids),executed=result.testsRun,passed=result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped)-len(result.expectedFailures)-len(result.unexpectedSuccesses),failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),expected_failures=len(result.expectedFailures),unexpected_successes=len(result.unexpectedSuccesses),unexecuted=len(ids)-result.testsRun,test_ids=ids,elapsed_seconds=time.monotonic()-start,not_added_to_full_unique_counts=True)
(out/'post-registration.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items() if k!='test_ids'})
assert result.wasSuccessful() and not any(report[k] for k in ('failures','errors','skipped','expected_failures','unexpected_successes','unexecuted'))
