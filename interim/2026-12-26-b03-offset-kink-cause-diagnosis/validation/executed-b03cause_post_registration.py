import json,pathlib,sys,time,unittest
root=pathlib.Path.cwd();out=root/'interim/2026-12-26-b03-offset-kink-cause-diagnosis/validation';patterns=['test_p6*.py','test_nc_physical_appendix.py','test_ci_python_partition.py','test_ci_rust_partition.py','test_b03_offset_kink_cause_diagnosis.py'];tests=[]
def flat(s):
 for x in s:
  if isinstance(x,unittest.TestSuite):yield from flat(x)
  else:yield x
for pattern in patterns:tests.extend(flat(unittest.TestLoader().discover(str(root/'tests'),pattern=pattern)))
ids=[t.id() for t in tests];assert len(ids)==len(set(ids));start=time.monotonic()
with (out/'post-registration-tests.log').open('w',encoding='utf-8') as log:result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.TestSuite(tests))
report=dict(success=result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped),tests_run=result.testsRun,failed=len(result.failures),errors=len(result.errors),skip=len(result.skipped),ids=ids,elapsed_seconds=time.monotonic()-start,not_added_to_unique_full_suite_counts=True)
(out/'post-registration.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='ids'}));sys.exit(0 if result.wasSuccessful() and result.testsRun==len(ids) and not result.skipped else 1)
