"""Run repository Python regressions while respecting the no-NC-retry request."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse,json,sys,time,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def run(output):
    assert Path.cwd().resolve()==ROOT, 'run from repository root'
    assert not output.exists(), 'new report required'
    sys.path.insert(0,str(ROOT/'tools'));sys.path.insert(0,str(ROOT/'tests'))
    def flatten(suite):
        for t in suite:
            if isinstance(t,unittest.TestSuite):yield from flatten(t)
            else:yield t
    tests=list(flatten(unittest.TestLoader().discover('tests',pattern='test_*.py')))
    prohibited=lambda t:t.id().startswith('test_low_preload_nc.') or t.id()=='test_low_preload_nc_followup.LowNCFollowup.test_real_separate_process_restarts_inside_every_phase'
    excluded=[t.id() for t in tests if prohibited(t)]
    selected=[t for t in tests if not prohibited(t)]
    start=time.monotonic();r=unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(selected))
    report=dict(discovered=len(tests),tests_run=r.testsRun,selected_successful=r.wasSuccessful(),
        failures=[(str(t),e) for t,e in r.failures],errors=[(str(t),e) for t,e in r.errors],skipped=[(str(t),e) for t,e in r.skipped],
        excluded=excluded,all_existing_python_tests_executed=not excluded,
        exclusion_reason='The excluded class/method actually launches the low-preload NC path, which this task forbids repeating.',
        elapsed_seconds=time.monotonic()-start)
    output.write_text(json.dumps(report,indent=2)+'\n');return r.wasSuccessful()
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();sys.exit(not run(a.output))
