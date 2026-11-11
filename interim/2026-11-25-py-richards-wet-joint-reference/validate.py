"""Explicit validation; preserves all failures and unique IDs. Does not alter CI."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import unittest
import numpy
import scipy
from reference import write
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def flatten(suite):
    for item in suite:
        if isinstance(item,unittest.TestSuite):yield from flatten(item)
        else:yield item

class Result(unittest.TextTestResult):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.ids=[]
    def startTest(self,test):self.ids.append(test.id());super().startTest(test)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--existing',action='store_true');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    fixed={str(path.relative_to(HERE)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(HERE.glob('*.py'))}
    write(a.output/'source-sha256-before.json',fixed)
    suite=unittest.TestLoader().discover(str(HERE),pattern='checks.py')
    ids=[t.id() for t in flatten(suite)]
    with (a.output/'new-tests.log').open('w',encoding='utf-8') as stream:
        result=unittest.TextTestRunner(stream=stream,resultclass=Result,verbosity=2).run(suite)
    report=dict(environment=dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,longdouble_mantissa_bits=numpy.finfo(numpy.longdouble).nmant),
        new=dict(discovered=len(ids),discovered_ids=ids,executed_ids=result.ids,tests_run=result.testsRun,
          failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),
          expected_failures=len(result.expectedFailures),unexpected_successes=len(result.unexpectedSuccesses)),existing=[])
    success=result.wasSuccessful() and result.testsRun==len(ids) and not(result.skipped or result.expectedFailures or result.unexpectedSuccesses)
    if a.existing:
        for group in ('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'):
            folder=a.output/group
            with (a.output/(group+'.driver.log')).open('w',encoding='utf-8') as log:
                proc=subprocess.run([sys.executable,str(ROOT/'tools/run_ci_python_tests.py'),'--group',group,'--output',str(folder)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            report['existing'].append(dict(group=group,exit_code=proc.returncode))
            success=success and proc.returncode==0
    after={str(path.relative_to(HERE)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(HERE.glob('*.py'))}
    write(a.output/'source-sha256-after.json',after)
    report.update(new_source_unchanged=fixed==after,software_success=bool(success and fixed==after),
                  elapsed_seconds=time.monotonic()-started,accepted=None,A04='未達',C04='未達')
    write(a.output/'validation.json',report)
    return 0 if report['software_success'] else 1

if __name__=='__main__':sys.exit(main())
