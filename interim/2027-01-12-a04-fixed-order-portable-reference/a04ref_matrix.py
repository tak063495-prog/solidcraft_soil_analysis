"""Explicit portable entry to the original a04cause_driver matrix (subprocesses included)."""
import argparse
from decimal import localcontext
import os
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key] = '1'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAUSE = ROOT / 'interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis'
sys.path.insert(0,str(CAUSE))
import a04cause_driver as driver
from a04ref_decimal import PRECISION, model_type
import numpy as np


def reroute(invoke):
    def wrapped(command,*args,**kwargs):
        if str(command[0]) == sys.executable:
            child = Path(command[1]).resolve()
            if child in (driver.OLD/'study.py',Path(driver.__file__)):
                command = [sys.executable,HERE/'a04ref_matrix.py','--force-binary64',
                    *(['--study'] if child.name=='study.py' else []),*command[2:]]
        return invoke(command,*args,**kwargs)
    return wrapped


def independent(folder,order):
    bound = driver.bind()
    soil = bound['soil_reference']
    # Library-local guard prevents accidental access to platform long double.
    class Binary64:
        longdouble = np.float64
        def __getattr__(self,name):
            return getattr(np,name)
    soil.np = soil.hyd.np = Binary64()
    mesh,ground,case = [driver.read(folder/(key+'.json')) for key in ('mesh','ground','case')]
    if order is not None:
        mesh['quadrature_order'] = order
    trace = []
    try:
        with localcontext() as ctx:
            ctx.prec = PRECISION
            model = model_type(soil.Model,trace)(mesh,ground,case)
            stream,sink = bound['study'].sink_file(folder/'independent-fields.jsonl')
            with stream:
                result = model.solve(sink)
    finally:
        driver.write(folder/'a04ref-line-search.json',dict(precision=PRECISION,
            numpy_state_mantissa=int(np.finfo(soil.np.longdouble).nmant),iterations=trace))
    result.update(diagnostic_hydraulic_order=mesh.get('quadrature_order',3),
                  input_only_no_native_solution_read=True,hydraulic_arithmetic='decimal40')
    driver.write(folder/'independent-result.json',result)


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--force-binary64',required=True,action='store_true')
    parser.add_argument('--study',action='store_true')
    options,rest = parser.parse_known_args()
    sys.argv[1:] = rest
    if options.study:
        study = driver.bind()['study']
        study.invoke = reroute(study.invoke)
        study.independent = lambda folder: independent(folder,None)
        study.main()
    else:
        driver.invoke = reroute(driver.invoke)
        driver.independent = independent
        driver.main()


if __name__=='__main__':
    main()
