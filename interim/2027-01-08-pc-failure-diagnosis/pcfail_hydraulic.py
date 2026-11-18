"""Explicit input-only A-04 precision observation. Original policy is unchanged."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'interim/2026-11-29-richards-wet-joint-a04'
HYD = ROOT/'interim/2026-11-22-wet-joint-transient-water/hydraulic_reference.py'


class Precision:
    """Model Windows NumPy longdouble==float64 without mutating NumPy globally."""
    def __init__(self, dtype):
        self.longdouble = dtype

    def __getattr__(self, key):
        return getattr(np, key)


def load(dtype):
    spec = importlib.util.spec_from_file_location('pcfail_a04_binding',
        ROOT/'interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis/a04cause_binding.py')
    binding = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(binding)
    bound = binding.bind()
    module = bound['soil_reference']
    module.pcfail_fixture = bound['study'].fixture
    module.np = Precision(dtype)
    module.hyd.np = Precision(dtype)
    return module


def run(output, steps=4, precision='host'):
    output.mkdir(parents=True, exist_ok=False)
    # Fixture import supplies inputs only, never native output.
    dtype = np.float64 if precision in ('binary64','state64-flux80') else np.longdouble
    flux_dtype = np.float64 if precision in ('binary64','state80-flux64') else np.longdouble
    module = load(dtype)
    module.hyd.np = Precision(flux_dtype)
    values = module.pcfail_fixture(steps, False)
    for key, value in zip(('mesh', 'ground', 'case'), values):
        (output/(key+'.json')).write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    model = module.Model(*values)
    iterations = []
    merit_line = next(i for i,line in enumerate(HYD.read_text(encoding='utf-8').splitlines(),1) if line.strip() == 'if merit < .01:')
    failed = {}
    snapshots = []

    def trace(frame, event, arg):
        if frame.f_code.co_filename != str(HYD) or frame.f_code.co_name != 'step':
            return None
        loc = frame.f_locals
        if event == 'line' and 'merit' in loc:
            # Record each iteration once, after merit has been computed.
            line = frame.f_lineno
            if line == merit_line:
                iterations.append(dict(iteration=loc['iteration'], dt=loc['dt'], period=loc['period']['name'],
                    merit=loc['merit'], residual_dtype=str(loc['r'].dtype), x_dtype=str(loc['x'].dtype)))
        if event == 'exception' and arg[0] is ArithmeticError and 'hydraulic line search' in str(arg[1]):
            for key in ('x','r','scale','jac','dx','h','trace','old_storage'):
                failed[key] = np.array(loc[key])
            base = loc['x']; evaluate = loc['evaluate']; dx = loc['dx']
            for cut in range(20):
                candidate = base+dx*2.**(-cut)
                rr = evaluate(candidate)[0]
                snapshots.append(dict(cut=cut, native_l2=float(np.linalg.norm(rr)),
                    changed_components=int(np.count_nonzero(candidate != base)),
                    scaled_residual=float(np.max(np.abs(rr)/(1e-16+1e-8*evaluate(candidate)[1])))))
            np.savez(output/'failed-state.npz', **failed)
        return trace

    rows = []
    previous = sys.gettrace()
    try:
        sys.settrace(trace)
        result = model.solve(lambda row, fields: rows.append(row))
        outcome = dict(execution_complete=True, accepted_rows=len(result['rows']))
    except ArithmeticError as exc:
        outcome = dict(execution_complete=False, error=repr(exc), accepted_rows=len(rows))
    finally:
        sys.settrace(previous)
    report = dict(schema='pcfail-hydraulic/1', passed=False, a04_complete=False, diagnostic_only=True,
        precision=precision, mantissa_bits=np.finfo(dtype).nmant, flux_mantissa_bits=np.finfo(flux_dtype).nmant, steps=steps,
        original_policy_unchanged=True, input_only=True, iterations=iterations,
        failed_candidates=snapshots, **outcome)
    (output/'observation.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', type=int, choices=(4,8,16), default=4)
    p.add_argument('--precision', choices=('host','binary64','state64-flux80','state80-flux64'), default='host')
    a=p.parse_args()
    print(json.dumps(run(a.output.resolve(), a.steps, a.precision)))
