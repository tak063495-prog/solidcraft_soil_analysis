"""Fresh standard solve-stages controls, distinct from strict same-binary restart checks."""
import argparse
import copy
import os
from pathlib import Path
from study import read,write,invoke,sha,exact,input_paths


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    binary=(a.run/'bin'/('solver.exe' if os.name=='nt' else 'solver')).resolve();rows={}
    for name in ('transfer_force','stress_free','stress_free-cycle'):
        folder=a.run/name;out=a.output/name;out.mkdir()
        result=invoke([binary,'solve-stages',*input_paths(folder),out/'run.json'],out/'solve.log')
        assert result['exit_code']==0 and not result['timed_out']
        raw=read(out/'run.json');expected=read(folder/'native/run.json');assert raw['complete']
        actual=copy.deepcopy(raw);reference=copy.deepcopy(expected)
        assert actual['provenance']['solver_binary_sha256']==sha(binary)
        # Different executable wrappers necessarily have different binary SHAs.
        # This cross-wrapper control is not a restart-equality test.
        for value in (actual,reference):value['provenance'].pop('solver_binary_sha256')
        assert exact(actual)==exact(reference)
        rows[name]=dict(fresh_standard_cli_solve=True,accepted=True,process=result,solver_sha256=sha(binary),
            inputs_sha256={p.name:sha(p) for p in input_paths(folder)},
            cross_wrapper_run_equal=True,cross_wrapper_only_exclusions=['elapsed_seconds','provenance.solver_binary_sha256'],
            not_counted_as_additional_independent_case=True)
    write(a.output/'verification.json',rows)


if __name__=='__main__':main()
