"""Check the single expected identity difference between two real executables."""
import argparse
import copy
from pathlib import Path
from b01_cases import read, write
from study import exact, sha


def run(study, controls, output):
    observer = next((study/'solver').iterdir())
    cli = next((controls/'solver').iterdir())
    rows = []
    for r in read(controls/'cli-controls.json')['records']:
        a = read(study/r['case']/'native/run.json')
        b = read(controls/r['case']/'run.json')
        assert a['provenance']['solver_binary_sha256'] == sha(observer)
        assert b['provenance']['solver_binary_sha256'] == sha(cli)
        original_equal = exact(a, b)
        identity_bound = copy.deepcopy(a)
        identity_bound['provenance']['solver_binary_sha256'] = b['provenance']['solver_binary_sha256']
        only_binary_differs = exact(identity_bound, b)
        rows.append(dict(case=r['case'], original_exact_except_elapsed=original_equal,
                         observer_sha256=sha(observer), production_cli_sha256=sha(cli),
                         all_other_run_fields_exact_except_elapsed=only_binary_differs,
                         original_runs_unchanged=True))
        assert not original_equal and only_binary_differs
    write(output, dict(passed=True, physics_acceptance=None, records=rows,
                       identity_difference='provenance.solver_binary_sha256 only; both actual file hashes verified',
                       restart_comparisons_exclude_only_elapsed_seconds=True,
                       scope='Cross-executable binding audit, not a restart comparison or a new solve'))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('study', 'controls', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    run(a.study, a.controls, a.output)
