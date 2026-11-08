"""Read saved evidence only. No solver or independent path is rerun."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import diagnose_cut_volume_history as cut
from probe_low_preload_nc import NAMES, deterministic, inspect_run, read, sha, snapshot, write
from probe_nc_midpoint import inputs, restart_observation


def measure(evidence):
    manifest = read(Path(__file__).parent/'archive-sha256.json')
    for name, digest in manifest['entries'].items():
        assert sha(evidence/name) == digest, ('evidence changed', name)
    rows = {}
    for key in ('nc-resumed', 'nc-retried'):
        folder = evidence/key; report = read(folder/'experiment.json'); job = folder/'job'
        values = [read(folder/'input'/(n+'.json')) for n in NAMES]
        assert values == list(inputs('midpoint')), '5 kPa midpoint input mismatch'
        binding = read(job/'job.json')['binding']
        for name, digest in binding['files_sha256'].items():
            if name == 'solver/soil-fem-3d.exe':
                assert digest == report['solver_sha256']
            else:
                assert sha(job/name) == digest, ('sealed input changed', name)
        _, native = snapshot(job)
        assert native == read(folder/'observed-run.json')
        assert inspect_run(*values, native) == report['observation']
        imported = read(job/'inputs/checkpoint.json')['run']
        restart = restart_observation(imported, native)
        restart['imported_checkpoint_sha256'] = sha(job/'inputs/checkpoint.json')
        assert restart == report['restart_observation']
        assert not report['passed'] and not native['complete']
        rows[key] = dict(observation=report['observation'], restart=restart, failure_reason=native['failure_reason'])
    a = read(evidence/'nc-resumed/observed-run.json'); b = read(evidence/'nc-retried/observed-run.json')
    assert deterministic(a) == deterministic(b), 'failed terminal restart differs'
    result = cut.measure(evidence/'cut-history')
    assert result == read(evidence/'cut-history/measurements.json'), 'cut-history remeasurement differs'
    return dict(schema='soil-fem-nc-resume-cut-history-remeasurement/1', passed=False,
                solver_rerun=False, independent_path_rerun=False, nc_saved_field_audits=rows,
                failed_terminal_restart_entire_run_exact_except_elapsed_seconds=True,
                cut_measurements_reproduced=True, history_converged=False,
                production_adopted=False, full_removal_reached=False, windows_ci='unverified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    write(args.output, measure(args.evidence.resolve()))
    print('Saved-field and history measurements reproduced; passed=false')
