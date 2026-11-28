"""Recompute measurements from newly extracted evidence; no new solve/CI claim."""
import argparse
import copy
import json
from pathlib import Path
import sys

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'Cargo.toml').exists())
sys.path.insert(0, str(ROOT/'tools'))
from diagnose_wet_joint_convergence import measure
from probe_low_preload_nc import inspect_run
from verify_wet_structures import read, sha, without_time
import fem_reference as fem
from verify_joint_water import audit


def remeasure(evidence):
    wet = evidence/'separated'
    restarts = read(wet/'restarts.json'); run = read(wet/'g256-m5e-06/stage-results.json')
    assert len(restarts) == 7
    for row in restarts:
        folder = wet/'g256-m5e-06'
        assert sha(folder/f"checkpoint-{row['stage']}.json") == row['checkpoint_sha256']
        assert without_time(read(folder/f"resumed-{row['stage']}.json")) == without_time(run)
    result = dict(schema='soil-fem-wet-history-and-nc-followup/1', passed=False,
                  starting_commit='a4f515636292be61069f5b2abc4279c1263d06be',
                  wet=measure(wet, 64, restarts), nc={})
    output_audits = read(wet/'output-audits.json')
    for key, expected_audit in output_audits.items():
        folder = wet/key
        assert audit(read(folder/'mesh.json'), fem.normalize_ground(read(folder/'ground.json')),
                     read(folder/'stage-results.json'), folder/'results') == expected_audit
        binding = read(folder/'results/export-verification.json')
        assert binding['result_sha256'] == sha(folder/'stage-results.json')
        assert binding['mesh_sha256'] == sha(folder/'mesh.json')
        assert binding['vtu_roundtrip'] == binding['joint_point_roundtrip'] == 'passed'
    result['wet_output_audits'] = output_audits
    for name in ('midpoint', 'backward'):
        folder = evidence/('nc-'+name); report = read(folder/'experiment.json')
        native = read(folder/'observed-run.json')
        values = [read(folder/'input'/(n+'.json')) for n in ('mesh', 'ground', 'stages')]
        assert report['observation'] == inspect_run(*values, native)
        assert not report['passed'] and not report['production_adopted']
        assert all(sha(folder/'input'/n) == h for n, h in report['input_sha256'].items())
        result['nc'][name] = copy.deepcopy(report)
        result['nc'][name]['termination_classification'] = ('wall_time_budget' if report['stop_requested']
            else 'completed' if report['complete'] else 'numerical_stop')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = remeasure(args.evidence)
    expected = Path(__file__).with_name('measurements.json')
    if expected.exists():
        assert report == read(expected), 'saved measurement mismatch'
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False)
    print('All seven-frame comparisons, saved whole-Run restarts, and NC stored-state audits reproduced; passed=false.')
