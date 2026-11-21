"""Authenticate a fresh extraction and remeasure; no new solve or resealing."""
import argparse
import json
from pathlib import Path
from slip_study import HERE, ROOT, sha, write, measure, inputs, load_module
from history_analysis import measure as history_measure


def remeasure(folder):
    from verify_wet_structures import without_time
    from verify_wet_joint_cycle import path_evidence
    from wet_joint_cycle_reference import serializable
    report = serializable(measure(folder))
    native = folder/'native-original'
    run = json.loads((native/'stage-results.json').read_text())
    checks = json.loads((native/'checks.json').read_text())
    assert checks['frames'] == path_evidence(run)
    prod = load_module('production_audit_remeasure', ROOT/'interim/2026-11-08-wet-joint-production-cut-volume/study.py')
    assert checks['accepted_gates'] == prod.native_gates(run)
    fields, audits = prod.fields_and_audit(*inputs()[:2], run, None)
    assert serializable(fields) == checks['fields'] and audits == checks['independent_force_audit']
    assert prod.holding_audit(*inputs()[:2], run) == checks['holding_audit']
    assert len(checks['restarts']) == 7
    for row in checks['restarts']:
        name = row['stage']
        checkpoint = native/f'checkpoint-{name}.json'
        assert sha(checkpoint) == row['checkpoint_sha256']
        resumed = json.loads((native/f'resumed-{name}.json').read_text())
        assert 0 < row['progress'] < 1
        assert without_time(resumed) == without_time(run)
    histories = history_measure(folder)
    assert report == json.loads((folder/'measurements.json').read_text())
    assert histories == json.loads((folder/'history-measurements.json').read_text())
    return dict(measurements_exact=True, history_measurements_exact=True,
        whole_native_run_restart_matches=7, checkpoint_sha256_matches=7,
        new_solves_executed=False, resealed_checkpoints=False,
        passed=False, history_converged=False)


def verify_archives(destination):
    manifest = json.loads((HERE/'archive-sha256.json').read_text())
    destination.mkdir(parents=True, exist_ok=False)
    import zipfile
    import hashlib
    seen = {}
    for row in manifest['archives']:
        path = HERE/row['file']
        assert sha(path) == row['sha256']
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                assert not Path(name).is_absolute() and '..' not in Path(name).parts
                assert name not in seen
                raw = archive.read(name)
                digest = hashlib.sha256(raw).hexdigest()
                assert manifest['entries'][name] == digest
                seen[name] = digest
                target = destination/name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(raw)
    assert seen == manifest['entries']
    result = remeasure(destination/'raw')
    result.update(archives_verified=len(manifest['archives']), entries_verified=len(seen))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extract', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    write(args.output, verify_archives(args.extract.resolve()))
