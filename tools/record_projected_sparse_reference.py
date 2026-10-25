"""Record bounded sparse-reference verification; never mark spatial/P6 complete."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import zipfile
import hashlib
import freeze_sources
from mcc_projected_volume_diagnostic import read, sha
import mcc_projected_sparse_reference as sparse


def verify_archive(folder):
    expected = read(folder/'source-snapshot.json')
    with zipfile.ZipFile(folder/'source-snapshot.zip') as archive:
        assert archive.testzip() is None
        assert set(archive.namelist()) == set(expected)
        for name, digest in expected.items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == digest, name
    return len(expected)


def record(root, output):
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    rate_folder = root/'runs/projected-sparse-rate-equivalence-20261025-v2'
    cycle_folder = root/'runs/projected-sparse-cycle-20261025-v1'
    inputs_folder = root/'runs/projected-mcc-spatial-inputs-20261025-v1'
    rate = read(rate_folder/'rate-equivalence.json')
    cycle = read(cycle_folder/'backend-equivalence.json')
    progress = read(cycle_folder/'progress.json')
    assert rate['passed'] and len(rate['cases']) == 12
    assert cycle['passed'] and len(cycle['full_field_checks']) == 38
    assert progress['status'] == 'completed' and all(p['exit_code'] == 0 for p in progress['phases'])
    for name in ('mcc_projected_sparse_reference.py', 'mcc_projected_volume_diagnostic.py',
                 'mcc_projected_cycle_diagnostic.py', 'mcc_rate_fem.py', 'mcc_rate_midpoint.py',
                 'fem_reference.py', 'verify_high_order.py', 'high_order_io.py'):
        assert sha(root/'tools'/name) == sha(cycle_folder/'reference-sources'/name), name
        assert sha(root/'tools'/name) == rate['source_sha256']['tools/'+name], name
    archives = {str(p): verify_archive(p) for p in (rate_folder, cycle_folder, inputs_folder)}
    for result in (read(cycle_folder/'initial.json'), read(cycle_folder/'cycle.json')):
        for name, digest in result['input_and_source_sha256'].items():
            assert sha(Path(name)) == digest, name
    command = [sys.executable, '-X', 'utf8', str(root/'tools/mcc_projected_sparse_run.py'), 'initial',
               *[str(cycle_folder/'inputs'/n) for n in ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json')],
               str(output/'wrapper-smoke.json'), '--fraction', '0.000001', '--steps', '2',
               '--dependency-lock', str(root/'requirements-reference.lock.txt')]
    smoke = subprocess.run(command, capture_output=True, text=True, encoding='utf-8')
    (output/'wrapper-smoke.log').write_text(smoke.stdout+smoke.stderr, encoding='utf-8')
    assert smoke.returncode == 0, smoke.stderr
    assert read(output/'wrapper-smoke.json')['load_fraction'] == 1e-6
    wrong = output/'wrong-environment'
    wrong.mkdir()
    (wrong/'requirements.lock.txt').write_text((root/'requirements.lock.txt').read_text(), encoding='utf-8')
    (wrong/'requirements-reference.lock.txt').write_text('-r requirements.lock.txt\nscipy==0.0.0\n', encoding='utf-8')
    command[command.index(str(output/'wrapper-smoke.json'))] = str(output/'must-not-exist.json')
    command[-1] = str(wrong/'requirements-reference.lock.txt')
    refusal = subprocess.run(command, capture_output=True, text=True, encoding='utf-8')
    (output/'wrong-environment.log').write_text(refusal.stdout+refusal.stderr, encoding='utf-8')
    assert refusal.returncode == 2 and 'versions differ' in refusal.stderr
    assert not (output/'must-not-exist.json').exists()
    build = read(root/'runs/projected-mcc-development-20261025-v1/build-sources/source-snapshot.json')
    rust = {k: v for k, v in build.items() if k.startswith(('src/', 'tests/', 'examples/')) or k in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml')}
    assert all(sha(root/k) == v for k, v in rust.items())
    evidence = [rate_folder/'rate-equivalence.json', cycle_folder/'backend-equivalence.json',
                cycle_folder/'initial.json', cycle_folder/'cycle.json', cycle_folder/'progress.json',
                inputs_folder/'preparation.json', inputs_folder/'area-measures-coarse.json',
                root/'runs/projected-sparse-reference-20261025-v1/install-report.json',
                output/'wrapper-smoke.json', output/'wrapper-smoke.log', output/'wrong-environment.log']
    report = dict(schema='soil-fem-sparse-reference-verification-record/1', passed=True,
                  goal_complete=False, p6_complete=False, production_backend_accepted=False,
                  spatial_convergence_accepted=False, rate_cases=12, nonlinear_field_checks=38,
                  source_archives_verified=archives, unchanged_rust_build_and_examples=len(rust),
                  scipy=sparse.scipy.__version__, numpy=sparse.np.__version__, python=sys.version,
                  current_wrapper_fractional_smoke_passed=True, wrong_dependency_version_refused=True,
                  evidence_sha256={str(p): sha(p) for p in evidence},
                  limitations='Reference backend verification and spatial inputs only; fine-mesh runs, integrator and spatial convergence remain pending')
    assert read(output/'source-snapshot.json') == freeze_sources.manifest()
    with (output/'verification.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(passed=True, rate_cases=12, nonlinear_field_checks=38, unchanged_rust_build_and_examples=len(rust))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record(Path(__file__).resolve().parents[1], args.output.resolve())
