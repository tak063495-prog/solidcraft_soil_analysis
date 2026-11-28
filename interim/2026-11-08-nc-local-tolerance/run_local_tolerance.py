"""Reproduce a bounded, unadopted local Newton tolerance experiment.

Always builds from the specified PR tip, never edits the production checkout or
reseals a checkpoint. Process success is not numerical acceptance.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

BASE = 'c4a2b39d7f420036a8903163a653a6ed943e718c'
ROOT = next(p for p in Path(__file__).resolve().parents
            if (p/'Cargo.toml').exists() and (p/'tools').is_dir())
sys.path.insert(0, str(ROOT/'tools'))
from probe_low_preload_nc import NAMES, expected_inputs, inspect_run, read, sha, write


def prepare(source, tolerance):
    source.mkdir(parents=True, exist_ok=False)
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE,
                                     'src', 'examples', 'Cargo.toml', 'Cargo.lock',
                                     'rust-toolchain.toml'], cwd=ROOT, text=True).splitlines()
    for name in names:
        path = source/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT))
    path = source/'src/mcc.rs'
    data = path.read_bytes()
    old = b'if r.norm() < 1e-11 {'
    assert data.count(old) == 1
    path.write_bytes(data.replace(old, ('if r.norm() < '+tolerance+' {').encode()))
    return {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*')) if p.is_file()}


def run(output, tolerance):
    # Input generation must still be the one used for this controlled study.
    for name in ('tools/probe_low_preload_nc.py', 'tools/make_free_contact.py',
                 'tools/fem_reference.py', 'examples/joint-column-tet10-mesh.json',
                 'examples/ground-joint-mcc.json', 'examples/stages-joint-explicit-normal-column.json'):
        committed = subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT)
        # Windows core.autocrlf may change only checkout line endings.
        assert (ROOT/name).read_bytes().replace(b'\r\n', b'\n') == committed.replace(
            b'\r\n', b'\n'), 'input generator/source changed: '+name
    output.mkdir(parents=True, exist_ok=False)
    source = output/'source'
    sources = prepare(source, tolerance)
    inputs = output/'input'; inputs.mkdir()
    values = expected_inputs(excavation_increment=1/128, minimum_increment=1/4194304,
                             dense_lu=True, tangent_perturbation=1e-9)
    for name, value in zip(NAMES, values):
        write(inputs/(name+'.json'), value)
    metadata = dict(base=BASE, local_newton_tolerance=tolerance,
                    unchanged_elastic_yield_branch_tolerance=1e-11,
                    relative_force_tolerance=1e-8, relative_linear_tolerance=1e-10,
                    source_sha256=sources, input_sha256={p.name: sha(p) for p in inputs.iterdir()},
                    passed=False, production_adopted=False)
    env = dict(os.environ)
    # Keep this build isolated even if the caller has a shared target override.
    env['CARGO_TARGET_DIR'] = str(source/'target')
    start = time.monotonic()
    with (output/'build.log').open('x') as log:
        build = subprocess.run(['cargo', 'build', '--locked', '--release', '--example', 'nc_newton_trace'],
                               cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
    metadata.update(build_exit_code=build.returncode, build_seconds=time.monotonic()-start)
    if build.returncode:
        write(output/'experiment.json', metadata)
        raise RuntimeError('build failed; see build.log')
    exe = source/'target/release/examples'/('nc_newton_trace.exe' if os.name == 'nt' else 'nc_newton_trace')
    metadata['executable_sha256'] = sha(exe)
    start = time.monotonic()
    with (output/'solver.log').open('x') as log:
        try:
            result = subprocess.run([str(exe), str(inputs), str(output/'result')], cwd=source,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=205)
            metadata['solver_exit_code'] = result.returncode
        except subprocess.TimeoutExpired:
            metadata['external_timeout_seconds'] = 205
    metadata['solver_seconds'] = time.monotonic()-start
    run_path = output/'result/run.json'
    if run_path.exists():
        native = read(run_path)
        metadata.update(complete=native['complete'], failure_reason=native['failure_reason'],
                        last_attempt=native['attempts'][-1],
                        observation=inspect_run(*values, native))
    write(output/'experiment.json', metadata)
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tolerance', choices=('1e-13', '1e-14'), default='1e-13')
    args = parser.parse_args()
    print(json.dumps(run(args.output.resolve(), args.tolerance), ensure_ascii=False, indent=2))
