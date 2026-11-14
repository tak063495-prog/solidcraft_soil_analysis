"""Build a separate read-only observer for rejected native dense LU systems.

The production source and solve result are unchanged. The isolated binary has
its own checkpoint identity: never use it to resume another binary's checkpoint.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HOOK = '''                if residual.is_finite() && residual > settings.relative_linear_tolerance {
                    crate::nc_linear_observer::save(a, rhs, &x, residual)?;
                }
'''
OBSERVER = '''use crate::Result;
use nalgebra::{DMatrix, DVector};
use std::io::Write;
pub(crate) fn save(a: &DMatrix<f64>, b: &DVector<f64>, x: &DVector<f64>, residual: f64) -> Result<()> {
    if let Some(path) = std::env::var_os("NC_LINEAR_SYSTEMS") {
        let matrix: Vec<Vec<f64>> = (0..a.nrows()).map(|i| a.row(i).iter().copied().collect()).collect();
        let mut file = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
        serde_json::to_writer(&mut file, &serde_json::json!({"matrix": matrix, "rhs": b.as_slice(),
            "solution": x.as_slice(), "native_true_relative_residual": residual,
            "rejected": true, "numerical_acceptance": false}))?;
        file.write_all(b"\\n")?;
    }
    Ok(())
}
'''


def instrument(source):
    marker = '                let residual = (a * &x - rhs).norm() / rhs.norm().max(1e-30);\n'
    if source.count(marker) != 1:
        raise ValueError('native dense solve observer insertion point differs')
    return source.replace(marker, marker + HOOK)


def run(output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    crate = output/'crate'
    crate.mkdir()
    source = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'src').rglob('*.rs')}
    shutil.copytree(ROOT/'src', crate/'src')
    for name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
        shutil.copyfile(ROOT/name, crate/name)
    p = crate/'src/stages.rs'
    p.write_text(instrument(p.read_text(encoding='utf-8')), encoding='utf-8')
    with (crate/'src/lib.rs').open('a', encoding='utf-8') as stream:
        stream.write('\nmod nc_linear_observer;\n')
    (crate/'src/nc_linear_observer.rs').write_text(OBSERVER, encoding='utf-8')
    (crate/'examples').mkdir()
    for name in ('nc_history_trace', 'nc_dense_lu_diagnostic'):
        shutil.copyfile(ROOT/f'examples/{name}.rs', crate/f'examples/{name}.rs')
    subprocess.run(['cargo', 'fmt', '--manifest-path', str(crate/'Cargo.toml')], check=True)
    subprocess.run(['cargo', 'build', '--locked', '--release', '--manifest-path', str(crate/'Cargo.toml'),
                    '--target-dir', str(output/'target'), '--examples'], check=True)
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != sha for p, sha in source.items()):
        raise ValueError('production source changed during observer preparation')
    suffix = '.exe' if sys.platform == 'win32' else ''
    binaries = {name: output/f'target/release/examples/{name}{suffix}'
                for name in ('nc_history_trace', 'nc_dense_lu_diagnostic')}
    record = dict(root_numerical_source_sha256=source, observer_only=True, production_modified=False,
                  equations_tolerances_and_solve_result_unchanged=True,
                  binaries={k: dict(path=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest())
                            for k, v in binaries.items()})
    (output/'observer.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    print(json.dumps(run(parser.parse_args().output), indent=2))
