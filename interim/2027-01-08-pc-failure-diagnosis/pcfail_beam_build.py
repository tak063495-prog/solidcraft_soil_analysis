"""Generate an isolated observation build; retain all original numerical code."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
ROOT=Path(__file__).resolve().parents[2]


def generate(output):
    output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(ROOT/'src', output/'src')
    for name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):
        shutil.copy2(ROOT/name, output/name)
    path=output/'src/beams.rs'; original=path.read_text(encoding='utf-8')
    marker='''        if !relative.is_finite() || relative > residual_tolerance {
            return Err(Error::Solve(
                "beam interior true residual exceeded tolerance".into(),
'''
    observation='''        if !relative.is_finite() || relative > residual_tolerance {
            // Explicit work-only observation, numerical result and error unchanged.
            if let Some(path) = std::env::var_os("PCFAIL_BEAM_TRACE") {
                if let Ok(file) = std::fs::OpenOptions::new().write(true).create_new(true).open(path) {
                    let value = serde_json::json!({
                        "factor": factor, "relative": relative, "tolerance": residual_tolerance,
                        "rows": c.kqq.nrows(), "column_major_kqq": c.kqq.as_slice(),
                        "rhs": rhs.as_slice(), "interior": interior.as_slice()
                    });
                    let _ = serde_json::to_writer(file, &value);
                }
            }
            return Err(Error::Solve(
                "beam interior true residual exceeded tolerance".into(),
'''
    if original.count(marker)!=1: raise ValueError('unreviewed source change')
    changed=original.replace(marker,observation)
    path.write_text(changed, encoding='utf-8', newline='\n')
    (output/'observation.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile='src/beams.rs',tofile='work-only/src/beams.rs')),encoding='utf-8')
    binding=dict(diagnostic_only=True, production_unchanged=True,
        original_sha256=hashlib.sha256(original.encode()).hexdigest(),
        observed_sha256=hashlib.sha256(changed.encode()).hexdigest())
    (output/'binding.json').write_text(json.dumps(binding,indent=2)+'\n',encoding='utf-8')
    return binding


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    generate(a.output.resolve())
    subprocess.run(['cargo','build','--locked','--release','--bin','soil-fem-3d','--manifest-path',str(a.output.resolve()/'Cargo.toml')],cwd=ROOT,check=True)
