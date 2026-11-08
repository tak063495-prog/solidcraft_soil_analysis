"""Prepare isolated saved-matrix / single-increment experiments from the verified head."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
BASE = 'eaeb160b68a1ec982634f45def65fc47bb279b45'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def prepare(out, candidate):
    out.mkdir(parents=True, exist_ok=False)
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE], cwd=ROOT, text=True).splitlines()
    hashes = {}
    for name in names:
        if name.startswith(('src/', 'examples/')) or name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
            data = subprocess.check_output(['git', 'show', f'{BASE}:{name}'], cwd=ROOT)
            dest = out/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            hashes[name] = hashlib.sha256(data).hexdigest()
    path = out/'src/stages.rs'
    s = path.read_text()
    old = '"internal": eval.internal, "target": target,'
    new = '''"joint_state": eval.joints.as_ref().map(|j| &j.state),
                    "joint_report": eval.joints.as_ref().map(|j| &j.report),
                    "joint_internal": eval.joints.as_ref().map(|j| &j.internal_kn),
                    "dof_rows": dofs.rows,
                    "internal": eval.internal, "target": target,'''
    assert s.count(old) == 1
    s = s.replace(old, new)
    old = '"alpha": alpha, "displacement_m": candidate, "states": trial.states,'
    new = '''"alpha": alpha, "displacement_m": candidate, "states": trial.states,
                                "joint_state": trial.joints.as_ref().map(|j| &j.state),
                                "joint_report": trial.joints.as_ref().map(|j| &j.report),'''
    assert s.count(old) == 1
    s = s.replace(old, new)
    if candidate:
        old = '''let x = a.clone().lu().solve(rhs).ok_or_else(|| {
                    Error::Solve("singular non-symmetric tangent; inspect active supports".into())
                })?;'''
        new = '''let (x, corrections, _, _) = crate::nc_linear_refinement::solve(
                    a, rhs, settings.relative_linear_tolerance,
                ).ok_or_else(|| Error::Solve("bounded compensated LU refinement failed both true residual gates".into()))?;'''
        assert s.count(old) == 1
        s = s.replace(old, new)
        s = s.replace('method: "dense_lu".into(),', 'method: "dense_lu_compensated_refinement_experiment".into(),')
        assert s.count('iterations: 1,') == 1
        s = s.replace('iterations: 1,', 'iterations: 1 + corrections,')
    s += '\n' + (HERE/'diagnostic_hook.rs').read_text()
    path.write_text(s)
    shutil.copyfile(HERE/'linear_refinement.rs', out/'src/nc_linear_refinement.rs')
    with (out/'src/lib.rs').open('a') as f:
        f.write('\n#[doc(hidden)]\npub mod nc_linear_refinement;\n')
    for name in ('nc_stop_probe.rs', 'linear_probe.rs'):
        shutil.copyfile(HERE/name, out/'examples'/name)
    (out/'tests').mkdir(exist_ok=True)
    shutil.copyfile(HERE/'linear_tests.rs', out/'tests/nc_linear_experiment.rs')
    (out/'original-source-sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    return hashes


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--candidate', action='store_true')
    a = p.parse_args()
    prepare(a.output.resolve(), a.candidate)
