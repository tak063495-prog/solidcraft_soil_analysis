"""Instrument only an isolated source copy of the existing dry probe."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLLOW = ROOT / 'interim/2026-12-31-a04-richards-wet-joint-cloud-followup'
sys.path.insert(0, str(FOLLOW))
from a04follow_build import replace
from a04follow_common import sha, write

# One evaluation on the solver thread. Collection starts only at open/.25.
CAPTURE = """
thread_local! {
    static ROWS: std::cell::RefCell<Option<Vec<serde_json::Value>>> = const { std::cell::RefCell::new(None) };
}
pub(crate) fn sample(row: impl FnOnce() -> serde_json::Value) {
    ROWS.with_borrow_mut(|rows| { if let Some(rows) = rows { rows.push(row()); } });
}
"""


def build(out):
    out.mkdir(parents=True, exist_ok=False)
    src = out / 'src'
    shutil.copytree(ROOT / 'src', src)
    for name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
        shutil.copyfile(ROOT / name, out / name)
    (src / 'bin').mkdir(exist_ok=True)
    shutil.copyfile(ROOT / 'examples/joint_a02_history_probe.rs', src / 'bin/a04ref_dry.rs')
    stages = src / 'stages.rs'
    code = stages.read_text()
    code = replace(code, '                        for i in 0..e.dofs.len() {\n                            internal[e.dofs[i]] += force[i];', '''
                        sample(|| serde_json::json!(["bulk", ei, qi, q.position, e.dofs, force.as_slice()]));
                        for i in 0..e.dofs.len() {
                            internal[e.dofs[i]] += force[i];''')
    call = '            let eval = self.evaluate(old, &u, &active, &dofs, settings, mode)?;'
    code = replace(code, call, '''
            let path = std::env::var_os("A04REF_DRY_TRACE")
                .filter(|_| stage == 2 && progress == 0.25);
            if path.is_some() { ROWS.with_borrow_mut(|r| *r = Some(Vec::new())); }
''' + call)
    mark = '            let relative = reduced.norm() / scale;'
    code = replace(code, mark, mark + '''
            if let Some(path) = path {
                use std::io::Write;
                let mut stream = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
                let joints = eval.joints.as_ref().unwrap();
                let rows = ROWS.with_borrow_mut(|r| r.take().unwrap());
                let record = serde_json::json!({"iteration":iteration,"u":u,"rows":rows,
                    "internal":eval.internal,"load":target,"relative":relative,"reduced":reduced.as_slice(),
                    "joint":joints.state.points,"old_joint":old.joints.as_ref().unwrap().points,
                    "partitions":&partitions[self.point_count..],"linear":linear_solves,"scale":scale});
                serde_json::to_writer(&mut stream, &record)?;
                writeln!(stream)?;
            }
''')
    stages.write_text(code + CAPTURE, newline='\n')
    joints = src / 'joints.rs'
    code = joints.read_text()
    mark = '                let q = Vector3::from(next.points[i].traction_kpa);'
    joints.write_text(replace(code, mark, mark + '''
                crate::stages::sample(|| serde_json::json!(["contact", i, p.info, p.coefficients,
                    p.coefficients.iter().map(|&(_, c)| [0,1,2].map(|a| p.info.weight_m2*c*q[a])).collect::<Vec<_>>()]));
'''), newline='\n')
    elements = src / 'elements.rs'
    code = replace(elements.read_text(), '    Ok(result)\n}\n\n#[cfg(test)]',
        '    if std::env::var_os("A04FOLLOW_Q").is_some() { result.reverse(); }\n    Ok(result)\n}\n\n#[cfg(test)]')
    elements.write_text(code, newline='\n')
    subprocess.run(['cargo', 'build', '--locked', '--release', '--bin', 'a04ref_dry'], cwd=out, check=True)
    write(out / 'binding.json', {p.relative_to(out).as_posix(): sha(p) for p in sorted(src.rglob('*.rs'))})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    build(p.parse_args().output.resolve())
