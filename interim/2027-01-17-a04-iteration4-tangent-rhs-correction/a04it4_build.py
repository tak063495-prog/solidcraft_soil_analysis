"""Build the 01-12 open/.25 probe in an isolated source tree."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / '2026-12-31-a04-richards-wet-joint-cloud-followup'))
from a04follow_build import replace
from a04follow_common import sha, write


def build(out):
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(ROOT / 'src', out / 'src')
    for name in ('Cargo.lock', 'Cargo.toml', 'rust-toolchain.toml'):
        shutil.copyfile(ROOT / name, out / name)
    binary = out / 'src/bin/a04it4_dry.rs'
    binary.parent.mkdir(exist_ok=True)
    shutil.copyfile(ROOT / 'examples/joint_a02_history_probe.rs', binary)
    elements = out / 'src/elements.rs'
    mark = '    Ok(result)\n}\n\n#[cfg(test)]'
    elements.write_text(replace(elements.read_text(), mark,
        '    if std::env::var_os("A04FOLLOW_Q").is_some() { result.reverse(); }\n' + mark), newline='\n')
    stages = out / 'src/stages.rs'
    source = stages.read_text()
    mark = '            let relative = reduced.norm() / scale;'
    source = replace(source, mark, mark + '''
            let trace = std::env::var_os("A04IT4_TRACE")
                .filter(|_| stage == 2 && progress == 0.25);
            let record = |value: serde_json::Value| -> Result<()> {
                use std::io::Write;
                let mut f = std::fs::OpenOptions::new().create(true).append(true)
                    .open(trace.as_ref().unwrap())?;
                serde_json::to_writer(&mut f, &value)?;
                writeln!(f)?;
                Ok(())
            };
            if trace.is_some() {
                record(serde_json::json!({"kind":"iterate","i":iteration,"u":u,"relative":relative,
                    "internal":eval.internal,"target":target,"old":old.joints.as_ref().unwrap().points,
                    "points":eval.joints.as_ref().unwrap().state.points,"counts":&partitions[self.point_count..],
                    "rhs":reduced.as_slice(),"scale":scale}))?;
            }''')
    mark = '            let mut best_candidate_relative = f64::INFINITY;'
    source = replace(source, mark, mark + '''
            if trace.is_some() && iteration == 4 {
                assert!(predictor_index.is_none());
                let Tangent::Dense(a) = &matrix else { panic!("expected dry dense LU") };
                record(serde_json::json!({"kind":"linear","i":iteration,"rows":dofs.rows,
                    "matrix":(0..a.nrows()).map(|i| a.row(i).iter().copied().collect::<Vec<_>>()).collect::<Vec<_>>(),
                    "rhs":reduced.as_slice(),"dx":correction.as_slice(),"report":linear_solves.last()}))?;
            }''')
    mark = '                        u = candidate;\n                        accepted = true;'
    source = replace(source, mark, '''
                        if trace.is_some() && iteration == 4 {
                            record(serde_json::json!({"kind":"accepted","i":iteration,"u":candidate,
                                "alpha":alpha,"relative":self.reduce(&r,&dofs).norm()/scale}))?;
                        }
''' + mark)
    stages.write_text(source, newline='\n')
    subprocess.run(['cargo','build','--locked','--release','--bin','a04it4_dry'], cwd=out, check=True)
    write(out / 'binding.json', {p.relative_to(out).as_posix():sha(p) for p in sorted((out/'src').rglob('*.rs'))})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    build(parser.parse_args().output.resolve())
