//! Explicit A-01 observation only. Same solver; no mutable state in callbacks.
use soil_fem_3d::{
    Result,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Observer, Progress},
    stages,
};
use std::{
    cell::RefCell,
    fs::{self, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};
struct Trace(RefCell<BufWriter<fs::File>>);
impl Trace {
    fn record(&self, value: &impl serde::Serialize) -> Result<()> {
        let mut f = self.0.borrow_mut();
        serde_json::to_writer(&mut *f, value)?;
        writeln!(f)?;
        f.flush()?;
        Ok(())
    }
}
impl Observer for Trace {
    fn on_progress(&self, p: &Progress) -> Result<()> {
        self.record(&serde_json::json!({"kind":"progress","event":p}))
    }
    fn diagnostic_enabled(&self, _: usize, _: f64) -> bool {
        true
    }
    fn a01_joint_cause_enabled(&self) -> bool {
        true
    }
    fn on_diagnostic(&self, v: &serde_json::Value) -> Result<()> {
        self.record(v)
    }
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert!(
        a.len() == 4 || a.len() == 5 && a[4] == "--unobserved",
        "MESH GROUND CASE NEW_FOLDER [--unobserved]"
    );
    let mesh: Mesh = serde_json::from_slice(&fs::read(&a[0])?)?;
    assert!(
        mesh.tetrahedra.len() <= 48 && mesh.nodes.len() * 3 <= 256,
        "small A-01 observation only"
    );
    let ground: Ground = serde_json::from_slice(&fs::read(&a[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&a[2])?)?;
    let out = Path::new(&a[3]);
    fs::create_dir(out)?;
    let trace = Trace(RefCell::new(BufWriter::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("trace.jsonl"))?,
    )));
    let run = if a.len() == 5 {
        stages::solve(&mesh, ground, case)?
    } else {
        stages::solve_observed(&mesh, ground, case, &trace)?
    };
    let mut f = BufWriter::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join("run.json"))?,
    );
    serde_json::to_writer(&mut f, &run)?;
    f.flush()?;
    Ok(())
}
