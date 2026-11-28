//! Observation only: capture converged trial partitions without changing a solve.
use soil_fem_3d::{
    Result,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Kind, Observer, Progress},
    stages,
};
use std::{
    cell::RefCell,
    fs::{self, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};
struct Trace {
    stream: RefCell<BufWriter<std::fs::File>>,
}
impl Observer for Trace {
    fn diagnostic_enabled(&self, _: usize, _: f64) -> bool {
        true
    }
    fn on_progress(&self, p: &Progress) -> Result<()> {
        if p.kind == Kind::Accepted {
            let mut out = self.stream.borrow_mut();
            serde_json::to_writer(
                &mut *out,
                &serde_json::json!({"kind":"accepted", "stage_index":p.stage_index,"stage":p.stage,"progress":p.progress}),
            )?;
            writeln!(out)?;
        }
        Ok(())
    }
    fn on_diagnostic(&self, v: &serde_json::Value) -> Result<()> {
        if v["kind"] == "iterate"
            && v["relative_force_residual"]
                .as_f64()
                .is_some_and(|r| r <= 1e-8)
        {
            let counts = v["partitions"].as_array().expect("diagnostic partitions");
            assert_eq!(counts.len(), 356, "fixed 324 bulk + 32 contact GP specimen");
            let mut out = self.stream.borrow_mut();
            serde_json::to_writer(
                &mut *out,
                &serde_json::json!({"kind":"converged_trial", "stage_index":v["stage"],"progress":v["progress"],"old_u":v["committed"]["displacement_m"],"u":v["displacement_m"],"contact_partition_counts": &counts[324..]}),
            )?;
            writeln!(out)?;
        }
        Ok(())
    }
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 4, "MESH GROUND CASE NEW_FOLDER");
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&args[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&args[2])?)?;
    assert_eq!(mesh.nodes.len(), 54);
    assert_eq!(mesh.tetrahedra.len(), 12);
    let folder = Path::new(&args[3]);
    fs::create_dir(folder)?;
    let trace = Trace {
        stream: RefCell::new(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(folder.join("trials.jsonl"))?,
        )),
    };
    let run = stages::solve_observed(&mesh, ground, case, &trace)?;
    trace.stream.borrow_mut().flush()?;
    let mut out = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(folder.join("run.json"))?,
    );
    serde_json::to_writer(&mut out, &run)?;
    out.flush()?;
    assert!(run.complete, "{:?}", run.failure_reason);
    println!("complete; observation only");
    Ok(())
}
