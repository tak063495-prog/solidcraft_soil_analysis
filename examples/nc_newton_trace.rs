//! Opt-in trace of the existing solver, not an independent solution.
use soil_fem_3d::{
    Error, Result,
    stage_control::{Observer, Progress},
    stages,
};
use std::{
    cell::RefCell,
    fs::{self, File, OpenOptions},
    io::Write,
    path::Path,
    time::Instant,
};
struct Trace {
    file: RefCell<File>,
    start: Instant,
}
impl Observer for Trace {
    fn on_progress(&self, p: &Progress) -> Result<()> {
        self.on_diagnostic(&serde_json::json!({"kind": "progress", "event": p}))
    }
    fn diagnostic_enabled(&self, stage: usize, progress: f64) -> bool {
        stage == 3 && (0.9512958..=0.951297).contains(&progress)
    }
    fn on_diagnostic(&self, value: &serde_json::Value) -> Result<()> {
        let mut f = self.file.borrow_mut();
        serde_json::to_writer(&mut *f, value)?;
        writeln!(f)?;
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(self.start.elapsed().as_secs() >= 180)
    }
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err(Error::Input(
            "nc_newton_trace INPUT_DIR NEW_OUTPUT_DIR".into(),
        ));
    }
    let input = Path::new(&args[1]);
    let out = Path::new(&args[2]);
    fs::create_dir(out)?;
    let observer = Trace {
        file: RefCell::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(out.join("trace.jsonl"))?,
        ),
        start: Instant::now(),
    };
    let run = stages::solve_observed(
        &serde_json::from_slice(&fs::read(input.join("mesh.json"))?)?,
        serde_json::from_slice(&fs::read(input.join("ground.json"))?)?,
        serde_json::from_slice(&fs::read(input.join("stages.json"))?)?,
        &observer,
    )?;
    fs::write(out.join("run.json"), serde_json::to_vec(&run)?)?;
    println!(
        "complete={} progress={} (diagnostic only)",
        run.complete, run.last_committed.stage_progress
    );
    Ok(())
}
