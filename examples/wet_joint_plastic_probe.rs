//! C-03 observation only: accepted constitutive partitions and real-file checkpoints.
use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Kind, Observer, Progress},
    stages,
};
use std::{
    cell::{Cell, RefCell},
    fs::{self, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};
fn write(path: &Path, value: &impl serde::Serialize) -> Result<()> {
    let mut f = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?);
    serde_json::to_writer(&mut f, value)?;
    f.flush()?;
    Ok(())
}
struct Trace {
    stream: RefCell<BufWriter<std::fs::File>>,
    accepted: Cell<usize>,
    diagnostics: RefCell<BufWriter<std::fs::File>>,
    last: RefCell<Option<serde_json::Value>>,
    diagnostic: bool,
}
impl Observer for Trace {
    fn on_progress(&self, p: &Progress) -> Result<()> {
        let mut f = self.stream.borrow_mut();
        serde_json::to_writer(&mut *f, &serde_json::json!({"kind":"progress","event":p}))?;
        writeln!(f)?;
        if p.kind == Kind::Accepted {
            self.accepted.set(self.accepted.get() + 1);
            if let Some(value) = self.last.borrow_mut().take() {
                let mut stream = self.diagnostics.borrow_mut();
                serde_json::to_writer(&mut *stream, &value)?;
                writeln!(stream)?;
                stream.flush()?;
            }
        }
        Ok(())
    }
    fn diagnostic_enabled(&self, _stage: usize, _progress: f64) -> bool {
        self.diagnostic
    }
    fn on_diagnostic(&self, value: &serde_json::Value) -> Result<()> {
        if value["kind"] == "iterate" {
            *self.last.borrow_mut() = Some(serde_json::json!({
                "stage": value["stage"], "progress": value["progress"],
                "partitions": value["partitions"], "strain_increment": value["strain_increment"],
                "relative_force_residual": value["relative_force_residual"],
            }));
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(self.accepted.get() >= 1)
    }
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert!(
        args.len() == 4 || args.len() == 6 || (args.len() == 5 && args[4] == "--whole"),
        "MESH GROUND CASE NEW_FOLDER [--resume CHECKPOINT]"
    );
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&args[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&args[2])?)?;
    assert!(mesh.tetrahedra.len() <= 768 || std::env::var_os("C03_USER_PC").is_some());
    let folder = Path::new(&args[3]);
    fs::create_dir(folder)?;
    if args.len() == 5 {
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete);
        write(&folder.join("run.json"), &run)?;
        return Ok(());
    }
    let trace = Trace {
        stream: RefCell::new(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(folder.join("progress.jsonl"))?,
        )),
        accepted: Cell::new(0),
        diagnostics: RefCell::new(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(folder.join("constitutive.jsonl"))?,
        )),
        last: RefCell::new(None),
        diagnostic: mesh.tetrahedra.len() <= 12 && mesh.nodes.len() * 3 <= 256,
    };
    let resumed = args.len() == 6;
    let mut run = if resumed {
        assert_eq!(args[4], "--resume");
        let cp: Checkpoint = serde_json::from_slice(&fs::read(&args[5])?)?;
        cp.resume_observed(&mesh, ground.clone(), case.clone(), &trace)?
    } else {
        stages::solve_observed(&mesh, ground.clone(), case.clone(), &trace)?
    };
    let mut stream = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(folder.join("accepted.jsonl"))?,
    );
    let mut count = 0;
    let steps = (1. / case.settings.initial_increment).round() as usize;
    loop {
        let actual = run.attempts.iter().filter(|a| a.accepted).count()
            + run
                .geostatic
                .as_ref()
                .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count());
        if trace.accepted.get() == 0 {
            assert!(run.complete);
            break;
        }
        assert_eq!(trace.accepted.get(), 1);
        count += 1;
        serde_json::to_writer(&mut stream, &run.last_committed)?;
        writeln!(stream)?;
        stream.flush()?;
        if run.complete {
            break;
        }
        assert_eq!(run.failure_reason.as_deref(), Some("checkpoint_requested"));
        let cp = Checkpoint::seal(&mesh, run)?;
        if !resumed && actual % steps == steps / 2 {
            write(&folder.join(format!("checkpoint-{actual:04}.json")), &cp)?;
        }
        trace.accepted.set(0);
        run = cp.resume_observed(&mesh, ground.clone(), case.clone(), &trace)?;
    }
    trace.stream.borrow_mut().flush()?;
    write(&folder.join("run.json"), &run)?;
    write(
        &folder.join("verification.json"),
        &serde_json::json!({"accepted_states_written":count,"complete":run.complete,"resumed_from_file":resumed,"process_id":std::process::id()}),
    )?;
    assert!(run.complete);
    println!("complete; wrote {count} accepted states");
    Ok(())
}
