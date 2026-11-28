//! Read-only diagnostic: checkpoint at every accepted increment, using public APIs.
//! No changes to the constitutive, contact, or nonlinear solution algorithms.
use serde::Serialize;
use soil_fem_3d::{
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stages::{self, Run, StageCase},
};
use std::{
    error::Error,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

type AnyResult<T> = std::result::Result<T, Box<dyn Error>>;
fn save(path: &Path, value: &impl Serialize) -> AnyResult<()> {
    let mut out = BufWriter::new(File::create_new(path)?);
    serde_json::to_writer(&mut out, value)?;
    out.flush()?;
    Ok(())
}
fn accepted(run: &Run) -> usize {
    run.attempts.iter().filter(|a| a.accepted).count()
        + run
            .geostatic
            .as_ref()
            .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count())
}
fn main() -> AnyResult<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        return Err(
            "usage: nc_history_trace trace|direct|resume MESH GROUND CASE|CHECKPOINT OUTPUT".into(),
        );
    }
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[2])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&args[3])?)?;
    let output = Path::new(&args[5]);
    fs::create_dir(output)?;
    if args[1] == "resume" {
        let cp: Checkpoint = serde_json::from_slice(&fs::read(&args[4])?)?;
        let case = cp.run.case.clone();
        let run = cp.resume(&mesh, ground, case, None)?;
        save(&output.join("run.json"), &run)?;
        return Ok(());
    }
    let case: StageCase = serde_json::from_slice(&fs::read(&args[4])?)?;
    if args[1] == "direct" {
        save(
            &output.join("run.json"),
            &stages::solve(&mesh, ground, case)?,
        )?;
        return Ok(());
    }
    if args[1] != "trace" {
        return Err("unknown diagnostic mode".into());
    }
    let mut run = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let mut stream = BufWriter::new(File::create_new(output.join("accepted.jsonl"))?);
    let mut seen = 0;
    let mut checkpoints = std::collections::BTreeSet::new();
    loop {
        let count = accepted(&run);
        if count > seen {
            if count != seen + 1 {
                return Err("accepted-increment trace has a gap".into());
            }
            serde_json::to_writer(
                &mut stream,
                &serde_json::json!({"accepted_count": count,
                "geostatic_complete": run.geostatic.as_ref().is_some_and(|g| g.complete),
                "state": run.last_committed}),
            )?;
            stream.write_all(b"\n")?;
            stream.flush()?;
            seen = count;
        }
        if run.failure_reason.as_deref() != Some("checkpoint_requested") {
            break;
        }
        if seen > 10000 {
            return Err("diagnostic trace budget exceeded; not a failure verdict".into());
        }
        let checkpoint = Checkpoint::seal(&mesh, run)?;
        let state = &checkpoint.run.last_committed;
        let phase = if checkpoint
            .run
            .geostatic
            .as_ref()
            .is_some_and(|g| !g.complete)
        {
            "geostatic".to_string()
        } else {
            format!("stage-{}", state.stage_index)
        };
        if state.stage_progress >= 0.5 && checkpoints.insert(phase.clone()) {
            save(
                &output.join(format!("checkpoint-{phase}.json")),
                &checkpoint,
            )?;
        }
        // Exercise exactly the serialized public checkpoint representation.
        let restored: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&checkpoint)?)?;
        run = restored.resume(&mesh, ground.clone(), case.clone(), Some(1))?;
    }
    save(&output.join("run.json"), &run)?;
    save(
        &output.join("trace-summary.json"),
        &serde_json::json!({
        "accepted_states_saved": seen, "complete": run.complete, "failure_reason": run.failure_reason,
        "numerical_acceptance": false, "all_steps_serialized_public_checkpoints": true,
        "checkpoints": checkpoints}),
    )?;
    Ok(())
}
