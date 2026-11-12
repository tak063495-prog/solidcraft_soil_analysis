// Standalone public-API evidence executable. Does not change production APIs.
use serde::de::DeserializeOwned;
use soil_fem_3d::{
    Result,
    checkpoint::{Checkpoint, hash_json},
    ground::Ground,
    mesh::Mesh,
    stages::{self, Run, StageCase},
};
use std::{
    fs::{self, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};
fn read<T: DeserializeOwned>(p: &str) -> Result<T> {
    Ok(serde_json::from_slice(&fs::read(p)?)?)
}
fn write(p: &Path, v: &impl serde::Serialize) -> Result<()> {
    let f = OpenOptions::new().write(true).create_new(true).open(p)?;
    Ok(serde_json::to_writer(BufWriter::new(f), v)?)
}
fn count(r: &Run) -> usize {
    r.attempts.iter().filter(|a| a.accepted).count()
        + r.geostatic
            .as_ref()
            .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count())
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    if a[1] == "rehash" {
        let mut cp: Checkpoint = read(&a[2])?;
        cp.run_sha256 = hash_json(&cp.run)?;
        return write(Path::new(&a[3]), &cp);
    }
    let mesh: Mesh = read(&a[2])?;
    let ground: Ground = read(&a[3])?;
    let case: StageCase = read(&a[4])?;
    let out = Path::new(&a[5]);
    fs::create_dir(out)?;
    if a[1] == "plain" {
        return write(&out.join("run.json"), &stages::solve(&mesh, ground, case)?);
    }
    let mut prior = 0;
    let mut run = if a[1] == "resume" {
        let cp: Checkpoint = read(&a[6])?;
        prior = count(&cp.run);
        cp.resume(&mesh, ground.clone(), case.clone(), Some(1))?
    } else {
        stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?
    };
    let mut history = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join("accepted.jsonl"))?,
    );
    loop {
        let n = count(&run);
        if n > prior {
            assert_eq!(n, prior + 1);
            let geo = run.attempts.iter().all(|a| !a.accepted);
            let state = if geo {
                &run.geostatic.as_ref().unwrap().last_committed
            } else {
                &run.last_committed
            };
            let progress = if geo {
                run.geostatic.as_ref().unwrap().progress
            } else {
                state.stage_progress
            };
            serde_json::to_writer(
                &mut history,
                &serde_json::json!({"accepted_index":n,
                "stage_index":if geo {0} else {state.stage_index}, "progress":progress, "state":state}),
            )?;
            history.write_all(b"\n")?;
            history.flush()?;
            prior = n;
        }
        if run.complete || run.failure_reason.as_deref() != Some("checkpoint_requested") {
            write(&out.join("run.json"), &run)?;
            break;
        }
        let cp = Checkpoint::seal(&mesh, run)?;
        if a[1] == "solve" {
            write(&out.join(format!("checkpoint-{n}.json")), &cp)?;
        }
        run = cp.resume(&mesh, ground.clone(), case.clone(), Some(1))?;
    }
    Ok(())
}
