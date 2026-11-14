//! Verification driver: observe committed states through the existing checkpoint API.
//! No new solver/contact path; each accepted boundary is sealed before resuming.
use soil_fem_3d::{Result, checkpoint::Checkpoint, ground::Ground, mesh::Mesh, stages};
use std::{
    fs::{self, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};
fn write(path: &Path, value: &impl serde::Serialize) -> Result<()> {
    let mut file = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?);
    serde_json::to_writer(&mut file, value)?;
    file.flush()?;
    Ok(())
}
fn strip_clock(v: &mut serde_json::Value) {
    match v {
        serde_json::Value::Object(m) => {
            m.remove("elapsed_seconds");
            for x in m.values_mut() {
                strip_clock(x);
            }
        }
        serde_json::Value::Array(a) => {
            for x in a {
                strip_clock(x);
            }
        }
        _ => {}
    }
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let user_pc = args.len() == 5 && args[4] == "--user-pc";
    assert!(
        args.len() == 4 || user_pc,
        "MESH GROUND STAGES NEW_FOLDER [--user-pc]"
    );
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&args[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&args[2])?)?;
    assert!(
        user_pc || mesh.tetrahedra.len() < 256,
        "large cases belong on the user PC"
    );
    let folder = Path::new(&args[3]);
    fs::create_dir(folder)?;
    let mut stream = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(folder.join("accepted.jsonl"))?,
    );
    let mut run = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let mut count = 0;
    loop {
        let actual = run.attempts.iter().filter(|a| a.accepted).count()
            + run
                .geostatic
                .as_ref()
                .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count());
        if actual == count {
            assert!(run.complete, "resume made no accepted progress");
            break;
        }
        assert_eq!(actual, count + 1, "missing accepted boundary");
        count = actual;
        serde_json::to_writer(&mut stream, &run.last_committed)?;
        writeln!(stream)?;
        stream.flush()?;
        if run.complete {
            break;
        }
        if run.failure_reason.as_deref() != Some("checkpoint_requested") {
            write(&folder.join("failed-run.json"), &run)?;
            return Err(soil_fem_3d::Error::Solve(format!(
                "incomplete: {:?}",
                run.failure_reason
            )));
        }
        let checkpoint = Checkpoint::seal(&mesh, run)?;
        run = checkpoint.resume(&mesh, ground.clone(), case.clone(), Some(1))?;
    }
    write(&folder.join("run.json"), &run)?;
    write(
        &folder.join("layout.json"),
        &serde_json::json!({"integration_points": run.integration_points, "joint_integration_points":run.joint_integration_points}),
    )?;
    // A separate uninterrupted solve authenticates the instrumentation at every endpoint.
    let whole = stages::solve(&mesh, ground, case)?;
    let mut a = serde_json::to_value(&run)?;
    let mut b = serde_json::to_value(&whole)?;
    strip_clock(&mut a);
    strip_clock(&mut b);
    if a != b {
        return Err(soil_fem_3d::Error::Solve(
            "checkpoint observer changed uninterrupted result".into(),
        ));
    }
    write(
        &folder.join("verification.json"),
        &serde_json::json!({"accepted_states":count,"uninterrupted_run_exact_except_elapsed":true,"complete":run.complete,"separate_process_restart_verified":false}),
    )?;
    println!("saved {count} accepted states; uninterrupted run matches");
    Ok(())
}
