use serde_json::{Value, json};
use soil_fem_3d::{
    Result,
    stage_job::{self, State, Status},
};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};
static COUNTER: AtomicU64 = AtomicU64::new(0);
fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
fn directory(tag: &str) -> PathBuf {
    let parent = root().join("runs");
    fs::create_dir_all(&parent).unwrap();
    parent.join(format!(
        "stage-job-test-{tag}-{}-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ))
}
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_soil-fem-3d"))
}
fn execute(exe: &Path, args: &[&Path], code: i32) -> Output {
    let output = Command::new(exe).args(args).output().unwrap();
    assert_eq!(
        output.status.code(),
        Some(code),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
fn cli(args: &[&Path], code: i32) -> Output {
    execute(&binary(), args, code)
}
fn cmd(s: &str) -> &Path {
    Path::new(s)
}
fn read(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn prepare(folder: &Path) {
    cli(
        &[
            cmd("prepare-stages-job"),
            &root().join("examples/hardening-column-tet10-mesh.json"),
            &root().join("examples/ground-mixed-mcc.json"),
            &root().join("examples/stages-mixed-mcc.json"),
            folder,
        ],
        0,
    );
}
fn status(folder: &Path) -> Status {
    serde_json::from_slice(&fs::read(folder.join("status.json")).unwrap()).unwrap()
}
fn no_time(mut v: Value) -> Value {
    fn visit(v: &mut Value) {
        match v {
            Value::Object(m) => {
                m.remove("elapsed_seconds");
                for v in m.values_mut() {
                    visit(v);
                }
            }
            Value::Array(a) => {
                for v in a {
                    visit(v);
                }
            }
            _ => (),
        }
    }
    visit(&mut v);
    v
}

fn prepare_durable_case(folder: &Path, interval: usize, case: &Path) {
    cli(
        &[
            cmd("prepare-durable-stages-job"),
            &root().join("examples/hardening-column-tet10-mesh.json"),
            &root().join("examples/ground-mixed-mcc.json"),
            case,
            Path::new(&interval.to_string()),
            folder,
        ],
        0,
    );
}

#[test]
fn durable_slots_and_snapshot_fork_preserve_all_non_time_result_fields() {
    let job = directory("durable-complete");
    prepare_durable_case(&job, 2, &root().join("examples/stages-mixed-mcc.json"));
    cli(&[cmd("run-stages-job"), &job], 0);
    assert_eq!(status(&job).state, State::Succeeded);
    let pointer = read(&job.join("recovery.json"));
    assert!(pointer["generation"].as_u64().unwrap() >= 2);
    assert!(job.join("recovery-checkpoint-0.json").is_file());
    assert!(job.join("recovery-checkpoint-1.json").is_file());
    assert!(
        !fs::read_dir(&job).unwrap().any(|q| q
            .unwrap()
            .path()
            .extension()
            .is_some_and(|e| e == "part"))
    );
    let control = directory("durable-control.json");
    cli(
        &[
            cmd("solve-stages"),
            &job.join("inputs/mesh.json"),
            &job.join("inputs/ground.json"),
            &job.join("inputs/case.json"),
            &control,
        ],
        0,
    );
    assert!(no_time(read(&control)) == no_time(read(&job.join("result.json"))));
    let fork = directory("durable-fork");
    cli(&[cmd("prepare-stages-recovery"), &job, &fork], 0);
    assert_eq!(
        read(&fork.join("job.json"))["binding"]["checkpoint_interval"],
        2
    );
    cli(&[cmd("run-stages-job"), &fork], 0);
    assert!(no_time(read(&control)) == no_time(read(&fork.join("result.json"))));
    assert_eq!(status(&job).state, State::Succeeded);
}

#[test]
fn durable_interval_and_recovery_tampering_are_rejected_without_claiming_an_output() {
    for interval in ["0", "4097"] {
        let destination = directory("invalid-durable-interval");
        cli(
            &[
                cmd("prepare-durable-stages-job"),
                &root().join("examples/hardening-column-tet10-mesh.json"),
                &root().join("examples/ground-mixed-mcc.json"),
                &root().join("examples/stages-mixed-mcc.json"),
                cmd(interval),
                &destination,
            ],
            1,
        );
        assert!(!destination.exists());
    }
    let job = directory("durable-tamper");
    prepare_durable_case(&job, 2, &root().join("examples/stages-mixed-mcc.json"));
    cli(&[cmd("run-stages-job"), &job], 0);
    let pointer_path = job.join("recovery.json");
    let bytes = fs::read(&pointer_path).unwrap();
    let mut pointer = read(&pointer_path);
    pointer["checkpoint"] = json!("../outside.json");
    fs::write(&pointer_path, serde_json::to_vec(&pointer).unwrap()).unwrap();
    let destination = directory("durable-path-rejected");
    cli(&[cmd("prepare-stages-recovery"), &job, &destination], 1);
    assert!(!destination.exists());
    fs::write(&pointer_path, &bytes).unwrap();
    let pointer = read(&pointer_path);
    let slot = job.join(pointer["checkpoint"].as_str().unwrap());
    let original = fs::read(&slot).unwrap();
    let mut changed = original.clone();
    changed.push(b' ');
    fs::write(&slot, changed).unwrap();
    let destination = directory("durable-checksum-rejected");
    cli(&[cmd("prepare-stages-recovery"), &job, &destination], 1);
    assert!(!destination.exists());
    fs::write(&slot, original).unwrap();
    let mut pointer = read(&pointer_path);
    pointer["job_id"] = json!("another-job");
    fs::write(&pointer_path, serde_json::to_vec(&pointer).unwrap()).unwrap();
    let destination = directory("durable-id-rejected");
    cli(&[cmd("prepare-stages-recovery"), &job, &destination], 1);
    assert!(!destination.exists());
    fs::write(pointer_path, bytes).unwrap();
}

#[test]
fn durable_checkpoint_survives_abrupt_owned_child_exit_and_recovers_exact_history() {
    use std::time::{Duration, Instant};
    let mut case = read(&root().join("examples/stages-mixed-mcc.json"));
    case["settings"]["initial_increment"] = json!(1. / 64.);
    let case_path = directory("durable-slow-case.json");
    fs::write(&case_path, serde_json::to_vec(&case).unwrap()).unwrap();
    let job = directory("durable-abrupt");
    prepare_durable_case(&job, 2, &case_path);
    let mut child = Command::new(job.join("solver/soil-fem-3d.exe"))
        .args([cmd("run-stages-job"), &job])
        .spawn()
        .unwrap();
    let started = Instant::now();
    loop {
        assert!(
            child.try_wait().unwrap().is_none(),
            "child finished before abrupt-exit observation"
        );
        if job.join("recovery.json").exists() {
            break;
        }
        assert!(
            started.elapsed() < Duration::from_secs(60),
            "no durable snapshot produced"
        );
        std::thread::sleep(Duration::from_millis(2));
    }
    child.kill().unwrap();
    child.wait().unwrap();
    assert_eq!(status(&job).state, State::Running);
    assert!(!job.join("checkpoint.json").exists());
    let recovered = directory("durable-after-abrupt");
    cli(&[cmd("prepare-stages-recovery"), &job, &recovered], 0);
    cli(&[cmd("run-stages-job"), &recovered], 0);
    let control = directory("durable-abrupt-control.json");
    cli(
        &[
            cmd("solve-stages"),
            &job.join("inputs/mesh.json"),
            &job.join("inputs/ground.json"),
            &job.join("inputs/case.json"),
            &control,
        ],
        0,
    );
    assert!(no_time(read(&control)) == no_time(read(&recovered.join("result.json"))));
}

#[test]
fn cli_pause_and_separate_process_resume_preserve_every_non_time_result_field() {
    let job = directory("pause");
    prepare(&job);
    assert_eq!(status(&job).state, State::Prepared);
    cli(&[cmd("stop-stages-job"), &job], 0);
    execute(
        &job.join("solver/soil-fem-3d.exe"),
        &[cmd("run-stages-job"), &job],
        2,
    );
    let saved = status(&job);
    assert_eq!(saved.state, State::Checkpointed);
    assert!(job.join("checkpoint.json").is_file() && !job.join("result.json").exists());
    let checkpoint = read(&job.join("checkpoint.json"));
    assert_eq!(
        checkpoint["run"]["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|a| a["accepted"] == true)
            .count(),
        1
    );
    let lines = fs::read_to_string(job.join("progress.jsonl")).unwrap();
    let events: Vec<Value> = lines
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    assert_eq!(events.len(), saved.progress_events);
    for (i, event) in events.iter().enumerate() {
        assert_eq!(event["sequence"], i + 1);
        assert_eq!(event["job_id"], saved.job_id);
    }
    assert!(events.iter().any(|e| e["progress"]["kind"] == "newton"
        && e["progress"]["relative_volume_residual"].is_number()));
    assert_eq!(events.last().unwrap()["progress"]["kind"], "accepted");
    let checkpoint_bytes = fs::read(job.join("checkpoint.json")).unwrap();
    cli(&[cmd("run-stages-job"), &job], 1);
    cli(&[cmd("stop-stages-job"), &job], 0);
    assert_eq!(status(&job).state, State::Checkpointed);
    assert!(checkpoint_bytes == fs::read(job.join("checkpoint.json")).unwrap());
    let resumed = directory("resume");
    execute(
        &job.join("solver/soil-fem-3d.exe"),
        &[cmd("prepare-stages-resume"), &job, &resumed],
        0,
    );
    assert_ne!(status(&resumed).job_id, saved.job_id);
    assert_eq!(
        read(&resumed.join("job.json"))["binding"]["parent_job_id"],
        saved.job_id
    );
    assert!(!resumed.join("stop-request.json").exists());
    execute(
        &resumed.join("solver/soil-fem-3d.exe"),
        &[cmd("run-stages-job"), &resumed],
        0,
    );
    assert_eq!(status(&resumed).state, State::Succeeded);
    let control = directory("control.json");
    cli(
        &[
            cmd("solve-stages"),
            &job.join("inputs/mesh.json"),
            &job.join("inputs/ground.json"),
            &job.join("inputs/case.json"),
            &control,
        ],
        0,
    );
    assert!(
        no_time(read(&control)) == no_time(read(&resumed.join("result.json"))),
        "CLI cooperative restart changed physical history/path or diagnostics"
    );
}

#[test]
fn a_stop_request_for_another_job_is_refused_and_cannot_pause_the_running_job() {
    let job = directory("wrong-request");
    prepare(&job);
    let wrong = json!({"schema":"soil-fem-stages-job-stop/1","job_id":"another-job"});
    fs::write(
        job.join("stop-request.json"),
        serde_json::to_vec(&wrong).unwrap(),
    )
    .unwrap();
    cli(&[cmd("stop-stages-job"), &job], 1);
    assert_eq!(read(&job.join("stop-request.json")), wrong);
    cli(&[cmd("run-stages-job"), &job], 0);
    assert_eq!(status(&job).state, State::Succeeded);
    assert!(!status(&job).stop_requested);
    assert!(!job.join("checkpoint.json").exists());
}

#[test]
fn immutable_input_changes_existing_output_and_second_execution_are_refused_before_solving() {
    let job = directory("changed-input");
    prepare(&job);
    let manifest = fs::read(job.join("job.json")).unwrap();
    let old = fs::read(job.join("inputs/case.json")).unwrap();
    let mut changed = old.clone();
    changed.push(b' ');
    fs::write(job.join("inputs/case.json"), &changed).unwrap();
    cli(&[cmd("run-stages-job"), &job], 1);
    assert_eq!(status(&job).state, State::Prepared);
    assert!(!job.join("run.claim").exists() && !job.join("progress.jsonl").exists());
    cli(
        &[
            cmd("prepare-stages-job"),
            &root().join("examples/hardening-column-tet10-mesh.json"),
            &root().join("examples/ground-mixed-mcc.json"),
            &root().join("examples/stages-mixed-mcc.json"),
            &job,
        ],
        1,
    );
    assert!(manifest == fs::read(job.join("job.json")).unwrap());
    assert!(changed == fs::read(job.join("inputs/case.json")).unwrap());
}

#[test]
fn a_different_running_solver_cannot_claim_a_valid_foreign_solver_snapshot() -> Result<()> {
    let job = directory("foreign-binary");
    // The library binds this test executable. The actual CLI is a different
    // valid binary, so the refusal exercises identity rather than corrupt JSON.
    stage_job::prepare(
        &root().join("examples/hardening-column-tet10-mesh.json"),
        &root().join("examples/ground-mixed-mcc.json"),
        &root().join("examples/stages-mixed-mcc.json"),
        &job,
    )?;
    let output = cli(&[cmd("run-stages-job"), &job], 1);
    assert!(String::from_utf8_lossy(&output.stderr).contains("frozen solver executable"));
    assert_eq!(status(&job).state, State::Prepared);
    assert!(!job.join("run.claim").exists());
    Ok(())
}

#[test]
fn invalid_material_and_numerical_failure_never_publish_success_or_a_false_checkpoint() {
    for (name, numerical) in [("invalid-material", false), ("numerical-failure", true)] {
        let source = directory(name);
        fs::create_dir(&source).unwrap();
        let mut ground = read(&root().join("examples/ground-mixed-mcc.json"));
        let mut case = read(&root().join("examples/stages-mixed-mcc.json"));
        if numerical {
            case["settings"]["maximum_newton_iterations"] = json!(1);
            case["settings"]["initial_increment"] = json!(1.);
            case["settings"]["minimum_increment"] = json!(1.);
        } else {
            ground["materials"][0]["mcc"]["lambda"] =
                ground["materials"][0]["mcc"]["kappa"].clone();
        }
        fs::write(
            source.join("ground.json"),
            serde_json::to_vec(&ground).unwrap(),
        )
        .unwrap();
        fs::write(source.join("case.json"), serde_json::to_vec(&case).unwrap()).unwrap();
        let job = source.join("job");
        cli(
            &[
                cmd("prepare-stages-job"),
                &root().join("examples/hardening-column-tet10-mesh.json"),
                &source.join("ground.json"),
                &source.join("case.json"),
                &job,
            ],
            0,
        );
        cli(&[cmd("run-stages-job"), &job], 1);
        let state = status(&job);
        assert_eq!(state.state, State::Failed);
        assert!(!job.join("result.json").exists() && !job.join("checkpoint.json").exists());
        if numerical {
            assert_eq!(state.artifact.as_deref(), Some("incomplete.json"));
            assert_eq!(read(&job.join("incomplete.json"))["complete"], false);
        } else {
            assert!(state.reason.is_some());
        }
    }
}
