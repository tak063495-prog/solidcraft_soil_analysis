use serde_json::{Value, json};
use soil_fem_3d::{
    Result, checkpoint, epoch_srm_gate::Specification, ground::Ground, mesh::Mesh,
    srm_bound::Problem, stages::StageCase,
};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

static COUNTER: AtomicU64 = AtomicU64::new(0);

fn cmd(name: &str) -> &Path {
    Path::new(name)
}
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_soil-fem-3d"))
}
fn execute(exe: &Path, args: &[&Path], expected: i32) -> Output {
    let out = Command::new(exe).args(args).output().unwrap();
    assert_eq!(
        out.status.code(),
        Some(expected),
        "{args:?}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}
fn cli(args: &[&Path], expected: i32) -> Output {
    execute(&binary(), args, expected)
}
fn read(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn write(path: &Path, value: &impl serde::Serialize) {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
}
fn directory(tag: &str) -> PathBuf {
    let parent = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("runs");
    fs::create_dir_all(&parent).unwrap();
    let path = parent.join(format!(
        "epoch SRM gate {tag} {}-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&path).unwrap();
    path
}

fn inputs(
    required: f64,
    attempts: usize,
    gates: Vec<usize>,
) -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let problem: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube": seed});
    raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    raw["materials"] = json!([{
        "id": 1, "name": "unconfined epoch gate", "unit_weight_kn_m3": 1e-10, "mohr_coulomb": problem.material
    }]);
    raw["bodies"] = json!([{"body": "TestCube", "material_id": 1}]);
    let ground = serde_json::from_value(raw)?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(problem.boundary);
    case.stages = serde_json::from_value(json!([
        {"kind": "boundary", "name": "load", "target_factor": 1.0, "initial_increment": 1.0},
        {"kind": "reset_display", "name": "reference"}
    ]))?;
    case.verification_reject_once_stages.clear();
    case.monitors.clear();
    case.monitor_sets.clear();
    let mut search: soil_fem_3d::srm_search::Specification =
        serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    search.through_stages = 0;
    search.required_factor = required;
    let spec = Specification {
        schema: "soil-fem-epoch-srm-stage-gate/1".into(),
        search: soil_fem_3d::epoch_srm_search::Specification {
            schema: "soil-fem-epoch-srm-search/1".into(),
            search,
            initial_increment: 0.25,
            minimum_increment: 1. / 4096.,
            maximum_attempts: attempts,
        },
        gate_after_stages: gates,
    };
    Ok((mesh, ground, case, spec))
}

fn fixture(tag: &str, required: f64, attempts: usize) -> PathBuf {
    let folder = directory(tag);
    let (mesh, ground, case, spec) = inputs(required, attempts, vec![1, 2]).unwrap();
    write(&folder.join("mesh.json"), &mesh);
    write(&folder.join("ground.json"), &ground);
    write(&folder.join("case.json"), &case);
    write(&folder.join("spec.json"), &spec);
    folder
}
fn prepare(folder: &Path, name: &str) -> PathBuf {
    let job = folder.join(name);
    cli(
        &[
            cmd("prepare-mc-epoch-srm-gate-job"),
            &folder.join("mesh.json"),
            &folder.join("ground.json"),
            &folder.join("case.json"),
            &folder.join("spec.json"),
            cmd("1"),
            &job,
        ],
        0,
    );
    job
}
fn worker(job: &Path) -> PathBuf {
    job.join("solver/soil-fem-3d.exe")
}
fn run(job: &Path, budget: Option<&str>, expected: i32) {
    let mut args = vec![cmd("run-mc-epoch-srm-gate-job"), job];
    if let Some(budget) = budget {
        args.push(cmd(budget));
    }
    execute(&worker(job), &args, expected);
}
fn scientific(value: &Value) -> Value {
    let mut value = value.clone();
    fn clean(value: &mut Value) {
        match value {
            Value::Object(map) => {
                map.remove("elapsed_seconds");
                map.remove("origin_run_sha256");
                map.remove("payload_sha256");
                for child in map.values_mut() {
                    clean(child);
                }
            }
            Value::Array(items) => {
                for child in items {
                    clean(child);
                }
            }
            _ => {}
        }
    }
    clean(&mut value);
    value
}
fn direct(folder: &Path, exe: &Path) -> PathBuf {
    let path = folder.join("direct-gate.json");
    execute(
        exe,
        &[
            cmd("solve-mc-epoch-srm-gate"),
            &folder.join("mesh.json"),
            &folder.join("ground.json"),
            &folder.join("case.json"),
            &folder.join("spec.json"),
            &path,
        ],
        0,
    );
    path
}

#[test]
fn paused_gate_recovery_matches_the_continuous_gate_and_the_unreduced_solve() {
    let folder = fixture("roundtrip", 1.1, 128);
    let job = prepare(&folder, "first");
    run(&job, Some("1"), 2);
    let info = execute(
        &worker(&job),
        &[cmd("verify-mc-epoch-srm-gate-job"), &job],
        0,
    );
    let info: Value = serde_json::from_slice(&info.stdout).unwrap();
    assert_eq!(info["complete"], false);
    assert!(info["generation"].as_u64().unwrap() > 1);
    let resumed = folder.join("recovered");
    execute(
        &worker(&job),
        &[cmd("prepare-mc-epoch-srm-gate-recovery"), &job, &resumed],
        0,
    );
    run(&resumed, None, 0);
    let report = direct(&folder, &worker(&resumed));
    assert_eq!(
        scientific(&read(&report)),
        scientific(&read(&resumed.join("report.json")))
    );
    let solved = folder.join("continuous.json");
    execute(
        &worker(&resumed),
        &[
            cmd("solve-stages"),
            &folder.join("mesh.json"),
            &folder.join("ground.json"),
            &folder.join("case.json"),
            &solved,
        ],
        0,
    );
    let gate = read(&resumed.join("report.json"));
    assert_eq!(gate["decision"], "passed");
    assert_eq!(
        scientific(&gate["construction"]),
        scientific(&read(&solved))
    );
    assert_eq!(gate["gates"][0]["effect"], "continue");
    assert_eq!(
        gate["continuum_stability_or_design_safety_certified"],
        false
    );
    run(&resumed, None, 1);
}

#[test]
fn exclusion_and_nonconvergence_remain_distinct_and_do_not_start_the_later_stage() {
    let excluded = fixture("excluded", 1.6, 128);
    let job = prepare(&excluded, "first");
    run(&job, None, 0);
    let report = read(&job.join("report.json"));
    assert_eq!(report["decision"], "stopped_excluded");
    assert_eq!(report["gates"].as_array().unwrap().len(), 1);
    assert_eq!(report["gates"][0]["effect"], "stop_excluded");
    assert_eq!(
        report["construction"]["case"]["stages"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(report["construction"]["complete"], true);
    assert!(report["construction"]["failure_reason"].is_null());
    assert_eq!(report["stopped_before_later_stages"], true);
    assert_eq!(
        report["nonconvergence_used_as_exclusion_or_collapse"],
        false
    );

    let unknown = fixture("unknown", 1.1, 1);
    let job = prepare(&unknown, "first");
    run(&job, None, 0);
    let report = read(&job.join("report.json"));
    assert_eq!(read(&job.join("status.json"))["state"], "completed");
    assert_eq!(report["decision"], "stopped_indeterminate");
    assert_eq!(report["gates"][0]["requirement"], "indeterminate");
    assert_eq!(report["gates"][0]["effect"], "stop_indeterminate");
    assert_ne!(report["decision"], "stopped_excluded");
    assert_eq!(
        report["construction"]["case"]["stages"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        report["nonconvergence_used_as_exclusion_or_collapse"],
        false
    );
}

#[test]
fn prepared_stop_invalid_budget_and_second_launch_keep_the_saved_checkpoint() {
    let folder = fixture("control", 1.1, 128);
    let untouched = prepare(&folder, "budget");
    let before = checkpoint::hash_file(&untouched.join("recovery.json")).unwrap();
    run(&untouched, Some("0"), 1);
    run(&untouched, Some("10001"), 1);
    assert_eq!(
        before,
        checkpoint::hash_file(&untouched.join("recovery.json")).unwrap()
    );
    assert!(!untouched.join("run.claim").exists());

    let job = prepare(&folder, "first");
    let pointer = fs::read(job.join("recovery.json")).unwrap();
    execute(&worker(&job), &[cmd("stop-mc-epoch-srm-gate-job"), &job], 0);
    run(&job, None, 2);
    assert_eq!(pointer, fs::read(job.join("recovery.json")).unwrap());
    run(&job, None, 1);
    let resumed = folder.join("resumed");
    execute(
        &worker(&job),
        &[cmd("prepare-mc-epoch-srm-gate-recovery"), &job, &resumed],
        0,
    );
    run(&resumed, None, 0);
    run(&resumed, None, 1);
}
