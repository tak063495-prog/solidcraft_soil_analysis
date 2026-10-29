use serde_json::{Value, json};
use soil_fem_3d::{
    checkpoint,
    ground::Ground,
    mc::{Flow, Parameters},
    mc_epoch_path::Specification,
    mesh::Mesh,
    soil,
    stages::StageCase,
};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    process::{Child, Command, Output, Stdio},
    sync::atomic::{AtomicU64, Ordering},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};
static COUNTER: AtomicU64 = AtomicU64::new(0);
fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_soil-fem-3d"))
}
fn cmd(s: &str) -> &Path {
    Path::new(s)
}
fn execute(exe: &Path, args: &[&Path], expected: i32) -> Output {
    let out = Command::new(exe).args(args).output().unwrap();
    assert_eq!(
        out.status.code(),
        Some(expected),
        "{:?}\n{}",
        args,
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
fn fixture(tag: &str, tension: bool, increment: f64) -> PathBuf {
    let parent = root().join("runs");
    fs::create_dir_all(&parent).unwrap();
    let path = parent.join(format!(
        "mc epoch CLI {tag} {}-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&path).unwrap();
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json")).unwrap();
    let mesh = mesh.elevate_straight().unwrap();
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json")).unwrap();
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json")).unwrap();
    for m in &mut ground.materials {
        m.unit_weight_kn_m3 = 1e-10;
        m.model = soil::Model::MohrCoulomb(Parameters {
            young_kpa: 30000.,
            poisson: 0.25,
            cohesion_kpa: if tension { 20. } else { 40. },
            friction_deg: if tension { 30. } else { 0. },
            dilation_deg: 0.,
            tension_cutoff_kpa: None,
            flow: Flow::NonAssociated,
        });
    }
    if tension {
        let constraints = (0..3)
            .map(|axis| {
                let mut direction = [0.; 3];
                direction[axis] = 1.;
                let mut gradient = [0.; 3];
                gradient[axis] = 1. / 6000.;
                soil_fem_3d::boundary::Constraint {
                    node_set: "all".into(),
                    direction,
                    gradient,
                    value_m: 0.,
                    reference_point_m: [0.; 3],
                }
            })
            .collect();
        case.boundary = Some(soil_fem_3d::boundary::Conditions {
            constraints,
            nodal_loads: vec![],
            surface_loads: vec![],
        });
        case.stages =
            serde_json::from_value(json!([{"kind":"boundary","name":"extend","target_factor":1.}]))
                .unwrap();
    } else {
        case.stages = serde_json::from_value(
            json!([{"kind":"pressure","name":"load","target_pressure_kpa":100.}]),
        )
        .unwrap();
    }
    let target_parameters = ground
        .materials
        .iter()
        .map(|m| {
            let soil::Model::MohrCoulomb(mut p) = m.model.clone() else {
                panic!("MC")
            };
            p.cohesion_kpa = if tension { 0. } else { 20. };
            p
        })
        .collect();
    let spec = Specification {
        schema: "soil-fem-mc-epoch-path/1".into(),
        target_parameters,
        initial_increment: increment,
        minimum_increment: 1. / 4096.,
        maximum_attempts: 4096,
    };
    write(&path.join("mesh.json"), &mesh);
    write(&path.join("ground.json"), &ground);
    write(&path.join("case.json"), &case);
    write(&path.join("spec.json"), &spec);
    cli(
        &[
            cmd("solve-stages"),
            &path.join("mesh.json"),
            &path.join("ground.json"),
            &path.join("case.json"),
            &path.join("origin.json"),
        ],
        0,
    );
    path
}
fn prepare(path: &Path, name: &str) -> PathBuf {
    let folder = path.join(name);
    cli(
        &[
            cmd("prepare-mc-epoch-job"),
            &path.join("mesh.json"),
            &path.join("origin.json"),
            &path.join("spec.json"),
            &folder,
        ],
        0,
    );
    folder
}
fn worker(job: &Path) -> PathBuf {
    job.join("solver/soil-fem-3d.exe")
}
fn run(job: &Path, budget: Option<&str>, expected: i32) {
    let mut a = vec![cmd("run-mc-epoch-job"), job];
    if let Some(b) = budget {
        a.push(cmd(b));
    }
    execute(&worker(job), &a, expected);
}
fn tree_hashes(folder: &Path) -> BTreeMap<PathBuf, String> {
    fn visit(path: &Path, base: &Path, out: &mut BTreeMap<PathBuf, String>) {
        for e in fs::read_dir(path).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                visit(&p, base, out);
            } else {
                out.insert(
                    p.strip_prefix(base).unwrap().to_path_buf(),
                    checkpoint::hash_file(&p).unwrap(),
                );
            }
        }
    }
    let mut out = BTreeMap::new();
    visit(folder, folder, &mut out);
    out
}
fn checkpoint_value(job: &Path) -> Value {
    read(&job.join("report.json"))["checkpoint"].clone()
}

#[test]
fn product_cli_pause_recovery_and_durable_checkpoint_match_continuous() {
    let path = fixture("plastic restart", false, 0.25);
    let continuous = prepare(&path, "continuous");
    run(&continuous, None, 0);
    let partial = prepare(&path, "paused plastic");
    run(&partial, Some("2"), 2);
    assert_eq!(read(&partial.join("status.json"))["state"], "paused");
    assert_eq!(
        checkpoint_value(&partial)["payload"]["continuation"]["progress"],
        0.5
    );
    let before = tree_hashes(&partial);
    let recovered = path.join("recovered plastic");
    cli(&[cmd("prepare-mc-epoch-recovery"), &partial, &recovered], 0);
    run(&recovered, None, 0);
    assert_eq!(checkpoint_value(&continuous), checkpoint_value(&recovered));
    assert_eq!(before, tree_hashes(&partial));
    execute(
        &worker(&recovered),
        &[cmd("verify-mc-epoch-job"), &recovered],
        0,
    );
    let report = read(&recovered.join("report.json"));
    assert_eq!(report["complete"], true);
    assert_eq!(report["stability_classification"], "not_evaluated");
    let mesh: Mesh = serde_json::from_slice(&fs::read(path.join("mesh.json")).unwrap()).unwrap();
    let displacement=report["checkpoint"]["payload"]["material_checkpoint"]["payload"]["committed"]["displacement_m"].as_array().unwrap();
    for (p, u) in mesh.nodes.iter().zip(displacement) {
        assert!((u[2].as_f64().unwrap() + 110. / 30000. * p[2]).abs() < 1e-9);
    }
}

#[test]
fn material_failure_rejected_prefix_recovers_without_inventing_stability_bounds() {
    let path = fixture("material failure", true, 1.);
    let full = prepare(&path, "incomplete continuous");
    run(&full, None, 1);
    let partial = prepare(&path, "after rejection");
    run(&partial, Some("1"), 2);
    let cp = checkpoint_value(&partial);
    let attempt = &cp["payload"]["continuation"]["attempts"][0];
    assert_eq!(attempt["failure_kind"], "material");
    assert_eq!(
        attempt["before_state_sha256"],
        attempt["after_state_sha256"]
    );
    assert_eq!(
        attempt["before_parameters_sha256"],
        attempt["after_parameters_sha256"]
    );
    let recovered = path.join("recovered rejection");
    cli(&[cmd("prepare-mc-epoch-recovery"), &partial, &recovered], 0);
    run(&recovered, None, 1);
    assert_eq!(checkpoint_value(&full), checkpoint_value(&recovered));
    let r = read(&recovered.join("report.json"));
    assert_eq!(r["complete"], false);
    assert_eq!(r["stop_reason"], "cutback_exhausted");
    assert_eq!(r["stability_classification"], "not_evaluated");
    let p = r["checkpoint"]["payload"]["continuation"]["progress"]
        .as_f64()
        .unwrap();
    let limit = 1. - 10. * 30_f64.to_radians().tan() / 20.;
    assert!(p <= limit + 1e-10 && limit - p < 2. / 4096.);
    execute(
        &worker(&recovered),
        &[cmd("verify-mc-epoch-job"), &recovered],
        0,
    );
}

#[test]
fn stop_before_run_one_shot_immutable_inputs_and_pointer_paths_are_guarded() {
    let path = fixture("guards", false, 0.25);
    let stopped = prepare(&path, "prepared stop");
    execute(&worker(&stopped), &[cmd("stop-mc-epoch-job"), &stopped], 0);
    run(&stopped, None, 2);
    assert_eq!(
        checkpoint_value(&stopped)["payload"]["continuation"]["attempts"]
            .as_array()
            .unwrap()
            .len(),
        0
    );
    let protected = tree_hashes(&stopped);
    run(&stopped, None, 1);
    assert_eq!(protected, tree_hashes(&stopped));
    let recovered = path.join("resume initial");
    cli(&[cmd("prepare-mc-epoch-recovery"), &stopped, &recovered], 0);
    run(&recovered, None, 0);
    let invalid = prepare(&path, "invalid budget");
    run(&invalid, Some("0"), 1);
    assert!(!invalid.join("run.claim").exists());
    let modified = prepare(&path, "modified input");
    fs::write(modified.join("inputs/specification.json"), b"changed").unwrap();
    run(&modified, None, 1);
    assert!(!modified.join("run.claim").exists());
    let escaped = prepare(&path, "bad pointer");
    let mut pointer = read(&escaped.join("recovery.json"));
    pointer["checkpoint"] = json!("../origin.json");
    write(&escaped.join("recovery.json"), &pointer);
    execute(
        &worker(&escaped),
        &[cmd("verify-mc-epoch-job"), &escaped],
        1,
    );
    run(&escaped, None, 1);
    assert!(!escaped.join("run.claim").exists());
    let existing = tree_hashes(&recovered);
    cli(&[cmd("prepare-mc-epoch-recovery"), &stopped, &recovered], 1);
    assert_eq!(existing, tree_hashes(&recovered));
}

struct OwnedChild(Child);
impl Drop for OwnedChild {
    fn drop(&mut self) {
        if self.0.try_wait().ok().flatten().is_none() {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
}
#[test]
fn owned_abrupt_exit_recovers_the_last_complete_plastic_epoch() {
    let path = fixture("owned abrupt loss", false, 1. / 128.);
    let full = prepare(&path, "continuous128");
    run(&full, None, 0);
    let abrupt = prepare(&path, "owned process");
    let mut child = OwnedChild(
        Command::new(worker(&abrupt))
            .args([cmd("run-mc-epoch-job"), &abrupt])
            .stdout(Stdio::from(
                fs::File::create(path.join("owned.stdout")).unwrap(),
            ))
            .stderr(Stdio::from(
                fs::File::create(path.join("owned.stderr")).unwrap(),
            ))
            .spawn()
            .unwrap(),
    );
    let start = Instant::now();
    loop {
        assert!(
            child.0.try_wait().unwrap().is_none(),
            "owned job ended before the planned crash point"
        );
        let pointer = read(&abrupt.join("recovery.json"));
        if pointer["generation"].as_u64().unwrap() >= 70 {
            break;
        }
        assert!(
            start.elapsed() < Duration::from_secs(60),
            "owned job did not reach its accepted plastic checkpoint"
        );
        std::thread::sleep(Duration::from_millis(5));
    }
    child.0.kill().unwrap();
    assert!(!child.0.wait().unwrap().success());
    assert_eq!(read(&abrupt.join("status.json"))["state"], "running");
    let output = execute(&worker(&abrupt), &[cmd("verify-mc-epoch-job"), &abrupt], 0);
    let info: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(info["progress"].as_f64().unwrap() > 0.5);
    let before = tree_hashes(&abrupt);
    let recovered = path.join("new process recovery");
    cli(&[cmd("prepare-mc-epoch-recovery"), &abrupt, &recovered], 0);
    run(&recovered, None, 0);
    assert_eq!(checkpoint_value(&full), checkpoint_value(&recovered));
    assert_eq!(before, tree_hashes(&abrupt));
}
