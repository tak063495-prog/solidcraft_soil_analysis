use serde_json::{Value, json};
use soil_fem_3d::{
    Result, checkpoint, epoch_srm_search::Specification, ground::Ground, mesh::Mesh, soil,
    srm_bound::Problem, stages::StageCase,
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
fn cmd(s: &str) -> &Path {
    Path::new(s)
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
fn read(q: &Path) -> Value {
    serde_json::from_slice(&fs::read(q).unwrap()).unwrap()
}
fn write(q: &Path, v: &impl serde::Serialize) {
    fs::write(q, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn directory(tag: &str) -> PathBuf {
    let parent = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("runs");
    fs::create_dir_all(&parent).unwrap();
    let q = parent.join(format!(
        "epoch SRM job {tag} {}-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&q).unwrap();
    q
}
fn inputs() -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let p: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
    raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    raw["materials"] = json!([{"id":1,"name":"unconfined epoch SRM","unit_weight_kn_m3":1e-10,"mohr_coulomb":p.material}]);
    raw["bodies"] = json!([{"body":"TestCube","material_id":1}]);
    let ground = serde_json::from_value(raw)?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(p.boundary);
    case.stages = serde_json::from_value(
        json!([{"kind":"boundary","name":"load","target_factor":1.,"initial_increment":1.}]),
    )?;
    case.verification_reject_once_stages.clear();
    case.monitors.clear();
    case.monitor_sets.clear();
    let search = serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    let spec = Specification {
        schema: "soil-fem-epoch-srm-search/1".into(),
        search,
        initial_increment: 0.25,
        minimum_increment: 1. / 4096.,
        maximum_attempts: 128,
    };
    Ok((mesh, ground, case, spec))
}
fn plastic_inputs(cut: bool) -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut g: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    for m in &mut g.materials {
        m.unit_weight_kn_m3 = 1e-10;
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            panic!("MC")
        };
        p.young_kpa = 100000.;
        p.poisson = 0.3;
        p.cohesion_kpa = 40.;
        p.friction_deg = 0.;
        p.dilation_deg = 0.;
        p.tension_cutoff_kpa = None;
    }
    let mut c: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    let mut steps = vec![
        json!({"kind":"pressure","name":"plastic_preload","target_pressure_kpa":200.}),
        json!({"kind":"pressure","name":"unload","target_pressure_kpa":150.}),
    ];
    if cut {
        steps.push(json!({"kind":"excavate","name":"cut_top","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.0}));
    }
    steps.push(json!({"kind":"reset_display","name":"reference"}));
    c.stages = serde_json::from_value(json!(steps))?;
    c.verification_reject_once_stages.clear();
    let (_, _, _, mut spec) = inputs()?;
    spec.search.material_ids = g.materials.iter().map(|m| m.id).collect();
    spec.search.required_factor = 2.;
    spec.search.maximum_factor = 3.;
    spec.search.through_stages = c.stages.len();
    Ok((mesh, g, c, spec))
}

fn fixture(tag: &str, plastic: bool, long: bool) -> PathBuf {
    let q = directory(tag);
    let (mesh, g, c, mut s) = if plastic {
        plastic_inputs(false).unwrap()
    } else {
        inputs().unwrap()
    };
    if long {
        s.initial_increment = 0.02;
        s.maximum_attempts = 512;
    }
    write(&q.join("mesh.json"), &mesh);
    write(&q.join("ground.json"), &g);
    write(&q.join("case.json"), &c);
    write(&q.join("spec.json"), &s);
    cli(
        &[
            cmd("solve-stages"),
            &q.join("mesh.json"),
            &q.join("ground.json"),
            &q.join("case.json"),
            &q.join("origin.json"),
        ],
        0,
    );
    assert_eq!(read(&q.join("origin.json"))["complete"], true);
    q
}
fn prepare(q: &Path, name: &str) -> PathBuf {
    let job = q.join(name);
    cli(
        &[
            cmd("prepare-mc-epoch-srm-job"),
            &q.join("mesh.json"),
            &q.join("origin.json"),
            &q.join("spec.json"),
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
    let mut args = vec![cmd("run-mc-epoch-srm-job"), job];
    if let Some(b) = budget {
        args.push(cmd(b));
    }
    execute(&worker(job), &args, expected);
}
fn recover(job: &Path, out: &Path) {
    execute(
        &worker(job),
        &[cmd("prepare-mc-epoch-srm-recovery"), job, out],
        0,
    );
}
fn verify(job: &Path, expected: i32) -> Value {
    let out = execute(
        &worker(job),
        &[cmd("verify-mc-epoch-srm-job"), job],
        expected,
    );
    if expected == 0 {
        serde_json::from_slice(&out.stdout).unwrap()
    } else {
        Value::Null
    }
}
fn expected(q: &Path) -> PathBuf {
    let path = q.join("direct report.json");
    cli(
        &[
            cmd("search-mc-epoch-srm"),
            &q.join("mesh.json"),
            &q.join("origin.json"),
            &q.join("spec.json"),
            &path,
        ],
        0,
    );
    path
}
fn hashes(q: &Path) -> BTreeMap<PathBuf, String> {
    fn visit(q: &Path, base: &Path, out: &mut BTreeMap<PathBuf, String>) {
        for entry in fs::read_dir(q).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(&path, base, out);
            } else {
                out.insert(
                    path.strip_prefix(base).unwrap().to_path_buf(),
                    checkpoint::hash_file(&path).unwrap(),
                );
            }
        }
    }
    let mut out = BTreeMap::new();
    visit(q, q, &mut out);
    out
}
fn same_bytes(a: &Path, b: &Path) {
    assert_eq!(fs::read(a).unwrap(), fs::read(b).unwrap());
}
fn spawn(job: &Path) -> Child {
    Command::new(worker(job))
        .arg("run-mc-epoch-srm-job")
        .arg(job)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap()
}
fn wait(mut child: Child) -> i32 {
    let deadline = Instant::now() + Duration::from_secs(45);
    loop {
        if let Some(code) = child.try_wait().unwrap() {
            return code.code().unwrap();
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("owned tiny job timeout");
        }
        std::thread::sleep(Duration::from_millis(2));
    }
}
fn wait_for_plastic_factor(child: &mut Child, job: &Path) {
    let deadline = Instant::now() + Duration::from_secs(30);
    loop {
        assert!(
            child.try_wait().unwrap().is_none(),
            "job ended before the actual in-progress second factor was observed"
        );
        if let Ok(raw) = fs::read(job.join("status.json"))
            && let Ok(s) = serde_json::from_slice::<Value>(&raw)
            && s["state"] == "running"
            && s["current_factor"] == 2.
            && s["recovery_generation"].as_u64().unwrap_or(0) > 55
            && s["progress_events"].as_u64().unwrap_or(0) > 0
        {
            return;
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("no second-factor accepted history observed");
        }
        std::thread::sleep(Duration::from_millis(2));
    }
}
#[test]
fn actual_cli_two_slot_pause_fork_and_complete_report_are_byte_exact_to_direct_search() {
    let q = fixture("roundtrip", false, false);
    let direct = expected(&q);
    let job = prepare(&q, "first");
    run(&job, Some("2"), 2);
    let info = verify(&job, 0);
    assert_eq!(info["complete"], false);
    assert!(info["generation"].as_u64().unwrap() > 1);
    assert!(
        job.join("recovery-checkpoint-0.json").exists()
            && job.join("recovery-checkpoint-1.json").exists()
    );
    let before = hashes(&job);
    let resumed = q.join("recovered");
    recover(&job, &resumed);
    run(&resumed, None, 0);
    assert_eq!(before, hashes(&job));
    same_bytes(&direct, &resumed.join("report.json"));
    let done = hashes(&resumed);
    run(&resumed, None, 1);
    assert_eq!(done, hashes(&resumed));
    let interval = read(&direct)["interval"].clone();
    assert!(interval["attained_equilibrium_factor"].as_f64().unwrap() < 4. / 3.);
    assert!(interval["independently_excluded_factor"].as_f64().unwrap() > 4. / 3.);
    assert!(interval["width"].as_f64().unwrap() <= 0.002);
}
#[test]
fn prepared_stop_and_existing_output_refusal_keep_the_exact_initial_checkpoint() {
    let q = fixture("early stop", false, false);
    let job = prepare(&q, "first");
    let pointer = fs::read(job.join("recovery.json")).unwrap();
    execute(&worker(&job), &[cmd("stop-mc-epoch-srm-job"), &job], 0);
    run(&job, None, 2);
    assert_eq!(pointer, fs::read(job.join("recovery.json")).unwrap());
    assert_eq!(verify(&job, 0)["generation"], 1);
    let before = hashes(&job);
    cli(
        &[
            cmd("prepare-mc-epoch-srm-job"),
            cmd("missing mesh"),
            cmd("missing origin"),
            cmd("missing spec"),
            cmd("1"),
            &job,
        ],
        1,
    );
    assert_eq!(before, hashes(&job));
    let resumed = q.join("resumed");
    recover(&job, &resumed);
    run(&resumed, None, 0);
}
#[test]
fn stop_during_actual_plastic_material_continuation_preserves_and_resumes_all_report_bytes() {
    let q = fixture("nonlinear cooperative stop", true, true);
    let direct = expected(&q);
    let job = prepare(&q, "first");
    let mut child = spawn(&job);
    wait_for_plastic_factor(&mut child, &job);
    execute(&worker(&job), &[cmd("stop-mc-epoch-srm-job"), &job], 0);
    assert_eq!(wait(child), 2);
    assert_eq!(read(&job.join("status.json"))["state"], "paused");
    assert_eq!(verify(&job, 0)["saved_material_epoch"], true);
    let before = hashes(&job);
    let resumed = q.join("resumed");
    recover(&job, &resumed);
    run(&resumed, None, 0);
    assert_eq!(before, hashes(&job));
    same_bytes(&direct, &resumed.join("report.json"));
}
#[test]
fn killed_owned_solver_recovers_only_its_last_fully_published_accepted_snapshot() {
    let q = fixture("actual process kill", true, true);
    let direct = expected(&q);
    let job = prepare(&q, "first");
    let mut child = spawn(&job);
    wait_for_plastic_factor(&mut child, &job);
    child.kill().unwrap();
    child.wait().unwrap();
    let info = verify(&job, 0);
    assert_eq!(info["complete"], false);
    assert_eq!(info["current_factor"], 2.);
    let pointer = read(&job.join("recovery.json"));
    let cp = read(&job.join(pointer["checkpoint"].as_str().unwrap()));
    let progress = cp["checkpoint"]["payload"]["inflight"]["payload"]["continuation"]["progress"]
        .as_f64()
        .unwrap();
    assert!(progress > 0. && progress < 1.);
    let before = hashes(&job);
    let resumed = q.join("resumed");
    recover(&job, &resumed);
    run(&resumed, None, 0);
    assert_eq!(before, hashes(&job));
    same_bytes(&direct, &resumed.join("report.json"));
}
#[test]
fn foreign_pointer_truncated_selected_slot_and_resealed_semantic_edit_are_refused() {
    let q = fixture("corruption", false, false);
    let foreign = prepare(&q, "foreign");
    let mut p = read(&foreign.join("recovery.json"));
    p["checkpoint"] = json!("../outside.json");
    write(&foreign.join("recovery.json"), &p);
    let before = hashes(&foreign);
    verify(&foreign, 1);
    assert_eq!(before, hashes(&foreign));
    let truncated = prepare(&q, "truncated");
    run(&truncated, Some("2"), 2);
    let p = read(&truncated.join("recovery.json"));
    assert!(
        truncated.join("recovery-checkpoint-0.json").exists()
            && truncated.join("recovery-checkpoint-1.json").exists()
    );
    fs::write(truncated.join(p["checkpoint"].as_str().unwrap()), b"{").unwrap();
    verify(&truncated, 1);
    let edited = prepare(&q, "resealed");
    let mut pointer = read(&edited.join("recovery.json"));
    let slot = edited.join(pointer["checkpoint"].as_str().unwrap());
    let mut saved = read(&slot);
    saved["checkpoint"]["payload"]["basis"]["physical_force_kn"][0] = json!(0.123);
    let mut typed: soil_fem_3d::epoch_srm_control::Snapshot =
        serde_json::from_value(saved["checkpoint"].clone()).unwrap();
    // Hash the typed serializer's exact payload field order. Value's sorted
    // object order would test only a stale seal, not semantic revalidation.
    use sha2::Digest;
    let compact = serde_json::to_string(&typed).unwrap();
    let start = compact.find("\"payload\":").unwrap() + "\"payload\":".len();
    let seal = format!(
        "{:x}",
        sha2::Sha256::digest(&compact.as_bytes()[start..compact.len() - 1])
    );
    let mut value = serde_json::to_value(&typed).unwrap();
    value["payload_sha256"] = json!(seal);
    typed = serde_json::from_value(value).unwrap();
    saved["checkpoint"] = serde_json::to_value(&typed).unwrap();
    write(&slot, &saved);
    pointer["checkpoint_sha256"] = json!(checkpoint::hash_file(&slot).unwrap());
    pointer["payload_sha256"] = json!(checkpoint::hash_json(&typed).unwrap());
    write(&edited.join("recovery.json"), &pointer);
    verify(&edited, 1);
}
#[test]
fn immutable_origin_and_wrong_job_stop_request_fail_without_replacing_the_saved_pointer() {
    let q = fixture("immutable bindings", false, false);
    let job = prepare(&q, "first");
    let mut origin = read(&job.join("inputs/origin.json"));
    origin["last_committed"]["applied_nodal_load_kn"][0] = json!(0.123);
    write(&job.join("inputs/origin.json"), &origin);
    verify(&job, 1);
    let other = prepare(&q, "wrong stop");
    let before = fs::read(other.join("recovery.json")).unwrap();
    write(
        &other.join("stop-request.json"),
        &json!({"schema":"soil-fem-epoch-srm-search-job-stop/1","job_id":"wrong"}),
    );
    run(&other, None, 1);
    assert_eq!(before, fs::read(other.join("recovery.json")).unwrap());
    verify(&other, 0);
}
#[test]
fn interrupted_part_files_are_ignored_and_invalid_budgets_do_not_claim_the_job() {
    let q = fixture("incomplete writes", false, false);
    let job = prepare(&q, "first");
    fs::write(job.join(".recovery-data-incomplete.part"), b"{").unwrap();
    let before = hashes(&job);
    run(&job, Some("0"), 1);
    run(&job, Some("10001"), 1);
    assert_eq!(before, hashes(&job));
    assert!(!job.join("run.claim").exists());
    verify(&job, 0);
    run(&job, None, 0);
    let direct = expected(&q);
    same_bytes(&direct, &job.join("report.json"));
}
#[test]
fn complete_job_can_be_forked_without_recomputing_or_changing_its_report() {
    let q = fixture("finished fork", false, false);
    let job = prepare(&q, "first");
    run(&job, None, 0);
    let before = hashes(&job);
    let fork = q.join("finished fork");
    recover(&job, &fork);
    run(&fork, None, 0);
    assert_eq!(before, hashes(&job));
    same_bytes(&job.join("report.json"), &fork.join("report.json"));
    assert_eq!(read(&fork.join("status.json"))["progress_events"], 0);
}
#[test]
fn numeric_attempt_limit_remains_unknown_even_when_the_disk_job_finishes() {
    let q = fixture("unknown", false, false);
    let mut s: Specification = serde_json::from_value(read(&q.join("spec.json"))).unwrap();
    s.maximum_attempts = 1;
    write(&q.join("spec.json"), &s);
    let job = prepare(&q, "first");
    run(&job, None, 0);
    let r = read(&job.join("report.json"));
    assert_eq!(read(&job.join("status.json"))["state"], "completed");
    assert_eq!(r["requirement"], "indeterminate");
    assert_eq!(r["stop_reason"], "indeterminate_trial");
    assert!(r["interval"].is_null());
    assert_eq!(r["continuum_stability_or_design_safety_certified"], false);
    assert!(
        r["trials"]
            .as_array()
            .unwrap()
            .iter()
            .any(|t| t["evidence"]["outcome"] == "indeterminate")
    );
}
#[test]
fn concurrent_launches_allow_exactly_one_worker_claim() {
    let q = fixture("one worker", false, false);
    let job = prepare(&q, "first");
    let a = spawn(&job);
    let b = spawn(&job);
    let mut codes = [wait(a), wait(b)];
    codes.sort();
    assert_eq!(codes, [0, 1]);
    verify(&job, 0);
    let direct = expected(&q);
    same_bytes(&direct, &job.join("report.json"));
}
