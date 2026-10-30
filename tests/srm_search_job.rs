use serde_json::{Value, json};
use soil_fem_3d::{
    Result, checkpoint, ground::Ground, mesh::Mesh, soil, srm_bound::Problem,
    srm_search::Specification, stages::StageCase,
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
fn write(path: &Path, v: &impl serde::Serialize) {
    fs::write(path, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn directory(tag: &str) -> PathBuf {
    let parent = root().join("runs");
    fs::create_dir_all(&parent).unwrap();
    let path = parent.join(format!(
        "SRM job CLI {tag} {}-{}-{}",
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
fn inputs(tag: &str, unknown: bool) -> PathBuf {
    let path = directory(tag);
    let (mesh, mut g, mut c, mut spec) = if unknown {
        plastic_column().unwrap()
    } else {
        fixture().unwrap()
    };
    spec.absolute_factor_tolerance = 0.04;
    if unknown {
        g.initialization
            .geostatic
            .as_mut()
            .unwrap()
            .verification_reject_once = false;
        c.settings.initial_increment = 1.;
        c.settings.minimum_increment = 1.;
        c.settings.maximum_newton_iterations = 3;
        spec.required_factor = 4.;
        spec.maximum_factor = 4.;
    } else {
        c.verification_reject_once_stages = vec!["load".into()];
    }
    write(&path.join("mesh.json"), &mesh);
    write(&path.join("ground.json"), &g);
    write(&path.join("case.json"), &c);
    write(&path.join("spec.json"), &spec);
    path
}
fn prepare(path: &Path, name: &str) -> PathBuf {
    let folder = path.join(name);
    cli(
        &[
            cmd("prepare-srm-search-job"),
            &path.join("mesh.json"),
            &path.join("ground.json"),
            &path.join("case.json"),
            &path.join("spec.json"),
            cmd("1"),
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
    let mut args = vec![cmd("run-srm-search-job"), job];
    if let Some(b) = budget {
        args.push(cmd(b));
    }
    execute(&worker(job), &args, expected);
}
fn recover(source: &Path, dest: &Path) {
    execute(
        &worker(source),
        &[cmd("prepare-srm-search-recovery"), source, dest],
        0,
    );
}
fn verify(job: &Path, expected: i32) -> Value {
    let out = execute(&worker(job), &[cmd("verify-srm-search-job"), job], expected);
    if expected == 0 {
        serde_json::from_slice(&out.stdout).unwrap()
    } else {
        Value::Null
    }
}
fn hashes(path: &Path) -> BTreeMap<PathBuf, String> {
    fn visit(path: &Path, base: &Path, out: &mut BTreeMap<PathBuf, String>) {
        for entry in fs::read_dir(path).unwrap() {
            let q = entry.unwrap().path();
            if q.is_dir() {
                visit(&q, base, out);
            } else {
                out.insert(
                    q.strip_prefix(base).unwrap().to_path_buf(),
                    checkpoint::hash_file(&q).unwrap(),
                );
            }
        }
    }
    let mut out = BTreeMap::new();
    visit(path, path, &mut out);
    out
}
fn scientific(mut v: Value) -> Value {
    fn clean(v: &mut Value) {
        match v {
            Value::Object(map) => {
                map.remove("elapsed_seconds");
                for v in map.values_mut() {
                    clean(v);
                }
            }
            Value::Array(values) => {
                for v in values {
                    clean(v);
                }
            }
            _ => {}
        }
    }
    clean(&mut v);
    v
}
fn wait(mut child: Child) -> i32 {
    let deadline = Instant::now() + Duration::from_secs(45);
    loop {
        if let Some(code) = child.try_wait().unwrap() {
            return code.code().unwrap();
        }
        if Instant::now() >= deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("owned tiny SRM test exceeded45seconds");
        }
        std::thread::sleep(Duration::from_millis(5));
    }
}

fn fixture() -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let p: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
    raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    raw["materials"] = json!([{"id":1,"name":"Unconfined proof-of-search",
        "unit_weight_kn_m3":1e-10,"mohr_coulomb":p.material}]);
    raw["bodies"] = json!([{"body":"TestCube","material_id":1}]);
    let ground: Ground = serde_json::from_value(raw)?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(p.boundary);
    case.stages = serde_json::from_value(
        json!([{"kind":"boundary","name":"load","target_factor":1.0,"initial_increment":1.0}]),
    )?;
    case.verification_reject_once_stages = vec![];
    case.monitor_sets.clear();
    case.monitors.clear();
    let spec: Specification = serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    Ok((mesh, ground, case, spec))
}

fn plastic_column() -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    for m in &mut ground.materials {
        m.unit_weight_kn_m3 = 1e-10;
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            panic!("MC fixture")
        };
        p.cohesion_kpa = 40.;
        p.friction_deg = 0.;
        p.dilation_deg = 0.;
        p.tension_cutoff_kpa = None;
    }
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.stages = serde_json::from_value(
        json!([{"kind":"pressure","name":"load","target_pressure_kpa":100.0}]),
    )?;
    case.verification_reject_once_stages = vec![];
    let (_, _, _, mut spec) = fixture()?;
    spec.material_ids = ground.materials.iter().map(|m| m.id).collect();
    spec.maximum_factor = 4.;
    spec.required_factor = 2.;
    Ok((mesh, ground, case, spec))
}

#[test]
fn actual_cli_alternating_snapshots_pause_recovery_and_one_shot_match_closed_form() {
    let path = inputs("closed form", false);
    let full = prepare(&path, "continuous");
    run(&full, None, 0);
    let partial = prepare(&path, "paused");
    run(&partial, Some("2"), 2);
    let info = verify(&partial, 0);
    assert_eq!(info["complete"], false);
    assert!(info["generation"].as_u64().unwrap() > 1);
    assert!(
        partial.join("recovery-checkpoint-0.json").is_file()
            && partial.join("recovery-checkpoint-1.json").is_file()
    );
    let before = hashes(&partial);
    let recovered = path.join("recovered");
    recover(&partial, &recovered);
    run(&recovered, None, 0);
    assert_eq!(before, hashes(&partial));
    let control = path.join("uninterrupted standalone.json");
    cli(
        &[
            cmd("search-mc-srm"),
            &path.join("mesh.json"),
            &path.join("ground.json"),
            &path.join("case.json"),
            &path.join("spec.json"),
            &control,
        ],
        0,
    );
    let expected = scientific(read(&control));
    let actual = read(&recovered.join("report.json"));
    assert_eq!(scientific(read(&full.join("report.json"))), expected);
    assert_eq!(scientific(actual.clone()), expected);
    assert!(
        actual["interval"]["attained_equilibrium_factor"]
            .as_f64()
            .unwrap()
            < 20. / 15.
    );
    assert!(
        actual["interval"]["independently_excluded_factor"]
            .as_f64()
            .unwrap()
            > 20. / 15.
    );
    assert_eq!(actual["requirement"], "equilibrium_requirement_met");
    let status = read(&recovered.join("status.json"));
    assert_eq!(status["state"], "completed");
    assert_eq!(status["requirement"], actual["requirement"]);
    assert_eq!(
        status["report_sha256"],
        checkpoint::hash_file(&recovered.join("report.json")).unwrap()
    );
    assert_eq!(verify(&recovered, 0)["complete"], true);
    let completed = hashes(&recovered);
    run(&recovered, None, 1);
    assert_eq!(completed, hashes(&recovered));
}

#[test]
fn real_nonlinear_failure_completes_an_indeterminate_search_without_an_upper() {
    let path = inputs("nonlinear unknown", true);
    let partial = prepare(&path, "before trial");
    run(&partial, Some("1"), 2);
    let recovered = path.join("resume unknown");
    recover(&partial, &recovered);
    run(&recovered, None, 0);
    let report = read(&recovered.join("report.json"));
    let status = read(&recovered.join("status.json"));
    assert_eq!(status["state"], "completed");
    assert_eq!(status["search_stop_reason"], "indeterminate_trial");
    assert_eq!(status["requirement"], "indeterminate");
    assert!(report["interval"].is_null());
    let failed = report["trials"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| v["factor"] == 4.)
        .unwrap();
    assert_eq!(failed["outcome"], "indeterminate");
    assert_eq!(failed["failure_kind"], "incomplete");
    let reason = failed["run"]["failure_reason"].as_str().unwrap();
    assert!(reason.contains("Newton iteration limit") && reason.contains("cutback exhausted"));
    assert!(
        failed["run"]["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .any(|a| a["accepted"] == false && a["committed_state_unchanged_on_rejection"] == true)
    );
    assert!(
        !report["trials"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v["outcome"] == "physically_excluded")
    );
    assert_eq!(
        report["continuum_stability_or_design_safety_certified"],
        false
    );
    let control = path.join("standalone unknown.json");
    cli(
        &[
            cmd("search-mc-srm"),
            &path.join("mesh.json"),
            &path.join("ground.json"),
            &path.join("case.json"),
            &path.join("spec.json"),
            &control,
        ],
        0,
    );
    assert_eq!(scientific(report), scientific(read(&control)));
}

#[test]
fn stop_before_run_and_stop_during_actual_mechanics_preserve_resumable_history() {
    let path = inputs("cooperative stop", false);
    let early = prepare(&path, "stop before run");
    execute(&worker(&early), &[cmd("stop-srm-search-job"), &early], 0);
    run(&early, None, 2);
    assert_eq!(verify(&early, 0)["generation"], 1);
    let mut c: StageCase = serde_json::from_value(read(&path.join("case.json"))).unwrap();
    c.verification_reject_once_stages.clear();
    c.stages=serde_json::from_value(Value::Array((1..=100).map(|i|json!({"kind":"boundary","name":format!("load{i}"),"target_factor":i as f64/100.,"initial_increment":1.})).collect())).unwrap();
    write(&path.join("case.json"), &c);
    let mut spec: Specification = serde_json::from_value(read(&path.join("spec.json"))).unwrap();
    spec.through_stages = 100;
    write(&path.join("spec.json"), &spec);
    let job = prepare(&path, "running mechanics");
    let mut child = Command::new(worker(&job))
        .arg("run-srm-search-job")
        .arg(&job)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(30);
    loop {
        assert!(
            child.try_wait().unwrap().is_none(),
            "actual job finished before cooperative stop could be tested"
        );
        if let Ok(raw) = fs::read(job.join("status.json"))
            && let Ok(s) = serde_json::from_slice::<Value>(&raw)
            && s["state"] == "running"
            && s["progress_events"].as_u64().unwrap_or(0) > 0
        {
            break;
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("no actual mechanical progress observed");
        }
        std::thread::sleep(Duration::from_millis(5));
    }
    execute(&worker(&job), &[cmd("stop-srm-search-job"), &job], 0);
    assert_eq!(wait(child), 2);
    assert_eq!(read(&job.join("status.json"))["state"], "paused");
    let info = verify(&job, 0);
    let pointer = read(&job.join("recovery.json"));
    let payload =
        read(&job.join(pointer["checkpoint"].as_str().unwrap()))["checkpoint"]["payload"].clone();
    assert!(
        info["accepted_mechanical_checkpoint"] == true || payload["original"]["complete"] == true
    );
    assert!(info["generation"].as_u64().unwrap() > 1);
    let original = hashes(&job);
    let restored = path.join("accepted resumed");
    recover(&job, &restored);
    run(&restored, Some("1"), 2);
    assert_eq!(original, hashes(&job));
    verify(&restored, 0);
}

#[test]
fn foreign_pointer_path_corruption_wrong_stop_and_immutable_changes_are_rejected() {
    let path = inputs("integrity", false);
    let a = prepare(&path, "a");
    let b = prepare(&path, "b");
    run(&a, Some("1"), 2);
    run(&b, Some("1"), 2);
    let pointer = fs::read(a.join("recovery.json")).unwrap();
    fs::copy(b.join("recovery.json"), a.join("recovery.json")).unwrap();
    verify(&a, 1);
    fs::write(a.join("recovery.json"), &pointer).unwrap();
    let mut bad = read(&a.join("recovery.json"));
    bad["checkpoint"] = json!("../outside.json");
    write(&a.join("recovery.json"), &bad);
    verify(&a, 1);
    fs::write(a.join("recovery.json"), &pointer).unwrap();
    let selected = a.join(
        read(&a.join("recovery.json"))["checkpoint"]
            .as_str()
            .unwrap(),
    );
    let old = fs::read(&selected).unwrap();
    fs::write(&selected, b"truncated snapshot").unwrap();
    verify(&a, 1);
    let refusal = path.join("must not create recovery");
    execute(
        &worker(&a),
        &[cmd("prepare-srm-search-recovery"), &a, &refusal],
        1,
    );
    assert!(!refusal.exists());
    fs::write(&selected, old).unwrap();
    verify(&a, 0);
    let changed = prepare(&path, "immutable");
    fs::write(changed.join("inputs/specification.json"), b"modified").unwrap();
    run(&changed, None, 1);
    assert!(!changed.join("run.claim").exists());
    let bad_budget = prepare(&path, "invalid budget");
    run(&bad_budget, Some("0"), 1);
    assert!(!bad_budget.join("run.claim").exists());
    let bad_solver = prepare(&path, "changed solver");
    use std::io::Write;
    fs::OpenOptions::new()
        .append(true)
        .open(worker(&bad_solver))
        .unwrap()
        .write_all(b"modified frozen executable")
        .unwrap();
    cli(&[cmd("run-srm-search-job"), &bad_solver], 1);
    assert!(!bad_solver.join("run.claim").exists());
    let stop = prepare(&path, "foreign stop");
    write(
        &stop.join("stop-request.json"),
        &json!({"schema":"soil-fem-srm-search-job-stop/1","job_id":"someone else"}),
    );
    run(&stop, None, 1);
    assert_eq!(read(&stop.join("status.json"))["state"], "failed");
    assert_eq!(verify(&stop, 0)["generation"], 1);
    let recovered = path.join("clean stop recovery");
    recover(&stop, &recovered);
    run(&recovered, Some("1"), 2);
    verify(&recovered, 0);
}

#[test]
fn snapshot_write_failure_keeps_previous_pointer_and_completed_report_can_be_recovered() {
    let path = inputs("write failure", false);
    let job = prepare(&path, "write fails");
    let before = fs::read(job.join("recovery.json")).unwrap();
    fs::create_dir(job.join("recovery-checkpoint-0.json")).unwrap();
    run(&job, None, 1);
    assert_eq!(fs::read(job.join("recovery.json")).unwrap(), before);
    assert_eq!(read(&job.join("status.json"))["state"], "failed");
    assert_eq!(verify(&job, 0)["generation"], 1);
    let restored = path.join("restored write");
    recover(&job, &restored);
    run(&restored, None, 0);
    let report = fs::read(restored.join("report.json")).unwrap();
    fs::write(
        restored.join("report.json"),
        b"damaged after final snapshot",
    )
    .unwrap();
    assert_eq!(verify(&restored, 0)["complete"], true);
    let finalcopy = path.join("regenerated report");
    recover(&restored, &finalcopy);
    run(&finalcopy, None, 0);
    assert_eq!(fs::read(finalcopy.join("report.json")).unwrap(), report);
    assert_eq!(
        read(&finalcopy.join("status.json"))["requirement"],
        "equilibrium_requirement_met"
    );
    assert_eq!(read(&finalcopy.join("status.json"))["progress_events"], 0);
}

struct BoundaryStop {
    queries: std::sync::atomic::AtomicUsize,
    events: std::sync::atomic::AtomicUsize,
    defer_queries: usize,
}
impl BoundaryStop {
    fn after(defer_queries: usize) -> Self {
        Self {
            queries: std::sync::atomic::AtomicUsize::new(0),
            events: std::sync::atomic::AtomicUsize::new(0),
            defer_queries,
        }
    }
}
impl soil_fem_3d::stage_control::Observer for BoundaryStop {
    fn on_progress(&self, _: &soil_fem_3d::stage_control::Progress) -> Result<()> {
        self.events.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(self.queries.fetch_add(1, Ordering::Relaxed) >= self.defer_queries)
    }
}

#[test]
fn stop_already_requested_at_completed_geostatic_checkpoint_keeps_every_saved_byte() -> Result<()> {
    use soil_fem_3d::srm_search_control::Snapshot;
    let (mesh, g, c, spec) = fixture()?;
    let saved = Snapshot::start(&mesh, &g, &c, &spec)?.advance(&mesh, &g, &c, &spec, 1, None)?;
    let raw = serde_json::to_value(&saved)?;
    assert_eq!(
        raw["payload"]["inflight"]["run"]["geostatic"]["complete"],
        true
    );
    assert!(saved.has_accepted_mechanical_checkpoint());
    let before = serde_json::to_vec(&saved)?;
    let stop = BoundaryStop::after(0);
    let after = saved.advance(&mesh, &g, &c, &spec, 1, Some(&stop))?;
    assert_eq!(serde_json::to_vec(&after)?, before);
    assert_eq!(stop.events.load(Ordering::Relaxed), 0);
    Ok(())
}

#[test]
fn fresh_stop_after_final_geostatic_acceptance_is_a_checkpoint_not_a_failed_origin() -> Result<()> {
    use soil_fem_3d::srm_search_control::Snapshot;
    let (mesh, g, c, spec) = fixture()?;
    let stop = BoundaryStop::after(1);
    let saved =
        Snapshot::start(&mesh, &g, &c, &spec)?.advance(&mesh, &g, &c, &spec, 100, Some(&stop))?;
    let raw = serde_json::to_value(&saved)?;
    assert!(saved.has_accepted_mechanical_checkpoint());
    assert_eq!(
        raw["payload"]["inflight"]["run"]["failure_reason"],
        "checkpoint_requested"
    );
    assert_eq!(
        raw["payload"]["inflight"]["run"]["geostatic"]["complete"],
        true
    );
    assert_eq!(
        raw["payload"]["inflight"]["run"]["last_committed"]["stage_index"],
        0
    );
    assert!(stop.events.load(Ordering::Relaxed) > 0);
    Ok(())
}

fn plastic_checkpoint() -> Result<(
    Mesh,
    Ground,
    StageCase,
    Specification,
    soil_fem_3d::srm_search_control::Snapshot,
)> {
    use soil_fem_3d::srm_search_control::Snapshot;
    let (mesh, mut g, mut c, mut spec) = plastic_column()?;
    g.initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    for m in &mut g.materials {
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            unreachable!()
        };
        p.cohesion_kpa = 20.;
    }
    c.settings.initial_increment = 0.25;
    spec.required_factor = 1.;
    spec.maximum_factor = 1.1;
    let mut saved = Snapshot::start(&mesh, &g, &c, &spec)?;
    for _ in 0..12 {
        saved = saved.advance(&mesh, &g, &c, &spec, 1, None)?;
        let raw = serde_json::to_value(&saved)?;
        let run = &raw["payload"]["inflight"]["run"];
        if run["last_committed"]["stage_index"] == 1
            && run["last_committed"]["stage_progress"]
                .as_f64()
                .is_some_and(|x| x < 1.)
            && run["last_committed"]["material_states"]
                .as_array()
                .is_some_and(|a| {
                    a.iter()
                        .any(|s| s["equivalent_plastic_strain"].as_f64().unwrap_or(0.) > 1e-9)
                })
        {
            return Ok((mesh, g, c, spec, saved));
        }
    }
    panic!("expected a real accepted plastic construction state before its endpoint")
}

#[test]
fn stop_before_resuming_plastic_construction_keeps_history_and_matches_uninterrupted_result()
-> Result<()> {
    let (mesh, g, c, spec, saved) = plastic_checkpoint()?;
    let bytes = serde_json::to_vec(&saved)?;
    let stop = BoundaryStop::after(0);
    let mut after = saved.advance(&mesh, &g, &c, &spec, 1, Some(&stop))?;
    assert_eq!(serde_json::to_vec(&after)?, bytes);
    assert_eq!(stop.events.load(Ordering::Relaxed), 0);
    for _ in 0..100 {
        if after.is_complete() {
            break;
        }
        after = after.advance(&mesh, &g, &c, &spec, 1, None)?;
    }
    assert!(after.is_complete());
    let expected = soil_fem_3d::srm_search::run(&mesh, g, c, spec)?;
    assert_eq!(
        scientific(serde_json::to_value(after.report().unwrap())?),
        scientific(serde_json::to_value(expected)?)
    );
    Ok(())
}

#[test]
fn stop_arriving_after_resume_precheck_never_rewinds_accepted_plastic_construction() -> Result<()> {
    let (mesh, g, c, spec, saved) = plastic_checkpoint()?;
    let before = serde_json::to_value(&saved)?;
    let before = &before["payload"]["inflight"]["run"];
    let stop = BoundaryStop::after(1);
    let after = saved.advance(&mesh, &g, &c, &spec, 1, Some(&stop))?;
    let raw = serde_json::to_value(after)?;
    let run = &raw["payload"]["inflight"]["run"];
    assert_eq!(run["failure_reason"], "checkpoint_requested");
    assert_eq!(run["last_committed"]["stage_index"], 1);
    assert!(
        run["last_committed"]["stage_progress"].as_f64().unwrap()
            > before["last_committed"]["stage_progress"].as_f64().unwrap()
    );
    assert!(
        run["attempts"].as_array().unwrap().len() > before["attempts"].as_array().unwrap().len()
    );
    let old = before["last_committed"]["material_states"]
        .as_array()
        .unwrap();
    let new = run["last_committed"]["material_states"].as_array().unwrap();
    assert!(
        old.iter()
            .any(|s| s["equivalent_plastic_strain"].as_f64().unwrap() > 1e-9)
    );
    for (old, new) in old.iter().zip(new) {
        assert!(
            new["equivalent_plastic_strain"].as_f64().unwrap()
                >= old["equivalent_plastic_strain"].as_f64().unwrap()
        );
    }
    assert!(stop.events.load(Ordering::Relaxed) > 0);
    Ok(())
}
