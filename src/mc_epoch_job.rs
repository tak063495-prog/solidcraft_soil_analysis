//! One-shot durable material-continuation jobs. Checkpoints always contain an
//! accepted mechanical state, including after rejected parameter trials.
use crate::{
    Error, Result,
    checkpoint::{self, Identity},
    mc_epoch_path::{self, Controller, Specification, StopReason},
    mesh::Mesh,
    stage_control::{Observer, Progress},
    stage_job::{copy_checked, read_small, replace_status},
    stages::Run,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    sync::{
        Mutex,
        atomic::{AtomicU64, Ordering},
    },
    time::{Instant, SystemTime, UNIX_EPOCH},
};
const MESH: &str = "inputs/mesh.json";
const ORIGIN: &str = "inputs/origin.json";
const SPEC: &str = "inputs/specification.json";
const SAVED: &str = "inputs/checkpoint.json";
const SOLVER: &str = "solver/soil-fem-3d.exe";
static TEMPORARIES: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parent {
    pub job_id: String,
    pub generation: u64,
    pub snapshot_sha256: String,
    pub checkpoint_payload_sha256: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Binding {
    pub identity: Identity,
    pub specification_sha256: String,
    pub files_sha256: BTreeMap<String, String>,
    pub nonce: String,
    pub parent: Option<Parent>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub schema: String,
    pub job_id: String,
    pub binding: Binding,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum State {
    Prepared,
    Running,
    Paused,
    Succeeded,
    Incomplete,
    Failed,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Status {
    pub schema: String,
    pub job_id: String,
    pub state: State,
    pub process_id: Option<u32>,
    pub parameter_progress: f64,
    pub trial_parameter_progress: Option<f64>,
    pub attempts: usize,
    pub accepted_epochs: usize,
    pub recovery_generation: u64,
    pub progress_events: u64,
    pub last_progress: Option<Progress>,
    pub report: Option<String>,
    pub reason: Option<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Succeeded,
    Paused,
    Incomplete,
    Failed,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub job_id: String,
    pub complete: bool,
    pub stability_classification: String,
    pub stop_reason: Option<String>,
    pub checkpoint: mc_epoch_path::Checkpoint,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Pointer {
    schema: String,
    job_id: String,
    generation: u64,
    checkpoint: String,
    checkpoint_sha256: String,
    payload_sha256: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Snapshot {
    schema: String,
    job_id: String,
    generation: u64,
    checkpoint: mc_epoch_path::Checkpoint,
}
#[derive(Serialize)]
struct SnapshotRef<'a> {
    schema: &'static str,
    job_id: &'a str,
    generation: u64,
    checkpoint: &'a mc_epoch_path::Checkpoint,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SnapshotInfo {
    pub schema: String,
    pub job_id: String,
    pub generation: u64,
    pub checkpoint_sha256: String,
    pub payload_sha256: String,
    pub progress: f64,
    pub attempts: usize,
    pub accepted_epochs: usize,
    pub complete: bool,
    pub stop_reason: Option<StopReason>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StopRequest {
    schema: String,
    job_id: String,
}

fn read<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T> {
    Ok(serde_json::from_reader(BufReader::with_capacity(
        1024 * 1024,
        File::open(path)?,
    ))?)
}
fn synced_new(path: &Path, value: &impl Serialize) -> Result<()> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut stream = BufWriter::with_capacity(1024 * 1024, file);
    serde_json::to_writer_pretty(&mut stream, value)?;
    stream.flush()?;
    stream.get_ref().sync_all()?;
    Ok(())
}
fn temporary(folder: &Path, label: &str) -> PathBuf {
    folder.join(format!(
        ".{label}-{}-{}.part",
        std::process::id(),
        TEMPORARIES.fetch_add(1, Ordering::Relaxed)
    ))
}
fn atomic_json(folder: &Path, name: &str, value: &impl Serialize) -> Result<()> {
    let temporary = temporary(folder, name);
    synced_new(&temporary, value)?;
    replace_status(&temporary, &folder.join(name))
}
fn copied(source: &Path, target: &Path) -> Result<String> {
    let digest = copy_checked(source, target)?;
    OpenOptions::new().write(true).open(target)?.sync_all()?;
    Ok(digest)
}
fn manifest(folder: &Path) -> Result<Manifest> {
    let m: Manifest = read_small(&folder.join("job.json"), 256 * 1024)?;
    if m.schema != "soil-fem-mc-epoch-job/1" || m.job_id != checkpoint::hash_json(&m.binding)? {
        return Err(Error::Input("MC epoch manifest identity mismatch".into()));
    }
    Ok(m)
}
fn verify_files(folder: &Path, m: &Manifest) -> Result<()> {
    let mut expected = BTreeSet::from([MESH, ORIGIN, SPEC, SOLVER]);
    if m.binding.parent.is_some() {
        expected.insert(SAVED);
    }
    if m.binding
        .files_sha256
        .keys()
        .map(String::as_str)
        .collect::<BTreeSet<_>>()
        != expected
    {
        return Err(Error::Input(
            "MC epoch immutable file inventory mismatch".into(),
        ));
    }
    for (name, hash) in &m.binding.files_sha256 {
        if checkpoint::hash_file(&folder.join(name))? != *hash {
            return Err(Error::Input(format!(
                "MC epoch immutable input differs: {name}"
            )));
        }
    }
    Ok(())
}
fn load(folder: &Path) -> Result<(Manifest, Mesh, Run, Specification)> {
    let m = manifest(folder)?;
    verify_files(folder, &m)?;
    let mesh: Mesh = read(&folder.join(MESH))?;
    let origin: Run = read(&folder.join(ORIGIN))?;
    let spec: Specification = read_small(&folder.join(SPEC), 1024 * 1024)?;
    if m.binding.identity != origin.provenance
        || m.binding.identity != Identity::current(&mesh, &origin.ground, &origin.case)?
        || m.binding.identity.solver_binary_sha256 != m.binding.files_sha256[SOLVER]
        || m.binding.specification_sha256 != checkpoint::hash_json(&spec)?
    {
        return Err(Error::Input(
            "MC epoch origin/canonical inputs/frozen executable mismatch".into(),
        ));
    }
    Ok((m, mesh, origin, spec))
}
fn save(
    folder: &Path,
    job_id: &str,
    generation: u64,
    cp: &mc_epoch_path::Checkpoint,
) -> Result<()> {
    if generation == 0 {
        return Err(Error::Control(
            "MC recovery generation must be positive".into(),
        ));
    }
    let name = format!("recovery-checkpoint-{}.json", generation % 2);
    let data = temporary(folder, "recovery-data");
    synced_new(
        &data,
        &SnapshotRef {
            schema: "soil-fem-mc-epoch-snapshot/1",
            job_id,
            generation,
            checkpoint: cp,
        },
    )?;
    let hash = checkpoint::hash_file(&data)?;
    replace_status(&data, &folder.join(&name))?;
    let pointer = Pointer {
        schema: "soil-fem-mc-epoch-recovery/1".into(),
        job_id: job_id.into(),
        generation,
        checkpoint: name,
        checkpoint_sha256: hash,
        payload_sha256: cp.payload_sha256.clone(),
    };
    atomic_json(folder, "recovery.json", &pointer)
}
fn selected(
    folder: &Path,
    m: &Manifest,
    mesh: &Mesh,
    origin: &Run,
    spec: &Specification,
) -> Result<(Pointer, mc_epoch_path::Checkpoint)> {
    for _ in 0..3 {
        let pointer: Pointer = read_small(&folder.join("recovery.json"), 8192)?;
        if pointer.schema != "soil-fem-mc-epoch-recovery/1"
            || pointer.job_id != m.job_id
            || pointer.generation == 0
            || pointer.checkpoint != format!("recovery-checkpoint-{}.json", pointer.generation % 2)
        {
            return Err(Error::Input(
                "MC recovery pointer identity/path mismatch".into(),
            ));
        }
        let result = (|| -> Result<mc_epoch_path::Checkpoint> {
            let path = folder.join(&pointer.checkpoint);
            if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
                return Err(Error::Input("MC snapshot file checksum mismatch".into()));
            }
            let s: Snapshot = read(&path)?;
            if s.schema != "soil-fem-mc-epoch-snapshot/1"
                || s.job_id != m.job_id
                || s.generation != pointer.generation
                || s.checkpoint.payload_sha256 != pointer.payload_sha256
                || checkpoint::hash_json(&s.checkpoint.payload.specification)?
                    != checkpoint::hash_json(spec)?
            {
                return Err(Error::Input("MC snapshot binding mismatch".into()));
            }
            Controller::resume(mesh, origin, &s.checkpoint)?;
            if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
                return Err(Error::Input("MC snapshot changed while read".into()));
            }
            Ok(s.checkpoint)
        })();
        match result {
            Ok(cp) => return Ok((pointer, cp)),
            Err(e) => {
                let now: Pointer = read_small(&folder.join("recovery.json"), 8192)?;
                if now == pointer {
                    return Err(e);
                }
            }
        }
    }
    Err(Error::Control(
        "MC snapshot advanced during all3 reads; retain the source and retry".into(),
    ))
}
fn prepare_impl(
    mesh_file: &Path,
    origin_file: &Path,
    spec_file: &Path,
    output: &Path,
    parent: Option<(Parent, &mc_epoch_path::Checkpoint)>,
) -> Result<Manifest> {
    if output.exists() {
        return Err(
            std::io::Error::new(std::io::ErrorKind::AlreadyExists, "MC job output exists").into(),
        );
    }
    let mesh: Mesh = read(mesh_file)?;
    let origin: Run = read(origin_file)?;
    let spec: Specification = read_small(spec_file, 1024 * 1024)?;
    let controller = if let Some((_, cp)) = &parent {
        if checkpoint::hash_json(&cp.payload.specification)? != checkpoint::hash_json(&spec)? {
            return Err(Error::Input("MC recovery specification differs".into()));
        }
        Controller::resume(&mesh, &origin, cp)?
    } else {
        Controller::new(&mesh, &origin, spec.clone())?
    };
    let initial = controller.checkpoint()?;
    fs::create_dir(output)?;
    fs::create_dir(output.join("inputs"))?;
    fs::create_dir(output.join("solver"))?;
    let mut files = BTreeMap::new();
    for (source, name) in [(mesh_file, MESH), (origin_file, ORIGIN), (spec_file, SPEC)] {
        files.insert(name.into(), copied(source, &output.join(name))?);
    }
    files.insert(
        SOLVER.into(),
        copied(&std::env::current_exe()?, &output.join(SOLVER))?,
    );
    if parent.is_some() {
        synced_new(&output.join(SAVED), &initial)?;
        files.insert(SAVED.into(), checkpoint::hash_file(&output.join(SAVED))?);
    }
    let binding = Binding {
        identity: origin.provenance.clone(),
        specification_sha256: checkpoint::hash_json(&spec)?,
        files_sha256: files,
        nonce: format!(
            "{}-{}-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|e| Error::Control(e.to_string()))?
                .as_nanos(),
            std::process::id(),
            TEMPORARIES.fetch_add(1, Ordering::Relaxed)
        ),
        parent: parent.map(|(p, _)| p),
    };
    let m = Manifest {
        schema: "soil-fem-mc-epoch-job/1".into(),
        job_id: checkpoint::hash_json(&binding)?,
        binding,
    };
    synced_new(&output.join("job.json"), &m)?;
    load(output)?;
    save(output, &m.job_id, 1, &initial)?;
    let c = controller.continuation();
    let status = Status {
        schema: "soil-fem-mc-epoch-job-status/1".into(),
        job_id: m.job_id.clone(),
        state: State::Prepared,
        process_id: None,
        parameter_progress: c.progress,
        trial_parameter_progress: None,
        attempts: c.attempts.len(),
        accepted_epochs: c.accepted.len(),
        recovery_generation: 1,
        progress_events: 0,
        last_progress: None,
        report: None,
        reason: None,
    };
    synced_new(&output.join("status.json"), &status)?;
    Ok(m)
}
pub fn prepare(mesh: &Path, origin: &Path, spec: &Path, output: &Path) -> Result<Manifest> {
    prepare_impl(mesh, origin, spec, output, None)
}
pub fn prepare_recovery(previous: &Path, output: &Path) -> Result<Manifest> {
    let folder = fs::canonicalize(previous)?;
    let (m, mesh, origin, spec) = load(&folder)?;
    let (pointer, cp) = selected(&folder, &m, &mesh, &origin, &spec)?;
    let parent = Parent {
        job_id: m.job_id,
        generation: pointer.generation,
        snapshot_sha256: pointer.checkpoint_sha256,
        checkpoint_payload_sha256: cp.payload_sha256.clone(),
    };
    prepare_impl(
        &folder.join(MESH),
        &folder.join(ORIGIN),
        &folder.join(SPEC),
        output,
        Some((parent, &cp)),
    )
}
pub fn status(folder: &Path) -> Result<Status> {
    let m = manifest(folder)?;
    let s: Status = read_small(&folder.join("status.json"), 256 * 1024)?;
    if s.schema != "soil-fem-mc-epoch-job-status/1" || s.job_id != m.job_id {
        return Err(Error::Input("MC job status identity mismatch".into()));
    }
    Ok(s)
}
fn requested(folder: &Path, job_id: &str) -> Result<bool> {
    let path = folder.join("stop-request.json");
    if !path.exists() {
        return Ok(false);
    }
    let r: StopRequest = read_small(&path, 4096)?;
    if r.schema != "soil-fem-mc-epoch-job-stop/1" || r.job_id != job_id {
        return Err(Error::Control(
            "MC stop request belongs to another job".into(),
        ));
    }
    Ok(true)
}
pub fn request_stop(folder: &Path) -> Result<bool> {
    let s = status(folder)?;
    if !matches!(s.state, State::Prepared | State::Running) {
        return Ok(false);
    }
    let path = folder.join("stop-request.json");
    if path.exists() {
        requested(folder, &s.job_id)?;
    } else {
        synced_new(
            &path,
            &StopRequest {
                schema: "soil-fem-mc-epoch-job-stop/1".into(),
                job_id: s.job_id,
            },
        )?;
    }
    Ok(true)
}
pub fn verify_snapshot(folder: &Path) -> Result<SnapshotInfo> {
    let (m, mesh, origin, spec) = load(folder)?;
    let (pointer, cp) = selected(folder, &m, &mesh, &origin, &spec)?;
    let c = &cp.payload.continuation;
    Ok(SnapshotInfo {
        schema: "soil-fem-mc-epoch-snapshot-verification/1".into(),
        job_id: m.job_id,
        generation: pointer.generation,
        checkpoint_sha256: pointer.checkpoint_sha256,
        payload_sha256: cp.payload_sha256.clone(),
        progress: c.progress,
        attempts: c.attempts.len(),
        accepted_epochs: c.accepted.len(),
        complete: c.progress == 1.,
        stop_reason: c.stop_reason,
    })
}

struct Trace {
    status: Status,
    events: BufWriter<File>,
    last_status: Instant,
}
struct FileObserver {
    folder: PathBuf,
    job_id: String,
    trace: Mutex<Trace>,
}
impl FileObserver {
    fn trial(&self, progress: f64) -> Result<()> {
        let mut t = self
            .trace
            .lock()
            .map_err(|_| Error::Control("MC progress lock poisoned".into()))?;
        t.status.trial_parameter_progress = Some(progress);
        Ok(())
    }
    fn update(&self, c: &mc_epoch_path::Continuation, generation: u64) -> Result<()> {
        let mut t = self
            .trace
            .lock()
            .map_err(|_| Error::Control("MC progress lock poisoned".into()))?;
        t.status.parameter_progress = c.progress;
        t.status.attempts = c.attempts.len();
        t.status.accepted_epochs = c.accepted.len();
        t.status.recovery_generation = generation;
        t.status.trial_parameter_progress = None;
        t.events.flush()?;
        atomic_json(&self.folder, "status.json", &t.status)
    }
    fn finish(&self, state: State, reason: Option<String>) -> Result<()> {
        let mut t = self
            .trace
            .lock()
            .map_err(|_| Error::Control("MC progress lock poisoned".into()))?;
        t.status.state = state;
        t.status.reason = reason;
        t.status.report = Some("report.json".into());
        t.events.flush()?;
        t.events.get_ref().sync_all()?;
        atomic_json(&self.folder, "status.json", &t.status)
    }
}
impl Observer for FileObserver {
    fn on_progress(&self, event: &Progress) -> Result<()> {
        let mut t = self
            .trace
            .lock()
            .map_err(|_| Error::Control("MC progress lock poisoned".into()))?;
        let mut event = event.clone();
        event.stage = "mc_parameter_epoch".into();
        event.progress = t
            .status
            .trial_parameter_progress
            .unwrap_or(t.status.parameter_progress);
        serde_json::to_writer(&mut t.events, &event)?;
        t.events.write_all(b"\n")?;
        t.status.progress_events = t
            .status
            .progress_events
            .checked_add(1)
            .ok_or_else(|| Error::Control("MC event counter overflow".into()))?;
        t.status.last_progress = Some(event);
        if t.last_status.elapsed().as_millis() >= 100 {
            t.events.flush()?;
            atomic_json(&self.folder, "status.json", &t.status)?;
            t.last_status = Instant::now();
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        requested(&self.folder, &self.job_id)
    }
}

pub fn run(folder: &Path, attempt_budget: Option<usize>) -> Result<Outcome> {
    if attempt_budget.is_some_and(|n| !(1..=10000).contains(&n)) {
        return Err(Error::Input("MC run attempt budget must be1..10000".into()));
    }
    let folder = fs::canonicalize(folder)?;
    let (m, mesh, origin, spec) = load(&folder)?;
    let mut state = status(&folder)?;
    if state.state != State::Prepared
        || ["run.claim", "report.json", "progress.jsonl"]
            .iter()
            .any(|n| folder.join(n).exists())
    {
        return Err(Error::Input(
            "MC job is one-shot; fork its accepted recovery snapshot into a new directory".into(),
        ));
    }
    let (pointer, cp) = selected(&folder, &m, &mesh, &origin, &spec)?;
    let initial = if m.binding.parent.is_some() {
        read::<mc_epoch_path::Checkpoint>(&folder.join(SAVED))?
    } else {
        Controller::new(&mesh, &origin, spec.clone())?.checkpoint()?
    };
    if pointer.generation != 1 || checkpoint::hash_json(&cp)? != checkpoint::hash_json(&initial)? {
        return Err(Error::Input(
            "prepared MC job must start from its exact initial snapshot".into(),
        ));
    }
    synced_new(
        &folder.join("run.claim"),
        &serde_json::json!({"job_id":m.job_id,"process_id":std::process::id()}),
    )?;
    state.state = State::Running;
    state.process_id = Some(std::process::id());
    atomic_json(&folder, "status.json", &state)?;
    let observer = FileObserver {
        folder: folder.clone(),
        job_id: m.job_id.clone(),
        trace: Mutex::new(Trace {
            status: state,
            events: BufWriter::new(
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(folder.join("progress.jsonl"))?,
            ),
            last_status: Instant::now(),
        }),
    };
    let result = (|| -> Result<Outcome> {
        let mut controller = Controller::resume(&mesh, &origin, &cp)?;
        controller.set_observer(&observer);
        let mut generation = pointer.generation;
        let start = controller.continuation().attempts.len();
        loop {
            let c = controller.continuation();
            let terminal = if controller.complete() {
                Some((Outcome::Succeeded, State::Succeeded, None))
            } else if let Some(stop) = c.stop_reason {
                let runtime = stop == StopReason::NonRecoverableTrial;
                let reason = match stop {
                    StopReason::CutbackExhausted => "cutback_exhausted",
                    StopReason::AttemptLimit => "attempt_limit",
                    StopReason::NonRecoverableTrial => "non_recoverable_trial",
                };
                Some((
                    if runtime {
                        Outcome::Failed
                    } else {
                        Outcome::Incomplete
                    },
                    if runtime {
                        State::Failed
                    } else {
                        State::Incomplete
                    },
                    Some(reason.into()),
                ))
            } else if observer.stop_requested()? {
                Some((
                    Outcome::Paused,
                    State::Paused,
                    Some("stop_requested".into()),
                ))
            } else if attempt_budget.is_some_and(|n| c.attempts.len() - start >= n) {
                Some((
                    Outcome::Paused,
                    State::Paused,
                    Some("attempt_budget".into()),
                ))
            } else {
                None
            };
            if let Some((outcome, state, reason)) = terminal {
                verify_files(&folder, &m)?;
                let report = Report {
                    schema: "soil-fem-mc-epoch-job-result/1".into(),
                    job_id: m.job_id.clone(),
                    complete: controller.complete(),
                    stability_classification: controller.stability_classification().into(),
                    stop_reason: reason.clone(),
                    checkpoint: controller.checkpoint()?,
                };
                synced_new(&folder.join("report.json"), &report)?;
                observer.finish(state, reason)?;
                return Ok(outcome);
            }
            observer.trial((c.progress + c.next_increment.min(1. - c.progress)).min(1.))?;
            controller.advance()?;
            verify_files(&folder, &m)?;
            if controller.continuation().stop_reason == Some(StopReason::NonRecoverableTrial) {
                // Runtime/control failures do not poison the recoverable path.
                // The failed report retains the diagnostic, while the durable
                // pointer stays at the preceding complete accepted epoch.
                observer.update(controller.continuation(), generation)?;
                continue;
            }
            generation = generation
                .checked_add(1)
                .ok_or_else(|| Error::Control("MC generation overflow".into()))?;
            save(&folder, &m.job_id, generation, &controller.checkpoint()?)?;
            observer.update(controller.continuation(), generation)?;
        }
    })();
    if let Err(error) = &result {
        let mut state = status(&folder)?;
        state.state = State::Failed;
        state.reason = Some(error.to_string());
        // A write failure never destroys or relabels the previous complete snapshot.
        let _ = atomic_json(&folder, "status.json", &state);
    }
    result
}
