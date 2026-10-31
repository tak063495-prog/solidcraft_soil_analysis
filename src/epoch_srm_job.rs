//! Immutable one-shot accepted-origin epoch SRM jobs. A completed search may still be indeterminate;
//! job completion is never a safety gate. Recovery selects only a verified,
//! fully written snapshot, including its same-factor accepted material history.
//! Each factor starts from the exact same completed unreduced construction Run.
use crate::{
    Error, Result,
    checkpoint::{self, Identity},
    epoch_srm_control::Snapshot as SearchSnapshot,
    epoch_srm_search::Specification,
    mesh::Mesh,
    srm_search::{Requirement, StopReason},
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
    pub material_attempts_per_save: usize,
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
    Completed,
    Failed,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Status {
    pub schema: String,
    pub job_id: String,
    pub state: State,
    pub process_id: Option<u32>,
    pub current_factor: Option<f64>,
    pub completed_trials: usize,
    pub saved_material_epoch: bool,
    pub recovery_generation: u64,
    pub progress_events: u64,
    pub last_progress: Option<Progress>,
    pub requirement: Option<Requirement>,
    pub search_stop_reason: Option<StopReason>,
    pub report: Option<String>,
    pub report_sha256: Option<String>,
    pub reason: Option<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Completed,
    Paused,
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
struct StoredSnapshot {
    schema: String,
    job_id: String,
    generation: u64,
    checkpoint: SearchSnapshot,
}
#[derive(Serialize)]
struct SnapshotRef<'a> {
    schema: &'static str,
    job_id: &'a str,
    generation: u64,
    checkpoint: &'a SearchSnapshot,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SnapshotInfo {
    pub schema: String,
    pub job_id: String,
    pub generation: u64,
    pub checkpoint_sha256: String,
    pub canonical_snapshot_sha256: String,
    pub current_factor: Option<f64>,
    pub completed_trials: usize,
    pub saved_material_epoch: bool,
    pub complete: bool,
    pub requirement: Option<Requirement>,
    pub stop_reason: Option<StopReason>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StopRequest {
    schema: String,
    job_id: String,
}
struct Context {
    manifest: Manifest,
    mesh: Mesh,
    origin: Run,
    spec: Specification,
}

fn read<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T> {
    Ok(serde_json::from_reader(BufReader::with_capacity(
        1024 * 1024,
        File::open(path)?,
    ))?)
}
fn synced_new(path: &Path, value: &impl Serialize) -> Result<()> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, file);
    serde_json::to_writer_pretty(&mut writer, value)?;
    writer.flush()?;
    writer.get_ref().sync_all()?;
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
    let temp = temporary(folder, name);
    synced_new(&temp, value)?;
    replace_status(&temp, &folder.join(name))
}
fn copied(source: &Path, target: &Path) -> Result<String> {
    let hash = copy_checked(source, target)?;
    OpenOptions::new().write(true).open(target)?.sync_all()?;
    Ok(hash)
}
fn manifest(folder: &Path) -> Result<Manifest> {
    let m: Manifest = read_small(&folder.join("job.json"), 256 * 1024)?;
    if m.schema != "soil-fem-epoch-srm-search-job/1"
        || m.job_id != checkpoint::hash_json(&m.binding)?
        || !(1..=10000).contains(&m.binding.material_attempts_per_save)
    {
        return Err(Error::Input(
            "SRM job manifest identity/interval mismatch".into(),
        ));
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
        return Err(Error::Input("SRM immutable file inventory mismatch".into()));
    }
    for (name, hash) in &m.binding.files_sha256 {
        if checkpoint::hash_file(&folder.join(name))? != *hash {
            return Err(Error::Input(format!("SRM immutable input differs: {name}")));
        }
    }
    Ok(())
}
fn load(folder: &Path) -> Result<Context> {
    let manifest = manifest(folder)?;
    verify_files(folder, &manifest)?;
    let mesh: Mesh = read(&folder.join(MESH))?;
    let origin: Run = read(&folder.join(ORIGIN))?;
    let spec: Specification = read_small(&folder.join(SPEC), 256 * 1024 * 1024)?;
    if manifest.binding.identity != origin.provenance
        || manifest.binding.identity != Identity::current(&mesh, &origin.ground, &origin.case)?
        || manifest.binding.identity.solver_binary_sha256 != manifest.binding.files_sha256[SOLVER]
        || manifest.binding.specification_sha256 != checkpoint::hash_json(&spec)?
    {
        return Err(Error::Input(
            "SRM canonical inputs/specification/frozen executable mismatch".into(),
        ));
    }
    Ok(Context {
        manifest,
        mesh,
        origin,
        spec,
    })
}
fn save(folder: &Path, job_id: &str, generation: u64, cp: &SearchSnapshot) -> Result<()> {
    if generation == 0 {
        return Err(Error::Control(
            "SRM recovery generation must be positive".into(),
        ));
    }
    let name = format!("recovery-checkpoint-{}.json", generation % 2);
    let temp = temporary(folder, "recovery-data");
    synced_new(
        &temp,
        &SnapshotRef {
            schema: "soil-fem-epoch-srm-search-snapshot/1",
            job_id,
            generation,
            checkpoint: cp,
        },
    )?;
    let hash = checkpoint::hash_file(&temp)?;
    replace_status(&temp, &folder.join(&name))?;
    atomic_json(
        folder,
        "recovery.json",
        &Pointer {
            schema: "soil-fem-epoch-srm-search-recovery/1".into(),
            job_id: job_id.into(),
            generation,
            checkpoint: name,
            checkpoint_sha256: hash,
            payload_sha256: checkpoint::hash_json(cp)?,
        },
    )
}
fn selected(folder: &Path, context: &Context) -> Result<(Pointer, SearchSnapshot)> {
    for _ in 0..3 {
        let pointer: Pointer = read_small(&folder.join("recovery.json"), 8192)?;
        if pointer.schema != "soil-fem-epoch-srm-search-recovery/1"
            || pointer.job_id != context.manifest.job_id
            || pointer.generation == 0
            || pointer.checkpoint != format!("recovery-checkpoint-{}.json", pointer.generation % 2)
        {
            return Err(Error::Input(
                "SRM recovery pointer identity/path mismatch".into(),
            ));
        }
        let result = (|| -> Result<SearchSnapshot> {
            let path = folder.join(&pointer.checkpoint);
            if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
                return Err(Error::Input("SRM snapshot file checksum mismatch".into()));
            }
            let saved: StoredSnapshot = read(&path)?;
            if saved.schema != "soil-fem-epoch-srm-search-snapshot/1"
                || saved.job_id != context.manifest.job_id
                || saved.generation != pointer.generation
                || checkpoint::hash_json(&saved.checkpoint)? != pointer.payload_sha256
            {
                return Err(Error::Input("SRM snapshot binding mismatch".into()));
            }
            saved
                .checkpoint
                .verify(&context.mesh, &context.origin, &context.spec)?;
            if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
                return Err(Error::Input("SRM snapshot changed while read".into()));
            }
            Ok(saved.checkpoint)
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
        "SRM snapshot advanced during all3 reads; retain source and retry".into(),
    ))
}
fn trial_count(cp: &SearchSnapshot) -> usize {
    cp.report()
        .map_or_else(|| cp.trials().len(), |r| r.trials.len())
}
fn next_factor(cp: &SearchSnapshot, spec: &Specification) -> Option<f64> {
    cp.current_factor(spec)
}
fn has_saved_epoch(cp: &SearchSnapshot) -> bool {
    cp.has_inflight_checkpoint() || cp.trials().iter().any(|t| t.epoch.is_some())
}

fn prepare_impl(
    inputs: [&Path; 3],
    interval: usize,
    output: &Path,
    parent: Option<(Parent, SearchSnapshot)>,
) -> Result<Manifest> {
    if output.exists() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "SRM job output exists",
        )
        .into());
    }
    if !(1..=10000).contains(&interval) {
        return Err(Error::Input(
            "epoch SRM material-attempt interval must be1..10000".into(),
        ));
    }
    let mesh: Mesh = read(inputs[0])?;
    let origin: Run = read(inputs[1])?;
    let spec: Specification = read_small(inputs[2], 256 * 1024 * 1024)?;
    let (parent, initial) = if let Some((parent, cp)) = parent {
        cp.verify(&mesh, &origin, &spec)?;
        if parent.generation == 0 || parent.checkpoint_payload_sha256 != checkpoint::hash_json(&cp)?
        {
            return Err(Error::Input(
                "SRM parent generation/snapshot binding differs".into(),
            ));
        }
        (Some(parent), cp)
    } else {
        (None, SearchSnapshot::new(&mesh, &origin, &spec)?)
    };
    fs::create_dir(output)?;
    fs::create_dir(output.join("inputs"))?;
    fs::create_dir(output.join("solver"))?;
    let mut files = BTreeMap::new();
    for (source, name) in inputs.into_iter().zip([MESH, ORIGIN, SPEC]) {
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
        material_attempts_per_save: interval,
        nonce: format!(
            "{}-{}-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|e| Error::Control(e.to_string()))?
                .as_nanos(),
            std::process::id(),
            TEMPORARIES.fetch_add(1, Ordering::Relaxed)
        ),
        parent,
    };
    let m = Manifest {
        schema: "soil-fem-epoch-srm-search-job/1".into(),
        job_id: checkpoint::hash_json(&binding)?,
        binding,
    };
    synced_new(&output.join("job.json"), &m)?;
    load(output)?;
    save(output, &m.job_id, 1, &initial)?;
    synced_new(
        &output.join("status.json"),
        &Status {
            schema: "soil-fem-epoch-srm-search-job-status/1".into(),
            job_id: m.job_id.clone(),
            state: State::Prepared,
            process_id: None,
            current_factor: next_factor(&initial, &spec),
            completed_trials: trial_count(&initial),
            saved_material_epoch: has_saved_epoch(&initial),
            recovery_generation: 1,
            progress_events: 0,
            last_progress: None,
            requirement: initial.report().map(|r| r.requirement),
            search_stop_reason: initial.report().map(|r| r.stop_reason),
            report: None,
            report_sha256: None,
            reason: None,
        },
    )?;
    Ok(m)
}
pub fn prepare(
    mesh: &Path,
    origin: &Path,
    spec: &Path,
    interval: usize,
    output: &Path,
) -> Result<Manifest> {
    prepare_impl([mesh, origin, spec], interval, output, None)
}
pub fn prepare_recovery(previous: &Path, output: &Path) -> Result<Manifest> {
    let folder = fs::canonicalize(previous)?;
    let context = load(&folder)?;
    let (pointer, cp) = selected(&folder, &context)?;
    let parent = Parent {
        job_id: context.manifest.job_id.clone(),
        generation: pointer.generation,
        snapshot_sha256: pointer.checkpoint_sha256,
        checkpoint_payload_sha256: pointer.payload_sha256,
    };
    prepare_impl(
        [&folder.join(MESH), &folder.join(ORIGIN), &folder.join(SPEC)],
        context.manifest.binding.material_attempts_per_save,
        output,
        Some((parent, cp)),
    )
}
pub fn status(folder: &Path) -> Result<Status> {
    let m = manifest(folder)?;
    let s: Status = read_small(&folder.join("status.json"), 256 * 1024)?;
    if s.schema != "soil-fem-epoch-srm-search-job-status/1" || s.job_id != m.job_id {
        return Err(Error::Input("SRM job status identity mismatch".into()));
    }
    Ok(s)
}
fn requested(folder: &Path, job_id: &str) -> Result<bool> {
    let path = folder.join("stop-request.json");
    if !path.exists() {
        return Ok(false);
    }
    let r: StopRequest = read_small(&path, 4096)?;
    if r.schema != "soil-fem-epoch-srm-search-job-stop/1" || r.job_id != job_id {
        return Err(Error::Control(
            "SRM stop request belongs to another job".into(),
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
        let temp = temporary(folder, "stop-request");
        synced_new(
            &temp,
            &StopRequest {
                schema: "soil-fem-epoch-srm-search-job-stop/1".into(),
                job_id: s.job_id.clone(),
            },
        )?;
        // Publish a complete request atomically without replacing another
        // request. Readers must never observe the newly created empty file.
        match fs::hard_link(&temp, &path) {
            Ok(()) => fs::remove_file(&temp)?,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                fs::remove_file(&temp)?;
                requested(folder, &s.job_id)?;
            }
            Err(e) => return Err(e.into()),
        }
    }
    Ok(true)
}
pub fn verify_snapshot(folder: &Path) -> Result<SnapshotInfo> {
    let c = load(folder)?;
    let (pointer, cp) = selected(folder, &c)?;
    Ok(SnapshotInfo {
        schema: "soil-fem-epoch-srm-search-snapshot-verification/1".into(),
        job_id: c.manifest.job_id,
        generation: pointer.generation,
        checkpoint_sha256: pointer.checkpoint_sha256,
        canonical_snapshot_sha256: pointer.payload_sha256,
        current_factor: next_factor(&cp, &c.spec),
        completed_trials: trial_count(&cp),
        saved_material_epoch: has_saved_epoch(&cp),
        complete: cp.is_complete(),
        requirement: cp.report().map(|r| r.requirement),
        stop_reason: cp.report().map(|r| r.stop_reason),
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
    fn update(&self, cp: &SearchSnapshot, spec: &Specification, generation: u64) -> Result<()> {
        let mut t = self
            .trace
            .lock()
            .map_err(|_| Error::Control("SRM progress lock poisoned".into()))?;
        t.status.current_factor = next_factor(cp, spec);
        t.status.completed_trials = trial_count(cp);
        t.status.saved_material_epoch = has_saved_epoch(cp);
        t.status.recovery_generation = generation;
        t.status.requirement = cp.report().map(|r| r.requirement);
        t.status.search_stop_reason = cp.report().map(|r| r.stop_reason);
        t.events.flush()?;
        atomic_json(&self.folder, "status.json", &t.status)
    }
    fn finish(&self, state: State, reason: Option<String>, report: Option<String>) -> Result<()> {
        let mut t = self
            .trace
            .lock()
            .map_err(|_| Error::Control("SRM progress lock poisoned".into()))?;
        t.status.state = state;
        t.status.reason = reason;
        t.status.report_sha256 = report
            .as_ref()
            .map(|n| checkpoint::hash_file(&self.folder.join(n)))
            .transpose()?;
        t.status.report = report;
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
            .map_err(|_| Error::Control("SRM progress lock poisoned".into()))?;
        let factor = t.status.current_factor;
        serde_json::to_writer(
            &mut t.events,
            &serde_json::json!({"strength_factor":factor,"mechanics":event}),
        )?;
        t.events.write_all(b"\n")?;
        t.status.progress_events = t
            .status
            .progress_events
            .checked_add(1)
            .ok_or_else(|| Error::Control("SRM event counter overflow".into()))?;
        t.status.last_progress = Some(event.clone());
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

pub fn run(folder: &Path, advance_budget: Option<usize>) -> Result<Outcome> {
    if advance_budget.is_some_and(|n| !(1..=10000).contains(&n)) {
        return Err(Error::Input("SRM advance budget must be1..10000".into()));
    }
    let folder = fs::canonicalize(folder)?;
    let c = load(&folder)?;
    let mut state = status(&folder)?;
    if state.state != State::Prepared
        || ["run.claim", "report.json", "progress.jsonl"]
            .iter()
            .any(|n| folder.join(n).exists())
    {
        return Err(Error::Input(
            "SRM job is one-shot; fork its verified recovery snapshot into a new directory".into(),
        ));
    }
    let (pointer, mut cp) = selected(&folder, &c)?;
    let initial = if c.manifest.binding.parent.is_some() {
        read::<SearchSnapshot>(&folder.join(SAVED))?
    } else {
        SearchSnapshot::new(&c.mesh, &c.origin, &c.spec)?
    };
    if pointer.generation != 1 || checkpoint::hash_json(&cp)? != checkpoint::hash_json(&initial)? {
        return Err(Error::Input(
            "prepared SRM job must start from its exact initial snapshot".into(),
        ));
    }
    synced_new(
        &folder.join("run.claim"),
        &serde_json::json!({"job_id":c.manifest.job_id,"process_id":std::process::id()}),
    )?;
    state.state = State::Running;
    state.process_id = Some(std::process::id());
    atomic_json(&folder, "status.json", &state)?;
    let observer = FileObserver {
        folder: folder.clone(),
        job_id: c.manifest.job_id.clone(),
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
        let mut generation = pointer.generation;
        let mut advances = 0;
        loop {
            if cp.is_complete() {
                verify_files(&folder, &c.manifest)?;
                synced_new(&folder.join("report.json"), cp.report().unwrap())?;
                observer.finish(State::Completed, None, Some("report.json".into()))?;
                return Ok(Outcome::Completed);
            }
            let reason = if observer.stop_requested()? {
                Some("stop_requested")
            } else if advance_budget.is_some_and(|n| advances >= n) {
                Some("advance_budget")
            } else {
                None
            };
            if let Some(reason) = reason {
                verify_files(&folder, &c.manifest)?;
                observer.finish(State::Paused, Some(reason.into()), None)?;
                return Ok(Outcome::Paused);
            }
            cp = cp.advance(
                &c.mesh,
                &c.origin,
                &c.spec,
                c.manifest.binding.material_attempts_per_save,
                Some(&observer),
            )?;
            verify_files(&folder, &c.manifest)?;
            generation = generation
                .checked_add(1)
                .ok_or_else(|| Error::Control("SRM generation overflow".into()))?;
            save(&folder, &c.manifest.job_id, generation, &cp)?;
            observer.update(&cp, &c.spec, generation)?;
            advances += 1;
        }
    })();
    if let Err(e) = &result {
        // Keep the preceding complete recovery pointer even on a write/control
        // failure. An uncommitted Newton trial is never made resumable.
        let _ = observer.finish(State::Failed, Some(e.to_string()), None);
    }
    result
}
