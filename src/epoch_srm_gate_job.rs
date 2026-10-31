//! One-shot jobs for the unreduced epoch Fs gate. Completion of the job is not a
//! safety certificate: read `decision`, each gate effect, and the search evidence.
//! Recovery publishes only a fully written snapshot. Later stages are absent
//! until the preceding Fs_req gate continues.
use crate::{
    Error, Result,
    checkpoint::{self, Identity},
    epoch_srm_gate::{Decision, GateEffect, Snapshot as GateSnapshot, Specification},
    ground::Ground,
    mesh::Mesh,
    srm_search::{Requirement, StopReason},
    stage_control::{Observer, Progress},
    stage_job::{copy_checked, read_small, replace_status},
    stages::StageCase,
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
const GROUND: &str = "inputs/ground.json";
const CASE: &str = "inputs/case.json";
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
    pub accepted_units_per_save: usize,
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
    pub phase: String,
    pub through_stages: Option<usize>,
    pub current_factor: Option<f64>,
    pub completed_gates: usize,
    pub recovery_generation: u64,
    pub progress_events: u64,
    pub last_progress: Option<Progress>,
    pub decision: Option<Decision>,
    pub gate_effect: Option<GateEffect>,
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
    checkpoint: GateSnapshot,
}
#[derive(Serialize)]
struct SnapshotRef<'a> {
    schema: &'static str,
    job_id: &'a str,
    generation: u64,
    checkpoint: &'a GateSnapshot,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SnapshotInfo {
    pub schema: String,
    pub job_id: String,
    pub generation: u64,
    pub checkpoint_sha256: String,
    pub canonical_snapshot_sha256: String,
    pub phase: String,
    pub through_stages: Option<usize>,
    pub current_factor: Option<f64>,
    pub completed_gates: usize,
    pub complete: bool,
    pub decision: Option<Decision>,
    pub gate_effect: Option<GateEffect>,
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
    ground: Ground,
    case: StageCase,
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
    let value: Manifest = read_small(&folder.join("job.json"), 256 * 1024)?;
    if value.schema != "soil-fem-epoch-srm-stage-gate-job/1"
        || value.job_id != checkpoint::hash_json(&value.binding)?
        || !(1..=10000).contains(&value.binding.accepted_units_per_save)
    {
        return Err(Error::Input(
            "epoch SRM gate job manifest identity/interval mismatch".into(),
        ));
    }
    Ok(value)
}
fn verify_files(folder: &Path, value: &Manifest) -> Result<()> {
    let mut expected = BTreeSet::from([MESH, GROUND, CASE, SPEC, SOLVER]);
    if value.binding.parent.is_some() {
        expected.insert(SAVED);
    }
    if value
        .binding
        .files_sha256
        .keys()
        .map(String::as_str)
        .collect::<BTreeSet<_>>()
        != expected
    {
        return Err(Error::Input(
            "epoch SRM gate immutable file inventory mismatch".into(),
        ));
    }
    for (name, hash) in &value.binding.files_sha256 {
        if checkpoint::hash_file(&folder.join(name))? != *hash {
            return Err(Error::Input(format!(
                "epoch SRM gate immutable input differs: {name}"
            )));
        }
    }
    Ok(())
}
fn load(folder: &Path) -> Result<Context> {
    let manifest = manifest(folder)?;
    verify_files(folder, &manifest)?;
    let mesh: Mesh = read(&folder.join(MESH))?;
    let ground: Ground = read(&folder.join(GROUND))?;
    let case: StageCase = read(&folder.join(CASE))?;
    let spec: Specification = read_small(&folder.join(SPEC), 256 * 1024 * 1024)?;
    if manifest.binding.identity != Identity::current(&mesh, &ground, &case)?
        || manifest.binding.identity.solver_binary_sha256 != manifest.binding.files_sha256[SOLVER]
        || manifest.binding.specification_sha256 != checkpoint::hash_json(&spec)?
    {
        return Err(Error::Input(
            "epoch SRM gate canonical inputs/specification/frozen executable mismatch".into(),
        ));
    }
    Ok(Context {
        manifest,
        mesh,
        ground,
        case,
        spec,
    })
}
fn save(folder: &Path, job_id: &str, generation: u64, checkpoint: &GateSnapshot) -> Result<()> {
    if generation == 0 {
        return Err(Error::Control(
            "epoch SRM gate recovery generation must be positive".into(),
        ));
    }
    let name = format!("recovery-checkpoint-{}.json", generation % 2);
    let temp = temporary(folder, "recovery-data");
    synced_new(
        &temp,
        &SnapshotRef {
            schema: "soil-fem-epoch-srm-stage-gate-snapshot/1",
            job_id,
            generation,
            checkpoint,
        },
    )?;
    let hash = checkpoint::hash_file(&temp)?;
    replace_status(&temp, &folder.join(&name))?;
    atomic_json(
        folder,
        "recovery.json",
        &Pointer {
            schema: "soil-fem-epoch-srm-stage-gate-recovery/1".into(),
            job_id: job_id.into(),
            generation,
            checkpoint: name,
            checkpoint_sha256: hash,
            payload_sha256: checkpoint::hash_json(checkpoint)?,
        },
    )
}
fn selected(folder: &Path, context: &Context) -> Result<(Pointer, GateSnapshot)> {
    for _ in 0..3 {
        let pointer: Pointer = read_small(&folder.join("recovery.json"), 8192)?;
        if pointer.schema != "soil-fem-epoch-srm-stage-gate-recovery/1"
            || pointer.job_id != context.manifest.job_id
            || pointer.generation == 0
            || pointer.checkpoint != format!("recovery-checkpoint-{}.json", pointer.generation % 2)
        {
            return Err(Error::Input(
                "epoch SRM gate recovery pointer identity/path mismatch".into(),
            ));
        }
        let result = (|| -> Result<GateSnapshot> {
            let path = folder.join(&pointer.checkpoint);
            if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
                return Err(Error::Input(
                    "epoch SRM gate snapshot file checksum mismatch".into(),
                ));
            }
            let saved: StoredSnapshot = read(&path)?;
            if saved.schema != "soil-fem-epoch-srm-stage-gate-snapshot/1"
                || saved.job_id != context.manifest.job_id
                || saved.generation != pointer.generation
                || checkpoint::hash_json(&saved.checkpoint)? != pointer.payload_sha256
            {
                return Err(Error::Input(
                    "epoch SRM gate snapshot binding mismatch".into(),
                ));
            }
            saved.checkpoint.verify(
                &context.mesh,
                &context.ground,
                &context.case,
                &context.spec,
            )?;
            if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
                return Err(Error::Input(
                    "epoch SRM gate snapshot changed while read".into(),
                ));
            }
            Ok(saved.checkpoint)
        })();
        match result {
            Ok(checkpoint) => return Ok((pointer, checkpoint)),
            Err(error) => {
                let now: Pointer = read_small(&folder.join("recovery.json"), 8192)?;
                if now == pointer {
                    return Err(error);
                }
            }
        }
    }
    Err(Error::Control(
        "epoch SRM gate snapshot advanced during all 3 reads; retain source and retry".into(),
    ))
}
fn effect_of(checkpoint: &GateSnapshot) -> Option<GateEffect> {
    checkpoint
        .report()
        .and_then(|report| report.gates.last().map(|gate| gate.effect))
}

fn prepare_impl(
    inputs: [&Path; 4],
    interval: usize,
    output: &Path,
    parent: Option<(Parent, GateSnapshot)>,
) -> Result<Manifest> {
    if output.exists() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "epoch SRM gate job output exists",
        )
        .into());
    }
    if !(1..=10000).contains(&interval) {
        return Err(Error::Input(
            "epoch SRM gate accepted-unit interval must be 1..10000".into(),
        ));
    }
    let mesh: Mesh = read(inputs[0])?;
    let ground: Ground = read(inputs[1])?;
    let case: StageCase = read(inputs[2])?;
    let spec: Specification = read_small(inputs[3], 256 * 1024 * 1024)?;
    let (parent, initial) = if let Some((parent, checkpoint)) = parent {
        checkpoint.verify(&mesh, &ground, &case, &spec)?;
        if parent.generation == 0
            || parent.checkpoint_payload_sha256 != checkpoint::hash_json(&checkpoint)?
        {
            return Err(Error::Input(
                "epoch SRM gate parent generation/snapshot binding differs".into(),
            ));
        }
        (Some(parent), checkpoint)
    } else {
        (None, GateSnapshot::new(&mesh, &ground, &case, &spec)?)
    };
    fs::create_dir(output)?;
    fs::create_dir(output.join("inputs"))?;
    fs::create_dir(output.join("solver"))?;
    let mut files = BTreeMap::new();
    for (source, name) in inputs.into_iter().zip([MESH, GROUND, CASE, SPEC]) {
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
        identity: Identity::current(&mesh, &ground, &case)?,
        specification_sha256: checkpoint::hash_json(&spec)?,
        files_sha256: files,
        accepted_units_per_save: interval,
        nonce: format!(
            "{}-{}-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|error| Error::Control(error.to_string()))?
                .as_nanos(),
            std::process::id(),
            TEMPORARIES.fetch_add(1, Ordering::Relaxed)
        ),
        parent,
    };
    let value = Manifest {
        schema: "soil-fem-epoch-srm-stage-gate-job/1".into(),
        job_id: checkpoint::hash_json(&binding)?,
        binding,
    };
    synced_new(&output.join("job.json"), &value)?;
    let context = load(output)?;
    save(output, &value.job_id, 1, &initial)?;
    synced_new(
        &output.join("status.json"),
        &Status {
            schema: "soil-fem-epoch-srm-stage-gate-job-status/1".into(),
            job_id: value.job_id.clone(),
            state: State::Prepared,
            process_id: None,
            phase: initial.phase_label().into(),
            through_stages: initial.through_stages(),
            current_factor: initial.current_factor(&context.spec),
            completed_gates: initial.completed_gates(),
            recovery_generation: 1,
            progress_events: 0,
            last_progress: None,
            decision: initial.decision(),
            gate_effect: effect_of(&initial),
            requirement: initial.latest_requirement(),
            search_stop_reason: initial.latest_stop_reason(),
            report: None,
            report_sha256: None,
            reason: None,
        },
    )?;
    Ok(value)
}
pub fn prepare(
    mesh: &Path,
    ground: &Path,
    case: &Path,
    spec: &Path,
    interval: usize,
    output: &Path,
) -> Result<Manifest> {
    prepare_impl([mesh, ground, case, spec], interval, output, None)
}
pub fn prepare_recovery(previous: &Path, output: &Path) -> Result<Manifest> {
    let folder = fs::canonicalize(previous)?;
    let context = load(&folder)?;
    let (pointer, checkpoint) = selected(&folder, &context)?;
    let parent = Parent {
        job_id: context.manifest.job_id.clone(),
        generation: pointer.generation,
        snapshot_sha256: pointer.checkpoint_sha256,
        checkpoint_payload_sha256: pointer.payload_sha256,
    };
    prepare_impl(
        [
            &folder.join(MESH),
            &folder.join(GROUND),
            &folder.join(CASE),
            &folder.join(SPEC),
        ],
        context.manifest.binding.accepted_units_per_save,
        output,
        Some((parent, checkpoint)),
    )
}
pub fn status(folder: &Path) -> Result<Status> {
    let value = manifest(folder)?;
    let state: Status = read_small(&folder.join("status.json"), 256 * 1024)?;
    if state.schema != "soil-fem-epoch-srm-stage-gate-job-status/1" || state.job_id != value.job_id
    {
        return Err(Error::Input(
            "epoch SRM gate job status identity mismatch".into(),
        ));
    }
    Ok(state)
}
fn requested(folder: &Path, job_id: &str) -> Result<bool> {
    let path = folder.join("stop-request.json");
    if !path.exists() {
        return Ok(false);
    }
    let request: StopRequest = read_small(&path, 4096)?;
    if request.schema != "soil-fem-epoch-srm-stage-gate-job-stop/1" || request.job_id != job_id {
        return Err(Error::Control(
            "epoch SRM gate stop request belongs to another job".into(),
        ));
    }
    Ok(true)
}
pub fn request_stop(folder: &Path) -> Result<bool> {
    let state = status(folder)?;
    if !matches!(state.state, State::Prepared | State::Running) {
        return Ok(false);
    }
    let path = folder.join("stop-request.json");
    if path.exists() {
        requested(folder, &state.job_id)?;
    } else {
        let temp = temporary(folder, "stop-request");
        synced_new(
            &temp,
            &StopRequest {
                schema: "soil-fem-epoch-srm-stage-gate-job-stop/1".into(),
                job_id: state.job_id.clone(),
            },
        )?;
        match fs::hard_link(&temp, &path) {
            Ok(()) => fs::remove_file(&temp)?,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
                fs::remove_file(&temp)?;
                requested(folder, &state.job_id)?;
            }
            Err(error) => return Err(error.into()),
        }
    }
    Ok(true)
}
pub fn verify_snapshot(folder: &Path) -> Result<SnapshotInfo> {
    let context = load(folder)?;
    let (pointer, checkpoint) = selected(folder, &context)?;
    Ok(SnapshotInfo {
        schema: "soil-fem-epoch-srm-stage-gate-snapshot-verification/1".into(),
        job_id: context.manifest.job_id,
        generation: pointer.generation,
        checkpoint_sha256: pointer.checkpoint_sha256,
        canonical_snapshot_sha256: pointer.payload_sha256,
        phase: checkpoint.phase_label().into(),
        through_stages: checkpoint.through_stages(),
        current_factor: checkpoint.current_factor(&context.spec),
        completed_gates: checkpoint.completed_gates(),
        complete: checkpoint.is_complete(),
        decision: checkpoint.decision(),
        gate_effect: effect_of(&checkpoint),
        requirement: checkpoint.latest_requirement(),
        stop_reason: checkpoint.latest_stop_reason(),
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
    fn update(
        &self,
        checkpoint: &GateSnapshot,
        spec: &Specification,
        generation: u64,
    ) -> Result<()> {
        let mut trace = self
            .trace
            .lock()
            .map_err(|_| Error::Control("epoch SRM gate progress lock poisoned".into()))?;
        trace.status.phase = checkpoint.phase_label().into();
        trace.status.through_stages = checkpoint.through_stages();
        trace.status.current_factor = checkpoint.current_factor(spec);
        trace.status.completed_gates = checkpoint.completed_gates();
        trace.status.recovery_generation = generation;
        trace.status.decision = checkpoint.decision();
        trace.status.gate_effect = effect_of(checkpoint);
        trace.status.requirement = checkpoint.latest_requirement();
        trace.status.search_stop_reason = checkpoint.latest_stop_reason();
        trace.events.flush()?;
        atomic_json(&self.folder, "status.json", &trace.status)
    }
    fn finish(&self, state: State, reason: Option<String>, report: Option<String>) -> Result<()> {
        let mut trace = self
            .trace
            .lock()
            .map_err(|_| Error::Control("epoch SRM gate progress lock poisoned".into()))?;
        trace.status.state = state;
        trace.status.reason = reason;
        trace.status.report_sha256 = report
            .as_ref()
            .map(|name| checkpoint::hash_file(&self.folder.join(name)))
            .transpose()?;
        trace.status.report = report;
        trace.events.flush()?;
        trace.events.get_ref().sync_all()?;
        atomic_json(&self.folder, "status.json", &trace.status)
    }
}
impl Observer for FileObserver {
    fn on_progress(&self, event: &Progress) -> Result<()> {
        let mut trace = self
            .trace
            .lock()
            .map_err(|_| Error::Control("epoch SRM gate progress lock poisoned".into()))?;
        let event_record = serde_json::json!({
            "phase": trace.status.phase,
            "through_stages": trace.status.through_stages,
            "strength_factor": trace.status.current_factor,
            "mechanics": event
        });
        serde_json::to_writer(&mut trace.events, &event_record)?;
        trace.events.write_all(b"\n")?;
        trace.status.progress_events = trace
            .status
            .progress_events
            .checked_add(1)
            .ok_or_else(|| Error::Control("epoch SRM gate event counter overflow".into()))?;
        trace.status.last_progress = Some(event.clone());
        if trace.last_status.elapsed().as_millis() >= 100 {
            trace.events.flush()?;
            atomic_json(&self.folder, "status.json", &trace.status)?;
            trace.last_status = Instant::now();
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        requested(&self.folder, &self.job_id)
    }
}

pub fn run(folder: &Path, advance_budget: Option<usize>) -> Result<Outcome> {
    if advance_budget.is_some_and(|count| !(1..=10000).contains(&count)) {
        return Err(Error::Input(
            "epoch SRM gate advance budget must be 1..10000".into(),
        ));
    }
    let folder = fs::canonicalize(folder)?;
    let context = load(&folder)?;
    let mut state = status(&folder)?;
    if state.state != State::Prepared
        || ["run.claim", "report.json", "progress.jsonl"]
            .iter()
            .any(|name| folder.join(name).exists())
    {
        return Err(Error::Input(
            "epoch SRM gate job is one-shot; fork its verified recovery snapshot into a new directory"
                .into(),
        ));
    }
    let (pointer, mut checkpoint) = selected(&folder, &context)?;
    let initial = if context.manifest.binding.parent.is_some() {
        read::<GateSnapshot>(&folder.join(SAVED))?
    } else {
        GateSnapshot::new(&context.mesh, &context.ground, &context.case, &context.spec)?
    };
    if pointer.generation != 1
        || checkpoint::hash_json(&checkpoint)? != checkpoint::hash_json(&initial)?
    {
        return Err(Error::Input(
            "prepared epoch SRM gate job must start from its exact initial snapshot".into(),
        ));
    }
    synced_new(
        &folder.join("run.claim"),
        &serde_json::json!({"job_id": context.manifest.job_id, "process_id": std::process::id()}),
    )?;
    state.state = State::Running;
    state.process_id = Some(std::process::id());
    atomic_json(&folder, "status.json", &state)?;
    let observer = FileObserver {
        folder: folder.clone(),
        job_id: context.manifest.job_id.clone(),
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
            if checkpoint.is_complete() {
                verify_files(&folder, &context.manifest)?;
                synced_new(
                    &folder.join("report.json"),
                    checkpoint.report().ok_or_else(|| {
                        Error::Control("completed epoch SRM gate omitted its report".into())
                    })?,
                )?;
                observer.finish(State::Completed, None, Some("report.json".into()))?;
                return Ok(Outcome::Completed);
            }
            let reason = if observer.stop_requested()? {
                Some("stop_requested")
            } else if advance_budget.is_some_and(|count| advances >= count) {
                Some("advance_budget")
            } else {
                None
            };
            if let Some(reason) = reason {
                verify_files(&folder, &context.manifest)?;
                observer.finish(State::Paused, Some(reason.into()), None)?;
                return Ok(Outcome::Paused);
            }
            checkpoint = checkpoint.advance(
                &context.mesh,
                &context.ground,
                &context.case,
                &context.spec,
                context.manifest.binding.accepted_units_per_save,
                Some(&observer),
            )?;
            verify_files(&folder, &context.manifest)?;
            generation = generation
                .checked_add(1)
                .ok_or_else(|| Error::Control("epoch SRM gate generation overflow".into()))?;
            save(&folder, &context.manifest.job_id, generation, &checkpoint)?;
            observer.update(&checkpoint, &context.spec, generation)?;
            advances += 1;
        }
    })();
    if let Err(error) = &result {
        let _ = observer.finish(State::Failed, Some(error.to_string()), None);
    }
    result
}
