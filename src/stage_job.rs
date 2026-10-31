//! File-based staged-mechanics jobs with immutable inputs and cooperative pause.
//! A status file is the last report, never proof of process liveness. No PID is
//! used to stop another process; only a matching request is observed by its job.
use crate::{
    Error, Result,
    checkpoint::{self, Checkpoint, Identity},
    ground::Ground,
    mesh::Mesh,
    stage_control::{Observer, Progress},
    stages::{self, StageCase},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    sync::{
        Mutex,
        atomic::{AtomicU64, Ordering},
    },
    time::{Instant, SystemTime, UNIX_EPOCH},
};

mod recovery;

const MESH: &str = "inputs/mesh.json";
const GROUND: &str = "inputs/ground.json";
const CASE: &str = "inputs/case.json";
const SAVED: &str = "inputs/checkpoint.json";
const SOLVER: &str = "solver/soil-fem-3d.exe";
static TOKENS: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Binding {
    pub identity: Identity,
    pub files_sha256: BTreeMap<String, String>,
    pub nonce: String,
    pub parent_job_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checkpoint_interval: Option<usize>,
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
    Checkpointed,
    Succeeded,
    Failed,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Status {
    pub schema: String,
    pub job_id: String,
    pub state: State,
    pub process_id: Option<u32>,
    pub progress_events: usize,
    pub last_progress: Option<Progress>,
    pub stop_requested: bool,
    pub artifact: Option<String>,
    pub artifact_sha256: Option<String>,
    pub reason: Option<String>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StopRequest {
    pub schema: String,
    pub job_id: String,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Succeeded,
    Checkpointed,
    Failed,
}

fn read<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T> {
    Ok(serde_json::from_reader(BufReader::new(File::open(path)?))?)
}
pub(crate) fn read_small<T: serde::de::DeserializeOwned>(path: &Path, limit: u64) -> Result<T> {
    let file = File::open(path)?;
    let mut data = Vec::new();
    file.take(limit + 1).read_to_end(&mut data)?;
    if data.len() as u64 > limit {
        return Err(Error::Input("job metadata exceeds its size limit".into()));
    }
    Ok(serde_json::from_slice(&data)?)
}
fn write_new(path: &Path, value: &impl Serialize) -> Result<()> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut stream = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut stream, value)?;
    stream.flush()?;
    Ok(())
}
fn atomic_status(folder: &Path, status: &Status, sequence: u64) -> Result<()> {
    let temporary = folder.join(format!(".status-{}-{sequence}.part", std::process::id()));
    write_new(&temporary, status)?;
    let target = folder.join("status.json");
    replace_status(&temporary, &target)
}
pub(crate) fn replace_status(temporary: &Path, target: &Path) -> Result<()> {
    // Windows observers/scanners may hold the destination without delete
    // sharing for longer than 40 ms. Back off for at most 630 ms, retrying
    // exactly this already complete owned temporary. Never reserialize a
    // progress event, select another temporary, or remove the old snapshot.
    let mut delays = [10, 20, 40, 80, 160, 320].into_iter();
    loop {
        match fs::rename(temporary, target) {
            Ok(()) => return Ok(()),
            Err(e) if cfg!(windows) && matches!(e.raw_os_error(), Some(5 | 32 | 33)) => {
                if let Some(delay) = delays.next() {
                    std::thread::sleep(std::time::Duration::from_millis(delay));
                } else {
                    return Err(e.into());
                }
            }
            Err(e) => return Err(e.into()),
        }
    }
}
pub(crate) fn copy_checked(source: &Path, target: &Path) -> Result<String> {
    let digest = checkpoint::hash_file(source)?;
    let mut incoming = File::open(source)?;
    let mut outgoing = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(target)?;
    std::io::copy(&mut incoming, &mut outgoing)?;
    outgoing.flush()?;
    drop(outgoing);
    // The byte hash is unchanged. Windows executes `*.exe` without a Unix mode;
    // other hosts must keep the source executable bits or the frozen solver cannot start.
    #[cfg(unix)]
    preserve_executable_bits(source, target)?;
    if checkpoint::hash_file(source)? != digest || checkpoint::hash_file(target)? != digest {
        return Err(Error::Input("job input changed while being copied".into()));
    }
    Ok(digest)
}
#[cfg(unix)]
fn preserve_executable_bits(source: &Path, target: &Path) -> Result<()> {
    use std::os::unix::fs::PermissionsExt;
    let executable = fs::metadata(source)?.permissions().mode() & 0o111;
    if executable == 0 {
        return Ok(());
    }
    let mut perms = fs::metadata(target)?.permissions();
    perms.set_mode((perms.mode() & 0o777) | executable);
    fs::set_permissions(target, perms)?;
    Ok(())
}
fn load(folder: &Path) -> Result<(Manifest, Mesh, Ground, StageCase)> {
    let manifest: Manifest = read_small(&folder.join("job.json"), 256 * 1024)?;
    if manifest.schema != "soil-fem-stages-job/1"
        || manifest.job_id != checkpoint::hash_json(&manifest.binding)?
    {
        return Err(Error::Input("job manifest identity mismatch".into()));
    }
    let mut expected = BTreeSet::from([MESH, GROUND, CASE, SOLVER]);
    if manifest.binding.parent_job_id.is_some() {
        expected.insert(SAVED);
    }
    if manifest
        .binding
        .files_sha256
        .keys()
        .map(String::as_str)
        .collect::<BTreeSet<_>>()
        != expected
    {
        return Err(Error::Input("job file inventory mismatch".into()));
    }
    for (name, digest) in &manifest.binding.files_sha256 {
        if checkpoint::hash_file(&folder.join(name))? != *digest {
            return Err(Error::Input(format!("immutable job input differs: {name}")));
        }
    }
    let mesh: Mesh = read(&folder.join(MESH))?;
    let ground: Ground = read(&folder.join(GROUND))?;
    let case: StageCase = read(&folder.join(CASE))?;
    validate_interval(manifest.binding.checkpoint_interval)?;
    let identity = &manifest.binding.identity;
    if identity.mesh_sha256 != checkpoint::hash_json(&mesh)?
        || identity.ground_sha256 != checkpoint::hash_json(&ground)?
        || identity.case_sha256 != checkpoint::hash_json(&case)?
        || identity.solver_binary_sha256 != manifest.binding.files_sha256[SOLVER]
    {
        return Err(Error::Input(
            "job canonical inputs differ from their identity".into(),
        ));
    }
    Ok((manifest, mesh, ground, case))
}
fn prepare_impl(
    mesh_file: &Path,
    ground_file: &Path,
    case_file: &Path,
    output: &Path,
    parent: Option<(&str, &Path)>,
    checkpoint_interval: Option<usize>,
) -> Result<Manifest> {
    validate_interval(checkpoint_interval)?;
    let mesh: Mesh = read(mesh_file)?;
    let ground: Ground = read(ground_file)?;
    let case: StageCase = read(case_file)?;
    let identity = Identity::current(&mesh, &ground, &case)?;
    if let Some((_, saved)) = parent {
        let checkpoint: Checkpoint = read(saved)?;
        checkpoint.verify(&mesh, &ground, &case)?;
    }
    fs::create_dir(output)?;
    fs::create_dir(output.join("inputs"))?;
    fs::create_dir(output.join("solver"))?;
    let mut files = BTreeMap::new();
    for (source, name) in [(mesh_file, MESH), (ground_file, GROUND), (case_file, CASE)] {
        files.insert(name.into(), copy_checked(source, &output.join(name))?);
    }
    files.insert(
        SOLVER.into(),
        copy_checked(&std::env::current_exe()?, &output.join(SOLVER))?,
    );
    if let Some((_, source)) = parent {
        files.insert(SAVED.into(), copy_checked(source, &output.join(SAVED))?);
    }
    let binding = Binding {
        identity,
        files_sha256: files,
        nonce: format!(
            "{}-{}-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|e| Error::Control(e.to_string()))?
                .as_nanos(),
            std::process::id(),
            TOKENS.fetch_add(1, Ordering::Relaxed)
        ),
        parent_job_id: parent.map(|(id, _)| id.into()),
        checkpoint_interval,
    };
    let manifest = Manifest {
        schema: "soil-fem-stages-job/1".into(),
        job_id: checkpoint::hash_json(&binding)?,
        binding,
    };
    write_new(&output.join("job.json"), &manifest)?;
    load(output)?;
    write_new(
        &output.join("status.json"),
        &Status {
            schema: "soil-fem-stages-job-status/1".into(),
            job_id: manifest.job_id.clone(),
            state: State::Prepared,
            process_id: None,
            progress_events: 0,
            last_progress: None,
            stop_requested: false,
            artifact: None,
            artifact_sha256: None,
            reason: None,
        },
    )?;
    Ok(manifest)
}
pub fn prepare(mesh: &Path, ground: &Path, case: &Path, output: &Path) -> Result<Manifest> {
    prepare_impl(mesh, ground, case, output, None, None)
}
pub fn prepare_resume(previous: &Path, output: &Path) -> Result<Manifest> {
    let (manifest, mesh, ground, case) = load(previous)?;
    if manifest.binding.identity != Identity::current(&mesh, &ground, &case)? {
        return Err(Error::Input(
            "resume the job with its frozen solver executable".into(),
        ));
    }
    let saved = previous.join("checkpoint.json");
    let checkpoint: Checkpoint = read(&saved)?;
    checkpoint.verify(&mesh, &ground, &case)?;
    prepare_impl(
        &previous.join(MESH),
        &previous.join(GROUND),
        &previous.join(CASE),
        output,
        Some((&manifest.job_id, &saved)),
        manifest.binding.checkpoint_interval,
    )
}
fn validate_interval(interval: Option<usize>) -> Result<()> {
    if interval.is_some_and(|n| !(1..=4096).contains(&n)) {
        return Err(Error::Input(
            "checkpoint interval must be1..4096 accepted increments".into(),
        ));
    }
    Ok(())
}
pub fn prepare_durable(
    mesh: &Path,
    ground: &Path,
    case: &Path,
    interval: usize,
    output: &Path,
) -> Result<Manifest> {
    prepare_impl(mesh, ground, case, output, None, Some(interval))
}
/// Fork a verified accepted snapshot into a new immutable job. This does not
/// stop its parent or infer liveness from a stored status/PID.
pub fn prepare_recovery(previous: &Path, output: &Path) -> Result<Manifest> {
    let (manifest, mesh, ground, case) = load(previous)?;
    if manifest.binding.identity != Identity::current(&mesh, &ground, &case)? {
        return Err(Error::Input(
            "recover with the frozen solver executable".into(),
        ));
    }
    let saved = recovery::saved_path(previous, &manifest, &mesh, &ground, &case)?;
    prepare_impl(
        &previous.join(MESH),
        &previous.join(GROUND),
        &previous.join(CASE),
        output,
        Some((&manifest.job_id, &saved)),
        manifest.binding.checkpoint_interval,
    )
}

pub fn status(folder: &Path) -> Result<Status> {
    let manifest: Manifest = read_small(&folder.join("job.json"), 256 * 1024)?;
    let status: Status = read_small(&folder.join("status.json"), 256 * 1024)?;
    if status.schema != "soil-fem-stages-job-status/1"
        || manifest.schema != "soil-fem-stages-job/1"
        || manifest.job_id != checkpoint::hash_json(&manifest.binding)?
        || status.job_id != manifest.job_id
    {
        return Err(Error::Input("job status identity mismatch".into()));
    }
    Ok(status)
}
pub fn request_stop(folder: &Path) -> Result<bool> {
    let state = status(folder)?;
    if matches!(
        state.state,
        State::Succeeded | State::Checkpointed | State::Failed
    ) {
        return Ok(false);
    }
    let request = StopRequest {
        schema: "soil-fem-stages-job-stop/1".into(),
        job_id: state.job_id,
    };
    let path = folder.join("stop-request.json");
    if path.exists() {
        let old: StopRequest = read_small(&path, 4096)?;
        if old.schema != request.schema || old.job_id != request.job_id {
            return Err(Error::Input(
                "existing stop request belongs to another job".into(),
            ));
        }
    } else {
        write_new(&path, &request)?;
    }
    Ok(true)
}

struct TraceState {
    status: Status,
    events: BufWriter<File>,
    sequence: u64,
    last_status: Instant,
    accepted_since_checkpoint: usize,
}
struct FileObserver {
    folder: PathBuf,
    job_id: String,
    trace: Mutex<TraceState>,
    checkpoint_interval: Option<usize>,
}
impl FileObserver {
    fn reset_checkpoint_count(&self) -> Result<()> {
        self.trace
            .lock()
            .map_err(|e| Error::Control(e.to_string()))?
            .accepted_since_checkpoint = 0;
        Ok(())
    }
    fn requested(&self) -> bool {
        read_small::<StopRequest>(&self.folder.join("stop-request.json"), 4096)
            .is_ok_and(|r| r.schema == "soil-fem-stages-job-stop/1" && r.job_id == self.job_id)
    }
    fn finish(&self, state: State, artifact: Option<&str>, reason: Option<String>) -> Result<()> {
        let mut trace = self
            .trace
            .lock()
            .map_err(|e| Error::Control(e.to_string()))?;
        trace.events.flush()?;
        trace.status.state = state;
        trace.status.stop_requested = self.requested();
        trace.status.reason = reason;
        trace.status.artifact = artifact.map(str::to_owned);
        trace.status.artifact_sha256 = artifact
            .map(|name| checkpoint::hash_file(&self.folder.join(name)))
            .transpose()?;
        trace.sequence += 1;
        atomic_status(&self.folder, &trace.status, trace.sequence)
    }
}
impl Observer for FileObserver {
    fn on_progress(&self, progress: &Progress) -> Result<()> {
        let mut trace = self
            .trace
            .lock()
            .map_err(|e| Error::Control(e.to_string()))?;
        let sequence = trace.status.progress_events + 1;
        serde_json::to_writer(
            &mut trace.events,
            &serde_json::json!({"schema":"soil-fem-stages-job-progress/1","job_id":self.job_id,"sequence":sequence,"progress":progress}),
        )?;
        trace.events.write_all(b"\n")?;
        trace.events.flush()?;
        if matches!(progress.kind, crate::stage_control::Kind::Accepted) {
            trace.accepted_since_checkpoint += 1;
        }
        trace.status.progress_events = sequence;
        trace.status.last_progress = Some(progress.clone());
        trace.status.stop_requested = self.requested();
        if matches!(
            progress.kind,
            crate::stage_control::Kind::Accepted | crate::stage_control::Kind::Rejected
        ) || trace.last_status.elapsed().as_millis() >= 200
        {
            trace.sequence += 1;
            atomic_status(&self.folder, &trace.status, trace.sequence)?;
            trace.last_status = Instant::now();
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        let accepted = self
            .trace
            .lock()
            .map_err(|e| Error::Control(e.to_string()))?
            .accepted_since_checkpoint;
        Ok(self.requested() || self.checkpoint_interval.is_some_and(|n| accepted >= n))
    }
}

pub fn run(folder: &Path) -> Result<Outcome> {
    let folder = fs::canonicalize(folder)?;
    let (manifest, mesh, ground, case) = load(&folder)?;
    if manifest.binding.identity != Identity::current(&mesh, &ground, &case)? {
        return Err(Error::Input(
            "run the job with its frozen solver executable".into(),
        ));
    }
    let mut state = status(&folder)?;
    if state.state != State::Prepared
        || [
            "result.json",
            "checkpoint.json",
            "incomplete.json",
            "progress.jsonl",
            "recovery.json",
            "recovery-checkpoint-0.json",
            "recovery-checkpoint-1.json",
        ]
        .iter()
        .any(|name| folder.join(name).exists())
    {
        return Err(Error::Input(
            "job is one-shot; resume a checkpoint into a new directory".into(),
        ));
    }
    write_new(
        &folder.join("run.claim"),
        &serde_json::json!({"job_id":manifest.job_id,"process_id":std::process::id()}),
    )?;
    state.state = State::Running;
    state.process_id = Some(std::process::id());
    atomic_status(&folder, &state, 0)?;
    let observer = FileObserver {
        job_id: manifest.job_id,
        folder: folder.clone(),
        checkpoint_interval: manifest.binding.checkpoint_interval,
        trace: Mutex::new(TraceState {
            status: state,
            events: BufWriter::new(
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(folder.join("progress.jsonl"))?,
            ),
            sequence: 0,
            last_status: Instant::now(),
            accepted_since_checkpoint: 0,
        }),
    };
    let result: Result<Outcome> = (|| {
        let mut continuation = if manifest.binding.parent_job_id.is_some() {
            Some(read::<Checkpoint>(&folder.join(SAVED))?)
        } else {
            None
        };
        let mut generation = 0_u64;
        loop {
            observer.reset_checkpoint_count()?;
            let result = if let Some(saved) = continuation.take() {
                saved.resume_observed(&mesh, ground.clone(), case.clone(), &observer)?
            } else {
                stages::solve_observed(&mesh, ground.clone(), case.clone(), &observer)?
            };
            load(&folder)?;
            if result.complete {
                write_new(&folder.join("result.json"), &result)?;
                observer.finish(State::Succeeded, Some("result.json"), None)?;
                return Ok(Outcome::Succeeded);
            }
            if result.failure_reason.as_deref() == Some("checkpoint_requested") {
                let saved = Checkpoint::seal(&mesh, result)?;
                if observer.requested() {
                    write_new(&folder.join("checkpoint.json"), &saved)?;
                    observer.finish(State::Checkpointed, Some("checkpoint.json"), None)?;
                    return Ok(Outcome::Checkpointed);
                }
                if manifest.binding.checkpoint_interval.is_none() {
                    return Err(Error::Control(
                        "unexpected internal checkpoint without a durable policy".into(),
                    ));
                }
                generation = generation
                    .checked_add(1)
                    .ok_or_else(|| Error::Control("recovery generation overflow".into()))?;
                recovery::save(&folder, &observer.job_id, generation, &saved)?;
                continuation = Some(saved);
            } else {
                let reason = result.failure_reason.clone();
                write_new(&folder.join("incomplete.json"), &result)?;
                observer.finish(State::Failed, Some("incomplete.json"), reason)?;
                return Ok(Outcome::Failed);
            }
        }
    })();
    if let Err(error) = &result {
        observer.finish(State::Failed, None, Some(error.to_string()))?;
    }
    result
}

#[cfg(all(test, windows))]
mod tests {
    use super::*;
    use std::{os::windows::fs::OpenOptionsExt, sync::Arc, time::Duration};

    fn snapshots() -> (PathBuf, PathBuf, PathBuf) {
        let directory = std::env::temp_dir().join(format!(
            "soil-fem-status-{}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            TOKENS.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&directory).unwrap();
        let target = directory.join("status.json");
        let temporary = directory.join(".owned-complete.part");
        fs::write(&target, b"{\"state\":\"old\"}").unwrap();
        fs::write(&temporary, b"{\"state\":\"new\",\"sequence\":7}").unwrap();
        (directory, temporary, target)
    }

    fn reader_without_delete_sharing(target: &Path) -> File {
        OpenOptions::new()
            .read(true)
            .share_mode(1 | 2)
            .open(target)
            .unwrap()
    }

    #[test]
    fn windows_status_replacement_waits_for_reader_without_rewriting_snapshot() {
        let (directory, temporary, target) = snapshots();
        let reader = reader_without_delete_sharing(&target);
        let error = fs::rename(&temporary, &target).unwrap_err();
        assert!(matches!(error.raw_os_error(), Some(5 | 32 | 33)));
        assert_eq!(fs::read(&target).unwrap(), b"{\"state\":\"old\"}");
        let expected = fs::read(&temporary).unwrap();
        let barrier = Arc::new(std::sync::Barrier::new(2));
        let held = barrier.clone();
        let release = std::thread::spawn(move || {
            held.wait();
            std::thread::sleep(Duration::from_millis(100));
            drop(reader);
        });
        barrier.wait();
        replace_status(&temporary, &target).unwrap();
        release.join().unwrap();
        assert_eq!(fs::read(&target).unwrap(), expected);
        assert!(!temporary.exists());
        fs::remove_file(target).unwrap();
        fs::remove_dir(directory).unwrap();
    }

    #[test]
    fn windows_status_replacement_failure_preserves_versions_and_owned_temporary() {
        let (directory, temporary, target) = snapshots();
        let old = fs::read(&target).unwrap();
        let next = fs::read(&temporary).unwrap();
        let reader = reader_without_delete_sharing(&target);
        let start = Instant::now();
        assert!(replace_status(&temporary, &target).is_err());
        assert!(start.elapsed() >= Duration::from_millis(600));
        assert_eq!(fs::read(&target).unwrap(), old);
        assert_eq!(fs::read(&temporary).unwrap(), next);
        drop(reader);
        replace_status(&temporary, &target).unwrap();
        assert_eq!(fs::read(&target).unwrap(), next);
        assert!(!temporary.exists());
        fs::remove_file(target).unwrap();
        fs::remove_dir(directory).unwrap();
    }
}
