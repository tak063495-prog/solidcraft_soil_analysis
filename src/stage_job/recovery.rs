//! Two alternating, fully flushed accepted-state snapshots and an atomic pointer.
//! A running job may be forked from this immutable accepted snapshot; status/PID
//! is deliberately not treated as evidence that its original process has died.
use super::{Manifest, read, read_small, replace_status};
use crate::{
    Error, Result,
    checkpoint::{self, Checkpoint},
    ground::Ground,
    mesh::Mesh,
    stages::StageCase,
};
use serde::{Deserialize, Serialize};
use std::{
    fs::OpenOptions,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Pointer {
    schema: String,
    job_id: String,
    generation: u64,
    checkpoint: String,
    checkpoint_sha256: String,
    run_sha256: String,
}

fn write_synced(path: &Path, value: &impl Serialize) -> Result<()> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut stream = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut stream, value)?;
    stream.flush()?;
    stream.get_ref().sync_all()?;
    Ok(())
}

pub(super) fn save(folder: &Path, job_id: &str, generation: u64, saved: &Checkpoint) -> Result<()> {
    if generation == 0 {
        return Err(Error::Control(
            "recovery generation must be positive".into(),
        ));
    }
    let name = format!("recovery-checkpoint-{}.json", generation % 2);
    let data_tmp = folder.join(format!(
        ".recovery-data-{}-{generation}.part",
        std::process::id()
    ));
    write_synced(&data_tmp, saved)?;
    let digest = checkpoint::hash_file(&data_tmp)?;
    replace_status(&data_tmp, &folder.join(&name))?;
    let pointer = Pointer {
        schema: "soil-fem-stages-recovery/1".into(),
        job_id: job_id.into(),
        generation,
        checkpoint: name,
        checkpoint_sha256: digest,
        run_sha256: saved.run_sha256.clone(),
    };
    let pointer_tmp = folder.join(format!(
        ".recovery-pointer-{}-{generation}.part",
        std::process::id()
    ));
    write_synced(&pointer_tmp, &pointer)?;
    replace_status(&pointer_tmp, &folder.join("recovery.json"))?;
    Ok(())
}

pub(super) fn saved_path(
    folder: &Path,
    manifest: &Manifest,
    mesh: &Mesh,
    ground: &Ground,
    case: &StageCase,
) -> Result<PathBuf> {
    let pointer: Pointer = read_small(&folder.join("recovery.json"), 8192)?;
    if manifest.binding.checkpoint_interval.is_none()
        || pointer.schema != "soil-fem-stages-recovery/1"
        || pointer.job_id != manifest.job_id
        || pointer.generation == 0
        || pointer.checkpoint != format!("recovery-checkpoint-{}.json", pointer.generation % 2)
    {
        return Err(Error::Input(
            "recovery pointer identity or path mismatch".into(),
        ));
    }
    let path = folder.join(&pointer.checkpoint);
    if checkpoint::hash_file(&path)? != pointer.checkpoint_sha256 {
        return Err(Error::Input(
            "recovery checkpoint file checksum mismatch".into(),
        ));
    }
    let saved: Checkpoint = read(&path)?;
    saved.verify(mesh, ground, case)?;
    if saved.run_sha256 != pointer.run_sha256
        || saved.identity != manifest.binding.identity
        || checkpoint::hash_file(&path)? != pointer.checkpoint_sha256
    {
        return Err(Error::Input(
            "recovery checkpoint changed or belongs to another job".into(),
        ));
    }
    Ok(path)
}
