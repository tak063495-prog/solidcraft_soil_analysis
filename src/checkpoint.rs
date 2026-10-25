//! Versioned checkpoints bind the entire accepted history and load path to
//! canonical mesh/ground/case data and the exact running solver binary.
use crate::{
    Error, Result,
    ground::Ground,
    mesh::Mesh,
    stages::{self, Run, StageCase},
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

struct HashWriter(Sha256);
impl Write for HashWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
pub fn hash_json(value: &impl Serialize) -> Result<String> {
    let mut writer = HashWriter(Sha256::new());
    serde_json::to_writer(&mut writer, value)?;
    Ok(format!("{:x}", writer.0.finalize()))
}
pub fn hash_file(path: &Path) -> Result<String> {
    let mut reader = File::open(path)?;
    let mut writer = HashWriter(Sha256::new());
    io::copy(&mut reader, &mut writer)?;
    Ok(format!("{:x}", writer.0.finalize()))
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Identity {
    pub mesh_sha256: String,
    pub ground_sha256: String,
    pub case_sha256: String,
    pub solver_binary_sha256: String,
    pub solver_version: String,
}
impl Identity {
    pub(crate) fn current(mesh: &Mesh, ground: &Ground, case: &StageCase) -> Result<Self> {
        Ok(Self {
            mesh_sha256: hash_json(mesh)?,
            ground_sha256: hash_json(ground)?,
            case_sha256: hash_json(case)?,
            solver_binary_sha256: hash_file(&std::env::current_exe()?)?,
            solver_version: env!("CARGO_PKG_VERSION").into(),
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Checkpoint {
    pub schema: String,
    pub identity: Identity,
    pub run_sha256: String,
    pub run: Run,
}
impl Checkpoint {
    pub fn seal(mesh: &Mesh, run: Run) -> Result<Self> {
        if run.complete
            || (run.continuation.is_none() && run.geostatic.is_none())
            || run.failure_reason.as_deref() != Some("checkpoint_requested")
        {
            return Err(Error::Input(
                "only a requested checkpoint at an accepted increment can be sealed".into(),
            ));
        }
        if run.continuation.is_none()
            && run
                .geostatic
                .as_ref()
                .is_none_or(|g| !g.attempts.iter().any(|a| a.accepted))
        {
            return Err(Error::Input(
                "geostatic checkpoint requires an accepted increment".into(),
            ));
        }
        let identity = Identity::current(mesh, &run.ground, &run.case)?;
        if run.provenance != identity {
            return Err(Error::Input(
                "checkpoint producer/input provenance mismatch".into(),
            ));
        }
        let run_sha256 = hash_json(&run)?;
        Ok(Self {
            schema: "soil-fem-checkpoint/1".into(),
            identity,
            run_sha256,
            run,
        })
    }
    pub fn verify(&self, mesh: &Mesh, ground: &Ground, case: &StageCase) -> Result<()> {
        if self.schema != "soil-fem-checkpoint/1"
            || self.identity != Identity::current(mesh, ground, case)?
        {
            return Err(Error::Input("checkpoint input/solver hash mismatch".into()));
        }
        if self.run_sha256 != hash_json(&self.run)? {
            return Err(Error::Input(
                "checkpoint history/load checksum mismatch".into(),
            ));
        }
        Ok(())
    }
    pub fn resume(
        self,
        mesh: &Mesh,
        ground: Ground,
        case: StageCase,
        stop_after: Option<usize>,
    ) -> Result<Run> {
        self.verify(mesh, &ground, &case)?;
        if stop_after == Some(0) {
            return Err(Error::Input("checkpoint interval must be positive".into()));
        }
        stages::resume(mesh, ground, case, self.run, stop_after)
    }
    pub fn resume_observed(
        self,
        mesh: &Mesh,
        ground: Ground,
        case: StageCase,
        observer: &dyn crate::stage_control::Observer,
    ) -> Result<Run> {
        self.verify(mesh, &ground, &case)?;
        stages::resume_observed(mesh, ground, case, self.run, observer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn standard_sha256_vector_and_canonical_json_are_stable() -> Result<()> {
        let mut writer = HashWriter(Sha256::new());
        writer.write_all(b"abc")?;
        assert_eq!(
            format!("{:x}", writer.0.finalize()),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        let a =
            serde_json::from_str::<std::collections::BTreeMap<String, f64>>(r#"{"b":2,"a":1}"#)?;
        let b = serde_json::from_str::<std::collections::BTreeMap<String, f64>>(
            r#"{"a":1.0,"b":2.0}"#,
        )?;
        assert_eq!(hash_json(&a)?, hash_json(&b)?);
        Ok(())
    }
}
