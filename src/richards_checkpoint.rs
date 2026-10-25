//! Hydraulic checkpoint binds all accepted heads, storage, surface decisions,
//! time-step history and runtime to exact inputs and the running executable.
use crate::{
    Error, Result, checkpoint,
    mesh::Mesh,
    richards::{self, Case, Report},
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Checkpoint {
    pub schema: String,
    pub payload_sha256: String,
    pub run: Report,
}
impl Checkpoint {
    pub fn seal(run: Report) -> Result<Self> {
        if run.complete
            || run.failure_reason.is_some()
            || run.runtime.is_none()
            || run.solver_binary_sha256 != checkpoint::hash_file(&std::env::current_exe()?)?
        {
            return Err(Error::Input("only an intermediate accepted hydraulic state can be checkpointed by its own executable".into()));
        }
        Ok(Self {
            schema: "soil-fem-richards-checkpoint/1".into(),
            payload_sha256: checkpoint::hash_json(&run)?,
            run,
        })
    }
    pub fn resume(
        &self,
        mesh: &Mesh,
        case: Case,
        active: &[bool],
        material_ids: Option<&[u32]>,
        accepted_steps: Option<usize>,
    ) -> Result<Report> {
        if self.schema != "soil-fem-richards-checkpoint/1"
            || self.payload_sha256 != checkpoint::hash_json(&self.run)?
        {
            return Err(Error::Input("hydraulic checkpoint payload changed".into()));
        }
        richards::resume(
            mesh,
            case,
            active,
            material_ids,
            self.run.clone(),
            accepted_steps,
        )
    }
}
