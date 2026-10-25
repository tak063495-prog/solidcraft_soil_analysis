//! Optional progress observation and cooperative pause for staged mechanics.
//! Observers cannot access trial/committed material state. A pause is checked
//! only after an accepted increment, including geostatic initialization.
use crate::{Result, linalg::LinearReport};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    AttemptStarted,
    Newton,
    LinearSolved,
    Accepted,
    Rejected,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Progress {
    pub schema: String,
    pub kind: Kind,
    /// Zero is geostatic initialization; construction indices start at one.
    pub stage_index: usize,
    pub stage: String,
    /// Actual attempted load/release fraction, not an overall time estimate.
    pub progress: f64,
    pub newton_iteration: Option<usize>,
    pub relative_force_residual: Option<f64>,
    pub relative_volume_residual: Option<f64>,
    pub relative_pressure_residual: Option<f64>,
    pub maximum_material_partitions: Option<usize>,
    pub linear_solve: Option<LinearReport>,
    pub reason: Option<String>,
}

/// Callback errors are terminal runtime-control errors, never numerical
/// failures that cause a smaller physical increment to be tried.
pub trait Observer {
    fn on_progress(&self, progress: &Progress) -> Result<()>;
    fn stop_requested(&self) -> Result<bool> {
        Ok(false)
    }
}
