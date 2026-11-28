//! Optional progress observation and cooperative pause for staged mechanics.
//! Optional diagnostics receive serialized copies, never mutable state. A pause is checked
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

/// Equilibrated coarse response used only by an explicitly selected C-01 continuation.
pub struct ArcSample<'a> {
    pub stage: usize,
    pub start: f64,
    pub end: f64,
    pub base: f64,
    pub old_u: &'a [[f64; 3]],
    pub trial_u: &'a [[f64; 3]],
    pub path_rejected: bool,
}

/// Callback errors are terminal runtime-control errors, never numerical
/// failures that cause a smaller physical increment to be tried.
pub trait Observer {
    fn on_progress(&self, progress: &Progress) -> Result<()>;
    /// Opt-in, bounded small-case diagnostics. Disabled for production observers.
    fn diagnostic_enabled(&self, _stage: usize, _progress: f64) -> bool {
        false
    }
    fn on_diagnostic(&self, _value: &serde_json::Value) -> Result<()> {
        Ok(())
    }
    /// Separate opt-in for A-01 cause records. Existing diagnostics are unchanged.
    fn a01_joint_cause_enabled(&self) -> bool {
        false
    }
    fn c01_initial_displacement(
        &self,
        _stage: usize,
        _progress: f64,
    ) -> Result<Option<Vec<[f64; 3]>>> {
        Ok(None)
    }
    fn c01_arc_retry(&self, _sample: ArcSample<'_>) -> Option<f64> {
        None
    }
    /// Called after acceptance, including geostatic initialization, never on trial states.
    fn on_accepted_state(
        &self,
        stage: usize,
        progress: f64,
        state: &crate::stages::Committed,
    ) -> Result<()> {
        self.on_accepted_displacement(stage, progress, &state.displacement_m)
    }
    /// Optional displacement stream at actual accepted load coordinates.
    fn on_accepted_displacement(
        &self,
        _stage: usize,
        _progress: f64,
        _u: &[[f64; 3]],
    ) -> Result<()> {
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(false)
    }
}
