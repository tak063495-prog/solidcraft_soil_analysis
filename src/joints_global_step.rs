//! Opt-in global contact-switch cutback. This never interpolates a state.
use crate::{Error, Result, joints};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GlobalStepControl {
    /// Resolve endpoint status changes down to base_increment / 2^bisections.
    /// This is a load-path resolution parameter, not a force tolerance.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub switch_bisections: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path_error: Option<PathErrorControl>,
}

fn is_zero(n: &u32) -> bool {
    *n == 0
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EventComparison {
    Observe,
    BisectRejections,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathErrorControl {
    /// Use signed-gap event localization and fit the proposal to the next base boundary.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub event_proposal: bool,
    /// Use the measured coarse/fine defect to propose a retry and clip growth at the next anchor.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub error_scaled_retry: bool,
    /// Bisect the local partition floor at previously slipped or currently open contact points.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub contact_history_bisection: bool,
    /// Record the first signed-gap zero. The bisection control varies only rejected steps.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub event_comparison: Option<EventComparison>,
    /// Use the generalized derivative of the open, sticking or sliding active branch.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub semismooth_contact: bool,
    /// Lower bound on local return-map subdivisions for every coarse and half trial.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub minimum_contact_substeps: Option<usize>,
    pub slip_scale_m: f64,
    pub dissipation_scale_kpa_m: f64,
}

impl GlobalStepControl {
    pub(crate) fn validate(&self, base: f64, minimum: f64) -> Result<()> {
        if let Some(c) = &self.path_error {
            if self.switch_bisections != 0
                || !c.slip_scale_m.is_finite()
                || c.slip_scale_m <= 0.
                || !c.dissipation_scale_kpa_m.is_finite()
                || c.dissipation_scale_kpa_m <= 0.
                || minimum > base * 0.5
            {
                return Err(Error::Input("JOINT path error control requires positive finite scales, no switch_bisections, and minimum increment no greater than half the base".into()));
            }
            if c.minimum_contact_substeps
                .is_some_and(|n| !(1..=50_000).contains(&n))
            {
                return Err(Error::Input(
                    "JOINT minimum_contact_substeps outside existing partition budget".into(),
                ));
            }
            if c.event_comparison.is_some() && !c.event_proposal {
                return Err(Error::Input(
                    "event_comparison requires event_proposal".into(),
                ));
            }
            if c.error_scaled_retry && (!c.event_proposal || c.event_comparison.is_some()) {
                return Err(Error::Input(
                    "error_scaled_retry requires unmodified event localization".into(),
                ));
            }
            if c.contact_history_bisection
                && (!c.event_proposal || c.minimum_contact_substeps != Some(4))
            {
                return Err(Error::Input(
                    "contact_history_bisection requires event4".into(),
                ));
            }
            return Ok(());
        }
        if !(1..=8).contains(&self.switch_bisections)
            || minimum > base / 2_f64.powi(self.switch_bisections as i32 + 1)
        {
            return Err(Error::Input(
                "JOINT global step control requires 1..=8 switch_bisections and a minimum increment no greater than half its switch resolution".into(),
            ));
        }
        Ok(())
    }

    pub(crate) fn next(progress: f64, increment: f64, base: f64) -> f64 {
        let mut boundary = ((progress / base).floor() + 1.) * base;
        // Only protect the arithmetic of a prescribed boundary (e.g. 3 * 0.1).
        // This is independent of all equilibrium and history acceptance gates.
        if boundary <= progress + 16. * f64::EPSILON {
            boundary += base;
        }
        (progress + increment).min(boundary).min(1.)
    }

    pub(crate) fn reject_switch(
        &self,
        old: &Option<joints::State>,
        new: &Option<joints::State>,
        actual_increment: f64,
        base: f64,
    ) -> bool {
        if self.path_error.is_some() {
            return false;
        }
        let changed = match (old, new) {
            (Some(a), Some(b)) => a
                .points
                .iter()
                .zip(&b.points)
                .any(|(a, b)| a.status != b.status),
            _ => false,
        };
        changed
            && actual_increment
                > (base / 2_f64.powi(self.switch_bisections as i32)) * (1. + 16. * f64::EPSILON)
    }
}
