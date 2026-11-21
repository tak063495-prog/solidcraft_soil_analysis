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

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathErrorControl {
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
