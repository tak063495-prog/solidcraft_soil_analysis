//! Single-valued van Genuchten retention and Mualem conductivity (l=1/2).
//! Heads in metres, theta per bulk volume. No hysteresis or conductivity floor.
use crate::{Error, Result};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VanGenuchten {
    pub residual_water_content: f64,
    pub saturated_water_content: f64,
    pub alpha_per_m: f64,
    pub n: f64,
    /// Saturated compressible storage [1/m], applied only for psi>=0.
    pub saturated_specific_storage_per_m: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub effective_saturation: f64,
    pub water_saturation: f64,
    pub water_content: f64,
    /// Kept separately even when theta rounds to theta_s near saturation.
    pub water_deficit: f64,
    /// Water content plus compressible storage, relative to theta_s.
    pub storage_offset: f64,
    pub capacity_per_m: f64,
    pub relative_conductivity: f64,
    /// Finite log value distinguishes true zero from floating point underflow.
    pub log_relative_conductivity: f64,
    pub d_log_relative_conductivity_per_m: f64,
}
fn log_one_plus_exp(x: f64) -> f64 {
    if x > 0. {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}
impl VanGenuchten {
    pub fn validate(&self) -> Result<()> {
        if ![
            self.residual_water_content,
            self.saturated_water_content,
            self.alpha_per_m,
            self.n,
            self.saturated_specific_storage_per_m,
        ]
        .iter()
        .all(|v| v.is_finite())
            || self.residual_water_content < 0.
            || self.saturated_water_content > 1.
            || self.saturated_water_content <= self.residual_water_content
            || self.alpha_per_m <= 0.
            || self.n <= 1. + 1e-8
            || self.n > 100.
            || self.saturated_specific_storage_per_m < 0.
        {
            return Err(Error::Input(
                "invalid van Genuchten parameters (1+1e-8<n<=100 required)".into(),
            ));
        }
        Ok(())
    }
    pub fn evaluate(&self, psi: f64) -> Result<State> {
        self.validate()?;
        if !psi.is_finite() {
            return Err(Error::Input("non-finite retention pressure head".into()));
        }
        if psi >= 0. {
            let offset = self.saturated_specific_storage_per_m * psi;
            if !offset.is_finite() {
                return Err(Error::Solve("saturated storage overflow".into()));
            }
            return Ok(State {
                effective_saturation: 1.,
                water_saturation: 1.,
                water_content: self.saturated_water_content,
                water_deficit: 0.,
                storage_offset: offset,
                capacity_per_m: self.saturated_specific_storage_per_m,
                relative_conductivity: 1.,
                log_relative_conductivity: 0.,
                d_log_relative_conductivity_per_m: 0.,
            });
        }
        let m = (self.n - 1.) / self.n;
        let log_x = self.alpha_per_m.ln() + (-psi).ln();
        let y = self.n * log_x;
        let l = log_one_plus_exp(y);
        let log_s = if log_x > 0. {
            -(self.n - 1.) * log_x - m * log_one_plus_exp(-y)
        } else {
            -m * l
        };
        let saturation = log_s.exp();
        let difference = self.saturated_water_content - self.residual_water_content;
        let deficit = difference * (-log_s.exp_m1());
        let log_a = -l;
        let log_b = if y > 0. { -log_one_plus_exp(-y) } else { y - l };
        // First-order expansion has relative truncation < exp(-36), roundoff scale.
        let log_inner = if log_a < -36. {
            m.ln() + log_a
        } else {
            (-((m * log_b).exp_m1())).ln()
        };
        let log_kr = 0.5 * log_s + 2. * log_inner;
        let kr = log_kr.exp();
        let log_ds = m.ln()
            + self.n.ln()
            + self.alpha_per_m.ln()
            + if log_x > 0. {
                -log_x - log_one_plus_exp(-y)
            } else {
                (self.n - 1.) * log_x - l
            };
        let log_ratio = (m - 1.) * log_b + log_a - log_inner;
        let derivative = (log_ds - std::f64::consts::LN_2).exp() + 2. * (log_ds + log_ratio).exp();
        let log_c = difference.ln()
            + m.ln()
            + self.n.ln()
            + self.alpha_per_m.ln()
            + if log_x > 0. {
                -self.n * log_x - (m + 1.) * log_one_plus_exp(-y)
            } else {
                (self.n - 1.) * log_x - (m + 1.) * l
            };
        let capacity = log_c.exp();
        let theta = if log_s < -std::f64::consts::LN_2 {
            self.residual_water_content + difference * saturation
        } else {
            self.saturated_water_content - deficit
        };
        if ![saturation, deficit, capacity, theta, kr, log_kr, derivative]
            .iter()
            .all(|v| v.is_finite())
            || !(0. ..=1.).contains(&kr)
            || !(0. ..=1.).contains(&saturation)
        {
            return Err(Error::Solve("retention arithmetic overflow".into()));
        }
        Ok(State {
            effective_saturation: saturation,
            water_saturation: theta / self.saturated_water_content,
            water_content: theta,
            water_deficit: deficit,
            storage_offset: -deficit,
            capacity_per_m: capacity,
            relative_conductivity: kr,
            log_relative_conductivity: log_kr,
            d_log_relative_conductivity_per_m: derivative,
        })
    }
}
