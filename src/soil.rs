//! Typed constitutive dispatch. MCC keeps its legacy JSON representation;
//! MC stores its own history and cannot masquerade as an MCC pressure state.
use crate::{Error, Result, elastic, hardening, mc, mcc};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone)]
pub enum Model {
    Mcc(mcc::Parameters),
    MohrCoulomb(mc::Parameters),
    InvariantDoubleHardening(hardening::Parameters),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum State {
    Mcc(mcc::State),
    MohrCoulomb(mc::State),
    InvariantDoubleHardening(hardening::State),
}
#[derive(Debug, Serialize)]
#[serde(untagged)]
pub enum InitialPointState {
    Mcc {
        mcc_state: mcc::State,
    },
    MohrCoulomb {
        mohr_coulomb_state: mc::State,
    },
    InvariantDoubleHardening {
        invariant_double_hardening_state: hardening::State,
    },
}
impl InitialPointState {
    pub fn into_state(self) -> State {
        match self {
            Self::Mcc { mcc_state } => State::Mcc(mcc_state),
            Self::MohrCoulomb { mohr_coulomb_state } => State::MohrCoulomb(mohr_coulomb_state),
            Self::InvariantDoubleHardening {
                invariant_double_hardening_state,
            } => State::InvariantDoubleHardening(invariant_double_hardening_state),
        }
    }
    pub fn as_mcc(&self) -> Option<&mcc::State> {
        match self {
            Self::Mcc { mcc_state } => Some(mcc_state),
            _ => None,
        }
    }
}
impl From<State> for InitialPointState {
    fn from(s: State) -> Self {
        match s {
            State::Mcc(mcc_state) => Self::Mcc { mcc_state },
            State::MohrCoulomb(mohr_coulomb_state) => Self::MohrCoulomb { mohr_coulomb_state },
            State::InvariantDoubleHardening(invariant_double_hardening_state) => {
                Self::InvariantDoubleHardening {
                    invariant_double_hardening_state,
                }
            }
        }
    }
}
impl State {
    pub fn stress(&self) -> &[f64; 6] {
        match self {
            Self::Mcc(s) => &s.effective_stress_compression_positive_kpa,
            Self::MohrCoulomb(s) => &s.effective_stress_compression_positive_kpa,
            Self::InvariantDoubleHardening(s) => &s.effective_stress_compression_positive_kpa,
        }
    }
    pub fn plastic_strain(&self) -> &[f64; 6] {
        match self {
            Self::Mcc(s) => &s.plastic_strain_compression_positive,
            Self::MohrCoulomb(s) => &s.plastic_strain_compression_positive,
            Self::InvariantDoubleHardening(s) => &s.plastic_strain_compression_positive,
        }
    }
    pub fn void_ratio(&self) -> f64 {
        match self {
            Self::Mcc(s) => s.void_ratio,
            Self::MohrCoulomb(s) => s.void_ratio,
            Self::InvariantDoubleHardening(s) => s.void_ratio,
        }
    }
    pub fn as_mcc(&self) -> Option<&mcc::State> {
        match self {
            Self::Mcc(s) => Some(s),
            _ => None,
        }
    }
}
pub struct Update {
    pub state: State,
    pub substeps: usize,
    pub plastic_substeps: usize,
}
fn mismatch() -> Error {
    Error::Input("constitutive state type does not match material model".into())
}
impl Model {
    pub fn as_mcc(&self) -> Option<&mcc::Parameters> {
        match self {
            Self::Mcc(p) => Some(p),
            _ => None,
        }
    }
    pub fn is_mcc(&self) -> bool {
        matches!(self, Self::Mcc(_))
    }
    pub fn validate(&self) -> Result<()> {
        match self {
            Self::Mcc(p) => p.validate(),
            Self::MohrCoulomb(p) => p.validate(),
            Self::InvariantDoubleHardening(p) => p.validate(),
        }
    }
    pub fn validate_state(&self, s: &State) -> Result<()> {
        match (self, s) {
            (Self::Mcc(p), State::Mcc(s)) => p.validate_state(s),
            (Self::MohrCoulomb(p), State::MohrCoulomb(s)) => p.validate_state(s),
            (Self::InvariantDoubleHardening(p), State::InvariantDoubleHardening(s)) => {
                p.validate_state(s)
            }
            _ => Err(mismatch()),
        }
    }
    pub fn normalized_yield(&self, s: &State) -> Result<f64> {
        match (self, s) {
            (Self::Mcc(p), State::Mcc(s)) => Ok(p.normalized_yield(s)),
            (Self::MohrCoulomb(p), State::MohrCoulomb(s)) => Ok(p.normalized_yield(s)),
            (Self::InvariantDoubleHardening(p), State::InvariantDoubleHardening(s)) => {
                Ok(p.normalized_yield(s))
            }
            _ => Err(mismatch()),
        }
    }
    pub fn moduli(&self, s: &State) -> Result<(f64, f64)> {
        match (self, s) {
            (Self::Mcc(p), State::Mcc(s)) => p.moduli(s),
            (Self::MohrCoulomb(p), State::MohrCoulomb(_)) => Ok(p.moduli()),
            (Self::InvariantDoubleHardening(p), State::InvariantDoubleHardening(s)) => p.moduli(s),
            _ => Err(mismatch()),
        }
    }
    pub fn integrate_partition(&self, s: &State, de: [f64; 6], steps: usize) -> Result<Update> {
        match (self, s) {
            (Self::InvariantDoubleHardening(p), State::InvariantDoubleHardening(s)) => {
                let u = p.integrate_partition(s, de, steps)?;
                Ok(Update {
                    state: State::InvariantDoubleHardening(u.state),
                    substeps: u.substeps,
                    plastic_substeps: u.plastic_substeps,
                })
            }
            (Self::Mcc(p), State::Mcc(s)) => {
                let u = p.integrate_partition(s, de, steps)?;
                Ok(Update {
                    state: State::Mcc(u.state),
                    substeps: u.substeps,
                    plastic_substeps: u.plastic_substeps,
                })
            }
            (Self::MohrCoulomb(p), State::MohrCoulomb(s)) => {
                let u = p.integrate_partition(s, de, steps)?;
                Ok(Update {
                    state: State::MohrCoulomb(u.state),
                    substeps: u.substeps,
                    plastic_substeps: u.plastic_substeps,
                })
            }
            _ => Err(mismatch()),
        }
    }
    pub fn integrate(&self, s: &State, de: [f64; 6], maximum: f64) -> Result<Update> {
        match (self, s) {
            (Self::InvariantDoubleHardening(p), State::InvariantDoubleHardening(s)) => {
                let u = p.integrate(s, de, maximum)?;
                Ok(Update {
                    state: State::InvariantDoubleHardening(u.state),
                    substeps: u.substeps,
                    plastic_substeps: u.plastic_substeps,
                })
            }
            (Self::Mcc(p), State::Mcc(s)) => {
                let u = p.integrate(s, de, maximum)?;
                Ok(Update {
                    state: State::Mcc(u.state),
                    substeps: u.substeps,
                    plastic_substeps: u.plastic_substeps,
                })
            }
            (Self::MohrCoulomb(p), State::MohrCoulomb(s)) => {
                let u = p.integrate(s, de, maximum)?;
                Ok(Update {
                    state: State::MohrCoulomb(u.state),
                    substeps: u.substeps,
                    plastic_substeps: u.plastic_substeps,
                })
            }
            _ => Err(mismatch()),
        }
    }
    pub fn stress_tangent_partition(
        &self,
        s: &State,
        de: [f64; 6],
        steps: usize,
        h: f64,
    ) -> Result<(Update, elastic::Tensor)> {
        match (self, s) {
            (Self::InvariantDoubleHardening(p), State::InvariantDoubleHardening(s)) => {
                let (u, d) = p.stress_tangent_partition(s, de, steps, h)?;
                Ok((
                    Update {
                        state: State::InvariantDoubleHardening(u.state),
                        substeps: u.substeps,
                        plastic_substeps: u.plastic_substeps,
                    },
                    d,
                ))
            }
            (Self::Mcc(p), State::Mcc(s)) => {
                let (u, d) = p.stress_tangent_partition(s, de, steps, h)?;
                Ok((
                    Update {
                        state: State::Mcc(u.state),
                        substeps: u.substeps,
                        plastic_substeps: u.plastic_substeps,
                    },
                    d,
                ))
            }
            (Self::MohrCoulomb(p), State::MohrCoulomb(s)) => {
                let (u, d) = p.stress_tangent_partition(s, de, steps, h)?;
                Ok((
                    Update {
                        state: State::MohrCoulomb(u.state),
                        substeps: u.substeps,
                        plastic_substeps: u.plastic_substeps,
                    },
                    d,
                ))
            }
            _ => Err(mismatch()),
        }
    }
}
