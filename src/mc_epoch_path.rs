//! Adaptive, post-construction MC parameter continuation. Lambda is a numerical
//! material path coordinate, never an inferred safety factor or collapse bound.
use crate::{
    Error, Result,
    checkpoint::hash_json,
    mc::Parameters,
    mesh::Mesh,
    stages::{self, McEpochCheckpoint, McParameterMechanics, Run},
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Specification {
    pub schema: String,
    pub target_parameters: Vec<Parameters>,
    pub initial_increment: f64,
    pub minimum_increment: f64,
    pub maximum_attempts: usize,
}
impl Specification {
    fn validate(&self, start: &[Parameters]) -> Result<()> {
        if self.schema != "soil-fem-mc-epoch-path/1"
            || !self.initial_increment.is_finite()
            || self.initial_increment <= 0.
            || self.initial_increment > 1.
            || !self.minimum_increment.is_finite()
            || self.minimum_increment <= 0.
            || self.minimum_increment > self.initial_increment
            || self.minimum_increment < 1e-8
            || !(1..=10000).contains(&self.maximum_attempts)
        {
            return Err(Error::Input("invalid MC parameter path controls".into()));
        }
        stages::validate_mc_epoch_parameters(start, &self.target_parameters)?;
        if start
            .iter()
            .zip(&self.target_parameters)
            .any(|(a, b)| a.tension_cutoff_kpa.is_some() != b.tension_cutoff_kpa.is_some())
        {
            return Err(Error::Input(
                "MC parameter path retains the Rankine surface family".into(),
            ));
        }
        Ok(())
    }
    /// c and tension interpolate linearly; tan(phi) and tan(psi) interpolate
    /// linearly. Endpoint values are exact. This policy is explicit, not an Fs law.
    fn parameters_at(&self, start: &[Parameters], lambda: f64) -> Result<Vec<Parameters>> {
        if !lambda.is_finite() || !(0. ..=1.).contains(&lambda) {
            return Err(Error::Input("MC path coordinate must be in[0,1]".into()));
        }
        if lambda == 0. {
            return Ok(start.to_vec());
        }
        if lambda == 1. {
            return Ok(self.target_parameters.clone());
        }
        start
            .iter()
            .zip(&self.target_parameters)
            .map(|(a, b)| {
                let mix = |x: f64, y: f64| (1. - lambda) * x + lambda * y;
                let angle = |x: f64, y: f64| {
                    mix(x.to_radians().tan(), y.to_radians().tan())
                        .atan()
                        .to_degrees()
                };
                let p = Parameters {
                    cohesion_kpa: mix(a.cohesion_kpa, b.cohesion_kpa),
                    friction_deg: angle(a.friction_deg, b.friction_deg),
                    dilation_deg: angle(a.dilation_deg, b.dilation_deg),
                    tension_cutoff_kpa: a
                        .tension_cutoff_kpa
                        .zip(b.tension_cutoff_kpa)
                        .map(|(x, y)| mix(x, y)),
                    ..a.clone()
                };
                p.validate()?;
                Ok(p)
            })
            .collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FailureKind {
    Input,
    Material,
    Solver,
    Linear,
    Control,
    Io,
    Json,
}
impl FailureKind {
    fn recoverable(self) -> bool {
        matches!(self, Self::Material | Self::Solver | Self::Linear)
    }
    fn of(error: &Error) -> Self {
        match error {
            Error::Input(_) => Self::Input,
            Error::Material(_) => Self::Material,
            Error::Solve(_) => Self::Solver,
            Error::Linear(_) => Self::Linear,
            Error::Control(_) => Self::Control,
            Error::Io(_) => Self::Io,
            Error::Json(_) => Self::Json,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StopReason {
    CutbackExhausted,
    AttemptLimit,
    NonRecoverableTrial,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Advance {
    Accepted,
    Rejected,
    Finished,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Attempt {
    pub from: f64,
    pub to: f64,
    pub increment: f64,
    pub accepted: bool,
    pub failure_kind: Option<FailureKind>,
    pub reason: Option<String>,
    pub before_state_sha256: String,
    pub after_state_sha256: String,
    pub before_parameters_sha256: String,
    pub after_parameters_sha256: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Accepted {
    pub progress: f64,
    pub state_sha256: String,
    pub parameters_sha256: String,
    pub newton_iterations: usize,
    pub linear_solves: usize,
    pub maximum_true_linear_residual: Option<f64>,
    pub relative_force_residual: f64,
    pub maximum_normalized_yield: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Continuation {
    pub progress: f64,
    pub next_increment: f64,
    pub attempts: Vec<Attempt>,
    pub accepted: Vec<Accepted>,
    pub stop_reason: Option<StopReason>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Checkpoint {
    pub schema: String,
    pub payload_sha256: String,
    pub payload: Payload,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Payload {
    pub specification: Specification,
    pub continuation: Continuation,
    pub material_checkpoint: McEpochCheckpoint,
}

pub struct Controller<'a> {
    mechanics: McParameterMechanics<'a>,
    origin: &'a Run,
    start: Vec<Parameters>,
    specification: Specification,
    continuation: Continuation,
}
impl<'a> Controller<'a> {
    pub fn new(mesh: &'a Mesh, origin: &'a Run, specification: Specification) -> Result<Self> {
        let mechanics = McParameterMechanics::new(mesh, origin)?;
        let start = mechanics.parameters().to_vec();
        specification.validate(&start)?;
        let continuation = Continuation {
            progress: 0.,
            next_increment: specification.initial_increment,
            attempts: vec![],
            accepted: vec![],
            stop_reason: None,
        };
        Ok(Self {
            mechanics,
            origin,
            start,
            specification,
            continuation,
        })
    }
    pub fn mechanics(&self) -> &McParameterMechanics<'a> {
        &self.mechanics
    }
    pub fn set_observer(&mut self, observer: &'a dyn crate::stage_control::Observer) {
        self.mechanics.set_observer(observer);
    }
    pub fn continuation(&self) -> &Continuation {
        &self.continuation
    }
    pub fn complete(&self) -> bool {
        self.continuation.progress == 1.
    }
    pub fn stability_classification(&self) -> &'static str {
        "not_evaluated"
    }
    pub fn advance(&mut self) -> Result<Advance> {
        if self.complete() || self.continuation.stop_reason.is_some() {
            return Ok(Advance::Finished);
        }
        let state = &mut self.continuation;
        if state.attempts.len() >= self.specification.maximum_attempts {
            state.stop_reason = Some(StopReason::AttemptLimit);
            return Ok(Advance::Finished);
        }
        let increment = state.next_increment.min(1. - state.progress);
        let to = (state.progress + increment).min(1.);
        let before_state = hash_json(self.mechanics.committed())?;
        let before_parameters = hash_json(&self.mechanics.parameters())?;
        let parameters = self.specification.parameters_at(&self.start, to)?;
        let result = self.mechanics.try_parameters(parameters);
        let after_state = hash_json(self.mechanics.committed())?;
        let after_parameters = hash_json(&self.mechanics.parameters())?;
        let mut attempt = Attempt {
            from: state.progress,
            to,
            increment,
            accepted: false,
            failure_kind: None,
            reason: None,
            before_state_sha256: before_state.clone(),
            before_parameters_sha256: before_parameters.clone(),
            after_state_sha256: after_state.clone(),
            after_parameters_sha256: after_parameters.clone(),
        };
        let outcome = match result {
            Ok(step) => {
                attempt.accepted = true;
                let maximum = step
                    .linear_solves
                    .iter()
                    .map(|r| r.true_relative_residual)
                    .reduce(f64::max);
                state.accepted.push(Accepted {
                    progress: to,
                    state_sha256: after_state,
                    parameters_sha256: after_parameters,
                    newton_iterations: step.iterations,
                    linear_solves: step.linear_solves.len(),
                    maximum_true_linear_residual: maximum,
                    relative_force_residual: step.frame.relative_force_residual,
                    maximum_normalized_yield: step.frame.maximum_normalized_yield,
                });
                state.progress = to;
                state.next_increment = increment;
                Advance::Accepted
            }
            Err(error) => {
                if before_state != after_state || before_parameters != after_parameters {
                    return Err(Error::Control(
                        "rejected MC path trial changed an accepted material epoch".into(),
                    ));
                }
                let kind = FailureKind::of(&error);
                attempt.failure_kind = Some(kind);
                attempt.reason = Some(error.to_string());
                state.next_increment = increment * 0.5;
                if !kind.recoverable() {
                    state.stop_reason = Some(StopReason::NonRecoverableTrial);
                } else if state.next_increment < self.specification.minimum_increment {
                    state.stop_reason = Some(StopReason::CutbackExhausted);
                }
                Advance::Rejected
            }
        };
        state.attempts.push(attempt);
        if state.progress < 1.
            && state.stop_reason.is_none()
            && state.attempts.len() >= self.specification.maximum_attempts
        {
            state.stop_reason = Some(StopReason::AttemptLimit);
        }
        Ok(outcome)
    }
    pub fn checkpoint(&self) -> Result<Checkpoint> {
        let payload = Payload {
            specification: self.specification.clone(),
            continuation: self.continuation.clone(),
            material_checkpoint: self.mechanics.checkpoint()?,
        };
        Ok(Checkpoint {
            schema: "soil-fem-mc-epoch-path-checkpoint/1".into(),
            payload_sha256: hash_json(&payload)?,
            payload,
        })
    }
    pub fn resume(mesh: &'a Mesh, origin: &'a Run, checkpoint: &Checkpoint) -> Result<Self> {
        if checkpoint.schema != "soil-fem-mc-epoch-path-checkpoint/1"
            || checkpoint.payload_sha256 != hash_json(&checkpoint.payload)?
        {
            return Err(Error::Input("MC path checkpoint seal mismatch".into()));
        }
        let payload = &checkpoint.payload;
        let mut controller = Self::new(mesh, origin, payload.specification.clone())?;
        let restored = McParameterMechanics::resume(mesh, origin, &payload.material_checkpoint)?;
        controller.validate_continuation(&payload.continuation, &restored)?;
        controller.mechanics = restored;
        controller.continuation = payload.continuation.clone();
        Ok(controller)
    }
    fn validate_continuation(
        &self,
        saved: &Continuation,
        restored: &McParameterMechanics<'_>,
    ) -> Result<()> {
        let invalid =
            || Error::Input("MC path checkpoint continuation/history chain mismatch".into());
        if !saved.progress.is_finite()
            || !(0. ..=1.).contains(&saved.progress)
            || !saved.next_increment.is_finite()
            || saved.next_increment <= 0.
            || saved.next_increment > self.specification.initial_increment
            || saved.attempts.len() > self.specification.maximum_attempts
        {
            return Err(invalid());
        }
        let mut progress = 0.;
        let mut increment = self.specification.initial_increment;
        let mut state = hash_json(&self.origin.last_committed)?;
        let mut parameters = hash_json(&self.start)?;
        let mut accepted = 0;
        let mut stopped = None;
        for (i, a) in saved.attempts.iter().enumerate() {
            if stopped.is_some() || progress == 1. {
                return Err(invalid());
            }
            let step = increment.min(1. - progress);
            let to = (progress + step).min(1.);
            if a.from != progress
                || a.to != to
                || a.increment != step
                || a.before_state_sha256 != state
                || a.before_parameters_sha256 != parameters
            {
                return Err(invalid());
            }
            if a.accepted {
                let s = saved.accepted.get(accepted).ok_or_else(invalid)?;
                let expected = hash_json(&self.specification.parameters_at(&self.start, to)?)?;
                if a.failure_kind.is_some()
                    || a.reason.is_some()
                    || s.progress != to
                    || a.after_parameters_sha256 != expected
                    || s.parameters_sha256 != expected
                    || s.state_sha256 != a.after_state_sha256
                    || s.newton_iterations > self.origin.case.settings.maximum_newton_iterations
                    || !s.relative_force_residual.is_finite()
                    || s.relative_force_residual < 0.
                    || s.relative_force_residual
                        > self.origin.case.settings.relative_force_tolerance
                    || !s.maximum_normalized_yield.is_finite()
                    || s.maximum_normalized_yield > 1e-8
                    || s.maximum_true_linear_residual.is_some_and(|r| {
                        !r.is_finite()
                            || r < 0.
                            || r > self.origin.case.settings.relative_linear_tolerance
                    })
                    || (s.linear_solves == 0) != s.maximum_true_linear_residual.is_none()
                {
                    return Err(invalid());
                }
                state = a.after_state_sha256.clone();
                parameters = expected;
                progress = to;
                increment = step;
                accepted += 1;
            } else {
                let kind = a.failure_kind.ok_or_else(invalid)?;
                if a.reason.as_ref().is_none_or(|r| r.is_empty())
                    || a.after_state_sha256 != state
                    || a.after_parameters_sha256 != parameters
                {
                    return Err(invalid());
                }
                increment = step * 0.5;
                if !kind.recoverable() {
                    stopped = Some(StopReason::NonRecoverableTrial);
                } else if increment < self.specification.minimum_increment {
                    stopped = Some(StopReason::CutbackExhausted);
                }
            }
            if progress < 1. && stopped.is_none() && i + 1 >= self.specification.maximum_attempts {
                stopped = Some(StopReason::AttemptLimit);
            }
        }
        if accepted != saved.accepted.len()
            || saved.progress != progress
            || saved.next_increment != increment
            || saved.stop_reason != stopped
            || hash_json(restored.committed())? != state
            || hash_json(&restored.parameters())? != parameters
            || hash_json(&self.specification.parameters_at(&self.start, progress)?)? != parameters
        {
            return Err(invalid());
        }
        Ok(())
    }
}
