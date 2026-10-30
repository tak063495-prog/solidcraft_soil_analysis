//! Accepted-increment SRM snapshots. Each strength factor starts from the same
//! original inputs; only a continuation of that SAME factor can reuse history.
use crate::{
    Error, Result,
    checkpoint::{self, Checkpoint, Identity},
    ground::Ground,
    mesh::Mesh,
    srm_search::{
        self, EndpointAudit, Outcome, Report, Requirement, Specification, StopReason, Trial,
    },
    stage_control::{Kind, Observer, Progress},
    stages::{self, Run, StageCase},
};
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Phase {
    Original,
    Searching,
    Finished,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Payload {
    phase: Phase,
    original: Option<Run>,
    basis: Option<EndpointAudit>,
    trials: Vec<Trial>,
    pending_trial: Option<Trial>,
    inflight: Option<Checkpoint>,
    report: Option<Report>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    schema: String,
    identity: Identity,
    specification_sha256: String,
    payload_sha256: String,
    payload: Payload,
}

struct BudgetObserver<'a> {
    accepted: AtomicUsize,
    limit: usize,
    external: Option<&'a dyn Observer>,
}
impl Observer for BudgetObserver<'_> {
    fn on_progress(&self, progress: &Progress) -> Result<()> {
        if let Some(external) = self.external {
            external.on_progress(progress)?;
        }
        if progress.kind == Kind::Accepted {
            self.accepted.fetch_add(1, Ordering::Relaxed);
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(self.accepted.load(Ordering::Relaxed) >= self.limit
            || self
                .external
                .map(|o| o.stop_requested())
                .transpose()?
                .unwrap_or(false))
    }
}

fn validate(mesh: &Mesh, ground: &Ground, case: &StageCase, spec: &Specification) -> Result<()> {
    mesh.validate()?;
    ground.validate()?;
    case.validate()?;
    spec.validate(ground, case)?;
    if spec
        .virtual_displacements_m
        .iter()
        .any(|v| v.len() != mesh.nodes.len())
    {
        return Err(Error::Input(
            "SRM virtual candidate node count differs from mesh".into(),
        ));
    }
    Ok(())
}

fn endpoint_case(case: &StageCase, spec: &Specification) -> StageCase {
    let mut prefix = case.clone();
    prefix.stages.truncate(spec.through_stages);
    prefix
        .verification_reject_once_stages
        .retain(|name| prefix.stages.iter().any(|stage| stage.name() == name));
    prefix
}

enum Action {
    Trial(f64),
    Finish(StopReason),
}

fn next(spec: &Specification, trials: &[Trial]) -> Action {
    for factor in [1., spec.required_factor, spec.maximum_factor] {
        if !trials.iter().any(|t| t.factor == factor) {
            return Action::Trial(factor);
        }
    }
    if trials
        .iter()
        .any(|t| t.failure_kind.as_deref() == Some("contradictory_evidence"))
    {
        return Action::Finish(StopReason::ContradictoryEvidence);
    }
    if trials.iter().any(|t| t.outcome == Outcome::Indeterminate) {
        return Action::Finish(StopReason::IndeterminateTrial);
    }
    let Some(bounds) = srm_search::interval(trials) else {
        return Action::Finish(StopReason::NoPhysicalUpper);
    };
    if bounds.width <= spec.absolute_factor_tolerance {
        Action::Finish(StopReason::FactorTolerance)
    } else if trials.len() >= spec.maximum_trials {
        Action::Finish(StopReason::TrialLimit)
    } else {
        Action::Trial(bounds.attained_equilibrium_factor + 0.5 * bounds.width)
    }
}

fn run_chunk(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    saved: Option<Checkpoint>,
    observer: &dyn Observer,
) -> Result<Run> {
    if let Some(saved) = saved {
        saved.resume_observed(mesh, ground, case, observer)
    } else {
        stages::solve_observed(mesh, ground, case, observer)
    }
}

impl Snapshot {
    pub fn start(
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        specification: &Specification,
    ) -> Result<Self> {
        validate(mesh, ground, case, specification)?;
        let mut snapshot = Self {
            schema: "soil-fem-srm-search-checkpoint/1".into(),
            identity: Identity::current(mesh, ground, case)?,
            specification_sha256: checkpoint::hash_json(specification)?,
            payload_sha256: String::new(),
            payload: Payload {
                phase: Phase::Original,
                original: None,
                basis: None,
                trials: vec![],
                pending_trial: None,
                inflight: None,
                report: None,
            },
        };
        snapshot.seal()?;
        Ok(snapshot)
    }

    fn seal(&mut self) -> Result<()> {
        self.payload_sha256 = checkpoint::hash_json(&self.payload)?;
        Ok(())
    }

    pub fn verify(
        &self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
    ) -> Result<()> {
        validate(mesh, ground, case, spec)?;
        if self.schema != "soil-fem-srm-search-checkpoint/1"
            || self.identity != Identity::current(mesh, ground, case)?
            || self.specification_sha256 != checkpoint::hash_json(spec)?
            || self.payload_sha256 != checkpoint::hash_json(&self.payload)?
        {
            return Err(Error::Input(
                "SRM snapshot input/solver/specification/history mismatch".into(),
            ));
        }
        let p = &self.payload;
        let valid = match p.phase {
            Phase::Original => {
                p.original.is_none()
                    && p.basis.is_none()
                    && p.trials.is_empty()
                    && p.pending_trial.is_none()
                    && p.report.is_none()
            }
            Phase::Searching => {
                p.original.is_some()
                    && p.basis.is_some()
                    && p.report.is_none()
                    && (p.inflight.is_none() || p.pending_trial.is_some())
            }
            Phase::Finished => {
                p.original.is_none()
                    && p.basis.is_none()
                    && p.trials.is_empty()
                    && p.pending_trial.is_none()
                    && p.inflight.is_none()
                    && p.report.is_some()
            }
        };
        if !valid {
            return Err(Error::Input("SRM snapshot phase inventory mismatch".into()));
        }
        Ok(())
    }

    pub fn is_complete(&self) -> bool {
        self.payload.phase == Phase::Finished
    }
    pub fn report(&self) -> Option<&Report> {
        self.payload.report.as_ref()
    }
    pub fn trials(&self) -> &[Trial] {
        &self.payload.trials
    }
    pub fn has_accepted_mechanical_checkpoint(&self) -> bool {
        self.payload.inflight.is_some()
    }
    pub fn current_factor(&self) -> Option<f64> {
        match self.payload.phase {
            Phase::Original => Some(1.),
            Phase::Searching => self.payload.pending_trial.as_ref().map(|t| t.factor),
            Phase::Finished => None,
        }
    }

    /// Consume an independently verified saved state and advance at most the
    /// declared number of ACCEPTED mechanical increments. On error, callers
    /// retain their last immutable saved snapshot; no failed trial is sealed.
    pub fn advance(
        mut self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
        accepted_budget: usize,
        observer: Option<&dyn Observer>,
    ) -> Result<Self> {
        self.verify(mesh, ground, case, spec)?;
        if accepted_budget == 0 || self.is_complete() {
            return Err(Error::Input(
                "SRM advance requires positive budget and unfinished snapshot".into(),
            ));
        }
        let prefix = endpoint_case(case, spec);
        let observer = BudgetObserver {
            accepted: AtomicUsize::new(0),
            limit: accepted_budget,
            external: observer,
        };
        if self.payload.phase == Phase::Original {
            let run = run_chunk(
                mesh,
                ground.clone(),
                prefix,
                self.payload.inflight.take(),
                &observer,
            )?;
            if run.failure_reason.as_deref() == Some("checkpoint_requested") {
                self.payload.inflight = Some(Checkpoint::seal(mesh, run)?);
            } else if run.complete && run.failure_reason.is_none() {
                self.payload.basis = Some(srm_search::audit_unreduced_endpoint(mesh, &run)?);
                self.payload.original = Some(run);
                self.payload.phase = Phase::Searching;
            } else {
                return Err(Error::Solve("unreduced construction prefix did not complete; no SRM origin or gate can be accepted".into()));
            }
            self.seal()?;
            return Ok(self);
        }
        if self.payload.pending_trial.is_none() {
            match next(spec, &self.payload.trials) {
                Action::Finish(reason) => {
                    self.finish(mesh, ground, case, spec, reason)?;
                    self.seal()?;
                    return Ok(self);
                }
                Action::Trial(factor) => {
                    let (trial, needs_solve) = srm_search::prepare_trial(
                        mesh,
                        ground,
                        spec,
                        self.payload.basis.as_ref().unwrap(),
                        factor,
                    )?;
                    if !needs_solve {
                        self.payload.trials.push(trial);
                        self.seal()?;
                        return Ok(self);
                    }
                    self.payload.pending_trial = Some(trial);
                }
            }
        }
        let trial = self.payload.pending_trial.take().unwrap();
        let adjusted = spec.replay(trial.factor).adjusted(ground, trial.factor)?;
        let result = run_chunk(
            mesh,
            adjusted,
            prefix,
            self.payload.inflight.take(),
            &observer,
        );
        match result {
            Ok(run) if run.failure_reason.as_deref() == Some("checkpoint_requested") => {
                self.payload.inflight = Some(Checkpoint::seal(mesh, run)?);
                self.payload.pending_trial = Some(trial);
            }
            Err(e @ (Error::Control(_) | Error::Io(_) | Error::Json(_) | Error::Input(_))) => {
                return Err(e);
            }
            result => {
                self.payload.trials.push(srm_search::finish_trial(
                    mesh,
                    self.payload.basis.as_ref().unwrap(),
                    trial,
                    result,
                )?);
            }
        }
        self.seal()?;
        Ok(self)
    }

    fn finish(
        &mut self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
        mut stop: StopReason,
    ) -> Result<()> {
        let trials = &self.payload.trials;
        if trials.iter().any(|a| {
            a.outcome == Outcome::EquilibriumAttained
                && trials
                    .iter()
                    .any(|b| b.outcome == Outcome::PhysicallyExcluded && a.factor >= b.factor)
        }) {
            stop = StopReason::ContradictoryEvidence;
        }
        let requirement = if stop == StopReason::ContradictoryEvidence {
            Requirement::Indeterminate
        } else {
            match trials
                .iter()
                .find(|t| t.factor == spec.required_factor)
                .map(|t| t.outcome)
            {
                Some(Outcome::EquilibriumAttained) => Requirement::EquilibriumRequirementMet,
                Some(Outcome::PhysicallyExcluded) => Requirement::IndependentlyExcluded,
                _ => Requirement::Indeterminate,
            }
        };
        let bounds = if stop == StopReason::ContradictoryEvidence {
            None
        } else {
            srm_search::interval(trials)
        };
        self.payload.report = Some(Report {
            schema: if spec.virtual_displacements_m.is_empty() {
                "soil-fem-srm-search-result/1"
            } else {
                "soil-fem-srm-search-result/2"
            }
            .into(),
            source_identity: self.identity.clone(),
            original_case_sha256: checkpoint::hash_json(case)?,
            original_ground_sha256: checkpoint::hash_json(ground)?,
            endpoint_case_sha256: checkpoint::hash_json(&endpoint_case(case, spec))?,
            specification: spec.clone(),
            original_unreduced_run: self.payload.original.take().unwrap(),
            original_endpoint_audit: self.payload.basis.take().unwrap(),
            trials: std::mem::take(&mut self.payload.trials),
            interval: bounds,
            requirement,
            stop_reason: stop,
            replay_origin: "original_inputs_before_initial_equilibration_for_every_factor".into(),
            unreduced_source_history_unchanged: self.identity
                == Identity::current(mesh, ground, case)?,
            continuum_stability_or_design_safety_certified: false,
        });
        self.payload.phase = Phase::Finished;
        Ok(())
    }
}
