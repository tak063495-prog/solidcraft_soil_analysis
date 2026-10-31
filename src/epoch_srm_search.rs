//! Strength-factor trials from one accepted unreduced construction origin.
//! Every factor owns a new material-path controller; trial histories never mix.
use crate::{
    Error, Result, checkpoint,
    mc_epoch_path::{self, Controller},
    mesh::Mesh,
    soil,
    srm_search::{self, EndpointAudit, Interval, Outcome, Requirement, StopReason},
    stages::{self, Run},
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Specification {
    pub schema: String,
    pub search: srm_search::Specification,
    pub initial_increment: f64,
    pub minimum_increment: f64,
    pub maximum_attempts: usize,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EpochTrial {
    /// Shared physical-exclusion/endpoint evidence. run is always None here:
    /// current material history is an explicit epoch, never a fabricated Run.
    pub evidence: srm_search::Trial,
    pub epoch: Option<mc_epoch_path::Checkpoint>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub source_identity: checkpoint::Identity,
    pub origin_run_sha256: String,
    pub origin_state_sha256: String,
    pub specification: Specification,
    pub original_unreduced_run: Run,
    pub original_endpoint_audit: EndpointAudit,
    pub trials: Vec<EpochTrial>,
    pub interval: Option<Interval>,
    pub requirement: Requirement,
    pub stop_reason: StopReason,
    pub trial_origin: String,
    pub unreduced_source_history_unchanged: bool,
    pub continuum_stability_or_design_safety_certified: bool,
}

impl Specification {
    pub(crate) fn path(&self, origin: &Run, factor: f64) -> Result<mc_epoch_path::Specification> {
        let adjusted = self
            .search
            .replay(factor)
            .adjusted(&origin.ground, factor)?;
        let target_parameters = adjusted
            .materials
            .into_iter()
            .map(|m| match m.model {
                soil::Model::MohrCoulomb(p) => Ok(p),
                _ => Err(Error::Input(
                    "epoch SRM never converts MCC/IDH to MC".into(),
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(mc_epoch_path::Specification {
            schema: "soil-fem-mc-epoch-path/1".into(),
            target_parameters,
            initial_increment: self.initial_increment,
            minimum_increment: self.minimum_increment,
            maximum_attempts: self.maximum_attempts,
        })
    }

    pub(crate) fn validate(&self, mesh: &Mesh, origin: &Run) -> Result<()> {
        if self.schema != "soil-fem-epoch-srm-search/1"
            || self.search.through_stages != origin.case.stages.len()
            || self
                .search
                .virtual_displacements_m
                .iter()
                .any(|v| v.len() != mesh.nodes.len())
        {
            return Err(Error::Input(
                "invalid epoch SRM schema, complete-origin stage count or virtual-field dimension"
                    .into(),
            ));
        }
        self.search.validate(&origin.ground, &origin.case)?;
        // Existing mechanics enforces exact producer identity, completed dry
        // all-MC origin, fixed construction data, parameter/path controls.
        for factor in [1., self.search.required_factor, self.search.maximum_factor] {
            let _validated = Controller::new(mesh, origin, self.path(origin, factor)?)?;
        }
        Ok(())
    }
}

fn evaluate(
    mesh: &Mesh,
    origin: &Run,
    specification: &Specification,
    basis: &EndpointAudit,
    factor: f64,
) -> Result<EpochTrial> {
    let (mut evidence, needs_solve) =
        srm_search::prepare_trial(mesh, &origin.ground, &specification.search, basis, factor)?;
    if !needs_solve {
        return Ok(EpochTrial {
            evidence,
            epoch: None,
        });
    }
    let mut controller = Controller::new(mesh, origin, specification.path(origin, factor)?)?;
    while !controller.complete() && controller.continuation().stop_reason.is_none() {
        controller.advance()?;
    }
    let saved = controller.checkpoint()?;
    // Validate all accepted current parameters/history and frozen metadata,
    // including after an incomplete path. A rejected Newton trial is absent.
    let checked = Controller::resume(mesh, origin, &saved)?;
    if controller.complete() {
        if !checked.complete() {
            return Err(Error::Control(
                "completed epoch path lost completion after validation".into(),
            ));
        }
        let audit =
            stages::srm_adapter::inspect_epoch(mesh, origin, &saved.payload.material_checkpoint)?;
        if audit.physical_force_sha256 != basis.physical_force_sha256
            || audit.active_elements != basis.active_elements
            || checkpoint::hash_json(&audit.boundary_constraints)?
                != checkpoint::hash_json(&basis.boundary_constraints)?
        {
            return Err(Error::Control(
                "epoch trial changed the bound construction endpoint".into(),
            ));
        }
        evidence.outcome = Outcome::EquilibriumAttained;
        evidence.endpoint_audit = Some(audit);
    } else {
        evidence.failure_kind = Some("epoch_path_incomplete".into());
        evidence.reason = Some(format!(
            "accepted parameter progress {}: {:?}",
            saved.payload.continuation.progress, saved.payload.continuation.stop_reason
        ));
    }
    Ok(EpochTrial {
        evidence,
        epoch: Some(saved),
    })
}

pub(crate) fn interval(trials: &[EpochTrial]) -> Option<Interval> {
    let (lo_index, lo) = trials
        .iter()
        .enumerate()
        .filter(|(_, t)| t.evidence.outcome == Outcome::EquilibriumAttained)
        .max_by(|(_, a), (_, b)| a.evidence.factor.total_cmp(&b.evidence.factor))?;
    let (hi_index, hi) = trials
        .iter()
        .enumerate()
        .filter(|(_, t)| t.evidence.outcome == Outcome::PhysicallyExcluded)
        .min_by(|(_, a), (_, b)| a.evidence.factor.total_cmp(&b.evidence.factor))?;
    (lo.evidence.factor < hi.evidence.factor).then_some(Interval {
        attained_equilibrium_factor: lo.evidence.factor,
        independently_excluded_factor: hi.evidence.factor,
        width: hi.evidence.factor - lo.evidence.factor,
        attained_trial_index: lo_index,
        excluded_trial_index: hi_index,
    })
}

/// Validate a saved epoch through the public search API. The internal
/// stages adapter stays private; no construction state is mutated.
pub fn audit_epoch_endpoint(
    mesh: &Mesh,
    origin: &Run,
    saved: &stages::McEpochCheckpoint,
) -> Result<EndpointAudit> {
    stages::srm_adapter::inspect_epoch(mesh, origin, saved)
}

pub fn run(mesh: &Mesh, origin: Run, specification: Specification) -> Result<Report> {
    mesh.validate()?;
    specification.validate(mesh, &origin)?;
    let origin_run_sha256 = checkpoint::hash_json(&origin)?;
    let origin_state_sha256 = checkpoint::hash_json(&origin.last_committed)?;
    let basis = srm_search::audit_unreduced_endpoint(mesh, &origin)?;
    let mut trials: Vec<EpochTrial> = Vec::new();
    for factor in [
        1.,
        specification.search.required_factor,
        specification.search.maximum_factor,
    ] {
        if trials.iter().any(|t| t.evidence.factor == factor) {
            continue;
        }
        trials.push(evaluate(mesh, &origin, &specification, &basis, factor)?);
    }
    let mut stop = StopReason::NoPhysicalUpper;
    if trials
        .iter()
        .any(|t| t.evidence.failure_kind.as_deref() == Some("contradictory_evidence"))
    {
        stop = StopReason::ContradictoryEvidence;
    } else if trials
        .iter()
        .any(|t| t.evidence.outcome == Outcome::Indeterminate)
    {
        stop = StopReason::IndeterminateTrial;
    } else {
        while let Some(bounds) = interval(&trials) {
            if bounds.width <= specification.search.absolute_factor_tolerance {
                stop = StopReason::FactorTolerance;
                break;
            }
            if trials.len() >= specification.search.maximum_trials {
                stop = StopReason::TrialLimit;
                break;
            }
            let factor = bounds.attained_equilibrium_factor + 0.5 * bounds.width;
            let trial = evaluate(mesh, &origin, &specification, &basis, factor)?;
            let unknown = trial.evidence.outcome == Outcome::Indeterminate;
            trials.push(trial);
            if unknown {
                stop = StopReason::IndeterminateTrial;
                break;
            }
        }
    }
    if trials.iter().any(|a| {
        a.evidence.outcome == Outcome::EquilibriumAttained
            && trials.iter().any(|b| {
                b.evidence.outcome == Outcome::PhysicallyExcluded
                    && a.evidence.factor >= b.evidence.factor
            })
    }) {
        stop = StopReason::ContradictoryEvidence;
    }
    let requirement = if stop == StopReason::ContradictoryEvidence {
        Requirement::Indeterminate
    } else {
        match trials
            .iter()
            .find(|t| t.evidence.factor == specification.search.required_factor)
            .map(|t| t.evidence.outcome)
        {
            Some(Outcome::EquilibriumAttained) => Requirement::EquilibriumRequirementMet,
            Some(Outcome::PhysicallyExcluded) => Requirement::IndependentlyExcluded,
            _ => Requirement::Indeterminate,
        }
    };
    let bounds = if stop == StopReason::ContradictoryEvidence {
        None
    } else {
        interval(&trials)
    };
    let unchanged = checkpoint::hash_json(&origin)? == origin_run_sha256
        && checkpoint::hash_json(&origin.last_committed)? == origin_state_sha256;
    if !unchanged {
        return Err(Error::Control(
            "epoch SRM changed its unreduced construction origin".into(),
        ));
    }
    Ok(Report {
        schema: "soil-fem-epoch-srm-search-result/1".into(),
        source_identity: origin.provenance.clone(),
        origin_run_sha256,
        origin_state_sha256,
        specification,
        original_unreduced_run: origin,
        original_endpoint_audit: basis,
        trials,
        interval: bounds,
        requirement,
        stop_reason: stop,
        trial_origin: "same_completed_unreduced_construction_Run_for_every_independent_factor"
            .into(),
        unreduced_source_history_unchanged: unchanged,
        continuum_stability_or_design_safety_certified: false,
    })
}
