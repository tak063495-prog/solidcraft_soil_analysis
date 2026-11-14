//! Read-only, evidence-based finite-element Fs interval reporting.
//! No Newton iterations, continuation advances or changes to production gates.
use crate::{
    Error, Result, checkpoint, epoch_srm_control, epoch_srm_search,
    ground::Ground,
    mesh::Mesh,
    srm_search::{self, Outcome, Requirement, StopReason, Trial},
    stages::StageCase,
};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct Assessment {
    pub schema: &'static str,
    pub exact_fs_for_this_case: bool,
    pub attained_equilibrium_factor: Option<f64>,
    pub independently_excluded_factor: Option<f64>,
    pub width: Option<f64>,
    pub absolute_factor_tolerance: f64,
    pub stop_reason: StopReason,
    pub source_report_sha256: String,
    pub continuum_stability_or_design_safety_certified: bool,
}
fn same(a: &impl Serialize, b: &impl Serialize) -> Result<()> {
    if checkpoint::hash_json(a)? != checkpoint::hash_json(b)? {
        return Err(Error::Input(
            "Fs evidence/input reconstruction mismatch".into(),
        ));
    }
    Ok(())
}
fn assessment(
    bounds: Option<&srm_search::Interval>,
    tolerance: f64,
    reason: StopReason,
    report: &impl Serialize,
) -> Result<Assessment> {
    let exact = reason == StopReason::FactorTolerance
        && bounds.is_some_and(|b| b.width > 0. && b.width <= tolerance);
    Ok(Assessment {
        schema: "soil-fem-exact-fs-assessment/1",
        exact_fs_for_this_case: exact,
        attained_equilibrium_factor: bounds.map(|b| b.attained_equilibrium_factor),
        independently_excluded_factor: bounds.map(|b| b.independently_excluded_factor),
        width: bounds.map(|b| b.width),
        absolute_factor_tolerance: tolerance,
        stop_reason: reason,
        source_report_sha256: checkpoint::hash_json(report)?,
        continuum_stability_or_design_safety_certified: false,
    })
}

/// Bind the separately supplied SPEC; widening the report's embedded tolerance
/// or replacing its path controls cannot make it pass.
pub fn epoch(
    mesh: &Mesh,
    origin: &crate::stages::Run,
    spec: &epoch_srm_search::Specification,
    report: &epoch_srm_search::Report,
) -> Result<Assessment> {
    same(origin, &report.original_unreduced_run)?;
    same(spec, &report.specification)?;
    epoch_srm_control::verify_report(mesh, report)?;
    assessment(
        report.interval.as_ref(),
        spec.search.absolute_factor_tolerance,
        report.stop_reason,
        report,
    )
}

fn next(
    spec: &srm_search::Specification,
    trials: &[Trial],
) -> std::result::Result<f64, StopReason> {
    for f in [1., spec.required_factor, spec.maximum_factor] {
        if !trials.iter().any(|t| t.factor == f) {
            return Ok(f);
        }
    }
    if trials
        .iter()
        .any(|t| t.failure_kind.as_deref() == Some("contradictory_evidence"))
    {
        return Err(StopReason::ContradictoryEvidence);
    }
    if trials.iter().any(|t| t.outcome == Outcome::Indeterminate) {
        return Err(StopReason::IndeterminateTrial);
    }
    let Some(b) = srm_search::interval(trials) else {
        return Err(StopReason::NoPhysicalUpper);
    };
    if b.width <= spec.absolute_factor_tolerance {
        return Err(StopReason::FactorTolerance);
    }
    if trials.len() >= spec.maximum_trials {
        return Err(StopReason::TrialLimit);
    }
    Ok(b.attained_equilibrium_factor + 0.5 * b.width)
}

pub fn fresh(
    mesh: &Mesh,
    ground: &Ground,
    case: &StageCase,
    spec: &srm_search::Specification,
    report: &srm_search::Report,
) -> Result<Assessment> {
    let _validated = crate::srm_search_control::Snapshot::start(mesh, ground, case, spec)?;
    same(spec, &report.specification)?;
    same(
        &checkpoint::Identity::current(mesh, ground, case)?,
        &report.source_identity,
    )?;
    same(
        &checkpoint::hash_json(ground)?,
        &report.original_ground_sha256,
    )?;
    same(&checkpoint::hash_json(case)?, &report.original_case_sha256)?;
    let mut prefix = case.clone();
    prefix.stages.truncate(spec.through_stages);
    prefix
        .verification_reject_once_stages
        .retain(|n| prefix.stages.iter().any(|s| s.name() == n));
    same(
        &checkpoint::hash_json(&prefix)?,
        &report.endpoint_case_sha256,
    )?;
    same(ground, &report.original_unreduced_run.ground)?;
    same(&prefix, &report.original_unreduced_run.case)?;
    let basis = srm_search::audit_unreduced_endpoint(mesh, &report.original_unreduced_run)?;
    same(&basis, &report.original_endpoint_audit)?;
    if !report.unreduced_source_history_unchanged
        || report.continuum_stability_or_design_safety_certified
        || report.replay_origin != "original_inputs_before_initial_equilibration_for_every_factor"
        || report.schema
            != if spec.virtual_displacements_m.is_empty() {
                "soil-fem-srm-search-result/1"
            } else {
                "soil-fem-srm-search-result/2"
            }
    {
        return Err(Error::Input(
            "Fs result schema/origin/claims changed".into(),
        ));
    }
    for (i, t) in report.trials.iter().enumerate() {
        if next(spec, &report.trials[..i]) != Ok(t.factor) {
            return Err(Error::Input(
                "Fs trial missing, reordered or added after termination".into(),
            ));
        }
        let (expected, solve) = srm_search::prepare_trial(mesh, ground, spec, &basis, t.factor)?;
        if solve {
            if let Some(run) = &t.run {
                let adjusted = spec.replay(t.factor).adjusted(ground, t.factor)?;
                same(&adjusted, &run.ground)?;
                same(&prefix, &run.case)?;
                let run = serde_json::from_slice(&serde_json::to_vec(run)?)?;
                let expected = srm_search::finish_trial(mesh, &basis, expected, Ok(run))?;
                same(&expected, t)?;
            } else {
                // An exception has no independently reproducible mechanics.
                // Preserve unknown, but never accept an equilibrium or upper.
                if t.outcome != Outcome::Indeterminate
                    || t.failure_kind.is_none()
                    || t.endpoint_audit.is_some()
                {
                    return Err(Error::Input(
                        "missing completed trial or forged exclusion".into(),
                    ));
                }
                same(&expected.physical_exclusion, &t.physical_exclusion)?;
                same(&expected.virtual_work, &t.virtual_work)?;
                same(&expected.adjusted_ground_sha256, &t.adjusted_ground_sha256)?;
                same(&expected.original_ground_sha256, &t.original_ground_sha256)?;
            }
        } else {
            same(&expected, t)?;
        }
    }
    let mut reason = next(spec, &report.trials)
        .err()
        .ok_or_else(|| Error::Input("unfinished Fs search".into()))?;
    if report.trials.iter().any(|a| {
        a.outcome == Outcome::EquilibriumAttained
            && report
                .trials
                .iter()
                .any(|b| b.outcome == Outcome::PhysicallyExcluded && a.factor >= b.factor)
    }) {
        reason = StopReason::ContradictoryEvidence;
    }
    let bounds = if reason == StopReason::ContradictoryEvidence {
        None
    } else {
        srm_search::interval(&report.trials)
    };
    let req = if reason == StopReason::ContradictoryEvidence {
        Requirement::Indeterminate
    } else {
        match report
            .trials
            .iter()
            .find(|t| t.factor == spec.required_factor)
            .map(|t| t.outcome)
        {
            Some(Outcome::EquilibriumAttained) => Requirement::EquilibriumRequirementMet,
            Some(Outcome::PhysicallyExcluded) => Requirement::IndependentlyExcluded,
            _ => Requirement::Indeterminate,
        }
    };
    same(&bounds, &report.interval)?;
    same(&req, &report.requirement)?;
    same(&reason, &report.stop_reason)?;
    assessment(
        bounds.as_ref(),
        spec.absolute_factor_tolerance,
        reason,
        report,
    )
}
