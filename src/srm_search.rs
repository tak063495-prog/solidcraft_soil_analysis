//! Fresh-history dry MC factor search with an independent physical exclusion.
//! An attained numerical equilibrium is never called a continuum stability proof.
pub use crate::stages::srm_adapter::Endpoint as EndpointAudit;
use crate::{
    Error, Result, checkpoint,
    ground::Ground,
    mc,
    mesh::Mesh,
    replay, soil,
    srm_bound::{self, Assessment, NodalProblem, WaterPolicy},
    srm_virtual_work,
    stages::{self, StageCase},
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Specification {
    pub schema: String,
    pub material_ids: Vec<u32>,
    pub psi_policy: replay::PsiPolicy,
    pub tension_policy: replay::TensionPolicy,
    pub maximum_factor: f64,
    pub required_factor: f64,
    pub absolute_factor_tolerance: f64,
    pub maximum_trials: usize,
    pub through_stages: usize,
    /// Explicit fields for independent verification, never preceding trial states.
    /// Schema2 is required when any candidate is supplied.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub virtual_displacements_m: Vec<Vec<[f64; 3]>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Outcome {
    EquilibriumAttained,
    PhysicallyExcluded,
    Indeterminate,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Trial {
    pub factor: f64,
    pub outcome: Outcome,
    pub failure_kind: Option<String>,
    pub reason: Option<String>,
    pub original_ground_sha256: String,
    pub adjusted_ground_sha256: String,
    pub endpoint_audit: Option<EndpointAudit>,
    pub physical_exclusion: Assessment,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub virtual_work: Vec<VirtualEvidence>,
    pub run: Option<stages::Run>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VirtualEvidence {
    pub candidate_index: usize,
    pub assessment: srm_virtual_work::Assessment,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Requirement {
    EquilibriumRequirementMet,
    IndependentlyExcluded,
    Indeterminate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StopReason {
    FactorTolerance,
    NoPhysicalUpper,
    IndeterminateTrial,
    TrialLimit,
    ContradictoryEvidence,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Interval {
    pub attained_equilibrium_factor: f64,
    pub independently_excluded_factor: f64,
    pub width: f64,
    pub attained_trial_index: usize,
    pub excluded_trial_index: usize,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub source_identity: checkpoint::Identity,
    pub original_case_sha256: String,
    pub original_ground_sha256: String,
    pub endpoint_case_sha256: String,
    pub specification: Specification,
    pub original_unreduced_run: stages::Run,
    pub original_endpoint_audit: EndpointAudit,
    pub trials: Vec<Trial>,
    pub interval: Option<Interval>,
    pub requirement: Requirement,
    pub stop_reason: StopReason,
    pub replay_origin: String,
    pub unreduced_source_history_unchanged: bool,
    pub continuum_stability_or_design_safety_certified: bool,
}

impl Specification {
    pub(crate) fn replay(&self, factor: f64) -> replay::Specification {
        replay::Specification {
            schema: "soil-fem-replay/1".into(),
            material_ids: self.material_ids.clone(),
            factors: vec![factor],
            psi_policy: self.psi_policy.clone(),
            tension_policy: self.tension_policy.clone(),
        }
    }
    pub fn validate(&self, ground: &Ground, case: &StageCase) -> Result<()> {
        if !["soil-fem-srm-search/1", "soil-fem-srm-search/2"].contains(&self.schema.as_str())
            || (self.schema == "soil-fem-srm-search/1" && !self.virtual_displacements_m.is_empty())
            || self.virtual_displacements_m.len() > 16
            || self
                .virtual_displacements_m
                .iter()
                .flatten()
                .flatten()
                .any(|x| !x.is_finite())
            || !self.maximum_factor.is_finite()
            || !(1. ..=100.).contains(&self.maximum_factor)
            || self.maximum_factor == 1.
            || !self.required_factor.is_finite()
            || !(1. ..=self.maximum_factor).contains(&self.required_factor)
            || !self.absolute_factor_tolerance.is_finite()
            || !(1e-6..=1.).contains(&self.absolute_factor_tolerance)
            || !(3..=50).contains(&self.maximum_trials)
            || self.through_stages > case.stages.len()
        {
            return Err(Error::Input(
                "invalid SRM search factors, controls or stage prefix".into(),
            ));
        }
        self.replay(1.).validate(ground)?;
        self.replay(self.maximum_factor)
            .adjusted(ground, self.maximum_factor)?;
        self.replay(self.required_factor)
            .adjusted(ground, self.required_factor)?;
        if ground.materials.iter().any(|m| {
            matches!(&m.model,
            soil::Model::MohrCoulomb(p) if p.flow == mc::Flow::Davis)
        }) || ground.initialization.water_level_m.is_some()
            || ground
                .initialization
                .geostatic
                .as_ref()
                .is_some_and(|g| !g.hydrostatic_face_sets.is_empty())
            || case.joint_system.is_some()
            || case.beam_system.is_some()
            || case.anchor_system.is_some()
            || case.stages[..self.through_stages].iter().any(|s| {
                !matches!(
                    s,
                    stages::Stage::Boundary { .. }
                        | stages::Stage::Pressure { .. }
                        | stages::Stage::Excavate { .. }
                        | stages::Stage::ResetDisplay { .. }
                )
            })
            || case.boundary.as_ref().is_none_or(|b| {
                b.constraints
                    .iter()
                    .any(|c| c.value_m != 0. || c.gradient != [0.; 3])
            })
        {
            return Err(Error::Input("this search requires dry MC, explicit zero constraints and pressure/boundary/excavation/display stages; Davis monotonic search, structures, material replacement and flow need separate policies".into()));
        }
        Ok(())
    }
}

fn exclusion(mesh: &Mesh, adjusted: &Ground, endpoint: &EndpointAudit) -> Result<Assessment> {
    let mut ids = mesh
        .material_ids
        .iter()
        .zip(&endpoint.active_elements)
        .filter_map(|(&id, &enabled)| enabled.then_some(id));
    let id = ids
        .next()
        .ok_or_else(|| Error::Input("empty active SRM endpoint".into()))?;
    if ids.any(|other| other != id) {
        return Ok(Assessment::NotApplicable {
            reason: "mean-stress witness requires one homogeneous active MC ID".into(),
        });
    }
    let material = adjusted
        .materials
        .iter()
        .find(|m| m.id == id)
        .ok_or_else(|| Error::Input("unknown active SRM material".into()))?;
    let soil::Model::MohrCoulomb(parameters) = &material.model else {
        return Err(Error::Input(
            "SRM never converts a hardening material to MC".into(),
        ));
    };
    srm_bound::analyze_nodal(
        mesh,
        &NodalProblem {
            schema: "soil-fem-mc-nodal-mean-stress-bound/1".into(),
            material_id: id,
            water_policy: WaterPolicy::Dry,
            material: parameters.clone(),
            constraints: endpoint.boundary_constraints.clone(),
            physical_nodal_force_kn: endpoint.physical_force_kn.clone(),
            active_elements: endpoint.active_elements.clone(),
        },
    )
}

fn failure_kind(error: &Error) -> &'static str {
    match error {
        Error::Input(_) => "input",
        Error::Material(_) => "material",
        Error::Solve(_) => "nonlinear",
        Error::Linear(_) => "linear",
        Error::Control(_) => "control",
        Error::Io(_) => "io",
        Error::Json(_) => "json",
    }
}

pub(crate) fn prepare_trial(
    mesh: &Mesh,
    ground: &Ground,
    spec: &Specification,
    basis: &EndpointAudit,
    factor: f64,
) -> Result<(Trial, bool)> {
    let adjusted = spec.replay(factor).adjusted(ground, factor)?;
    let assessment = exclusion(mesh, &adjusted, basis)?;
    let mut excluded = matches!(&assessment, Assessment::Evaluated {certificate}
        if certificate.mean_stress_excludes_admissible_discrete_equilibrium);
    let mut virtual_work = Vec::new();
    if !excluded {
        let materials = adjusted
            .materials
            .iter()
            .map(|m| {
                let soil::Model::MohrCoulomb(p) = &m.model else {
                    unreachable!("validated all-MC")
                };
                srm_virtual_work::Material {
                    id: m.id,
                    parameters: p.clone(),
                }
            })
            .collect::<Vec<_>>();
        for (candidate_index, velocity) in spec.virtual_displacements_m.iter().enumerate() {
            let problem = srm_virtual_work::Problem {
                schema: "soil-fem-mc-virtual-work-bound/1".into(),
                water_policy: WaterPolicy::Dry,
                materials: materials.clone(),
                constraints: basis.boundary_constraints.clone(),
                physical_nodal_force_kn: basis.physical_force_kn.clone(),
                active_elements: basis.active_elements.clone(),
                virtual_displacement_m: velocity.clone(),
            };
            let result = srm_virtual_work::analyze(mesh, &problem)?;
            excluded = matches!(&result, srm_virtual_work::Assessment::Evaluated {certificate}
                if certificate.excludes_admissible_discrete_equilibrium);
            virtual_work.push(VirtualEvidence {
                candidate_index,
                assessment: result,
            });
            if excluded {
                break;
            }
        }
    }
    let mut trial = Trial {
        factor,
        outcome: Outcome::Indeterminate,
        failure_kind: None,
        reason: None,
        original_ground_sha256: checkpoint::hash_json(ground)?,
        adjusted_ground_sha256: checkpoint::hash_json(&adjusted)?,
        endpoint_audit: None,
        physical_exclusion: assessment,
        virtual_work,
        run: None,
    };
    if excluded {
        if factor == 1. {
            trial.reason = Some("independent exclusion at the already attained unreduced factor; contradictory evidence requires investigation".into());
            trial.failure_kind = Some("contradictory_evidence".into());
        } else {
            trial.outcome = Outcome::PhysicallyExcluded;
        }
        return Ok((trial, false));
    }
    Ok((trial, true))
}

pub(crate) fn finish_trial(
    mesh: &Mesh,
    basis: &EndpointAudit,
    mut trial: Trial,
    result: Result<stages::Run>,
) -> Result<Trial> {
    match result {
        Ok(run) => {
            if run.complete && run.failure_reason.is_none() {
                match stages::srm_adapter::inspect(mesh, &run) {
                    Ok(audit)
                        if audit.physical_force_sha256 == basis.physical_force_sha256
                            && audit.active_elements == basis.active_elements
                            && checkpoint::hash_json(&audit.boundary_constraints)?
                                == checkpoint::hash_json(&basis.boundary_constraints)? =>
                    {
                        trial.outcome = Outcome::EquilibriumAttained;
                        trial.endpoint_audit = Some(audit);
                    }
                    Ok(_) => {
                        trial.failure_kind = Some("endpoint_binding".into());
                        trial.reason = Some(
                            "trial changed the reconstructed endpoint load, domain or constraints"
                                .into(),
                        );
                    }
                    Err(e) => {
                        trial.failure_kind = Some(failure_kind(&e).into());
                        trial.reason = Some(e.to_string());
                    }
                }
            } else {
                trial.failure_kind = Some("incomplete".into());
                trial.reason = run.failure_reason.clone();
            }
            trial.run = Some(run);
        }
        Err(e) => {
            trial.failure_kind = Some(failure_kind(&e).into());
            trial.reason = Some(e.to_string());
        }
    }
    Ok(trial)
}

fn evaluate(
    mesh: &Mesh,
    ground: &Ground,
    case: &StageCase,
    spec: &Specification,
    basis: &EndpointAudit,
    factor: f64,
) -> Result<Trial> {
    let (trial, solve_required) = prepare_trial(mesh, ground, spec, basis, factor)?;
    if !solve_required {
        return Ok(trial);
    }
    // Fresh Ground and StageCase on EVERY call. No preceding Run is input here.
    let adjusted = spec.replay(factor).adjusted(ground, factor)?;
    finish_trial(
        mesh,
        basis,
        trial,
        stages::solve(mesh, adjusted, case.clone()),
    )
}

pub(crate) fn interval(trials: &[Trial]) -> Option<Interval> {
    let (lo_index, lo) = trials
        .iter()
        .enumerate()
        .filter(|(_, t)| t.outcome == Outcome::EquilibriumAttained)
        .max_by(|(_, a), (_, b)| a.factor.total_cmp(&b.factor))?;
    let (hi_index, hi) = trials
        .iter()
        .enumerate()
        .filter(|(_, t)| t.outcome == Outcome::PhysicallyExcluded)
        .min_by(|(_, a), (_, b)| a.factor.total_cmp(&b.factor))?;
    (lo.factor < hi.factor).then_some(Interval {
        attained_equilibrium_factor: lo.factor,
        independently_excluded_factor: hi.factor,
        width: hi.factor - lo.factor,
        attained_trial_index: lo_index,
        excluded_trial_index: hi_index,
    })
}

pub fn run(
    mesh: &Mesh,
    ground: Ground,
    original_case: StageCase,
    spec: Specification,
) -> Result<Report> {
    mesh.validate()?;
    ground.validate()?;
    original_case.validate()?;
    spec.validate(&ground, &original_case)?;
    if spec
        .virtual_displacements_m
        .iter()
        .any(|v| v.len() != mesh.nodes.len())
    {
        return Err(Error::Input(
            "SRM virtual candidate node count differs from mesh".into(),
        ));
    }
    let original_case_sha256 = checkpoint::hash_json(&original_case)?;
    let original_ground_sha256 = checkpoint::hash_json(&ground)?;
    let source_identity = checkpoint::Identity::current(mesh, &ground, &original_case)?;
    let mut case = original_case.clone();
    case.stages.truncate(spec.through_stages);
    case.verification_reject_once_stages
        .retain(|name| case.stages.iter().any(|stage| stage.name() == name));
    let endpoint_case_sha256 = checkpoint::hash_json(&case)?;
    let origin = stages::solve(mesh, ground.clone(), case.clone())?;
    if !origin.complete || origin.failure_reason.is_some() {
        return Err(Error::Solve(
            "unreduced construction prefix did not complete; no SRM origin or gate can be accepted"
                .into(),
        ));
    }
    let origin_sha = checkpoint::hash_json(&origin)?;
    let basis = stages::srm_adapter::inspect(mesh, &origin)?;
    let mut trials = Vec::new();
    let mut stop = StopReason::NoPhysicalUpper;
    for factor in [1., spec.required_factor, spec.maximum_factor] {
        if trials.iter().any(|t: &Trial| t.factor == factor) {
            continue;
        }
        trials.push(evaluate(mesh, &ground, &case, &spec, &basis, factor)?);
    }
    if trials
        .iter()
        .any(|t| t.failure_kind.as_deref() == Some("contradictory_evidence"))
    {
        stop = StopReason::ContradictoryEvidence;
    } else if trials.iter().any(|t| t.outcome == Outcome::Indeterminate) {
        stop = StopReason::IndeterminateTrial;
    } else {
        while let Some(bounds) = interval(&trials) {
            if bounds.width <= spec.absolute_factor_tolerance {
                stop = StopReason::FactorTolerance;
                break;
            }
            if trials.len() >= spec.maximum_trials {
                stop = StopReason::TrialLimit;
                break;
            }
            let midpoint = bounds.attained_equilibrium_factor + 0.5 * bounds.width;
            let t = evaluate(mesh, &ground, &case, &spec, &basis, midpoint)?;
            let indeterminate = t.outcome == Outcome::Indeterminate;
            trials.push(t);
            if indeterminate {
                stop = StopReason::IndeterminateTrial;
                break;
            }
        }
    }
    let conflicting = trials.iter().any(|a| {
        a.outcome == Outcome::EquilibriumAttained
            && trials
                .iter()
                .any(|b| b.outcome == Outcome::PhysicallyExcluded && a.factor >= b.factor)
    });
    if conflicting {
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
        interval(&trials)
    };
    let unchanged = checkpoint::hash_json(&ground)? == original_ground_sha256
        && checkpoint::hash_json(&original_case)? == original_case_sha256
        && checkpoint::hash_json(&origin)? == origin_sha;
    if !unchanged {
        return Err(Error::Solve(
            "SRM changed its immutable unreduced inputs or history".into(),
        ));
    }
    Ok(Report {
        schema: if spec.virtual_displacements_m.is_empty() {
            "soil-fem-srm-search-result/1".into()
        } else {
            "soil-fem-srm-search-result/2".into()
        },
        source_identity,
        original_case_sha256,
        original_ground_sha256,
        endpoint_case_sha256,
        specification: spec,
        original_unreduced_run: origin,
        original_endpoint_audit: basis,
        trials,
        interval: bounds,
        requirement,
        stop_reason: stop,
        replay_origin: "original_inputs_before_initial_equilibration_for_every_factor".into(),
        unreduced_source_history_unchanged: unchanged,
        continuum_stability_or_design_safety_certified: false,
    })
}

/// Reconstruct and audit a supplied original Run; this never accepts an external
/// exclusion certificate and is useful for checking saved construction origins.
pub fn audit_unreduced_endpoint(mesh: &Mesh, run: &stages::Run) -> Result<EndpointAudit> {
    stages::srm_adapter::inspect(mesh, run)
}
