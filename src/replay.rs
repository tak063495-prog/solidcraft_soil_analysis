//! Fresh constitutive/stage replay for every MC strength factor. This module
//! deliberately does not infer collapse or an Fs interval from solver status.
use crate::{
    Error, Result,
    checkpoint::{self, Identity},
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{self, Run, StageCase},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PsiPolicy {
    Keep,
    Zero,
    Cap,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TensionPolicy {
    Keep,
    DivideByFactor,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Specification {
    pub schema: String,
    pub material_ids: Vec<u32>,
    pub factors: Vec<f64>,
    pub psi_policy: PsiPolicy,
    pub tension_policy: TensionPolicy,
}
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Outcome {
    Completed,
    Incomplete,
    RejectedInput,
    MaterialError,
    SolverError,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Trial {
    pub factor: f64,
    pub outcome: Outcome,
    pub reason: Option<String>,
    pub adjusted_ground: Option<Ground>,
    pub adjusted_ground_sha256: Option<String>,
    pub run: Option<Run>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub source_identity: Identity,
    pub source_ground: Ground,
    pub case: StageCase,
    pub specification: Specification,
    pub specification_sha256: String,
    pub replay_origin: String,
    pub stability_classification: String,
    pub trials: Vec<Trial>,
}

impl Specification {
    pub fn validate(&self, ground: &Ground) -> Result<()> {
        if self.schema != "soil-fem-replay/1"
            || self.material_ids.is_empty()
            || self.material_ids.iter().collect::<BTreeSet<_>>().len() != self.material_ids.len()
            || self.factors.is_empty()
            || self.factors.len() > 50
            || self
                .factors
                .iter()
                .any(|f| !f.is_finite() || !(1. ..=100.).contains(f))
            || self
                .factors
                .iter()
                .map(|f| f.to_bits())
                .collect::<BTreeSet<_>>()
                .len()
                != self.factors.len()
        {
            return Err(Error::Input("replay requires unique material IDs, unique finite factors in [1,100], at most 50 trials".into()));
        }
        if ground
            .materials
            .iter()
            .any(|m| !matches!(m.model, soil::Model::MohrCoulomb(_)))
        {
            return Err(Error::Input("c-phi replay currently requires all-MC ground; MCC is not converted or strength-reduced".into()));
        }
        if self
            .material_ids
            .iter()
            .any(|id| !ground.materials.iter().any(|m| m.id == *id))
        {
            return Err(Error::Input("unknown replay material ID".into()));
        }
        Ok(())
    }
    pub fn adjusted(&self, source: &Ground, factor: f64) -> Result<Ground> {
        self.validate(source)?;
        if !self.factors.contains(&factor) {
            return Err(Error::Input("factor is not in replay specification".into()));
        }
        let mut ground = source.clone();
        for material in &mut ground.materials {
            if !self.material_ids.contains(&material.id) {
                continue;
            }
            let soil::Model::MohrCoulomb(p) = &mut material.model else {
                unreachable!()
            };
            p.cohesion_kpa /= factor;
            if factor != 1. {
                p.friction_deg = (p.friction_deg.to_radians().tan() / factor)
                    .atan()
                    .to_degrees();
            }
            p.dilation_deg = match self.psi_policy {
                PsiPolicy::Keep => p.dilation_deg,
                PsiPolicy::Zero => 0.,
                PsiPolicy::Cap => p.dilation_deg.min(p.friction_deg),
            };
            if matches!(self.tension_policy, TensionPolicy::DivideByFactor) {
                p.tension_cutoff_kpa = p.tension_cutoff_kpa.map(|t| t / factor);
            }
        }
        ground.validate()?;
        Ok(ground)
    }
}

pub fn run(
    mesh: &Mesh,
    source: Ground,
    case: StageCase,
    specification: Specification,
) -> Result<Report> {
    mesh.validate()?;
    source.validate()?;
    case.validate()?;
    specification.validate(&source)?;
    let source_identity = Identity::current(mesh, &source, &case)?;
    let specification_sha256 = checkpoint::hash_json(&specification)?;
    let mut trials = Vec::new();
    for &factor in &specification.factors {
        // Never use any preceding factor's displacement, material state or load path.
        let adjusted = specification.adjusted(&source, factor);
        let (adjusted_ground, result) = match adjusted {
            Ok(ground) => {
                let result = stages::solve(mesh, ground.clone(), case.clone());
                (Some(ground), result)
            }
            Err(e) => (None, Err(e)),
        };
        let (outcome, reason, run) = match result {
            Ok(run) => {
                let outcome = if run.complete {
                    Outcome::Completed
                } else {
                    Outcome::Incomplete
                };
                (outcome, run.failure_reason.clone(), Some(run))
            }
            Err(e) => {
                let outcome = match &e {
                    Error::Input(_) => Outcome::RejectedInput,
                    Error::Material(_) => Outcome::MaterialError,
                    Error::Solve(_) => Outcome::SolverError,
                    _ => return Err(e),
                };
                (outcome, Some(e.to_string()), None)
            }
        };
        let adjusted_ground_sha256 = adjusted_ground
            .as_ref()
            .map(checkpoint::hash_json)
            .transpose()?;
        trials.push(Trial {
            factor,
            outcome,
            reason,
            adjusted_ground,
            adjusted_ground_sha256,
            run,
        });
    }
    Ok(Report {
        schema: "soil-fem-replay-result/1".into(),
        source_identity,
        source_ground: source,
        case,
        specification,
        specification_sha256,
        replay_origin: "original_inputs_before_initial_equilibration".into(),
        stability_classification: "not_evaluated".into(),
        trials,
    })
}

/// A construction endpoint for independent reduced-strength replays. This is
/// not reduction applied after an unreduced endpoint, nor an Fs classification.
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StageGateReport {
    pub schema: String,
    pub source_identity: Identity,
    pub original_case: StageCase,
    pub requested_through_construction_stages: usize,
    pub strength_applied_from: String,
    pub stability_classification: String,
    pub replay: Report,
}

pub fn run_through_stage(
    mesh: &Mesh,
    source: Ground,
    case: StageCase,
    specification: Specification,
    through_stages: usize,
) -> Result<StageGateReport> {
    mesh.validate()?;
    source.validate()?;
    case.validate()?;
    specification.validate(&source)?;
    if through_stages > case.stages.len() {
        return Err(Error::Input(
            "replay stage gate exceeds the original construction stages".into(),
        ));
    }
    let source_identity = Identity::current(mesh, &source, &case)?;
    let mut prefix = case.clone();
    prefix.stages.truncate(through_stages);
    let names: BTreeSet<_> = prefix
        .stages
        .iter()
        .map(|stage| stage.name().to_owned())
        .collect();
    // A discarded suffix's diagnostic injection is unexecuted, rather than an
    // unknown stage in the effective prefix. Keep the original case in full.
    prefix
        .verification_reject_once_stages
        .retain(|name| names.contains(name));
    prefix.validate()?;
    let replay = run(mesh, source, prefix, specification)?;
    Ok(StageGateReport {
        schema: "soil-fem-replay-stage-gate/1".into(),
        source_identity,
        original_case: case,
        requested_through_construction_stages: through_stages,
        strength_applied_from: "original_inputs_before_initial_equilibration".into(),
        stability_classification: "not_evaluated".into(),
        replay,
    })
}
