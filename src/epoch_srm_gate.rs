//! Construction-stage Fs gate driven by an unreduced accepted epoch.
//! Later stages start only when Fs_req attained numerical equilibrium and the
//! search finished without a trial limit or contradictory evidence.
//! Nonconvergence stays indeterminate: it is neither an exclusion nor collapse.
use crate::{
    Error, Result,
    checkpoint::{self, Identity},
    epoch_srm_control::Snapshot as SearchSnapshot,
    epoch_srm_search::{Report as SearchReport, Specification as SearchSpecification},
    ground::Ground,
    mesh::Mesh,
    srm_search::{Requirement, StopReason},
    stage_control::Observer,
    stages::{self, Run, StageCase},
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Specification {
    pub schema: String,
    pub search: SearchSpecification,
    /// Completed stage counts, strictly increasing, each in `1..=stages.len()`.
    /// `search.search.through_stages` must be 0; each gate assigns its own prefix.
    pub gate_after_stages: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GateEffect {
    Continue,
    StopExcluded,
    StopIndeterminate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Decision {
    Passed,
    StoppedExcluded,
    StoppedIndeterminate,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GateRecord {
    pub through_stages: usize,
    pub requirement: Requirement,
    pub stop_reason: StopReason,
    pub effect: GateEffect,
    pub search: SearchReport,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub source_identity: Identity,
    pub specification: Specification,
    pub gates: Vec<GateRecord>,
    pub decision: Decision,
    pub construction: Run,
    pub construction_stage_count: usize,
    pub stopped_before_later_stages: bool,
    pub nonconvergence_used_as_exclusion_or_collapse: bool,
    pub continuum_stability_or_design_safety_certified: bool,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CompletedGate {
    through_stages: usize,
    origin: Run,
    search: SearchSnapshot,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ConstructPhase {
    through_stages: usize,
    checkpoint: Option<Run>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SearchPhase {
    through_stages: usize,
    origin: Run,
    search: SearchSnapshot,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct FinishPhase {
    checkpoint: Option<Run>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Phase {
    Construct(Box<ConstructPhase>),
    Search(Box<SearchPhase>),
    Finish(Box<FinishPhase>),
    Done,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Payload {
    phase: Phase,
    records: Vec<CompletedGate>,
    construction: Option<Run>,
    report: Option<Report>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    schema: String,
    identity: Identity,
    case_stage_count: usize,
    specification_sha256: String,
    payload_sha256: String,
    payload: Payload,
}

fn duplicate<T: Serialize + for<'de> Deserialize<'de>>(value: &T) -> Result<T> {
    Ok(serde_json::from_slice(&serde_json::to_vec(value)?)?)
}

fn prefix_case(case: &StageCase, through_stages: usize) -> Result<StageCase> {
    if through_stages == 0 || through_stages > case.stages.len() {
        return Err(Error::Input(
            "epoch SRM gate prefix is outside the construction".into(),
        ));
    }
    if through_stages == case.stages.len() {
        return Ok(case.clone());
    }
    let mut prefix = case.clone();
    prefix.stages.truncate(through_stages);
    prefix
        .verification_reject_once_stages
        .retain(|name| prefix.stages.iter().any(|stage| stage.name() == name));
    Ok(prefix)
}

fn epoch_spec(template: &SearchSpecification, through_stages: usize) -> SearchSpecification {
    let mut spec = template.clone();
    spec.search.through_stages = through_stages;
    spec
}

/// Fs_req authorizes later stages only when that factor attained equilibrium and
/// the search itself reached a finished non-contradictory reason. An independent
/// exclusion of Fs_req stops construction. A trial limit or nonconvergence does
/// not authorize continuation and is not relabeled as exclusion or collapse.
pub fn effect(report: &SearchReport) -> GateEffect {
    match report.requirement {
        Requirement::Indeterminate => GateEffect::StopIndeterminate,
        Requirement::IndependentlyExcluded => GateEffect::StopExcluded,
        Requirement::EquilibriumRequirementMet => match report.stop_reason {
            StopReason::FactorTolerance | StopReason::NoPhysicalUpper => GateEffect::Continue,
            StopReason::TrialLimit
            | StopReason::IndeterminateTrial
            | StopReason::ContradictoryEvidence => GateEffect::StopIndeterminate,
        },
    }
}

fn validate_template(
    mesh: &Mesh,
    ground: &Ground,
    case: &StageCase,
    spec: &Specification,
) -> Result<()> {
    let stages_ok = !spec.gate_after_stages.is_empty()
        && spec.gate_after_stages.len() <= case.stages.len()
        && spec
            .gate_after_stages
            .windows(2)
            .all(|pair| pair[0] < pair[1])
        && spec
            .gate_after_stages
            .iter()
            .all(|stage| (1..=case.stages.len()).contains(stage));
    let path = &spec.search;
    let controls_ok = path.initial_increment.is_finite()
        && path.initial_increment > 0.
        && path.initial_increment <= 1.
        && path.minimum_increment.is_finite()
        && path.minimum_increment >= 1e-8
        && path.minimum_increment <= path.initial_increment
        && (1..=10000).contains(&path.maximum_attempts)
        && path
            .search
            .virtual_displacements_m
            .iter()
            .all(|field| field.len() == mesh.nodes.len());
    if spec.schema != "soil-fem-epoch-srm-stage-gate/1"
        || path.schema != "soil-fem-epoch-srm-search/1"
        || path.search.through_stages != 0
        || !stages_ok
        || !controls_ok
    {
        return Err(Error::Input(
            "invalid epoch SRM stage-gate schema, stage list or path controls".into(),
        ));
    }
    let mut probe = spec.search.search.clone();
    probe.through_stages = *spec
        .gate_after_stages
        .last()
        .ok_or_else(|| Error::Input("epoch SRM stage gate has no stage".into()))?;
    probe.validate(ground, case)
}

fn accepted_checkpoint(run: &Run) -> bool {
    !run.complete && run.failure_reason.as_deref() == Some("checkpoint_requested")
}

fn construction_step(
    mesh: &Mesh,
    ground: &Ground,
    case: StageCase,
    checkpoint: Option<Run>,
    observer: Option<&dyn Observer>,
) -> Result<Run> {
    let ground = ground.clone();
    let run = match (checkpoint, observer) {
        (None, None) => stages::solve_until(mesh, ground, case, 1)?,
        (None, Some(observer)) => stages::solve_until_observed(mesh, ground, case, 1, observer)?,
        (Some(saved), None) => stages::resume(mesh, ground, case, saved, Some(1))?,
        (Some(saved), Some(observer)) => {
            stages::resume_limited_observed(mesh, ground, case, saved, Some(1), observer)?
        }
    };
    if run.complete || accepted_checkpoint(&run) {
        Ok(run)
    } else {
        Err(Error::Solve(run.failure_reason.unwrap_or_else(|| {
            "unreduced construction increment failed; no Fs gate was evaluated".into()
        })))
    }
}

impl Snapshot {
    pub fn new(
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
    ) -> Result<Self> {
        mesh.validate()?;
        ground.validate()?;
        case.validate()?;
        validate_template(mesh, ground, case, spec)?;
        let through_stages = spec.gate_after_stages[0];
        let mut value = Self {
            schema: "soil-fem-epoch-srm-stage-gate-checkpoint/1".into(),
            identity: Identity::current(mesh, ground, case)?,
            case_stage_count: case.stages.len(),
            specification_sha256: checkpoint::hash_json(spec)?,
            payload_sha256: String::new(),
            payload: Payload {
                phase: Phase::Construct(Box::new(ConstructPhase {
                    through_stages,
                    checkpoint: None,
                })),
                records: Vec::new(),
                construction: None,
                report: None,
            },
        };
        value.seal()?;
        Ok(value)
    }

    fn seal(&mut self) -> Result<()> {
        self.payload_sha256 = checkpoint::hash_json(&self.payload)?;
        Ok(())
    }

    pub fn is_complete(&self) -> bool {
        matches!(self.payload.phase, Phase::Done)
    }

    pub fn report(&self) -> Option<&Report> {
        self.payload.report.as_ref()
    }

    pub fn phase_label(&self) -> &'static str {
        match self.payload.phase {
            Phase::Construct(_) => "construct",
            Phase::Search(_) => "search",
            Phase::Finish(_) => "finish",
            Phase::Done => "done",
        }
    }

    pub fn through_stages(&self) -> Option<usize> {
        match &self.payload.phase {
            Phase::Construct(phase) => Some(phase.through_stages),
            Phase::Search(phase) => Some(phase.through_stages),
            Phase::Finish(_) | Phase::Done => None,
        }
    }

    pub fn completed_gates(&self) -> usize {
        self.payload.records.len()
    }

    pub fn current_factor(&self, spec: &Specification) -> Option<f64> {
        let Phase::Search(phase) = &self.payload.phase else {
            return None;
        };
        phase
            .search
            .current_factor(&epoch_spec(&spec.search, phase.through_stages))
    }

    pub fn decision(&self) -> Option<Decision> {
        self.report().map(|report| report.decision)
    }

    pub fn latest_requirement(&self) -> Option<Requirement> {
        self.report()
            .and_then(|report| report.gates.last().map(|gate| gate.requirement))
            .or_else(|| {
                self.payload
                    .records
                    .last()
                    .and_then(|gate| gate.search.report().map(|report| report.requirement))
            })
    }

    pub fn latest_stop_reason(&self) -> Option<StopReason> {
        self.report()
            .and_then(|report| report.gates.last().map(|gate| gate.stop_reason))
            .or_else(|| {
                self.payload
                    .records
                    .last()
                    .and_then(|gate| gate.search.report().map(|report| report.stop_reason))
            })
    }

    fn assemble(&self, spec: &Specification) -> Result<Report> {
        let gates =
            self.payload
                .records
                .iter()
                .map(|gate| {
                    let search = gate.search.report().ok_or_else(|| {
                        Error::Control("epoch SRM gate record is unfinished".into())
                    })?;
                    Ok(GateRecord {
                        through_stages: gate.through_stages,
                        requirement: search.requirement,
                        stop_reason: search.stop_reason,
                        effect: effect(search),
                        search: duplicate(search)?,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
        let construction =
            self.payload.construction.as_ref().ok_or_else(|| {
                Error::Control("epoch SRM gate has no construction result".into())
            })?;
        let decision = match gates.last().map(|gate| gate.effect) {
            Some(GateEffect::StopExcluded) => Decision::StoppedExcluded,
            Some(GateEffect::StopIndeterminate) => Decision::StoppedIndeterminate,
            Some(GateEffect::Continue) => Decision::Passed,
            None => {
                return Err(Error::Control(
                    "epoch SRM gate finished without a stage decision".into(),
                ));
            }
        };
        if decision == Decision::Passed
            && gates.iter().any(|gate| gate.effect != GateEffect::Continue)
        {
            return Err(Error::Control(
                "passed epoch SRM gate contains a stopping stage".into(),
            ));
        }
        Ok(Report {
            schema: "soil-fem-epoch-srm-stage-gate-result/1".into(),
            source_identity: self.identity.clone(),
            specification: spec.clone(),
            gates,
            decision,
            construction: duplicate(construction)?,
            construction_stage_count: construction.case.stages.len(),
            stopped_before_later_stages: decision != Decision::Passed
                && construction.case.stages.len() < self.case_stage_count,
            nonconvergence_used_as_exclusion_or_collapse: false,
            continuum_stability_or_design_safety_certified: false,
        })
    }

    pub fn verify(
        &self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
    ) -> Result<()> {
        mesh.validate()?;
        ground.validate()?;
        case.validate()?;
        validate_template(mesh, ground, case, spec)?;
        if self.schema != "soil-fem-epoch-srm-stage-gate-checkpoint/1"
            || self.identity != Identity::current(mesh, ground, case)?
            || self.case_stage_count != case.stages.len()
            || self.specification_sha256 != checkpoint::hash_json(spec)?
            || self.payload_sha256 != checkpoint::hash_json(&self.payload)?
        {
            return Err(Error::Input(
                "epoch SRM gate snapshot identity/specification/seal mismatch".into(),
            ));
        }
        self.verify_records(mesh, spec)?;
        match &self.payload.phase {
            Phase::Done => self.verify_done(mesh, ground, case, spec),
            Phase::Construct(phase) => {
                self.verify_open(spec, phase.through_stages)?;
                if let Some(saved) = &phase.checkpoint {
                    self.verify_checkpoint(
                        mesh,
                        ground,
                        &prefix_case(case, phase.through_stages)?,
                        saved,
                    )?;
                }
                Ok(())
            }
            Phase::Search(phase) => {
                self.verify_open(spec, phase.through_stages)?;
                self.verify_origin(mesh, ground, case, phase.through_stages, &phase.origin)?;
                phase.search.verify(
                    mesh,
                    &phase.origin,
                    &epoch_spec(&spec.search, phase.through_stages),
                )?;
                if phase.search.is_complete() {
                    return Err(Error::Input(
                        "completed epoch search remains in the active gate phase".into(),
                    ));
                }
                Ok(())
            }
            Phase::Finish(phase) => {
                if self.payload.records.len() != spec.gate_after_stages.len()
                    || self.payload.records.iter().any(|gate| {
                        gate.search
                            .report()
                            .is_none_or(|report| effect(report) != GateEffect::Continue)
                    })
                    || self
                        .payload
                        .records
                        .last()
                        .is_none_or(|gate| gate.through_stages >= case.stages.len())
                {
                    return Err(Error::Input(
                        "epoch SRM finish phase does not follow passed intermediate gates".into(),
                    ));
                }
                if let Some(saved) = &phase.checkpoint {
                    self.verify_checkpoint(mesh, ground, case, saved)?;
                }
                Ok(())
            }
        }
    }

    fn verify_open(&self, spec: &Specification, through_stages: usize) -> Result<()> {
        if self.payload.report.is_some() || self.payload.construction.is_some() {
            return Err(Error::Input(
                "unfinished epoch SRM gate already stores a final result".into(),
            ));
        }
        if spec
            .gate_after_stages
            .get(self.payload.records.len())
            .copied()
            != Some(through_stages)
        {
            return Err(Error::Input("epoch SRM gate stage sequence changed".into()));
        }
        if self.payload.records.last().is_some_and(|gate| {
            gate.search
                .report()
                .is_some_and(|report| effect(report) != GateEffect::Continue)
        }) {
            return Err(Error::Input(
                "epoch SRM gate continued after a stopping decision".into(),
            ));
        }
        Ok(())
    }

    fn verify_records(&self, mesh: &Mesh, spec: &Specification) -> Result<()> {
        if self.payload.records.len() > spec.gate_after_stages.len() {
            return Err(Error::Input("epoch SRM gate has unexpected stages".into()));
        }
        for (index, gate) in self.payload.records.iter().enumerate() {
            if gate.through_stages != spec.gate_after_stages[index] {
                return Err(Error::Input("epoch SRM gate stage sequence changed".into()));
            }
            gate.search.verify(
                mesh,
                &gate.origin,
                &epoch_spec(&spec.search, gate.through_stages),
            )?;
            let report = gate.search.report().ok_or_else(|| {
                Error::Input("epoch SRM gate record omitted its search report".into())
            })?;
            if !gate.search.is_complete()
                || gate.origin.case.stages.len() != gate.through_stages
                || checkpoint::hash_json(&gate.origin)? != report.origin_run_sha256
            {
                return Err(Error::Input(
                    "epoch SRM gate record origin/report binding changed".into(),
                ));
            }
            if effect(report) != GateEffect::Continue && index + 1 != self.payload.records.len() {
                return Err(Error::Input(
                    "epoch SRM gate continued after a stopping decision".into(),
                ));
            }
        }
        Ok(())
    }

    fn verify_origin(
        &self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        through_stages: usize,
        origin: &Run,
    ) -> Result<()> {
        let prefix = prefix_case(case, through_stages)?;
        if !origin.complete
            || origin.failure_reason.is_some()
            || origin.provenance != Identity::current(mesh, ground, &prefix)?
            || checkpoint::hash_json(&origin.case)? != checkpoint::hash_json(&prefix)?
            || checkpoint::hash_json(&origin.ground)? != checkpoint::hash_json(ground)?
        {
            return Err(Error::Input(
                "epoch SRM gate origin is not the completed unreduced prefix".into(),
            ));
        }
        Ok(())
    }

    fn verify_checkpoint(
        &self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        saved: &Run,
    ) -> Result<()> {
        if !accepted_checkpoint(saved)
            || saved.provenance != Identity::current(mesh, ground, case)?
            || checkpoint::hash_json(&saved.case)? != checkpoint::hash_json(case)?
            || checkpoint::hash_json(&saved.ground)? != checkpoint::hash_json(ground)?
        {
            return Err(Error::Input(
                "epoch SRM gate construction checkpoint binding changed".into(),
            ));
        }
        Ok(())
    }

    fn verify_done(
        &self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
    ) -> Result<()> {
        let report = self
            .payload
            .report
            .as_ref()
            .ok_or_else(|| Error::Input("finished epoch SRM gate has no report".into()))?;
        let construction =
            self.payload.construction.as_ref().ok_or_else(|| {
                Error::Input("finished epoch SRM gate has no construction".into())
            })?;
        let rebuilt = self.assemble(spec)?;
        if checkpoint::hash_json(report)? != checkpoint::hash_json(&rebuilt)?
            || checkpoint::hash_json(&report.construction)? != checkpoint::hash_json(construction)?
            || report.nonconvergence_used_as_exclusion_or_collapse
            || report.continuum_stability_or_design_safety_certified
        {
            return Err(Error::Input(
                "epoch SRM gate report does not match its saved decisions".into(),
            ));
        }
        let last =
            self.payload.records.last().ok_or_else(|| {
                Error::Input("finished epoch SRM gate has no stage decision".into())
            })?;
        let last_effect =
            effect(last.search.report().ok_or_else(|| {
                Error::Input("finished epoch SRM gate record has no search".into())
            })?);
        match report.decision {
            Decision::Passed => {
                if last_effect != GateEffect::Continue
                    || self.payload.records.len() != spec.gate_after_stages.len()
                    || !construction.complete
                    || construction.failure_reason.is_some()
                    || construction.case.stages.len() != case.stages.len()
                    || construction.provenance != Identity::current(mesh, ground, case)?
                    || report.stopped_before_later_stages
                {
                    return Err(Error::Input(
                        "passed epoch SRM gate construction is not the full unreduced run".into(),
                    ));
                }
            }
            Decision::StoppedExcluded | Decision::StoppedIndeterminate => {
                let expected = if report.decision == Decision::StoppedExcluded {
                    GateEffect::StopExcluded
                } else {
                    GateEffect::StopIndeterminate
                };
                if last_effect != expected
                    || !construction.complete
                    || construction.failure_reason.is_some()
                    || construction.case.stages.len() != last.through_stages
                    || checkpoint::hash_json(construction)? != checkpoint::hash_json(&last.origin)?
                    || report.stopped_before_later_stages
                        != (last.through_stages < case.stages.len())
                {
                    return Err(Error::Input(
                        "stopped epoch SRM gate does not retain its unreduced prefix".into(),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Advance by accepted construction increments, or by one material-attempt
    /// boundary inside the active strength trial. A rejected material trial does
    /// not replace the accepted state. On error the caller keeps its previous snapshot.
    pub fn advance(
        mut self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
        units: usize,
        observer: Option<&dyn Observer>,
    ) -> Result<Self> {
        self.verify(mesh, ground, case, spec)?;
        if units == 0 || units > 10000 || self.is_complete() {
            return Err(Error::Input(
                "epoch SRM gate advance needs 1..10000 units and an unfinished snapshot".into(),
            ));
        }
        for _ in 0..units {
            if self.is_complete()
                || observer
                    .map(|item| item.stop_requested())
                    .transpose()?
                    .unwrap_or(false)
            {
                break;
            }
            self.one_unit(mesh, ground, case, spec, observer)?;
        }
        self.seal()?;
        Ok(self)
    }

    fn one_unit(
        &mut self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
        observer: Option<&dyn Observer>,
    ) -> Result<()> {
        let phase = std::mem::replace(&mut self.payload.phase, Phase::Done);
        let phase = match phase {
            Phase::Construct(phase) => {
                self.construct_unit(mesh, ground, case, spec, *phase, observer)?
            }
            Phase::Search(phase) => self.search_unit(mesh, spec, *phase, observer)?,
            Phase::Finish(phase) => self.finish_unit(mesh, ground, case, spec, *phase, observer)?,
            Phase::Done => {
                return Err(Error::Input(
                    "epoch SRM gate advance continued after its decision".into(),
                ));
            }
        };
        self.payload.phase = phase;
        Ok(())
    }

    fn construct_unit(
        &mut self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
        phase: ConstructPhase,
        observer: Option<&dyn Observer>,
    ) -> Result<Phase> {
        let prefix = prefix_case(case, phase.through_stages)?;
        let run = construction_step(mesh, ground, prefix, phase.checkpoint, observer)?;
        if !run.complete {
            return Ok(Phase::Construct(Box::new(ConstructPhase {
                through_stages: phase.through_stages,
                checkpoint: Some(run),
            })));
        }
        let search =
            SearchSnapshot::new(mesh, &run, &epoch_spec(&spec.search, phase.through_stages))?;
        Ok(Phase::Search(Box::new(SearchPhase {
            through_stages: phase.through_stages,
            origin: run,
            search,
        })))
    }

    fn search_unit(
        &mut self,
        mesh: &Mesh,
        spec: &Specification,
        phase: SearchPhase,
        observer: Option<&dyn Observer>,
    ) -> Result<Phase> {
        let SearchPhase {
            through_stages,
            origin,
            search,
        } = phase;
        let search = search.advance(
            mesh,
            &origin,
            &epoch_spec(&spec.search, through_stages),
            1,
            observer,
        )?;
        if !search.is_complete() {
            return Ok(Phase::Search(Box::new(SearchPhase {
                through_stages,
                origin,
                search,
            })));
        }
        let gate_effect =
            effect(search.report().ok_or_else(|| {
                Error::Control("completed epoch search omitted its report".into())
            })?);
        self.payload.records.push(CompletedGate {
            through_stages,
            origin: duplicate(&origin)?,
            search,
        });
        match gate_effect {
            GateEffect::Continue => {
                if let Some(next) = spec
                    .gate_after_stages
                    .iter()
                    .copied()
                    .find(|stage| *stage > through_stages)
                {
                    Ok(Phase::Construct(Box::new(ConstructPhase {
                        through_stages: next,
                        checkpoint: None,
                    })))
                } else if through_stages < self.case_stage_count {
                    Ok(Phase::Finish(Box::new(FinishPhase { checkpoint: None })))
                } else {
                    self.finish_with(spec, origin)
                }
            }
            GateEffect::StopExcluded | GateEffect::StopIndeterminate => {
                self.finish_with(spec, origin)
            }
        }
    }

    fn finish_unit(
        &mut self,
        mesh: &Mesh,
        ground: &Ground,
        case: &StageCase,
        spec: &Specification,
        phase: FinishPhase,
        observer: Option<&dyn Observer>,
    ) -> Result<Phase> {
        let run = construction_step(mesh, ground, case.clone(), phase.checkpoint, observer)?;
        if !run.complete {
            return Ok(Phase::Finish(Box::new(FinishPhase {
                checkpoint: Some(run),
            })));
        }
        self.finish_with(spec, run)
    }

    fn finish_with(&mut self, spec: &Specification, construction: Run) -> Result<Phase> {
        self.payload.construction = Some(construction);
        self.payload.report = Some(self.assemble(spec)?);
        Ok(Phase::Done)
    }
}

pub fn run(mesh: &Mesh, ground: &Ground, case: &StageCase, spec: &Specification) -> Result<Report> {
    let mut snapshot = Snapshot::new(mesh, ground, case, spec)?;
    for _ in 0..2_000 {
        if snapshot.is_complete() {
            return duplicate(snapshot.report().ok_or_else(|| {
                Error::Control("completed epoch SRM gate omitted its report".into())
            })?);
        }
        snapshot = snapshot.advance(mesh, ground, case, spec, 64, None)?;
    }
    Err(Error::Control(
        "epoch SRM gate exceeded its accepted-unit limit".into(),
    ))
}
