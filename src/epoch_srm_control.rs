//! Accepted material-attempt boundaries for an immutable-origin SRM search.
//! No unfinished Newton trial is serialized and no factor reuses another factor.
use crate::{
    Error, Result,
    checkpoint::{self, Identity},
    epoch_srm_search::{self as search, EpochTrial, Report, Specification},
    mc_epoch_path::{self, Controller},
    mesh::Mesh,
    srm_search::{self, EndpointAudit, Interval, Outcome, Requirement, StopReason, Trial},
    stage_control::{Observer, Progress},
    stages::Run,
};
use serde::{Deserialize, Serialize};
use std::cell::RefCell;

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Payload {
    basis: Option<EndpointAudit>,
    trials: Vec<EpochTrial>,
    pending: Option<Trial>,
    inflight: Option<mc_epoch_path::Checkpoint>,
    report: Option<Report>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    schema: String,
    identity: Identity,
    origin_run_sha256: String,
    origin_state_sha256: String,
    specification_sha256: String,
    payload_sha256: String,
    payload: Payload,
}

// Forward progress errors as terminal control failures. A stop request waits
// for the current material attempt, then preserves the last accepted epoch.
struct ProgressOnly<'a> {
    external: Option<&'a dyn Observer>,
    failed: RefCell<Option<String>>,
}
impl Observer for ProgressOnly<'_> {
    fn on_progress(&self, event: &Progress) -> Result<()> {
        if let Some(observer) = self.external
            && let Err(error) = observer.on_progress(event)
        {
            *self.failed.borrow_mut() = Some(error.to_string());
            return Err(error);
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(false)
    }
}

enum Action {
    Trial(f64),
    Finish(StopReason),
}
fn next(spec: &Specification, trials: &[EpochTrial]) -> Action {
    for factor in [1., spec.search.required_factor, spec.search.maximum_factor] {
        if !trials.iter().any(|t| t.evidence.factor == factor) {
            return Action::Trial(factor);
        }
    }
    if trials
        .iter()
        .any(|t| t.evidence.failure_kind.as_deref() == Some("contradictory_evidence"))
    {
        return Action::Finish(StopReason::ContradictoryEvidence);
    }
    if trials
        .iter()
        .any(|t| t.evidence.outcome == Outcome::Indeterminate)
    {
        return Action::Finish(StopReason::IndeterminateTrial);
    }
    let Some(bounds) = search::interval(trials) else {
        return Action::Finish(StopReason::NoPhysicalUpper);
    };
    if bounds.width <= spec.search.absolute_factor_tolerance {
        return Action::Finish(StopReason::FactorTolerance);
    }
    if trials.len() >= spec.search.maximum_trials {
        return Action::Finish(StopReason::TrialLimit);
    }
    Action::Trial(bounds.attained_equilibrium_factor + 0.5 * bounds.width)
}

fn decision(
    spec: &Specification,
    trials: &[EpochTrial],
) -> Result<(StopReason, Requirement, Option<Interval>)> {
    let Action::Finish(mut reason) = next(spec, trials) else {
        return Err(Error::Input(
            "epoch SRM report was finalized before its next trial".into(),
        ));
    };
    if trials.iter().any(|a| {
        a.evidence.outcome == Outcome::EquilibriumAttained
            && trials.iter().any(|b| {
                b.evidence.outcome == Outcome::PhysicallyExcluded
                    && a.evidence.factor >= b.evidence.factor
            })
    }) {
        reason = StopReason::ContradictoryEvidence;
    }
    let requirement = if reason == StopReason::ContradictoryEvidence {
        Requirement::Indeterminate
    } else {
        match trials
            .iter()
            .find(|t| t.evidence.factor == spec.search.required_factor)
            .map(|t| t.evidence.outcome)
        {
            Some(Outcome::EquilibriumAttained) => Requirement::EquilibriumRequirementMet,
            Some(Outcome::PhysicallyExcluded) => Requirement::IndependentlyExcluded,
            _ => Requirement::Indeterminate,
        }
    };
    let bounds = if reason == StopReason::ContradictoryEvidence {
        None
    } else {
        search::interval(trials)
    };
    Ok((reason, requirement, bounds))
}

fn finish_trial(
    mesh: &Mesh,
    origin: &Run,
    mut evidence: Trial,
    saved: &mc_epoch_path::Checkpoint,
) -> Result<Trial> {
    let checked = Controller::resume(mesh, origin, saved)?;
    if checked.complete() {
        evidence.endpoint_audit = Some(search::audit_epoch_endpoint(
            mesh,
            origin,
            &saved.payload.material_checkpoint,
        )?);
        evidence.outcome = Outcome::EquilibriumAttained;
    } else {
        evidence.failure_kind = Some("epoch_path_incomplete".into());
        evidence.reason = Some(format!(
            "accepted parameter progress {}: {:?}",
            saved.payload.continuation.progress, saved.payload.continuation.stop_reason
        ));
    }
    Ok(evidence)
}

fn verify_trial(
    mesh: &Mesh,
    origin: &Run,
    spec: &Specification,
    basis: &EndpointAudit,
    trial: &EpochTrial,
) -> Result<()> {
    let (expected, needs_solve) = srm_search::prepare_trial(
        mesh,
        &origin.ground,
        &spec.search,
        basis,
        trial.evidence.factor,
    )?;
    let expected = if needs_solve {
        let saved = trial
            .epoch
            .as_ref()
            .ok_or_else(|| Error::Input("epoch trial omitted its material checkpoint".into()))?;
        if checkpoint::hash_json(&saved.payload.specification)?
            != checkpoint::hash_json(&spec.path(origin, trial.evidence.factor)?)?
        {
            return Err(Error::Input(
                "saved epoch belongs to a different strength factor/path".into(),
            ));
        }
        finish_trial(mesh, origin, expected, saved)?
    } else {
        if trial.epoch.is_some() {
            return Err(Error::Input(
                "physically excluded trial has an unexpected material epoch".into(),
            ));
        }
        expected
    };
    if checkpoint::hash_json(&expected)? != checkpoint::hash_json(&trial.evidence)? {
        return Err(Error::Input(
            "epoch trial evidence differs from independently reconstructed endpoint/physical work"
                .into(),
        ));
    }
    Ok(())
}

fn verify_trials(
    mesh: &Mesh,
    origin: &Run,
    spec: &Specification,
    basis: &EndpointAudit,
    trials: &[EpochTrial],
) -> Result<()> {
    for (index, trial) in trials.iter().enumerate() {
        let Action::Trial(factor) = next(spec, &trials[..index]) else {
            return Err(Error::Input(
                "epoch SRM has trials after its terminal decision".into(),
            ));
        };
        if trial.evidence.factor != factor {
            return Err(Error::Input(
                "epoch SRM trial sequence/factor changed".into(),
            ));
        }
        verify_trial(mesh, origin, spec, basis, trial)?;
    }
    Ok(())
}

impl Snapshot {
    pub fn new(mesh: &Mesh, origin: &Run, spec: &Specification) -> Result<Self> {
        mesh.validate()?;
        spec.validate(mesh, origin)?;
        let basis = srm_search::audit_unreduced_endpoint(mesh, origin)?;
        let mut value = Self {
            schema: "soil-fem-epoch-srm-checkpoint/1".into(),
            identity: origin.provenance.clone(),
            origin_run_sha256: checkpoint::hash_json(origin)?,
            origin_state_sha256: checkpoint::hash_json(&origin.last_committed)?,
            specification_sha256: checkpoint::hash_json(spec)?,
            payload_sha256: String::new(),
            payload: Payload {
                basis: Some(basis),
                trials: vec![],
                pending: None,
                inflight: None,
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
        self.payload.report.is_some()
    }
    pub fn report(&self) -> Option<&Report> {
        self.payload.report.as_ref()
    }
    pub fn trials(&self) -> &[EpochTrial] {
        self.payload
            .report
            .as_ref()
            .map_or(&self.payload.trials, |r| &r.trials)
    }
    pub fn current_factor(&self, spec: &Specification) -> Option<f64> {
        if self.is_complete() {
            return None;
        }
        if let Some(pending) = &self.payload.pending {
            return Some(pending.factor);
        }
        match next(spec, &self.payload.trials) {
            Action::Trial(f) => Some(f),
            Action::Finish(_) => None,
        }
    }
    pub fn accepted_parameter_progress(&self) -> Option<f64> {
        self.payload
            .inflight
            .as_ref()
            .map(|s| s.payload.continuation.progress)
    }
    pub fn has_inflight_checkpoint(&self) -> bool {
        self.payload.inflight.is_some()
    }

    pub fn verify(&self, mesh: &Mesh, origin: &Run, spec: &Specification) -> Result<()> {
        mesh.validate()?;
        spec.validate(mesh, origin)?;
        if self.schema != "soil-fem-epoch-srm-checkpoint/1"
            || self.identity != origin.provenance
            || self.origin_run_sha256 != checkpoint::hash_json(origin)?
            || self.origin_state_sha256 != checkpoint::hash_json(&origin.last_committed)?
            || self.specification_sha256 != checkpoint::hash_json(spec)?
            || self.payload_sha256 != checkpoint::hash_json(&self.payload)?
        {
            return Err(Error::Input(
                "epoch SRM snapshot origin/input/solver/specification/seal mismatch".into(),
            ));
        }
        let expected_basis = srm_search::audit_unreduced_endpoint(mesh, origin)?;
        let p = &self.payload;
        if let Some(report) = &p.report {
            if p.basis.is_some()
                || !p.trials.is_empty()
                || p.pending.is_some()
                || p.inflight.is_some()
            {
                return Err(Error::Input(
                    "finished epoch SRM has unfinished phase data".into(),
                ));
            }
            verify_trials(mesh, origin, spec, &expected_basis, &report.trials)?;
            let (reason, requirement, bounds) = decision(spec, &report.trials)?;
            if report.schema != "soil-fem-epoch-srm-search-result/1"
                || report.source_identity != self.identity
                || report.origin_run_sha256 != self.origin_run_sha256
                || report.origin_state_sha256 != self.origin_state_sha256
                || checkpoint::hash_json(&report.specification)? != self.specification_sha256
                || checkpoint::hash_json(&report.original_unreduced_run)? != self.origin_run_sha256
                || checkpoint::hash_json(&report.original_endpoint_audit)?
                    != checkpoint::hash_json(&expected_basis)?
                || report.stop_reason != reason
                || report.requirement != requirement
                || checkpoint::hash_json(&report.interval)? != checkpoint::hash_json(&bounds)?
                || !report.unreduced_source_history_unchanged
                || report.continuum_stability_or_design_safety_certified
                || report.trial_origin
                    != "same_completed_unreduced_construction_Run_for_every_independent_factor"
            {
                return Err(Error::Input(
                    "epoch SRM final report/evidence/requirement mismatch".into(),
                ));
            }
            return Ok(());
        }
        let basis = p
            .basis
            .as_ref()
            .ok_or_else(|| Error::Input("epoch SRM origin audit missing".into()))?;
        if checkpoint::hash_json(basis)? != checkpoint::hash_json(&expected_basis)? {
            return Err(Error::Input(
                "epoch SRM saved origin endpoint audit changed".into(),
            ));
        }
        verify_trials(mesh, origin, spec, basis, &p.trials)?;
        match (&p.pending, &p.inflight) {
            (None, None) => {}
            (Some(pending), Some(saved)) => {
                let Action::Trial(factor) = next(spec, &p.trials) else {
                    return Err(Error::Input(
                        "epoch SRM pending trial after terminal decision".into(),
                    ));
                };
                let (expected, needs_solve) =
                    srm_search::prepare_trial(mesh, &origin.ground, &spec.search, basis, factor)?;
                if pending.factor != factor
                    || !needs_solve
                    || checkpoint::hash_json(pending)? != checkpoint::hash_json(&expected)?
                    || checkpoint::hash_json(&saved.payload.specification)?
                        != checkpoint::hash_json(&spec.path(origin, factor)?)?
                {
                    return Err(Error::Input(
                        "epoch SRM pending factor/path binding changed".into(),
                    ));
                }
                let checked = Controller::resume(mesh, origin, saved)?;
                if checked.complete() || checked.continuation().stop_reason.is_some() {
                    return Err(Error::Input(
                        "terminal material path left in pending epoch SRM slot".into(),
                    ));
                }
            }
            _ => {
                return Err(Error::Input(
                    "epoch SRM pending/checkpoint inventory mismatch".into(),
                ));
            }
        }
        Ok(())
    }

    /// At most material_attempt_budget attempts are run. Rejected attempts may
    /// update the next increment but never replace the accepted state/parameters.
    /// On error the caller must keep its last immutable saved snapshot.
    pub fn advance(
        mut self,
        mesh: &Mesh,
        origin: &Run,
        spec: &Specification,
        material_attempt_budget: usize,
        observer: Option<&dyn Observer>,
    ) -> Result<Self> {
        self.verify(mesh, origin, spec)?;
        if material_attempt_budget == 0 || material_attempt_budget > 10000 || self.is_complete() {
            return Err(Error::Input(
                "epoch SRM advance needs1..10000 attempts and an unfinished snapshot".into(),
            ));
        }
        if observer
            .map(|o| o.stop_requested())
            .transpose()?
            .unwrap_or(false)
        {
            return Ok(self);
        }
        if self.payload.pending.is_none() {
            match next(spec, &self.payload.trials) {
                Action::Finish(_) => {
                    let (reason, requirement, bounds) = decision(spec, &self.payload.trials)?;
                    let clone: Run = serde_json::from_slice(&serde_json::to_vec(origin)?)?;
                    self.payload.report = Some(Report {
                        schema: "soil-fem-epoch-srm-search-result/1".into(),
                        source_identity: self.identity.clone(),
                        origin_run_sha256: self.origin_run_sha256.clone(),
                        origin_state_sha256: self.origin_state_sha256.clone(),
                        specification: spec.clone(),
                        original_unreduced_run: clone,
                        original_endpoint_audit: self.payload.basis.take().unwrap(),
                        trials: std::mem::take(&mut self.payload.trials),
                        interval: bounds,
                        requirement,
                        stop_reason: reason,
                        trial_origin:
                            "same_completed_unreduced_construction_Run_for_every_independent_factor"
                                .into(),
                        unreduced_source_history_unchanged: true,
                        continuum_stability_or_design_safety_certified: false,
                    });
                    self.seal()?;
                    return Ok(self);
                }
                Action::Trial(factor) => {
                    let (trial, needs_solve) = srm_search::prepare_trial(
                        mesh,
                        &origin.ground,
                        &spec.search,
                        self.payload.basis.as_ref().unwrap(),
                        factor,
                    )?;
                    if !needs_solve {
                        self.payload.trials.push(EpochTrial {
                            evidence: trial,
                            epoch: None,
                        });
                        self.seal()?;
                        return Ok(self);
                    }
                    self.payload.pending = Some(trial);
                }
            }
        }
        let trial = self.payload.pending.take().unwrap();
        let progress = ProgressOnly {
            external: observer,
            failed: RefCell::new(None),
        };
        let mut controller = match self.payload.inflight.take() {
            Some(saved) => Controller::resume(mesh, origin, &saved)?,
            None => Controller::new(mesh, origin, spec.path(origin, trial.factor)?)?,
        };
        controller.set_observer(&progress);
        for _ in 0..material_attempt_budget {
            if controller.complete()
                || controller.continuation().stop_reason.is_some()
                || observer
                    .map(|o| o.stop_requested())
                    .transpose()?
                    .unwrap_or(false)
            {
                break;
            }
            controller.advance()?;
            if let Some(error) = progress.failed.borrow().as_ref() {
                return Err(Error::Control(format!(
                    "epoch SRM observer failed; retain previous snapshot: {error}"
                )));
            }
        }
        let saved = controller.checkpoint()?;
        if controller.complete() || controller.continuation().stop_reason.is_some() {
            let evidence = finish_trial(mesh, origin, trial, &saved)?;
            self.payload.trials.push(EpochTrial {
                evidence,
                epoch: Some(saved),
            });
        } else {
            self.payload.pending = Some(trial);
            self.payload.inflight = Some(saved);
        }
        self.seal()?;
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ground::Ground,
        soil,
        srm_bound::Problem,
        stages::{self, StageCase},
    };
    use serde_json::json;
    use std::cell::Cell;

    fn inputs() -> Result<(Mesh, Ground, StageCase, Specification)> {
        let mesh: Mesh =
            serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
        let mesh = mesh.elevate_straight()?;
        let p: Problem =
            serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
        let mut raw: serde_json::Value =
            serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
        let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
        raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
        raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
        raw["materials"] = json!([{"id":1,"name":"unconfined epoch SRM","unit_weight_kn_m3":1e-10,"mohr_coulomb":p.material}]);
        raw["bodies"] = json!([{"body":"TestCube","material_id":1}]);
        let ground = serde_json::from_value(raw)?;
        let mut case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
        case.boundary = Some(p.boundary);
        case.stages = serde_json::from_value(
            json!([{"kind":"boundary","name":"load","target_factor":1.,"initial_increment":1.}]),
        )?;
        case.verification_reject_once_stages.clear();
        case.monitors.clear();
        case.monitor_sets.clear();
        let search = serde_json::from_str(include_str!("../examples/srm-search.json"))?;
        let spec = Specification {
            schema: "soil-fem-epoch-srm-search/1".into(),
            search,
            initial_increment: 0.25,
            minimum_increment: 1. / 4096.,
            maximum_attempts: 128,
        };
        Ok((mesh, ground, case, spec))
    }
    fn plastic_origin(cut: bool) -> Result<(Mesh, Run, Specification)> {
        let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
        let mesh = mesh.elevate_straight()?;
        let mut g: Ground =
            serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
        for m in &mut g.materials {
            m.unit_weight_kn_m3 = 1e-10;
            let soil::Model::MohrCoulomb(p) = &mut m.model else {
                panic!("MC")
            };
            p.young_kpa = 100000.;
            p.poisson = 0.3;
            p.cohesion_kpa = 40.;
            p.friction_deg = 0.;
            p.dilation_deg = 0.;
            p.tension_cutoff_kpa = None;
        }
        let mut c: StageCase =
            serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
        let mut steps = vec![
            json!({"kind":"pressure","name":"plastic_preload","target_pressure_kpa":200.}),
            json!({"kind":"pressure","name":"unload","target_pressure_kpa":150.}),
        ];
        if cut {
            steps.push(json!({"kind":"excavate","name":"cut_top","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.0}));
        }
        steps.push(json!({"kind":"reset_display","name":"reference"}));
        c.stages = serde_json::from_value(json!(steps))?;
        c.verification_reject_once_stages.clear();
        let (_, _, _, mut spec) = inputs()?;
        spec.search.material_ids = g.materials.iter().map(|m| m.id).collect();
        spec.search.required_factor = 2.;
        spec.search.maximum_factor = 3.;
        spec.search.through_stages = c.stages.len();
        let origin = stages::solve(&mesh, g, c)?;
        assert!(origin.complete, "{:?}", origin.failure_reason);
        Ok((mesh, origin, spec))
    }

    fn restore(s: &Snapshot) -> Result<Snapshot> {
        let bytes = serde_json::to_vec(s)?;
        let value = serde_json::from_slice(&bytes)?;
        Ok(value)
    }
    fn copy_run(r: &Run) -> Result<Run> {
        Ok(serde_json::from_slice(&serde_json::to_vec(r)?)?)
    }
    fn complete(
        mut s: Snapshot,
        mesh: &Mesh,
        origin: &Run,
        spec: &Specification,
    ) -> Result<Snapshot> {
        for _ in 0..1024 {
            s.verify(mesh, origin, spec)?;
            if s.is_complete() {
                return Ok(s);
            }
            s = restore(&s)?.advance(mesh, origin, spec, 1, None)?;
        }
        Err(Error::Control(
            "tiny epoch SRM test exceeded1024 snapshot steps".into(),
        ))
    }
    fn unit() -> Result<(Mesh, Run, Specification)> {
        let (mesh, g, c, spec) = inputs()?;
        let origin = stages::solve(&mesh, g, c)?;
        Ok((mesh, origin, spec))
    }
    fn at_required(mesh: &Mesh, origin: &Run, spec: &Specification) -> Result<Snapshot> {
        let mut s = Snapshot::new(mesh, origin, spec)?;
        for _ in 0..1024 {
            if s.current_factor(spec) == Some(spec.search.required_factor) {
                return Ok(s);
            }
            s = s.advance(mesh, origin, spec, 1, None)?;
        }
        Err(Error::Control("required factor not reached".into()))
    }
    fn same_report(s: &Snapshot, expected: &Report) -> Result<()> {
        assert_eq!(
            serde_json::to_vec(s.report().unwrap())?,
            serde_json::to_vec(expected)?
        );
        Ok(())
    }

    #[test]
    fn all_attempt_roundtrips_match_the_complete_direct_report_bytes() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let expected = search::run(&mesh, copy_run(&origin)?, spec.clone())?;
        let finished = complete(Snapshot::new(&mesh, &origin, &spec)?, &mesh, &origin, &spec)?;
        same_report(&finished, &expected)?;
        assert_eq!(finished.current_factor(&spec), None);
        assert_eq!(checkpoint::hash_json(&origin)?, finished.origin_run_sha256);
        Ok(())
    }
    struct AlreadyStop;
    impl Observer for AlreadyStop {
        fn on_progress(&self, _: &Progress) -> Result<()> {
            panic!("no mechanical attempt should start");
        }
        fn stop_requested(&self) -> Result<bool> {
            Ok(true)
        }
    }
    #[test]
    fn stop_before_an_attempt_keeps_the_whole_snapshot_byte_exact() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let s = Snapshot::new(&mesh, &origin, &spec)?;
        let before = serde_json::to_vec(&s)?;
        let stopped = s.advance(&mesh, &origin, &spec, 32, Some(&AlreadyStop))?;
        assert_eq!(before, serde_json::to_vec(&stopped)?);
        assert!(!stopped.has_inflight_checkpoint());
        Ok(())
    }
    struct StopFromProgress {
        stop: Cell<bool>,
        calls: Cell<usize>,
    }
    impl Observer for StopFromProgress {
        fn on_progress(&self, _: &Progress) -> Result<()> {
            self.calls.set(self.calls.get() + 1);
            self.stop.set(true);
            Ok(())
        }
        fn stop_requested(&self) -> Result<bool> {
            Ok(self.stop.get())
        }
    }
    #[test]
    fn stop_from_a_nonlinear_callback_finishes_only_the_current_attempt_and_resumes() -> Result<()>
    {
        let (mesh, origin, spec) = plastic_origin(false)?;
        let expected = search::run(&mesh, copy_run(&origin)?, spec.clone())?;
        let ready = at_required(&mesh, &origin, &spec)?;
        let observer = StopFromProgress {
            stop: Cell::new(false),
            calls: Cell::new(0),
        };
        let stopped = ready.advance(&mesh, &origin, &spec, 128, Some(&observer))?;
        assert!(observer.calls.get() > 0);
        let saved = stopped
            .payload
            .inflight
            .as_ref()
            .expect("stop must preserve an unfinished parameter path");
        assert!(
            saved.payload.continuation.progress > 0. && saved.payload.continuation.progress < 1.
        );
        assert!(saved.payload.continuation.stop_reason.is_none());
        assert!(saved.payload.continuation.attempts.last().unwrap().accepted);
        stopped.verify(&mesh, &origin, &spec)?;
        let finished = complete(restore(&stopped)?, &mesh, &origin, &spec)?;
        same_report(&finished, &expected)
    }
    struct BrokenProgress;
    impl Observer for BrokenProgress {
        fn on_progress(&self, _: &Progress) -> Result<()> {
            Err(Error::Control("progress sink failed".into()))
        }
    }
    #[test]
    fn callback_error_is_terminal_control_and_previous_saved_snapshot_remains_usable() -> Result<()>
    {
        let (mesh, origin, spec) = plastic_origin(false)?;
        let ready = at_required(&mesh, &origin, &spec)?;
        let bytes = serde_json::to_vec(&ready)?;
        let outcome = restore(&ready)?.advance(&mesh, &origin, &spec, 128, Some(&BrokenProgress));
        assert!(matches!(outcome, Err(Error::Control(_))));
        assert_eq!(bytes, serde_json::to_vec(&ready)?);
        let expected = search::run(&mesh, copy_run(&origin)?, spec.clone())?;
        let finished = complete(restore(&ready)?, &mesh, &origin, &spec)?;
        same_report(&finished, &expected)
    }
    #[test]
    fn altered_origin_specification_and_producer_identity_are_rejected() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let s = Snapshot::new(&mesh, &origin, &spec)?;
        let mut changed = spec.clone();
        changed.maximum_attempts += 1;
        assert!(s.verify(&mesh, &origin, &changed).is_err());
        let mut changed = copy_run(&origin)?;
        changed.case.boundary.as_mut().unwrap().constraints[0].value_m = 1e-4;
        assert!(s.verify(&mesh, &changed, &spec).is_err());
        let mut changed = restore(&s)?;
        changed.identity.solver_binary_sha256 = "0".repeat(64);
        changed.seal()?;
        assert!(changed.verify(&mesh, &origin, &spec).is_err());
        Ok(())
    }
    #[test]
    fn resealed_pending_factor_and_target_path_cannot_cross_factor_history() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let s = Snapshot::new(&mesh, &origin, &spec)?.advance(&mesh, &origin, &spec, 1, None)?;
        let mut changed = restore(&s)?;
        changed.payload.pending.as_mut().unwrap().factor = 1.1;
        changed.seal()?;
        assert!(changed.verify(&mesh, &origin, &spec).is_err());
        let mut changed = restore(&s)?;
        let saved = changed.payload.inflight.as_mut().unwrap();
        saved.payload.specification.target_parameters[0].cohesion_kpa *= 0.99;
        saved.payload_sha256 = checkpoint::hash_json(&saved.payload)?;
        changed.seal()?;
        assert!(changed.verify(&mesh, &origin, &spec).is_err());
        Ok(())
    }
    #[test]
    fn resealed_frozen_load_metadata_fails_even_when_all_outer_hashes_are_updated() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let mut s =
            Snapshot::new(&mesh, &origin, &spec)?.advance(&mesh, &origin, &spec, 1, None)?;
        let saved = s.payload.inflight.as_mut().unwrap();
        let inner = &mut saved.payload.material_checkpoint;
        inner.payload.committed.applied_nodal_load_kn[0] += 1e-12;
        inner.payload_sha256 = checkpoint::hash_json(&inner.payload)?;
        saved.payload_sha256 = checkpoint::hash_json(&saved.payload)?;
        s.seal()?;
        assert!(s.verify(&mesh, &origin, &spec).is_err());
        Ok(())
    }
    #[test]
    fn resealed_finished_requirement_interval_and_safety_claims_are_recomputed() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let finished = complete(Snapshot::new(&mesh, &origin, &spec)?, &mesh, &origin, &spec)?;
        let mut bad = restore(&finished)?;
        bad.payload.report.as_mut().unwrap().requirement = Requirement::IndependentlyExcluded;
        bad.seal()?;
        assert!(bad.verify(&mesh, &origin, &spec).is_err());
        let mut bad = restore(&finished)?;
        bad.payload
            .report
            .as_mut()
            .unwrap()
            .interval
            .as_mut()
            .unwrap()
            .width = 0.;
        bad.seal()?;
        assert!(bad.verify(&mesh, &origin, &spec).is_err());
        let mut bad = restore(&finished)?;
        bad.payload
            .report
            .as_mut()
            .unwrap()
            .continuum_stability_or_design_safety_certified = true;
        bad.seal()?;
        assert!(bad.verify(&mesh, &origin, &spec).is_err());
        Ok(())
    }
    #[test]
    fn distributed_support_limit_has_real_rejections_accepted_rollback_and_unknown_result()
    -> Result<()> {
        let (mut mesh, g, mut c, mut spec) = inputs()?;
        let base: Vec<_> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, p)| (p[2].abs() < 1e-12).then_some(i))
            .collect();
        assert!(base.len() > 3);
        mesh.node_sets.insert("whole_base".into(), base);
        let boundary = c.boundary.as_mut().unwrap();
        boundary.constraints.retain(|q| q.direction[2] == 0.);
        boundary.constraints.push(crate::boundary::Constraint {
            node_set: "whole_base".into(),
            direction: [0., 0., 1.],
            value_m: 0.,
            gradient: [0.; 3],
            reference_point_m: [0.; 3],
        });
        spec.search.required_factor = 2.;
        spec.search.maximum_factor = 2.;
        spec.minimum_increment = 1. / 64.;
        spec.maximum_attempts = 32;
        let origin = stages::solve(&mesh, g, c)?;
        assert!(origin.complete);
        let mut s = Snapshot::new(&mesh, &origin, &spec)?;
        let mut rejected = 0;
        for _ in 0..1024 {
            if s.is_complete() {
                break;
            }
            let before = s.payload.inflight.as_ref().map(|v| {
                (
                    v.payload.continuation.attempts.len(),
                    checkpoint::hash_json(&v.payload.material_checkpoint.payload.committed)
                        .unwrap(),
                    checkpoint::hash_json(&v.payload.material_checkpoint.payload.parameters)
                        .unwrap(),
                )
            });
            s = restore(&s)?.advance(&mesh, &origin, &spec, 1, None)?;
            let after = s
                .payload
                .inflight
                .as_ref()
                .or_else(|| s.trials().last().and_then(|t| t.epoch.as_ref()));
            if let (Some((count, state, parameters)), Some(cp)) = (before, after)
                && cp.payload.continuation.attempts.len() > count
                && !cp.payload.continuation.attempts.last().unwrap().accepted
            {
                rejected += 1;
                assert_eq!(
                    state,
                    checkpoint::hash_json(&cp.payload.material_checkpoint.payload.committed)?
                );
                assert_eq!(
                    parameters,
                    checkpoint::hash_json(&cp.payload.material_checkpoint.payload.parameters)?
                );
            }
            s.verify(&mesh, &origin, &spec)?;
        }
        assert!(s.is_complete() && rejected > 0);
        let report = s.report().unwrap();
        assert_eq!(report.requirement, Requirement::Indeterminate);
        assert_eq!(report.stop_reason, StopReason::IndeterminateTrial);
        assert!(report.interval.is_none());
        let t = report
            .trials
            .iter()
            .find(|t| t.evidence.factor == 2.)
            .unwrap();
        assert_eq!(t.evidence.outcome, Outcome::Indeterminate);
        assert!(t.epoch.as_ref().unwrap().payload.continuation.progress < 1.);
        Ok(())
    }
    #[test]
    fn post_cut_dormant_history_and_final_report_match_uninterrupted_direct_search() -> Result<()> {
        let (mesh, origin, spec) = plastic_origin(true)?;
        let expected = search::run(&mesh, copy_run(&origin)?, spec.clone())?;
        let dormant: Vec<_> = origin
            .integration_points
            .iter()
            .enumerate()
            .filter_map(|(i, p)| {
                (!origin.last_committed.active_elements[p.element_index]).then_some(i)
            })
            .collect();
        assert!(!dormant.is_empty());
        let mut s = Snapshot::new(&mesh, &origin, &spec)?;
        for _ in 0..1024 {
            if s.is_complete() {
                break;
            }
            s = restore(&s)?.advance(&mesh, &origin, &spec, 1, None)?;
            for cp in s
                .payload
                .inflight
                .iter()
                .chain(s.trials().iter().filter_map(|t| t.epoch.as_ref()))
            {
                for &i in &dormant {
                    assert_eq!(
                        checkpoint::hash_json(
                            &cp.payload
                                .material_checkpoint
                                .payload
                                .committed
                                .material_states[i]
                        )?,
                        checkpoint::hash_json(&origin.last_committed.material_states[i])?
                    );
                }
            }
        }
        s.verify(&mesh, &origin, &spec)?;
        same_report(&s, &expected)
    }
    #[test]
    fn invalid_budget_and_advance_after_finished_are_rejected() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let s = Snapshot::new(&mesh, &origin, &spec)?;
        assert!(
            restore(&s)?
                .advance(&mesh, &origin, &spec, 0, None)
                .is_err()
        );
        assert!(
            restore(&s)?
                .advance(&mesh, &origin, &spec, 10001, None)
                .is_err()
        );
        let done = complete(s, &mesh, &origin, &spec)?;
        assert!(done.advance(&mesh, &origin, &spec, 1, None).is_err());
        Ok(())
    }
    #[test]
    fn resealed_completed_trial_and_origin_audit_evidence_are_semantically_checked() -> Result<()> {
        let (mesh, origin, spec) = unit()?;
        let start = Snapshot::new(&mesh, &origin, &spec)?;
        let mut changed = restore(&start)?;
        changed.payload.basis.as_mut().unwrap().physical_force_kn[0] += 1e-12;
        changed.seal()?;
        assert!(changed.verify(&mesh, &origin, &spec).is_err());
        let complete = complete(start, &mesh, &origin, &spec)?;
        let mut changed = restore(&complete)?;
        changed.payload.report.as_mut().unwrap().trials[0]
            .evidence
            .adjusted_ground_sha256 = "0".repeat(64);
        changed.seal()?;
        assert!(changed.verify(&mesh, &origin, &spec).is_err());
        Ok(())
    }
}
