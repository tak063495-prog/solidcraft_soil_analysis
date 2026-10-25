use serde_json::Value;
use soil_fem_3d::{
    Error, Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Kind, Observer, Progress},
    stages::{self, StageCase},
};
use std::cell::{Cell, RefCell};

#[derive(Clone, Copy)]
enum Trigger {
    Never,
    Accepted(usize),
    Rejected,
    Newton,
    Error,
}
struct Trace {
    events: RefCell<Vec<Progress>>,
    trigger: Trigger,
    accepted: Cell<usize>,
    stop: Cell<bool>,
}
impl Trace {
    fn new(trigger: Trigger) -> Self {
        Self {
            events: RefCell::new(vec![]),
            trigger,
            accepted: Cell::new(0),
            stop: Cell::new(false),
        }
    }
}
impl Observer for Trace {
    fn on_progress(&self, event: &Progress) -> Result<()> {
        self.events.borrow_mut().push(event.clone());
        if event.kind == Kind::Accepted {
            self.accepted.set(self.accepted.get() + 1);
        }
        match self.trigger {
            Trigger::Accepted(n) if self.accepted.get() >= n => self.stop.set(true),
            Trigger::Rejected if event.kind == Kind::Rejected => self.stop.set(true),
            Trigger::Newton if event.kind == Kind::Newton => self.stop.set(true),
            Trigger::Error if event.kind == Kind::Newton => {
                return Err(Error::Input("diagnostic sink unavailable".into()));
            }
            _ => (),
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(self.stop.get())
    }
}
fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(include_str!("../examples/hardening-column-tet10-mesh.json"))?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-mixed-mcc.json"))?;
    for layer in &mut ground.layers {
        layer.k0_x = 1.;
        layer.k0_y = 1.;
        layer.consolidation =
            soil_fem_3d::ground::Consolidation::MccPcOverMeanEffective { ratio: 1. };
    }
    Ok((
        mesh,
        ground,
        serde_json::from_str(include_str!("../examples/stages-mixed-mcc.json"))?,
    ))
}
fn normalized(mut value: Value) -> Value {
    fn visit(value: &mut Value) {
        match value {
            Value::Object(map) => {
                map.remove("elapsed_seconds");
                for child in map.values_mut() {
                    visit(child);
                }
            }
            Value::Array(items) => {
                for child in items {
                    visit(child);
                }
            }
            _ => (),
        }
    }
    visit(&mut value);
    value
}
fn same(a: &impl serde::Serialize, b: &impl serde::Serialize) -> Result<()> {
    assert!(
        normalized(serde_json::to_value(a)?) == normalized(serde_json::to_value(b)?),
        "physical states/path or non-time diagnostics changed"
    );
    Ok(())
}

#[test]
fn observation_preserves_every_result_field_and_reports_actual_newton_linear_and_attempt_values()
-> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let control = stages::solve(&mesh, ground.clone(), case.clone())?;
    let trace = Trace::new(Trigger::Never);
    let observed = stages::solve_observed(&mesh, ground, case, &trace)?;
    assert!(observed.complete);
    same(&control, &observed)?;
    let events = trace.events.borrow();
    assert!(events.iter().any(|e| e.kind == Kind::Newton
        && e.relative_volume_residual.is_some()
        && e.relative_pressure_residual.is_some()));
    assert!(events.iter().any(|e| {
        e.kind == Kind::LinearSolved
            && e.linear_solve
                .as_ref()
                .is_some_and(|r| r.method == "csr_gmres_ilu0_volume_first_rcm_up")
    }));
    let accepted: Vec<_> = events.iter().filter(|e| e.kind == Kind::Accepted).collect();
    assert_eq!(accepted.len(), observed.attempts.len());
    for (event, attempt) in accepted.iter().zip(&observed.attempts) {
        assert_eq!(event.stage, attempt.stage);
        assert_eq!(event.progress.to_bits(), attempt.progress_to.to_bits());
        assert_eq!(event.newton_iteration, attempt.newton_iterations);
        assert_eq!(
            event.relative_force_residual,
            attempt.relative_force_residual
        );
    }
    Ok(())
}

#[test]
fn requests_from_newton_or_rejected_trial_pause_only_after_acceptance_and_resume_exactly()
-> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.verification_reject_once_stages = vec!["load".into()];
    let control = stages::solve(&mesh, ground.clone(), case.clone())?;
    for trigger in [Trigger::Newton, Trigger::Rejected, Trigger::Accepted(2)] {
        let trace = Trace::new(trigger);
        let partial = stages::solve_observed(&mesh, ground.clone(), case.clone(), &trace)?;
        assert!(
            !partial.complete && partial.failure_reason.as_deref() == Some("checkpoint_requested")
        );
        assert!(
            partial
                .attempts
                .iter()
                .any(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
        );
        assert!(
            partial
                .last_committed
                .material_states
                .iter()
                .any(|s| s.plastic_strain().iter().any(|x| x.abs() > 1e-8))
        );
        let count = partial.attempts.iter().filter(|a| a.accepted).count();
        assert_eq!(
            count,
            if matches!(trigger, Trigger::Accepted(2)) {
                2
            } else {
                1
            }
        );
        let expected = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        same(&partial, &expected)?;
        let resumed_trace = Trace::new(Trigger::Never);
        let resumed = Checkpoint::seal(&mesh, partial)?.resume_observed(
            &mesh,
            ground.clone(),
            case.clone(),
            &resumed_trace,
        )?;
        assert!(resumed.complete);
        same(&control, &resumed)?;
        assert_eq!(
            resumed_trace.accepted.get(),
            control.attempts.iter().filter(|a| a.accepted).count() - count
        );
    }
    Ok(())
}

#[test]
fn observer_errors_remain_runtime_errors_and_never_trigger_physical_cutback() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let trace = Trace::new(Trigger::Error);
    let result = stages::solve_observed(&mesh, ground, case, &trace);
    assert!(
        matches!(result,Err(Error::Control(reason)) if reason.contains("diagnostic sink unavailable"))
    );
    assert_eq!(
        trace
            .events
            .borrow()
            .iter()
            .filter(|e| e.kind == Kind::AttemptStarted)
            .count(),
        1
    );
    assert!(
        !trace
            .events
            .borrow()
            .iter()
            .any(|e| e.kind == Kind::Rejected)
    );
    Ok(())
}

#[test]
fn geostatic_and_initialization_endpoint_pause_keep_reference_and_all_history_on_resume()
-> Result<()> {
    let (mut mesh, _, _) = fixture()?;
    mesh.face_sets.insert(
        "top".into(),
        soil_fem_3d::boundary::surface_owners(&mesh)
            .into_iter()
            .filter(|(face, owners)| {
                owners.len() == 1 && face.iter().all(|&n| mesh.nodes[n][2] == 1.)
            })
            .map(|(face, _)| face)
            .collect(),
    );
    let ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.schema = "soil-fem-mixed-mcc-stages/1".into();
    let control = stages::solve(&mesh, ground.clone(), case.clone())?;
    let total = control
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    for n in [1, total] {
        let trace = Trace::new(Trigger::Accepted(n));
        let partial = stages::solve_observed(&mesh, ground.clone(), case.clone(), &trace)?;
        assert!(partial.attempts.is_empty());
        assert_eq!(partial.geostatic.as_ref().unwrap().complete, n == total);
        assert_eq!(
            partial.failure_reason.as_deref(),
            Some("checkpoint_requested")
        );
        let expected = stages::solve_until(&mesh, ground.clone(), case.clone(), n)?;
        same(&partial, &expected)?;
        let resumed = Checkpoint::seal(&mesh, partial)?.resume_observed(
            &mesh,
            ground.clone(),
            case.clone(),
            &Trace::new(Trigger::Never),
        )?;
        same(&control, &resumed)?;
    }
    Ok(())
}
