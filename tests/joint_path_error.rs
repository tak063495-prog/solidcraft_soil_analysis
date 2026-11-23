use soil_fem_3d::{
    Result,
    checkpoint::{Checkpoint, hash_json},
    ground::Ground,
    joints::GlobalStepControl,
    mesh::Mesh,
    stages::{self, Stage, StageCase},
};

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(include_str!(
        "../examples/joint-water-column-tet10-mesh.json"
    ))?;
    let ground = serde_json::from_str(include_str!("../examples/ground-joint-water-mc.json"))?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-water-cycle.json"))?;
    case.settings.initial_increment = 1. / 16.;
    case.verification_reject_once_stages.clear();
    case.joint_system.as_mut().unwrap().maximum_increment_m = 5e-5;
    for stage in &mut case.stages {
        if let Stage::Boundary {
            initial_increment, ..
        } = stage
        {
            *initial_increment = Some(1. / 16.);
        }
    }
    Ok((mesh, ground, case))
}

#[test]
fn path_error_option_is_explicit_and_separate_from_legacy_switch_control() -> Result<()> {
    let (_, _, mut case) = fixture()?;
    let legacy = hash_json(&case)?;
    let config = |value| serde_json::from_value::<GlobalStepControl>(value);
    let valid =
        serde_json::json!({"path_error":{"slip_scale_m":1e-4,"dissipation_scale_kpa_m":1e-2}});
    case.joint_system.as_mut().unwrap().global_step_control = Some(config(valid.clone())?);
    case.validate()?;
    assert_ne!(hash_json(&case)?, legacy);
    assert_eq!(
        serde_json::to_value(
            case.joint_system
                .as_ref()
                .unwrap()
                .global_step_control
                .as_ref()
                .unwrap()
        )?,
        valid
    );
    for value in [
        serde_json::json!({"switch_bisections":4,"path_error":{"slip_scale_m":1e-4,"dissipation_scale_kpa_m":1e-2}}),
        serde_json::json!({"path_error":{"slip_scale_m":0.,"dissipation_scale_kpa_m":1e-2}}),
        serde_json::json!({"path_error":{"slip_scale_m":1e-4,"dissipation_scale_kpa_m":-1.}}),
    ] {
        case.joint_system.as_mut().unwrap().global_step_control = Some(config(value)?);
        assert!(case.validate().is_err());
    }
    Ok(())
}

#[test]
fn global_error_rejections_leave_history_intact_and_checkpoint_replays_exactly() -> Result<()> {
    check_path_replay(false)
}

#[test]
fn event_proposals_recheck_all_gates_and_checkpoint_replays_exactly() -> Result<()> {
    check_path_replay(true)
}

fn check_path_replay(event_proposal: bool) -> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.settings.minimum_increment = 1. / (16. * 65536.);
    case.joint_system.as_mut().unwrap().global_step_control = Some(serde_json::from_value(
        serde_json::json!({"path_error":{"slip_scale_m":1e-4,"dissipation_scale_kpa_m":1e-2,"event_proposal":event_proposal}}),
    )?);
    if event_proposal {
        case.joint_system
            .as_mut()
            .unwrap()
            .global_step_control
            .as_mut()
            .unwrap()
            .path_error
            .as_mut()
            .unwrap()
            .minimum_contact_substeps = Some(4);
    }
    let whole = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(whole.complete, "{:?}", whole.failure_reason);
    let attempts = whole
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .chain(&whole.attempts);
    let mut rejects = 0;
    for a in attempts {
        let audit = a
            .joint_path_error
            .as_ref()
            .expect("every converged trial has a global estimate");
        assert!(
            audit.midpoint_progress > a.progress_from && audit.midpoint_progress < a.progress_to
        );
        assert!(audit.midpoint_relative_force_residual <= case.settings.relative_force_tolerance);
        assert!(audit.fine_relative_force_residual <= case.settings.relative_force_tolerance);
        assert!(
            audit
                .midpoint_linear_solves
                .iter()
                .chain(&audit.fine_linear_solves)
                .all(|l| l.true_relative_residual <= case.settings.relative_linear_tolerance)
        );
        assert_eq!(audit.proposed_increment.is_some(), event_proposal);
        if a.accepted {
            assert!(audit.maximum_slip_error_m <= audit.slip_budget_m);
            assert!(audit.maximum_dissipation_error_kpa_m <= audit.dissipation_budget_kpa_m);
            assert!(!audit.reject);
        } else {
            assert!(a.committed_state_unchanged_on_rejection);
            assert_eq!(a.reason, "joint_global_path_error_cutback");
            rejects += 1;
        }
    }
    assert!(rejects > 0);
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 35)?;
    let checkpoint = Checkpoint::seal(&mesh, partial)?;
    let checkpoint: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&checkpoint)?)?;
    let resumed = checkpoint.resume(&mesh, ground, case, None)?;
    fn strip_clock(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Object(m) => {
                m.remove("elapsed_seconds");
                for x in m.values_mut() {
                    strip_clock(x);
                }
            }
            serde_json::Value::Array(a) => {
                for x in a {
                    strip_clock(x);
                }
            }
            _ => {}
        }
    }
    let mut a = serde_json::to_value(whole)?;
    let mut b = serde_json::to_value(resumed)?;
    strip_clock(&mut a);
    strip_clock(&mut b);
    assert_eq!(a, b);
    Ok(())
}

#[test]
fn diagnostic_gap_bisection_reuses_the_error_budget_and_reports_the_crossing() -> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.settings.minimum_increment = 1. / (16. * 65536.);
    let control = serde_json::json!({"path_error":{
        "slip_scale_m":1e-4, "dissipation_scale_kpa_m":1e-2,
        "event_proposal":true, "minimum_contact_substeps":4,
        "event_comparison":"bisect_rejections"
    }});
    case.joint_system.as_mut().unwrap().global_step_control =
        Some(serde_json::from_value(control)?);
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete);
    let mut rejected = 0;
    let mut crossings = 0;
    for attempt in &run.attempts {
        let audit = attempt.joint_path_error.as_ref().unwrap();
        if let Some((position, gp)) = audit.first_gap_zero {
            assert!((attempt.progress_from..=attempt.progress_to).contains(&position));
            assert!(gp < audit.old_status.len());
            crossings += 1;
        }
        if audit.reject {
            assert!(!attempt.accepted && attempt.committed_state_unchanged_on_rejection);
            assert_eq!(audit.proposed_increment, Some(audit.actual_increment / 2.));
            rejected += 1;
        }
    }
    assert!(crossings > 0 && rejected > 0);
    Ok(())
}

#[test]
fn proportional_stage_budget_counts_spatial_cells_and_preserves_legacy_identity() -> Result<()> {
    use soil_fem_3d::stages::AttemptBudgetPolicy;
    let (mesh, _, case) = fixture()?;
    let legacy = serde_json::to_value(&case.settings)?;
    assert!(legacy.get("attempt_budget_policy").is_none());
    assert_eq!(
        case.settings.attempt_budget_policy.stage_limit(&mesh)?,
        2000
    );
    let mut settings = legacy.clone();
    settings["attempt_budget_policy"] = "fixed".into();
    assert_eq!(
        serde_json::to_value(serde_json::from_value::<stages::Settings>(
            settings.clone()
        )?)?,
        legacy
    );
    settings["attempt_budget_policy"] = "mesh_proportional".into();
    let opt: stages::Settings = serde_json::from_value(settings.clone())?;
    assert_eq!(opt.attempt_budget_policy.stage_limit(&mesh)?, 2065);
    assert_eq!(serde_json::to_value(opt)?, settings);
    // A separate single-cell mesh exercises the ceiling on a non-divisible product.
    let mut cell = mesh;
    cell.tetrahedra.truncate(1);
    cell.tetrahedra10.truncate(1);
    assert_eq!(
        AttemptBudgetPolicy::MeshProportional.stage_limit(&cell)?,
        173
    );
    Ok(())
}
