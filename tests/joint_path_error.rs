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
    let (mesh, ground, mut case) = fixture()?;
    case.settings.minimum_increment = 1. / (16. * 65536.);
    case.joint_system.as_mut().unwrap().global_step_control = Some(serde_json::from_value(
        serde_json::json!({"path_error":{"slip_scale_m":1e-4,"dissipation_scale_kpa_m":1e-2}}),
    )?);
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
        if a.accepted {
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
