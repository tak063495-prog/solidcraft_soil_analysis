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
fn absent_and_null_control_preserve_the_pinned_legacy_case_hash() -> Result<()> {
    let (_, _, case) = fixture()?;
    // Recorded by the unmodified 9338204 executable, before adding this option.
    let expected = "ef1073793309d4a0a4b9d7d3f0c9ca290ac38734ce7b11022433eb0a0290e7e5";
    assert_eq!(hash_json(&case)?, expected);
    let mut json = serde_json::to_value(&case)?;
    assert!(json["joint_system"].get("global_step_control").is_none());
    json["joint_system"]["global_step_control"] = serde_json::Value::Null;
    let decoded: StageCase = serde_json::from_value(json)?;
    assert_eq!(hash_json(&decoded)?, expected);
    Ok(())
}

#[test]
fn invalid_resolution_and_unsupported_stage_are_refused() -> Result<()> {
    let (_, _, mut case) = fixture()?;
    for depth in [0, 9, u32::MAX] {
        case.joint_system.as_mut().unwrap().global_step_control = Some(GlobalStepControl {
            switch_bisections: depth,
        });
        assert!(case.validate().is_err());
    }
    case.joint_system.as_mut().unwrap().global_step_control = Some(GlobalStepControl {
        switch_bisections: 2,
    });
    case.validate()?;
    case.settings.minimum_increment = case.settings.initial_increment;
    assert!(case.validate().is_err());
    case.settings.minimum_increment = 1. / 1024.;
    case.stages[0] = Stage::ResetDisplay {
        name: "unsupported".into(),
    };
    assert!(case.validate().is_err());
    Ok(())
}

#[test]
fn adaptive_global_equilibria_keep_anchors_reject_immutably_and_resume_exactly() -> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.joint_system.as_mut().unwrap().global_step_control = Some(GlobalStepControl {
        switch_bisections: 2,
    });
    let whole = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(whole.complete, "{:?}", whole.failure_reason);
    let mut phases: Vec<Vec<_>> = vec![whole.geostatic.as_ref().unwrap().attempts.iter().collect()];
    for stage in &case.stages {
        let Stage::Boundary { name, .. } = stage else {
            unreachable!()
        };
        phases.push(whole.attempts.iter().filter(|a| &a.stage == name).collect());
    }
    let mut cuts = 0;
    for attempts in phases {
        for i in 1..=16 {
            assert!(
                attempts
                    .iter()
                    .any(|a| a.accepted && a.progress_to == i as f64 / 16.)
            );
        }
        for a in attempts {
            if a.accepted {
                assert!(
                    a.relative_force_residual.unwrap() <= case.settings.relative_force_tolerance
                );
            } else {
                assert!(a.committed_state_unchanged_on_rejection);
                cuts += usize::from(a.reason == "joint_global_switch_cutback");
            }
        }
    }
    assert!(cuts > 0);
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 35)?;
    let checkpoint = Checkpoint::seal(&mesh, partial)?;
    let checkpoint: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&checkpoint)?)?;
    let resumed = checkpoint.resume(&mesh, ground, case, None)?;
    fn strip_clock(value: &mut serde_json::Value) {
        match value {
            serde_json::Value::Object(map) => {
                map.remove("elapsed_seconds");
                for child in map.values_mut() {
                    strip_clock(child);
                }
            }
            serde_json::Value::Array(values) => {
                for child in values {
                    strip_clock(child);
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
