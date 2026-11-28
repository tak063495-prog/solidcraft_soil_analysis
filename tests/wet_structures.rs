//! Small integration regressions; independent 3D validation is in Python.
use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    Ok((
        serde_json::from_str(include_str!(
            "../examples/joint-water-column-tet10-mesh.json"
        ))?,
        serde_json::from_str(include_str!("../examples/ground-wet-plate-beam.json"))?,
        serde_json::from_str(include_str!("../examples/stages-wet-plate-beam.json"))?,
    ))
}

#[test]
fn wet_joint_plate_beam_keeps_physical_pressure_separate_through_flow_and_removal() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.schema, "soil-fem-stage-result/14");
    assert_eq!(run.frames.len(), 8);
    // Seam is at z=0.5. Independent two-layer resistance gives H=2+2/3.
    let pressures = [
        5.,
        5.,
        7.5,
        7.5,
        10. * (2. + 2. / 3. - 0.5),
        10. * (2. + 2. / 3. - 0.5),
        5.,
        5.,
    ];
    for (frame, expected) in run.frames.iter().zip(pressures) {
        assert!(frame.relative_force_residual <= 1e-8);
        let joint = frame.joint_results.as_ref().unwrap();
        for p in &joint.points {
            let water = p.hydraulic.as_ref().unwrap();
            assert!((water.physical_pore_pressure_kpa - expected).abs() < 1e-8);
            assert!(p.normal_pressure_kpa > 0.);
            assert!(
                (water.total_normal_pressure_kpa - p.normal_pressure_kpa - expected).abs() < 1e-8
            );
            assert!((water.friction_limit_kpa - 0.3 * p.normal_pressure_kpa).abs() < 1e-8);
        }
        assert!(
            frame
                .state
                .material_states
                .iter()
                .all(|s| serde_json::to_value(s).unwrap()["equivalent_plastic_strain"] == 0.)
        );
    }
    assert!(
        run.frames[3]
            .beam_results
            .as_ref()
            .unwrap()
            .plates
            .iter()
            .any(|p| p.bending_energy_kn_m > 0.)
    );
    let last = run.last_committed.beams.as_ref().unwrap();
    assert!(last.active_plates.iter().all(|v| !v));
    assert!(last.active_members.iter().all(|v| !v));
    assert_eq!(run.flow_audits.len(), 1);
    assert!(
        run.attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    Ok(())
}

#[test]
fn wet_flow_transition_with_active_plate_and_beam_resumes_exactly() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete);
    let initial = full
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let mut count = initial;
    for attempt in &full.attempts {
        if attempt.accepted {
            count += 1;
            if attempt.stage == "flow" {
                break;
            }
        }
    }
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
    assert!(!partial.complete);
    assert!(
        partial
            .last_committed
            .beams
            .as_ref()
            .unwrap()
            .active_plates
            .iter()
            .all(|v| *v)
    );
    assert!(partial.last_committed.joint_water.is_some());
    let checkpoint = Checkpoint::seal(&mesh, partial)?;
    let resumed = checkpoint.resume(&mesh, ground, case, None)?;
    assert!(resumed.complete);
    assert_eq!(
        serde_json::to_value(&full.last_committed)?,
        serde_json::to_value(&resumed.last_committed)?
    );
    assert_eq!(
        serde_json::to_value(&full.frames)?,
        serde_json::to_value(&resumed.frames)?
    );
    Ok(())
}
