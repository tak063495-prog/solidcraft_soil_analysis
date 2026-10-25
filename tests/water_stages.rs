use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool, mcc: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(if quadratic {
        include_str!("../examples/hardening-column-tet10-mesh.json")
    } else {
        include_str!("../examples/mc-column-mesh.json")
    })?;
    let ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-water-column-mcc.json")
    } else {
        include_str!("../examples/ground-water-column-mc.json")
    })?;
    let case = serde_json::from_str(include_str!("../examples/stages-water-column.json"))?;
    Ok((mesh, ground, case))
}

#[test]
fn hydrostatic_drawdown_flooding_and_buoyant_excavation_follow_column_theory() -> Result<()> {
    for quadratic in [false, true] {
        let (mut mesh, ground, case) = fixture(quadratic, false)?;
        if quadratic {
            mesh.quadrature_order = Some(5);
        }
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(run.schema, "soil-fem-stage-result/11");
        let modulus = 12000.;
        let gamma = 10.;
        for (i, frame) in run.frames.iter().enumerate() {
            assert!(frame.relative_force_residual < 1e-8);
            let h = frame.state.hydrostatic_water_level_m.unwrap();
            for (n, p) in mesh.nodes.iter().enumerate() {
                if i == run.frames.len() - 1 && p[2] > 0.5 + 1e-12 {
                    continue;
                }
                let z = p[2];
                let wet = h.max(0.).min(z);
                let integral = z - z * z / 2. - h * wet + wet * wet / 2. + (h - 1.).max(0.) * z;
                let expected = if i == run.frames.len() - 1 {
                    5. * z / modulus
                } else {
                    -gamma * integral / modulus
                };
                let error = (frame.state.displacement_m[n][2] - expected).abs();
                // A water horizon crossing an element is integrated at the
                // declared volume material rule; its kink has quadrature error.
                let tolerance = if frame.name == "drawdown" {
                    if quadratic { 2e-6 } else { 2e-5 }
                } else {
                    2e-11
                };
                assert!(
                    error < tolerance,
                    "quadratic={quadratic}, stage={}, z={z}, error={error:e}",
                    frame.name
                );
            }
            assert!(
                frame
                    .state
                    .material_states
                    .iter()
                    .all(|s| s.plastic_strain().iter().all(|e| e.abs() < 1e-12))
            );
        }
        assert!(run.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_water_changed
            && a.committed_state_unchanged_on_rejection));
    }
    Ok(())
}

#[test]
fn plastic_water_history_and_partial_excavation_resume_exactly_and_reject_tampering() -> Result<()>
{
    let (mesh, ground, case) = fixture(true, true)?;
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.attempts.iter().any(|a| !a.accepted
        && a.rejected_trial_water_changed
        && a.rejected_trial_plastic_changed));
    let mut clean = case.clone();
    clean.verification_reject_once_stages.clear();
    let mut value = serde_json::to_value(&clean.stages)?;
    for s in value.as_array_mut().unwrap() {
        let name = s["name"].as_str().unwrap();
        let increment = run
            .attempts
            .iter()
            .find(|a| a.accepted && a.stage == name)
            .unwrap()
            .progress_to;
        s["initial_increment"] = serde_json::json!(increment);
    }
    clean.stages = serde_json::from_value(value)?;
    let direct = stages::solve(&mesh, ground.clone(), clean)?;
    assert!(direct.complete);
    assert_eq!(
        serde_json::to_value(&direct.frames)?,
        serde_json::to_value(&run.frames)?
    );
    let mut count = 0;
    for stage in &run.case.stages {
        let stage = serde_json::to_value(stage)?;
        let name = stage["name"].as_str().unwrap();
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count + 1)?;
        let cp = Checkpoint::seal(&mesh, partial)?;
        let restored: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&cp)?)?;
        let resumed = restored.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert!(resumed.complete);
        assert_eq!(
            serde_json::to_value(&resumed.frames)?,
            serde_json::to_value(&run.frames)?
        );
        count += run
            .attempts
            .iter()
            .filter(|a| a.accepted && a.stage == name)
            .count();
    }
    let mut partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    partial.last_committed.hydrostatic_water_level_m = Some(0.1234);
    assert!(
        Checkpoint::seal(&mesh, partial)?
            .resume(&mesh, ground, case, None)
            .is_err()
    );
    let before = &run.frames[run.frames.len() - 2].state;
    let last = &run.last_committed;
    for (index, p) in run.integration_points.iter().enumerate() {
        if !last.active_elements[p.element_index] {
            assert_eq!(
                serde_json::to_value(&before.material_states[index])?,
                serde_json::to_value(&last.material_states[index])?
            );
        }
    }
    Ok(())
}

#[test]
fn intermediate_water_head_rebuilds_actual_pressure_and_does_not_interpolate_endpoint_loads()
-> Result<()> {
    let (mesh, ground, case) = fixture(true, false)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    let count = full
        .attempts
        .iter()
        .filter(|a| a.accepted && a.stage == "drawdown")
        .count();
    let partial = stages::solve_until(&mesh, ground, case, count + 1)?;
    let s = &partial.last_committed;
    assert_eq!(s.stage_index, 2);
    let h = s.hydrostatic_water_level_m.unwrap();
    assert!(h > 0.8 && h < 1.);
    // h=0.835 has no external water above the one-metre column. A linear
    // blend towards the h=1.3 endpoint would incorrectly load its top.
    let vertical: f64 = s
        .applied_nodal_load_kn
        .as_chunks::<3>()
        .0
        .iter()
        .map(|f| f[2])
        .sum();
    assert!((vertical + 120.).abs() < 1e-10);
    assert!(
        full.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_water_changed
            && !a.rejected_trial_load_changed)
    );
    Ok(())
}
