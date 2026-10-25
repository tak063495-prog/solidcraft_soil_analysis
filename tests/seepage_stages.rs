use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool, mcc: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(if quadratic {
        include_str!("../examples/seepage-column-tet10-mesh.json")
    } else {
        include_str!("../examples/seepage-column-tet4-mesh.json")
    })?;
    let ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-seepage-column-mcc.json")
    } else {
        include_str!("../examples/ground-water-column-mc.json")
    })?;
    let case = serde_json::from_str(include_str!("../examples/stages-seepage-column.json"))?;
    Ok((mesh, ground, case))
}
#[test]
fn upward_downward_flow_excavation_and_dewatering_match_column_theory() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, case) = fixture(quadratic, false)?;
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(run.schema, "soil-fem-stage-result/12");
        assert_eq!(run.flow_audits.len(), 3);
        for (index, frame) in run.frames.iter().enumerate() {
            let active_nodes: std::collections::BTreeSet<_> = frame
                .state
                .active_elements
                .iter()
                .enumerate()
                .filter(|(_, v)| **v)
                .flat_map(|(i, _)| mesh.cell_nodes(i))
                .collect();
            for n in active_nodes {
                let z = mesh.nodes[n][2];
                let expected = match index {
                    0 | 3 => 0.,
                    1 => 4. * (z - z * z / 2.) / 12000.,
                    2 => -4. * (z - z * z / 2.) / 12000.,
                    4 => 5. * z / 12000.,
                    5 => (z + 4. * z * z) / 12000.,
                    _ => unreachable!(),
                };
                let tolerance = 2e-11;
                assert!(
                    (frame.state.displacement_m[n][2] - expected).abs() < tolerance,
                    "quadratic={quadratic} frame={} z={z} actual={} expected={expected}",
                    frame.name,
                    frame.state.displacement_m[n][2]
                );
            }
            assert!(frame.relative_force_residual < 1e-8);
        }
        assert!(run.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_water_changed
            && !a.rejected_trial_load_changed));
        for solution in &run.flow_audits {
            assert!(
                solution
                    .solution
                    .balances
                    .iter()
                    .all(|b| b.relative_nodal_balance_error < 1e-10)
            );
        }
    }
    Ok(())
}
#[test]
fn mcc_seepage_cutback_and_each_phase_checkpoint_restore_pressure_and_plastic_history() -> Result<()>
{
    let (mesh, ground, case) = fixture(true, true)?;
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.attempts.iter().any(|a| !a.accepted
        && a.rejected_trial_plastic_changed
        && a.rejected_trial_water_changed));
    let mut clean = case.clone();
    clean.verification_reject_once_stages.clear();
    let mut values = serde_json::to_value(&clean.stages)?;
    for stage in values.as_array_mut().unwrap() {
        let first = run
            .attempts
            .iter()
            .find(|a| a.accepted && a.stage == stage["name"].as_str().unwrap())
            .unwrap();
        stage["initial_increment"] = serde_json::json!(first.progress_to);
    }
    clean.stages = serde_json::from_value(values)?;
    let direct = stages::solve(&mesh, ground.clone(), clean)?;
    assert!(direct.complete);
    assert_eq!(
        serde_json::to_value(&run.frames)?,
        serde_json::to_value(&direct.frames)?
    );
    let mut count = 0;
    for stage in &case.stages {
        let value = serde_json::to_value(stage)?;
        let name = value["name"].as_str().unwrap();
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count + 1)?;
        let cp = Checkpoint::seal(&mesh, partial)?;
        let saved: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&cp)?)?;
        let restored = saved.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert!(restored.complete);
        assert_eq!(
            serde_json::to_value(&restored.frames)?,
            serde_json::to_value(&run.frames)?
        );
        for (a, b) in restored.flow_audits.iter().zip(&run.flow_audits) {
            assert!(soil_fem_3d::darcy::same_solution(&a.solution, &b.solution)?);
        }
        count += run
            .attempts
            .iter()
            .filter(|a| a.accepted && a.stage == name)
            .count();
    }
    for modify_audit in [false, true] {
        let mut partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
        if modify_audit {
            partial.flow_audits[0].solution.field.head_m[0] += 0.1;
        } else {
            partial.last_committed.pore_pressure_kpa.as_mut().unwrap()[0] += 0.1;
        }
        assert!(
            Checkpoint::seal(&mesh, partial)?
                .resume(&mesh, ground.clone(), case.clone(), None)
                .is_err()
        );
    }
    let stop = run
        .attempts
        .iter()
        .filter(|a| a.accepted && a.stage == "upward_flow")
        .count()
        + 1;
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), stop)?;
    let mut altered = serde_json::to_value(partial)?;
    altered["continuation"]["origin"]["pore_pressure_kpa"] = serde_json::json!([1.]);
    altered["frames"][1]["state"]["pore_pressure_kpa"] = serde_json::json!([1.]);
    let malformed: stages::Run = serde_json::from_value(altered)?;
    assert!(
        Checkpoint::seal(&mesh, malformed)?
            .resume(&mesh, ground.clone(), case.clone(), None)
            .is_err()
    );
    let before = &run.frames[4].state;
    let end = &run.last_committed;
    for (i, p) in run.integration_points.iter().enumerate() {
        if !end.active_elements[p.element_index] {
            assert_eq!(
                before.pore_pressure_kpa.as_ref().unwrap()[i],
                end.pore_pressure_kpa.as_ref().unwrap()[i]
            );
            assert_eq!(
                serde_json::to_value(&before.material_states[i])?,
                serde_json::to_value(&end.material_states[i])?
            );
        }
    }
    Ok(())
}
#[test]
fn unsaturated_candidate_does_not_change_the_committed_mechanical_state() -> Result<()> {
    let (mesh, ground, mut case) = fixture(true, false)?;
    let stages::Stage::SteadySeepage { flow, .. } = &mut case.stages[0] else {
        unreachable!()
    };
    for bc in &mut flow.heads {
        bc.head_m = 0.;
    }
    let run = stages::solve(&mesh, ground, case)?;
    assert!(!run.complete);
    assert_eq!(run.frames.len(), 1);
    assert!(!run.flow_audits[0].solution.complete);
    assert_eq!(
        serde_json::to_value(&run.last_committed)?,
        serde_json::to_value(&run.frames[0].state)?
    );
    Ok(())
}
