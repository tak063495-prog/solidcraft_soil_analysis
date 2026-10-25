use serde_json::json;
use soil_fem_3d::{
    Result,
    checkpoint::{self, Checkpoint},
    geostatic::Seed,
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{self, Stage, StageCase},
};

fn fixture(mcc: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let top = soil_fem_3d::boundary::surface_owners(&mesh)
        .into_iter()
        .filter(|(f, o)| o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][2] == 1.))
        .map(|(f, _)| f)
        .collect();
    mesh.face_sets.insert("top".into(), top);
    let ground: Ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-geostatic-mcc-column.json")
    } else {
        include_str!("../examples/ground-geostatic-mc-column.json")
    })?;
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    Ok((mesh.elevate_straight()?, ground, case))
}

#[test]
fn display_reset_preserves_every_physical_history_for_mc_and_mcc() -> Result<()> {
    for mcc in [false, true] {
        let (mesh, ground, mut case) = fixture(mcc)?;
        case.stages.insert(
            1,
            Stage::ResetDisplay {
                name: "zero_display".into(),
            },
        );
        case.verification_reject_once_stages = vec!["zero_display".into()];
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete, "{:?}", full.failure_reason);
        assert_eq!(full.schema, "soil-fem-stage-result/6");
        let before = &full.frames[1].state;
        let after = &full.frames[2].state;
        assert_eq!(before.displacement_m, after.displacement_m);
        assert_eq!(
            serde_json::to_vec(&before.material_states)?,
            serde_json::to_vec(&after.material_states)?
        );
        assert_eq!(before.applied_nodal_load_kn, after.applied_nodal_load_kn);
        assert_eq!(
            after.display_reference_displacement_m,
            before.displacement_m
        );
        assert_eq!(
            full.frames[2].monitors["surface"].mean_displacement_z_m,
            Some(0.)
        );
        assert!(
            before
                .material_states
                .iter()
                .any(|s| s.plastic_strain().iter().any(|x| x.abs() > 1e-8))
        );
        assert!(full.attempts.iter().any(|a| a.stage == "zero_display"
            && !a.accepted
            && a.committed_state_unchanged_on_rejection));
        let mut plain = case.clone();
        plain.stages.remove(1);
        plain.verification_reject_once_stages.clear();
        let plain = stages::solve(&mesh, ground.clone(), plain)?;
        assert_eq!(
            full.last_committed.displacement_m,
            plain.last_committed.displacement_m
        );
        assert_eq!(
            serde_json::to_vec(&full.last_committed.material_states)?,
            serde_json::to_vec(&plain.last_committed.material_states)?
        );
        let count = full
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .count()
            + full
                .attempts
                .iter()
                .take_while(|a| a.stage != "zero_display")
                .filter(|a| a.accepted)
                .count()
            + 1;
        let cp = Checkpoint::seal(
            &mesh,
            stages::solve_until(&mesh, ground.clone(), case.clone(), count)?,
        )?;
        assert_eq!(cp.run.last_committed.stage_progress, 1.);
        let cloned: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&cp)?)?;
        let resumed = cloned.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert_eq!(
            serde_json::to_vec(&full.frames)?,
            serde_json::to_vec(&resumed.frames)?
        );
        let mut bad = cp;
        bad.run.last_committed.display_reference_displacement_m[0][2] += 0.1;
        bad.run_sha256 = checkpoint::hash_json(&bad.run)?;
        assert!(bad.resume(&mesh, ground, case, None).is_err());
    }
    Ok(())
}

#[test]
fn explicit_state_reset_reequilibrates_without_resetting_constitutive_displacement() -> Result<()> {
    let (mesh, mut ground, mut case) = fixture(false)?;
    for m in &mut ground.materials {
        if let soil::Model::MohrCoulomb(p) = &mut m.model {
            p.cohesion_kpa = 1000.;
        }
    }
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    let Seed::BodyStates { states } = &ground.initialization.geostatic.as_ref().unwrap().seed
    else {
        unreachable!()
    };
    case.stages = vec![Stage::ResetState {
        name: "replace_state".into(),
        bodies: vec!["LowerClay".into(), "UpperClay".into()],
        states: states.clone(),
        initial_increment: Some(0.25),
    }];
    case.verification_reject_once_stages = vec!["replace_state".into()];
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    let origin = &full.frames[0].state;
    let after = &full.frames[1].state;
    let audit = &full.transition_audits[0];
    assert_eq!(audit.prepared_state.displacement_m, origin.displacement_m);
    assert!(
        audit
            .prepared_state
            .material_states
            .iter()
            .all(|s| s.stress().iter().all(|x| *x == 0.))
    );
    for (a, b) in after
        .displacement_m
        .iter()
        .flatten()
        .zip(origin.displacement_m.iter().flatten())
    {
        assert!((a - 2. * b).abs() < 1e-9);
    }
    for (a, b) in after.material_states.iter().zip(&origin.material_states) {
        for (s, t) in a.stress().iter().zip(b.stress()) {
            assert!((s - t).abs() < 1e-6);
        }
    }
    assert!(after.display_reference_displacement_m.is_empty());
    let count = full
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count()
        + 1;
    let cp = Checkpoint::seal(
        &mesh,
        stages::solve_until(&mesh, ground.clone(), case.clone(), count)?,
    )?;
    assert!(
        cp.run
            .last_committed
            .transition_release_load_kn
            .iter()
            .any(|x| x.abs() > 1e-4)
    );
    let resumed = cp.resume(&mesh, ground, case, None)?;
    assert_eq!(
        serde_json::to_vec(&full.last_committed)?,
        serde_json::to_vec(&resumed.last_committed)?
    );
    assert_eq!(
        serde_json::to_vec(&full.transition_audits)?,
        serde_json::to_vec(&resumed.transition_audits)?
    );
    Ok(())
}

#[test]
fn born_nodes_get_their_new_display_origin_after_reset() -> Result<()> {
    let (mesh, ground, _) = fixture(false)?;
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    let Seed::BodyStates { states } = &ground.initialization.geostatic.as_ref().unwrap().seed
    else {
        unreachable!()
    };
    case["initially_inactive_bodies"] = json!(["UpperClay"]);
    case["stages"] = json!([
        {"kind":"birth","name":"fill","placement":"cad_new_nodes","bodies":["UpperClay"],"states":{"UpperClay":states["UpperClay"]}},
        {"kind":"reset_display","name":"display_epoch"},
        {"kind":"excavate","name":"remove","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.},
        {"kind":"birth","name":"refill","placement":"cad_new_nodes","bodies":["UpperClay"],"states":{"UpperClay":states["UpperClay"]}}
    ]);
    let run = stages::solve(&mesh, ground, serde_json::from_value(case)?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for &n in &mesh.node_sets["top"] {
        assert_eq!(
            run.last_committed.display_reference_displacement_m[n],
            [0.; 3]
        );
    }
    Ok(())
}

#[test]
fn mcc_reset_uses_explicit_admissible_pc_and_void_ratio_then_forms_new_history() -> Result<()> {
    let (mesh, ground, mut case) = fixture(true)?;
    let Seed::BodyStates { states } = &ground.initialization.geostatic.as_ref().unwrap().seed
    else {
        unreachable!()
    };
    case.stages = vec![Stage::ResetState {
        name: "mcc_state".into(),
        bodies: vec!["LowerClay".into(), "UpperClay".into()],
        states: states.clone(),
        initial_increment: Some(0.25),
    }];
    case.verification_reject_once_stages = vec!["mcc_state".into()];
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let audit = &run.transition_audits[0];
    assert!(
        audit.prepared_state.material_states.iter().all(|s| s
            .as_mcc()
            .unwrap()
            .preconsolidation_pressure_kpa
            == 140.)
    );
    assert!(
        run.last_committed.material_states.iter().any(|s| s
            .as_mcc()
            .unwrap()
            .preconsolidation_pressure_kpa
            > 140.01)
    );
    assert_eq!(
        audit.prepared_state.displacement_m,
        run.frames[0].state.displacement_m
    );
    assert!(
        run.last_committed
            .display_reference_displacement_m
            .is_empty()
    );
    Ok(())
}
