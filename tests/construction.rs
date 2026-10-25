use serde_json::json;
use soil_fem_3d::{
    Result,
    checkpoint::{self, Checkpoint},
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};

fn fixture(plastic: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let top = soil_fem_3d::boundary::surface_owners(&mesh)
        .into_iter()
        .filter(|(f, o)| o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][2] == 1.))
        .map(|(f, _)| f)
        .collect();
    mesh.face_sets.insert("top".into(), top);
    let mesh = mesh.elevate_straight()?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    for (i, m) in ground.materials.iter_mut().enumerate() {
        if let soil::Model::MohrCoulomb(p) = &mut m.model {
            p.cohesion_kpa = if plastic && i == 1 { 2. } else { 1000. };
        }
    }
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    let mut replacement = ground.materials[1].clone();
    replacement.id = 3;
    replacement.name = "stiffer fill".into();
    if let soil::Model::MohrCoulomb(p) = &mut replacement.model {
        p.young_kpa *= 2.;
    }
    ground.materials.push(replacement);
    ground.construction_material_ids = vec![3];
    let mut case = serde_json::to_value(serde_json::from_str::<StageCase>(include_str!(
        "../examples/stages-geostatic-column.json"
    ))?)?;
    case["initially_inactive_bodies"] = json!(["UpperClay"]);
    let seed = match &ground.initialization.geostatic.as_ref().unwrap().seed {
        soil_fem_3d::geostatic::Seed::BodyStates { states } => {
            serde_json::to_value(&states["UpperClay"])?
        }
        _ => unreachable!(),
    };
    case["stages"] = json!([
        {"kind":"birth","placement":"cad_new_nodes","name":"place_fill","bodies":["UpperClay"],"states":{"UpperClay":seed}},
        {"kind":"material_set","name":"improve_fill","bodies":["UpperClay"],"material_id":3,"history":{"mode":"preserve"}},
        {"kind":"pressure","name":"load","target_pressure_kpa":10.},
        {"kind":"pressure","name":"unload","target_pressure_kpa":0.}
    ]);
    case["verification_reject_once_stages"] = json!(["place_fill", "improve_fill"]);
    Ok((mesh, ground, serde_json::from_value(case)?))
}

#[test]
fn stress_free_fill_and_material_change_match_oedometer_theory_and_restore_rejected_trials()
-> Result<()> {
    let (mesh, ground, case) = fixture(false)?;
    let run = stages::solve(&mesh, ground.clone(), case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.schema, "soil-fem-stage-result/5");
    assert_eq!(run.transition_audits.len(), 2);
    let born = &run.transition_audits[0];
    assert!(born.zero_release_relative_force_residual < 1e-12);
    assert_eq!(
        born.prepared_state.displacement_m,
        run.frames[0].state.displacement_m
    );
    let before = &run.frames[0].state;
    assert!(before.active_elements.iter().any(|a| !a));
    for p in &run.integration_points {
        let ei = p.element_index;
        let state = &born.prepared_state.material_states[run
            .integration_points
            .iter()
            .position(|a| a.element_index == ei && a.local_index == p.local_index)
            .unwrap()];
        if born.selected_elements.contains(&ei) {
            assert!(state.stress().iter().all(|x| *x == 0.));
            assert!(state.plastic_strain().iter().all(|x| *x == 0.));
        }
    }
    assert!(run.frames[1].state.active_elements.iter().all(|a| *a));
    assert_eq!(
        run.frames[1].state.displacement_m,
        run.frames[2].state.displacement_m
    );
    assert_eq!(
        serde_json::to_vec(&run.frames[1].state.material_states)?,
        serde_json::to_vec(&run.frames[2].state.material_states)?
    );
    let modulus = |id| {
        let m = ground.materials.iter().find(|m| m.id == id).unwrap();
        if let soil::Model::MohrCoulomb(p) = &m.model {
            p.young_kpa * (1. - p.poisson) / ((1. + p.poisson) * (1. - 2. * p.poisson))
        } else {
            unreachable!()
        }
    };
    let top = &mesh.node_sets["top"];
    let mean = |frame: usize| {
        top.iter()
            .map(|&n| run.frames[frame].state.displacement_m[n][2])
            .sum::<f64>()
            / top.len() as f64
    };
    let expected_birth = -25. / modulus(1) - 12.5 / modulus(2);
    assert!(
        (mean(1) - expected_birth).abs() < 1e-9,
        "{} vs {}",
        mean(1),
        expected_birth
    );
    assert!((mean(3) - mean(2) + 5. / modulus(1) + 5. / modulus(3)).abs() < 1e-9);
    assert!((mean(4) - mean(2)).abs() < 1e-9);
    for frame in &run.frames {
        assert!(frame.relative_force_residual < 1e-8);
    }
    assert!(
        run.attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    assert!(run.attempts.iter().any(|a| !a.accepted
        && a.rejected_trial_active_set_changed
        && a.rejected_trial_material_changed));
    Ok(())
}

#[test]
fn plastic_birth_history_material_ids_and_holding_load_resume_exactly() -> Result<()> {
    let (mesh, ground, case) = fixture(true)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    assert!(full.attempts.iter().any(|a| !a.accepted
        && a.rejected_trial_plastic_changed
        && a.committed_state_unchanged_on_rejection));
    let geostatic_count = full
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let birth_count = full
        .attempts
        .iter()
        .filter(|a| a.accepted && a.stage == "place_fill")
        .count();
    for count in [geostatic_count + 1, geostatic_count + birth_count + 1] {
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        assert!(!partial.complete);
        let checkpoint = Checkpoint::seal(&mesh, partial)?;
        let resumed = checkpoint.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert!(resumed.complete);
        assert_eq!(
            serde_json::to_vec(&full.last_committed)?,
            serde_json::to_vec(&resumed.last_committed)?
        );
        assert_eq!(
            serde_json::to_vec(&full.frames)?,
            serde_json::to_vec(&resumed.frames)?
        );
        assert_eq!(
            serde_json::to_vec(&full.transition_audits)?,
            serde_json::to_vec(&resumed.transition_audits)?
        );
    }
    Ok(())
}

#[test]
fn inconsistent_replacement_born_active_wet_or_aged_material_is_rejected() -> Result<()> {
    let (mesh, ground, case) = fixture(false)?;
    let mut missing = case.clone();
    if let stages::Stage::Birth { states, .. } = &mut missing.stages[0] {
        states.clear();
    }
    assert!(stages::solve(&mesh, ground.clone(), missing).is_err());
    let mut active = case.clone();
    active.initially_inactive_bodies.clear();
    assert!(stages::solve(&mesh, ground.clone(), active).is_err());
    let mut aged = case.clone();
    if let stages::Stage::Birth { states, .. } = &mut aged.stages[0]
        && let soil::State::MohrCoulomb(s) = states.get_mut("UpperClay").unwrap()
    {
        s.equivalent_plastic_strain = 0.1;
    }
    assert!(stages::solve(&mesh, ground.clone(), aged).is_err());
    let mut changed = case.clone();
    if let stages::Stage::MaterialSet { material_id, .. } = &mut changed.stages[1] {
        *material_id = 99;
    }
    assert!(stages::solve(&mesh, ground.clone(), changed).is_err());
    let mut wet = ground.clone();
    wet.initialization.water_level_m = Some(2.);
    assert!(
        stages::solve(&mesh, wet, case.clone())
            .unwrap_err()
            .to_string()
            .contains("dry")
    );
    let mut cp = Checkpoint::seal(
        &mesh,
        stages::solve_until(&mesh, ground.clone(), case.clone(), 2)?,
    )?;
    cp.run.last_committed.element_material_ids[0] = 99;
    cp.run_sha256 = checkpoint::hash_json(&cp.run)?;
    assert!(cp.resume(&mesh, ground, case, None).is_err());
    Ok(())
}

#[test]
fn material_density_change_rebuilds_gravity_and_matches_layered_column_theory() -> Result<()> {
    let (mesh, mut ground, case) = fixture(false)?;
    ground.materials[2].unit_weight_kn_m3 = 200.;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let difference = mesh.node_sets["top"]
        .iter()
        .map(|&n| {
            run.frames[2].state.displacement_m[n][2] - run.frames[1].state.displacement_m[n][2]
        })
        .sum::<f64>()
        / mesh.node_sets["top"].len() as f64;
    let constrained = |e: f64| e * 0.9 / (1.1 * 0.8);
    let expected = -25. / constrained(10000.) - 12.5 / constrained(40000.);
    assert!(
        (difference - expected).abs() < 1e-9,
        "{difference} != {expected}"
    );
    assert!(
        run.transition_audits[1]
            .initial_holding_load_kn
            .iter()
            .any(|x| x.abs() > 1e-4)
    );
    Ok(())
}

#[test]
fn mcc_foundation_accepts_mc_fill_and_rejects_implicit_history_conversion() -> Result<()> {
    let (mesh, mut ground, case) = fixture(false)?;
    let template: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    ground.materials[0].model = template.materials[0].model.clone();
    if let (
        soil_fem_3d::geostatic::Seed::BodyStates { states },
        soil_fem_3d::geostatic::Seed::BodyStates { states: other },
    ) = (
        &mut ground.initialization.geostatic.as_mut().unwrap().seed,
        &template.initialization.geostatic.as_ref().unwrap().seed,
    ) {
        states.insert("LowerClay".into(), other["LowerClay"].clone());
    }
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.frames
            .last()
            .unwrap()
            .state
            .material_states
            .iter()
            .any(|s| s.as_mcc().is_some())
    );
    let mut bad = case;
    bad.stages.insert(0,serde_json::from_value(json!({"kind":"material_set","name":"bad_conversion","bodies":["LowerClay"],"material_id":2,"history":{"mode":"preserve"}}))?);
    assert!(
        stages::solve(&mesh, ground, bad)
            .unwrap_err()
            .to_string()
            .contains("type")
    );
    Ok(())
}

#[test]
fn rebirth_placement_distinguishes_cad_nodes_from_archived_excavated_displacements() -> Result<()> {
    let (mesh, ground, mut case) = fixture(false)?;
    case.stages.push(serde_json::from_value(json!({"kind":"excavate","name":"remove","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.}))?);
    let mut birth = case.stages[0].clone();
    if let stages::Stage::Birth { name, .. } = &mut birth {
        *name = "rebirth".into();
    }
    case.stages.push(birth);
    let cad = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(cad.complete, "{:?}", cad.failure_reason);
    if let stages::Stage::Birth { placement, .. } = case.stages.last_mut().unwrap() {
        *placement = stages::BirthPlacement::CurrentDisplacement;
    }
    let archived = stages::solve(&mesh, ground, case)?;
    assert!(archived.complete, "{:?}", archived.failure_reason);
    let cad_origin = &cad.transition_audits.last().unwrap().prepared_state;
    let old_origin = &archived.transition_audits.last().unwrap().prepared_state;
    for &n in &mesh.node_sets["top"] {
        assert_eq!(cad_origin.displacement_m[n], [0.; 3]);
        assert!(old_origin.displacement_m[n][2] < -1e-4);
        let difference =
            archived.last_committed.displacement_m[n][2] - cad.last_committed.displacement_m[n][2];
        assert!((difference - old_origin.displacement_m[n][2]).abs() < 1e-9);
    }
    assert!(
        (cad.frames.last().unwrap().monitors["surface"]
            .stage_increment_z_m
            .unwrap()
            - archived.frames.last().unwrap().monitors["surface"]
                .stage_increment_z_m
                .unwrap())
        .abs()
            < 1e-9
    );
    Ok(())
}
