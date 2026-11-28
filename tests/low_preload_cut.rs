use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    Ok((
        serde_json::from_str(include_str!("fixtures/low-preload-nc/mesh.json"))?,
        serde_json::from_str(include_str!("fixtures/low-preload-nc/ground.json"))?,
        serde_json::from_str(include_str!("fixtures/low-preload-nc/stages.json"))?,
    ))
}

fn deterministic(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::Object(map) => {
            map.remove("elapsed_seconds");
            for child in map.values_mut() {
                deterministic(child);
            }
        }
        serde_json::Value::Array(items) => items.iter_mut().for_each(deterministic),
        _ => {}
    }
}

#[test]
fn low_preload_nc_partial_cut_freezes_removed_plastic_and_contact_history() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let run = stages::solve_until(&mesh, ground, case, 15)?;
    assert!(!run.complete);
    assert_eq!(run.last_committed.stage_index, 3);
    assert!(run.last_committed.stage_progress > 0. && run.last_committed.stage_progress < 1.);
    let unload = &run
        .frames
        .iter()
        .find(|f| f.name == "unload")
        .unwrap()
        .state;
    assert_eq!(unload.joints, run.last_committed.joints);
    let mut removed = 0;
    let mut plastic = 0;
    for (i, point) in run.integration_points.iter().enumerate() {
        let state = &run.last_committed.material_states[i];
        if !run.last_committed.active_elements[point.element_index] {
            assert_eq!(
                serde_json::to_value(state)?,
                serde_json::to_value(&unload.material_states[i])?
            );
            removed += 1;
        }
        let value = serde_json::to_value(state)?;
        if value["plastic_strain_compression_positive"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v.as_f64().unwrap() != 0.)
        {
            plastic += 1;
        }
    }
    assert_eq!(removed, 162);
    assert_eq!(plastic, 324);
    assert!(
        run.last_committed
            .excavation_release_load_kn
            .iter()
            .any(|v| *v != 0.)
    );
    assert!(run.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
    assert!(
        run.attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    Ok(())
}

#[test]
fn low_preload_nc_cut_checkpoint_resumes_exact_prefix_and_rejects_preload_change() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 15)?;
    let saved = Checkpoint::seal(&mesh, partial)?;
    let mut wrong = serde_json::to_value(&case)?;
    wrong["boundary"]["surface_loads"][0]["traction_kpa"][2] = (-100.).into();
    assert!(
        saved
            .verify(&mesh, &ground, &serde_json::from_value(wrong)?)
            .is_err()
    );
    let resumed = saved.resume(&mesh, ground.clone(), case.clone(), Some(1))?;
    let direct = stages::solve_until(&mesh, ground, case, 16)?;
    let mut a = serde_json::to_value(resumed)?;
    let mut b = serde_json::to_value(direct)?;
    deterministic(&mut a);
    deterministic(&mut b);
    assert_eq!(a, b);
    Ok(())
}
