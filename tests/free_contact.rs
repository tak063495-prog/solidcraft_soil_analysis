use nalgebra::{Rotation3, Vector3};
use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(mcc: bool, angle: f64) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/joint-free-mesh.json"))?;
    let ground = if mcc {
        include_str!("../examples/ground-joint-mcc.json")
    } else {
        include_str!("../examples/ground-free-contact-mc.json")
    };
    let mut ground: serde_json::Value = serde_json::from_str(ground)?;
    ground["initialization"]["geostatic"]["verification_reject_once"] = true.into();
    if mcc {
        ground["initialization"]["geostatic"]["base_loads"]["surface_loads"]
            .as_array_mut()
            .unwrap()
            .push(
                serde_json::json!({"kind":"traction","face_set":"top","traction_kpa":[1.,0.5,0.]}),
            );
    }
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-free-contact.json"))?;
    let r = Rotation3::from_axis_angle(&Vector3::y_axis(), angle.to_radians());
    for p in &mut mesh.nodes {
        *p = (r * Vector3::from(*p)).into();
    }
    for region in &mut mesh.regions {
        let mut bounds = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
        for (ei, id) in mesh.region_ids.iter().enumerate() {
            if *id == region.id {
                let nodes = if mesh.tetrahedra10.is_empty() {
                    mesh.tetrahedra[ei].to_vec()
                } else {
                    mesh.tetrahedra10[ei].to_vec()
                };
                for n in nodes {
                    for (a, &coordinate) in mesh.nodes[n].iter().enumerate() {
                        bounds[0][a] = bounds[0][a].min(coordinate);
                        bounds[1][a] = bounds[1][a].max(coordinate);
                    }
                }
            }
        }
        region.bounds_m = bounds;
    }
    Ok((mesh, serde_json::from_value(ground)?, case))
}

#[test]
fn free_upper_block_closes_without_physical_pins_and_restores_predictor_history() -> Result<()> {
    for (mcc, angle) in [(false, 0.), (false, 10.), (true, 0.)] {
        let (mesh, ground, case) = fixture(mcc, angle)?;
        assert!(
            case.boundary
                .as_ref()
                .unwrap()
                .constraints
                .iter()
                .all(|c| c.node_set == "bottom")
        );
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        eprintln!(
            "cold block mcc={mcc}, angle={angle}, complete={}, geostatic attempts={}",
            full.complete,
            full.geostatic.as_ref().map_or(0, |g| g.attempts.len())
        );
        assert!(
            full.complete,
            "mcc={mcc}, angle={angle}: {:?}",
            full.failure_reason
        );
        assert!(
            full.frames
                .iter()
                .all(|f| f.joint_support.as_ref().unwrap().unrestrained_modes == 0)
        );
        let geo = full.geostatic.as_ref().unwrap();
        assert!(
            geo.attempts
                .iter()
                .any(|a| a.contact_predictors.iter().any(|p| p.activation_predictor))
        );
        for p in geo.attempts.iter().flat_map(|a| &a.contact_predictors) {
            if let Some(after) = p.full_relative_force_after {
                assert!(after < p.full_relative_force_before);
            }
            assert!(
                p.temporary_linear_solve.true_relative_residual
                    < case.settings.relative_linear_tolerance
            );
        }
        assert!(geo.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_structure_changed
            && !a.contact_predictors.is_empty()
            && a.committed_state_unchanged_on_rejection));
        let mean_x = mesh.node_sets["top"]
            .iter()
            .map(|&n| full.frames[0].state.displacement_m[n][0])
            .sum::<f64>()
            / mesh.node_sets["top"].len() as f64;
        assert!(mean_x.abs() > 1e-7, "upper block was artificially clamped");
        assert!(full.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
        assert!(
            full.frames[2].state.joints == full.frames[3].state.joints,
            "removed contact history changed"
        );
        let mut clean = case.clone();
        clean.settings.initial_increment = geo
            .attempts
            .iter()
            .find(|a| a.accepted)
            .unwrap()
            .progress_to;
        clean.verification_reject_once_stages.clear();
        let mut stages = serde_json::to_value(&clean.stages)?;
        for s in stages.as_array_mut().unwrap() {
            let name = s["name"].as_str().unwrap();
            s["initial_increment"] = serde_json::json!(
                full.attempts
                    .iter()
                    .find(|a| a.accepted && a.stage == name)
                    .unwrap()
                    .progress_to
            );
        }
        clean.stages = serde_json::from_value(stages)?;
        let mut clean_ground = serde_json::to_value(&ground)?;
        clean_ground["initialization"]["geostatic"]["verification_reject_once"] = false.into();
        let fresh = stages::solve(&mesh, serde_json::from_value(clean_ground)?, clean)?;
        assert!(fresh.complete, "{:?}", fresh.failure_reason);
        assert!(
            fresh
                .geostatic
                .as_ref()
                .unwrap()
                .attempts
                .iter()
                .filter(|a| a.accepted)
                .map(|a| a.progress_to)
                .eq(geo
                    .attempts
                    .iter()
                    .filter(|a| a.accepted)
                    .map(|a| a.progress_to)),
            "geostatic accepted paths differ"
        );
        assert!(
            fresh
                .attempts
                .iter()
                .filter(|a| a.accepted)
                .map(|a| (&a.stage, a.progress_to))
                .eq(full
                    .attempts
                    .iter()
                    .filter(|a| a.accepted)
                    .map(|a| (&a.stage, a.progress_to))),
            "construction accepted paths differ"
        );
        assert!(
            serde_json::to_value(&fresh.last_committed)?
                == serde_json::to_value(&full.last_committed)?,
            "mcc={mcc} angle={angle}: rejected trial changed committed history"
        );
        assert!(
            serde_json::to_value(&fresh.frames)? == serde_json::to_value(&full.frames)?,
            "mcc={mcc} angle={angle}: clean frames differ"
        );
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
        let saved: Checkpoint =
            serde_json::from_slice(&serde_json::to_vec(&Checkpoint::seal(&mesh, partial)?)?)?;
        let resumed = saved.resume(&mesh, ground, case, None)?;
        assert!(
            serde_json::to_value(&resumed.frames)? == serde_json::to_value(&full.frames)?,
            "mcc={mcc} angle={angle}: resumed frames differ"
        );
    }
    Ok(())
}

#[test]
fn free_frictionless_or_overcapacity_blocks_never_commit_a_false_equilibrium() -> Result<()> {
    for (mu, force) in [
        (0., [0., 0., 0.]),
        (0.3, [5., 0., 0.]),
        (0.3, [3., 0., 0.]),
        (0.3, [0., 0., 15.]),
    ] {
        let (mesh, ground, case) = fixture(false, 0.)?;
        let mut g = serde_json::to_value(ground)?;
        g["initialization"]["geostatic"]["base_loads"] = serde_json::json!({"surface_loads":[{"kind":"traction","face_set":"top","traction_kpa":force}]});
        let mut c = serde_json::to_value(case)?;
        c["joint_system"]["joints"][0]["material"]["friction_coefficient"] = mu.into();
        let run = stages::solve(
            &mesh,
            serde_json::from_value(g)?,
            serde_json::from_value(c)?,
        )?;
        assert!(!run.complete, "mu={mu}, force={force:?}");
        assert!(run.frames.is_empty());
        let geo = run.geostatic.as_ref().unwrap();
        assert!(
            geo.attempts
                .iter()
                .all(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
        );
        assert_eq!(
            serde_json::to_value(&run.last_committed)?,
            serde_json::to_value(&geo.seed)?
        );
        if mu == 0. {
            assert!(
                geo.attempts
                    .iter()
                    .any(|a| a.reason.contains("JOINT rigid mechanism"))
            );
        }
    }
    Ok(())
}
