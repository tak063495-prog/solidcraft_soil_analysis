use soil_fem_3d::{
    Result,
    checkpoint::{self, Checkpoint},
    geostatic::{self, PointSeed, Seed},
    ground::Ground,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(mcc: bool, quadratic: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let top = soil_fem_3d::boundary::surface_owners(&mesh)
        .into_iter()
        .filter(|(face, owners)| owners.len() == 1 && face.iter().all(|&i| mesh.nodes[i][2] == 1.))
        .map(|(face, _)| face)
        .collect();
    mesh.face_sets.insert("top".into(), top);
    let ground: Ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-geostatic-mcc-column.json")
    } else {
        include_str!("../examples/ground-geostatic-mc-column.json")
    })?;
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    Ok((
        if quadratic {
            mesh.elevate_straight()?
        } else {
            mesh
        },
        ground,
        case,
    ))
}

#[test]
fn mc_and_mcc_form_plastic_geostatic_history_and_keep_a_separate_construction_reference()
-> Result<()> {
    for mcc in [false, true] {
        let (mesh, ground, case) = fixture(mcc, true)?;
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(run.schema, "soil-fem-stage-result/4");
        let g = run.geostatic.as_ref().unwrap();
        assert!(g.complete && g.final_relative_physical_force_residual.unwrap() < 1e-9);
        assert!(
            g.last_committed
                .material_states
                .iter()
                .any(|p| p.plastic_strain().iter().any(|v| v.abs() > 1e-7))
        );
        if mcc {
            assert!(
                g.last_committed.material_states.iter().any(|p| p
                    .as_mcc()
                    .unwrap()
                    .preconsolidation_pressure_kpa
                    > 140.01)
            );
        }
        assert!(
            run.reference_displacement_m
                .iter()
                .flatten()
                .any(|v| v.abs() > 1e-5)
        );
        assert!(
            run.frames[0]
                .maximum_absolute_engineering_strain_from_k0
                .is_none()
        );
        assert!(
            run.frames[0]
                .maximum_absolute_engineering_strain_from_seed
                .unwrap()
                > 1e-5
        );
        assert_eq!(
            run.frames[0].maximum_absolute_engineering_strain_from_initial_equilibrium,
            Some(0.)
        );
        assert_eq!(
            run.frames[0].monitors["surface"].mean_displacement_z_m,
            Some(0.)
        );
        let last = &run.last_committed;
        let nodes = &mesh.node_sets["top"];
        let expected = nodes
            .iter()
            .map(|&n| last.displacement_m[n][2] - run.reference_displacement_m[n][2])
            .sum::<f64>()
            / nodes.len() as f64;
        assert!(
            (run.frames.last().unwrap().monitors["surface"]
                .mean_displacement_z_m
                .unwrap()
                - expected)
                .abs()
                < 1e-14
        );
        assert!(g.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_plastic_changed
            && a.committed_state_unchanged_on_rejection));
    }
    Ok(())
}
#[test]
fn geostatic_partial_boundary_and_construction_checkpoints_restore_identically() -> Result<()> {
    let (mesh, ground, case) = fixture(true, false)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    let accepted_geo = full
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    for count in [1, accepted_geo, accepted_geo + 1] {
        let saved = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        assert!(!saved.complete);
        let sealed = Checkpoint::seal(&mesh, saved)?;
        let encoded = serde_json::to_vec(&sealed)?;
        let cp: Checkpoint = serde_json::from_slice(&encoded)?;
        if count == 1 {
            let mut tampered: Checkpoint = serde_json::from_slice(&encoded)?;
            tampered
                .run
                .geostatic
                .as_mut()
                .unwrap()
                .initial_holding_load_kn[0] += 1.;
            tampered.run_sha256 = checkpoint::hash_json(&tampered.run)?;
            assert!(
                tampered
                    .resume(&mesh, ground.clone(), case.clone(), None)
                    .is_err()
            );
        }
        let resumed = cp.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert!(resumed.complete, "{:?}", resumed.failure_reason);
        assert_eq!(
            serde_json::to_vec(&full.last_committed)?,
            serde_json::to_vec(&resumed.last_committed)?
        );
        assert_eq!(
            serde_json::to_vec(&full.frames)?,
            serde_json::to_vec(&resumed.frames)?
        );
        let normalize = |mut value: serde_json::Value| {
            for attempt in value["attempts"].as_array_mut().unwrap() {
                for linear in attempt["linear_solves"].as_array_mut().unwrap() {
                    linear["elapsed_seconds"] = serde_json::json!(0.);
                }
            }
            value
        };
        // Wall-clock timings are observations, not deterministic histories.
        assert!(
            normalize(serde_json::to_value(&full.geostatic)?)
                == normalize(serde_json::to_value(&resumed.geostatic)?),
            "geostatic state/path differs at checkpoint {count}"
        );
    }
    Ok(())
}
#[test]
fn integration_point_seed_import_checks_mesh_position_material_and_body_mapping() -> Result<()> {
    let (mesh, mut ground, _) = fixture(true, true)?;
    let states = geostatic::states(&mesh, &ground)?;
    let layout = mesh.integration_layout()?;
    let points = layout
        .iter()
        .zip(states.clone())
        .map(|(q, state)| PointSeed {
            element_index: q.element_index,
            local_index: q.local_index,
            position_m: q.position_m,
            material_id: mesh.material_ids[q.element_index],
            state,
        })
        .collect();
    ground.initialization.geostatic.as_mut().unwrap().seed = Seed::IntegrationPoints {
        mesh_sha256: checkpoint::hash_json(&mesh)?,
        points,
    };
    assert_eq!(
        serde_json::to_vec(&states)?,
        serde_json::to_vec(&geostatic::states(&mesh, &ground)?)?
    );
    let mut bad = mesh.clone();
    bad.crs = "CHANGED".into();
    assert!(geostatic::states(&bad, &ground).is_err());
    if let Seed::IntegrationPoints { points, .. } =
        &mut ground.initialization.geostatic.as_mut().unwrap().seed
    {
        points[0].position_m[0] += 0.01;
    }
    assert!(geostatic::states(&mesh, &ground).is_err());
    let (_, mut bad, _) = fixture(false, false)?;
    bad.bodies[0].body = "Wrong".into();
    assert!(geostatic::states(&mesh, &bad).is_err());
    Ok(())
}
#[test]
fn excessive_initial_strain_preserves_seed_and_never_creates_an_equilibrated_frame() -> Result<()> {
    let (mesh, mut ground, case) = fixture(false, false)?;
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .maximum_initial_engineering_strain = 1e-8;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(!run.complete && run.frames.is_empty());
    let g = run.geostatic.unwrap();
    assert!(!g.complete && g.progress == 0.);
    assert!(g.failure_reason.unwrap().contains("strain"));
    assert_eq!(
        serde_json::to_vec(&g.seed)?,
        serde_json::to_vec(&g.last_committed)?
    );
    assert!(
        g.attempts
            .iter()
            .all(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    Ok(())
}

#[test]
fn saturated_general_column_matches_effective_weight_and_rejects_internal_water_faces() -> Result<()>
{
    let (mesh, mut ground, mut case) = fixture(false, true)?;
    ground.initialization.water_level_m = Some(1.);
    let cfg = ground.initialization.geostatic.as_mut().unwrap();
    cfg.verification_reject_once = false;
    cfg.hydrostatic_face_sets = vec!["top".into()];
    for m in &mut ground.materials {
        if let soil_fem_3d::soil::Model::MohrCoulomb(p) = &mut m.model {
            p.cohesion_kpa = 1000.;
        }
    }
    case.stages.clear();
    let run = stages::solve(&mesh, ground.clone(), case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for (p, state) in run
        .integration_points
        .iter()
        .zip(&run.last_committed.material_states)
    {
        let expected =
            (100. - ground.initialization.water_unit_weight_kn_m3) * (1. - p.position_m[2]);
        assert!(
            (state.stress()[2] - expected).abs() < 1e-5,
            "computed {:?}, expected {expected}",
            state.stress()
        );
    }
    let mut bad = mesh;
    let internal = soil_fem_3d::boundary::surface_owners(&bad)
        .into_iter()
        .find(|(_, o)| o.len() == 2)
        .unwrap()
        .0;
    bad.face_sets.insert("internal".into(), vec![internal]);
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .hydrostatic_face_sets = vec!["internal".into()];
    assert!(geostatic::states(&bad, &ground).is_err());
    Ok(())
}

#[test]
fn sloping_exposed_surface_allows_excavation_without_a_horizontal_top_facet() -> Result<()> {
    let (mut mesh, mut ground, mut case) = fixture(false, false)?;
    for p in &mut mesh.nodes {
        p[2] *= 1. - 0.25 * p[0];
    }
    let layout = mesh.integration_layout()?;
    for region in &mut mesh.regions {
        region.volume_m3 = layout
            .iter()
            .filter(|q| mesh.region_ids[q.element_index] == region.id)
            .map(|q| q.weight_m3)
            .sum();
        let nodes: Vec<_> = mesh
            .tetrahedra
            .iter()
            .enumerate()
            .filter(|(i, _)| mesh.region_ids[*i] == region.id)
            .flat_map(|(_, t)| t.iter().map(|&n| mesh.nodes[n]))
            .collect();
        for a in 0..3 {
            region.bounds_m[0][a] = nodes.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min);
            region.bounds_m[1][a] = nodes.iter().map(|p| p[a]).fold(f64::NEG_INFINITY, f64::max);
        }
    }
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    for m in &mut ground.materials {
        if let soil_fem_3d::soil::Model::MohrCoulomb(p) = &mut m.model {
            p.cohesion_kpa = 1000.;
        }
    }
    case.stages = vec![soil_fem_3d::stages::Stage::Excavate {
        name: "sloping_excavation".into(),
        bodies: vec!["UpperClay".into()],
        remaining_surface_pressure_kpa: 0.,
        initial_increment: None,
    }];
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.last_committed.active_elements.iter().any(|a| !a));
    assert!(
        run.last_committed
            .excavation_release_load_kn
            .iter()
            .all(|x| *x == 0.)
    );
    assert!(run.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
    Ok(())
}
