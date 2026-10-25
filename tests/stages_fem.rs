use soil_fem_3d::{
    Result,
    ground::{self, Ground},
    mesh::{Mesh, Region},
    stages::{self, Kinematics, Stage, StageCase},
};
use std::collections::BTreeMap;

fn column() -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh = Mesh {
        schema: "soil-fem-mesh/2".into(),
        units: "m-kN-kPa".into(),
        axis: "z_up".into(),
        crs: "LOCAL".into(),
        origin_m: [0.; 3],
        nodes: vec![],
        tetrahedra: vec![],
        tetrahedra10: vec![],
        quadrature_order: None,
        face_sets: BTreeMap::new(),
        material_ids: vec![],
        node_sets: BTreeMap::new(),
        regions: vec![],
        region_ids: vec![],
        excavation_sets: BTreeMap::new(),
    };
    for z in 0..3 {
        for y in 0..2 {
            for x in 0..2 {
                mesh.nodes.push([x as f64, y as f64, 5. * z as f64]);
            }
        }
    }
    for z in 0usize..2 {
        let c = [
            4 * z,
            4 * z + 1,
            4 * z + 3,
            4 * z + 2,
            4 * z + 4,
            4 * z + 5,
            4 * z + 7,
            4 * z + 6,
        ];
        for t in [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ] {
            mesh.tetrahedra.push(t.map(|i| c[i]));
            mesh.material_ids.push(z as u32 + 1);
            mesh.region_ids.push(z as u32 + 1);
        }
        mesh.regions.push(Region {
            id: z as u32 + 1,
            body: if z == 0 { "LowerClay" } else { "UpperClay" }.into(),
            material_id: z as u32 + 1,
            bounds_m: [[0., 0., 5. * z as f64], [1., 1., 5. * (z + 1) as f64]],
            volume_m3: 5.,
            source_sha256: "a".repeat(64),
        });
    }
    let ground = serde_json::from_str(include_str!("../examples/ground.json"))?;
    let mut case: StageCase = serde_json::from_str(include_str!("../examples/stages.json"))?;
    case.settings.kinematics = Kinematics::Oedometer;
    Ok((mesh, ground, case))
}

#[test]
fn sparse_assembly_sums_tied_plane_dofs_and_matches_dense_oedometer() -> Result<()> {
    use soil_fem_3d::linalg::{GmresOptions, LinearSolver};
    let (mesh, ground, mut case) = column()?;
    let dense = stages::solve(&mesh, ground.clone(), case.clone())?;
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions::default(),
    };
    let sparse = stages::solve(&mesh, ground, case)?;
    assert!(
        dense.complete && sparse.complete,
        "{:?} {:?}",
        dense.failure_reason,
        sparse.failure_reason
    );
    for (a, b) in dense
        .last_committed
        .displacement_m
        .iter()
        .flatten()
        .zip(sparse.last_committed.displacement_m.iter().flatten())
    {
        assert!((a - b).abs() < 1e-9);
    }
    assert!(
        sparse
            .frames
            .iter()
            .all(|f| f.relative_force_residual < 1e-8)
    );
    Ok(())
}

#[test]
fn no_load_does_not_reload_initial_gravity() -> Result<()> {
    let (mesh, ground, mut case) = column()?;
    case.stages = vec![Stage::Pressure {
        name: "no_change".into(),
        target_pressure_kpa: 0.,
        initial_increment: None,
    }];
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.last_committed
            .displacement_m
            .iter()
            .flatten()
            .all(|u| *u == 0.)
    );
    assert_eq!(
        serde_json::to_string(&run.frames[0].state.material_states)?,
        serde_json::to_string(&run.last_committed.material_states)?
    );
    Ok(())
}

#[test]
fn tiny_elastic_pressure_matches_constrained_modulus_theory() -> Result<()> {
    let (mesh, mut ground, mut case) = column()?;
    ground.initialization.water_level_m = None;
    ground.initialization.surface_surcharge_kpa = 100.;
    for layer in &mut ground.layers {
        layer.consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 100. };
    }
    let initial = ground::initialize(&mesh, ground.clone())?;
    let mut c = [0.; 2];
    for point in &initial.point_states {
        let m = &ground.materials[(point.material_id - 1) as usize];
        let p = point
            .constitutive_state
            .as_mcc()
            .unwrap()
            .effective_stress_compression_positive_kpa[..3]
            .iter()
            .sum::<f64>()
            / 3.;
        let k = (1. + point.constitutive_state.as_mcc().unwrap().void_ratio) * p
            / m.model.as_mcc().unwrap().kappa;
        c[(point.material_id - 1) as usize] += 3. * k * (1. - m.model.as_mcc().unwrap().poisson)
            / (1. + m.model.as_mcc().unwrap().poisson)
            / 6.;
    }
    let dq = 1e-4;
    case.stages = vec![Stage::Pressure {
        name: "elastic_increment".into(),
        target_pressure_kpa: 100. + dq,
        initial_increment: Some(1.),
    }];
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let actual = run.frames[1].monitors["ground_surface"]
        .mean_displacement_z_m
        .unwrap_or(0.);
    let expected = -dq * (5. / c[0] + 5. / c[1]);
    println!(
        "elastic theory: u={actual:.12e} m; expected={expected:.12e} m; relative difference={:.3e}",
        (actual / expected - 1.).abs()
    );
    assert!(
        (actual / expected - 1.).abs() < 1e-3,
        "actual={actual:e}, theory={expected:e}"
    );
    Ok(())
}

#[test]
fn preload_unload_and_excavation_preserve_plasticity_and_heave() -> Result<()> {
    let (mesh, ground, case) = column()?;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.frames.len(), 4);
    assert!(
        run.frames[1].monitors["ground_surface"]
            .mean_displacement_z_m
            .unwrap_or(0.)
            < 0.
    );
    assert!(
        run.frames[2].monitors["ground_surface"]
            .stage_increment_z_m
            .unwrap_or(0.)
            > 0.
    );
    assert!(
        run.frames[3].monitors["excavation_floor"]
            .stage_increment_z_m
            .unwrap_or(0.)
            > 0.
    );
    assert_eq!(
        run.frames[3].monitors["ground_surface"].active_node_count,
        0
    );
    assert!(
        run.frames[2]
            .state
            .material_states
            .iter()
            .any(|s| { s.plastic_strain().iter().any(|v| v.abs() > 1e-5) })
    );
    for i in 6..12 {
        assert!(!run.last_committed.active_elements[i]);
        assert_eq!(
            serde_json::to_string(&run.frames[2].state.material_states[i])?,
            serde_json::to_string(&run.last_committed.material_states[i])?
        );
    }
    for frame in &run.frames {
        assert!(frame.relative_force_residual < 1e-8);
        assert!(frame.maximum_normalized_yield <= 1e-8);
        assert!(frame.force_balance_error_kn.iter().all(|f| f.abs() < 1e-5));
    }
    Ok(())
}

#[test]
fn rejected_plastic_and_excavation_trials_restore_entire_committed_history() -> Result<()> {
    let (mesh, ground, mut injected) = column()?;
    injected.settings.kinematics = Kinematics::SideRollers;
    injected.verification_reject_once_stages =
        vec!["preload".into(), "excavate_upper_layer".into()];
    let mut clean = injected.clone();
    clean.verification_reject_once_stages.clear();
    for s in &mut clean.stages {
        match s {
            Stage::Pressure {
                name,
                initial_increment,
                ..
            } if name == "preload" => *initial_increment = Some(0.125),
            Stage::Excavate {
                initial_increment, ..
            } => *initial_increment = Some(0.125),
            _ => (),
        }
    }
    let a = stages::solve(&mesh, ground.clone(), injected)?;
    let b = stages::solve(&mesh, ground, clean)?;
    assert!(
        a.complete && b.complete,
        "{:?} {:?}",
        a.failure_reason,
        b.failure_reason
    );
    let rejects: Vec<_> = a.attempts.iter().filter(|a| !a.accepted).collect();
    assert_eq!(rejects.len(), 2);
    assert!(
        rejects
            .iter()
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    assert!(rejects.iter().any(|a| a.rejected_trial_material_changed));
    assert!(rejects.iter().any(|a| a.rejected_trial_active_set_changed));
    assert!(rejects.iter().all(|a| a.rejected_trial_load_changed));
    assert_eq!(
        serde_json::to_string(&a.last_committed)?,
        serde_json::to_string(&b.last_committed)?
    );
    Ok(())
}

#[test]
fn real_nonconvergence_cutback_and_exhaustion_are_not_success() -> Result<()> {
    let (mesh, ground, mut case) = column()?;
    case.settings.maximum_newton_iterations = 6;
    case.settings.initial_increment = 1.;
    case.settings.minimum_increment = 1. / 512.;
    let recovered = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(recovered.complete, "{:?}", recovered.failure_reason);
    assert!(
        recovered
            .attempts
            .iter()
            .any(|a| !a.accepted && a.reason.contains("Newton iteration limit"))
    );
    assert!(
        recovered
            .attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    case.stages = vec![case.stages[2].clone()];
    case.settings.maximum_newton_iterations = 1;
    case.settings.initial_increment = 0.25;
    case.settings.minimum_increment = 0.25;
    let failed = stages::solve(&mesh, ground, case)?;
    assert!(!failed.complete);
    assert_eq!(failed.frames.len(), 1);
    assert!(
        failed
            .failure_reason
            .as_deref()
            .unwrap_or("")
            .contains("cutback exhausted")
    );
    assert_eq!(
        serde_json::to_string(&failed.frames[0].state)?,
        serde_json::to_string(&failed.last_committed)?
    );
    Ok(())
}
