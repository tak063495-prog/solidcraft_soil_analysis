use soil_fem_3d::{
    Result,
    boundary::{Conditions, Constraint},
    checkpoint::Checkpoint,
    ground::{self, Ground},
    linalg::{GmresOptions, LinearSolver},
    mcc,
    mesh::Mesh,
    stages::{self, Kinematics, Stage, StageCase},
};

fn constraint(set: &str, axis: usize, value: f64) -> Constraint {
    let mut direction = [0.; 3];
    direction[axis] = 1.;
    Constraint {
        node_set: set.into(),
        direction,
        value_m: value,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    }
}

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    let mesh: Mesh =
        serde_json::from_str(include_str!("../examples/hardening-column-tet10-mesh.json"))?;
    let mut ground: Ground = serde_json::from_str(include_str!("../examples/ground.json"))?;
    ground.initialization.water_level_m = None;
    ground.initialization.surface_surcharge_kpa = 100.;
    for material in &mut ground.materials {
        material.unit_weight_kn_m3 = 1e-8;
        if let soil_fem_3d::soil::Model::Mcc(p) = &mut material.model {
            *p = mcc::Parameters {
                critical_state_m: 1.2,
                lambda: 0.2,
                kappa: 0.04,
                poisson: 0.3,
                reference_pressure_kpa: 1.,
                specific_volume_ncl_at_reference: 3.5,
            };
        }
    }
    for layer in &mut ground.layers {
        layer.k0_x = 1.;
        layer.k0_y = 1.;
        layer.consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 1. };
    }
    let mut case: StageCase = serde_json::from_str(include_str!("../examples/stages.json"))?;
    case.schema = "soil-fem-mixed-mcc-stages/1".into();
    case.settings.kinematics = Kinematics::Explicit;
    case.settings.maximum_free_dofs = 2000;
    case.settings.maximum_newton_iterations = 40;
    case.settings.initial_increment = 1.;
    case.settings.minimum_increment = 1. / 256.;
    case.settings.relative_force_tolerance = 1e-9;
    case.settings.relative_linear_tolerance = 1e-11;
    case.settings.material_maximum_substep = 2e-4;
    case.boundary = Some(Conditions {
        constraints: vec![
            constraint("all", 0, 0.),
            constraint("all", 1, 0.),
            constraint("bottom", 2, 0.),
            constraint("top", 2, -0.005),
        ],
        ..Conditions::default()
    });
    case.monitors.clear();
    case.monitor_sets.insert("top".into(), "top".into());
    case.stages = vec![
        Stage::Boundary {
            name: "load".into(),
            target_factor: 1.,
            initial_increment: None,
        },
        Stage::Boundary {
            name: "unload".into(),
            target_factor: 0.8,
            initial_increment: None,
        },
    ];
    Ok((mesh, ground, case))
}

#[test]
fn plastic_affine_patch_satisfies_material_and_volume_equations() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let run = stages::solve(&mesh, ground.clone(), case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.schema, "soil-fem-mixed-mcc-stage-result/1");
    assert_eq!(run.attempts.len(), 2);
    let loaded = &run.frames[1].state;
    let mixed = loaded.mixed_mcc.as_ref().unwrap();
    let maximum_patch_error = loaded
        .displacement_m
        .iter()
        .zip(&mesh.nodes)
        .map(|(u, x)| (u[2] + 0.005 * x[2]).abs())
        .fold(0., f64::max);
    println!(
        "patch maximum displacement error={maximum_patch_error:e}; partitions={:?}",
        mixed
            .last_increment_partitions
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
    );
    for (u, x) in loaded.displacement_m.iter().zip(&mesh.nodes) {
        assert!((u[2] + 0.005 * x[2]).abs() < 1e-10);
    }
    for (i, s) in loaded.material_states.iter().enumerate() {
        let old = run.frames[0].state.material_states[i].as_mcc().unwrap();
        let point = &run.integration_points[i];
        let material = ground.materials[(mesh.material_ids[point.element_index] - 1) as usize]
            .model
            .as_mcc()
            .unwrap();
        let direct = material.integrate_partition(
            old,
            [0., 0., 0.005, 0., 0., 0.],
            mixed.last_increment_partitions[i],
        )?;
        let actual = s.as_mcc().unwrap();
        for (a, b) in actual
            .effective_stress_compression_positive_kpa
            .iter()
            .zip(direct.state.effective_stress_compression_positive_kpa)
        {
            assert!((a - b).abs() < 1e-6, "stress {a} vs {b}");
        }
        assert!(
            (actual.preconsolidation_pressure_kpa - direct.state.preconsolidation_pressure_kpa)
                .abs()
                < 1e-6
        );
        assert!((mixed.cumulative_volume_strain_compression_positive[i] - 0.005).abs() < 1e-10);
        assert!(actual.plastic_strain_compression_positive[2] > 0.);
    }
    for (a, b) in loaded
        .material_states
        .iter()
        .zip(&run.last_committed.material_states)
    {
        assert_eq!(
            a.as_mcc().unwrap().preconsolidation_pressure_kpa,
            b.as_mcc().unwrap().preconsolidation_pressure_kpa
        );
        assert_eq!(a.plastic_strain(), b.plastic_strain());
    }
    assert!(run.frames.iter().all(|f| f.relative_force_residual < 1e-9));
    assert!(mixed.relative_volume_equation_residual < 1e-9);
    assert!(mixed.relative_pressure_consistency_residual < 1e-9);
    Ok(())
}

#[test]
fn nearly_incompressible_elastic_pressure_matches_dense_sparse_and_column_theory() -> Result<()> {
    let (mesh, mut ground, mut case) = fixture()?;
    for material in &mut ground.materials {
        if let soil_fem_3d::soil::Model::Mcc(p) = &mut material.model {
            p.poisson = 0.4999;
        }
    }
    for layer in &mut ground.layers {
        layer.consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 10. };
    }
    case.boundary.as_mut().unwrap().constraints.pop();
    case.stages = vec![Stage::Pressure {
        name: "tiny_elastic".into(),
        target_pressure_kpa: 100.01,
        initial_increment: None,
    }];
    let dense = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(dense.complete, "{:?}", dense.failure_reason);
    let initial = dense.frames[0].state.material_states[0].as_mcc().unwrap();
    let p = ground.materials[0].model.as_mcc().unwrap();
    let bulk = p.moduli(initial)?.0;
    let modulus = 3. * bulk * (1. - p.poisson) / (1. + p.poisson);
    let top = dense.last_committed.displacement_m[mesh.node_sets["top"][0]][2];
    assert!((top / (-0.01 / modulus) - 1.).abs() < 1e-5, "{top}");
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions::default(),
    };
    let sparse = stages::solve(&mesh, ground, case)?;
    assert!(sparse.complete, "{:?}", sparse.failure_reason);
    for (a, b) in dense
        .last_committed
        .displacement_m
        .iter()
        .flatten()
        .zip(sparse.last_committed.displacement_m.iter().flatten())
    {
        assert!((a - b).abs() < 1e-10);
    }
    Ok(())
}

#[test]
fn rejection_and_checkpoint_retain_primary_pressure_volume_and_plastic_history() -> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.verification_reject_once_stages = vec!["load".into()];
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let rejected = run.attempts.iter().find(|a| !a.accepted).unwrap();
    assert!(
        rejected.committed_state_unchanged_on_rejection && rejected.rejected_trial_plastic_changed
    );
    let mut control = case.clone();
    control.verification_reject_once_stages.clear();
    if let Stage::Boundary {
        initial_increment, ..
    } = &mut control.stages[0]
    {
        *initial_increment = Some(0.5);
    }
    let control = stages::solve(&mesh, ground.clone(), control)?;
    assert!(control.complete, "{:?}", control.failure_reason);
    assert!(
        serde_json::to_vec(&run.last_committed)? == serde_json::to_vec(&control.last_committed)?,
        "rejected trial changed the control path"
    );
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let checkpoint = Checkpoint::seal(&mesh, partial)?;
    let resumed = checkpoint.resume(&mesh, ground.clone(), case.clone(), None)?;
    assert!(
        serde_json::to_vec(&run.last_committed)? == serde_json::to_vec(&resumed.last_committed)?,
        "restart changed committed state"
    );
    let mut partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    partial
        .last_committed
        .mixed_mcc
        .as_mut()
        .unwrap()
        .pressure_dofs[0]
        .compression_positive_kpa += 1.;
    // Re-seal a maliciously changed value: semantic replay/constraint checks,
    // rather than just the outer checksum, must reject it.
    let changed = Checkpoint::seal(&mesh, partial)?;
    assert!(changed.resume(&mesh, ground, case, None).is_err());
    Ok(())
}

#[test]
fn excavation_drops_inactive_pressure_dofs_and_freezes_removed_history() -> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.stages.push(Stage::Excavate {
        name: "cut".into(),
        bodies: vec!["UpperClay".into()],
        remaining_surface_pressure_kpa: 100.,
        initial_increment: Some(0.5),
    });
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let before = &run.frames[2].state;
    let after = &run.last_committed;
    assert_eq!(after.active_elements.iter().filter(|a| **a).count(), 6);
    assert_eq!(before.mixed_mcc.as_ref().unwrap().pressure_dofs.len(), 16);
    assert_eq!(after.mixed_mcc.as_ref().unwrap().pressure_dofs.len(), 8);
    for (i, q) in run.integration_points.iter().enumerate() {
        if !after.active_elements[q.element_index] {
            assert_eq!(
                serde_json::to_vec(&before.material_states[i])?,
                serde_json::to_vec(&after.material_states[i])?
            );
            assert_eq!(
                before
                    .mixed_mcc
                    .as_ref()
                    .unwrap()
                    .cumulative_volume_strain_compression_positive[i],
                after
                    .mixed_mcc
                    .as_ref()
                    .unwrap()
                    .cumulative_volume_strain_compression_positive[i]
            );
        }
    }
    assert!(run.frames.last().unwrap().relative_force_residual < 1e-9);
    Ok(())
}

#[test]
fn plastic_geostatic_equilibration_and_restart_use_the_mixed_backend() -> Result<()> {
    let (mut mesh, _, _) = fixture()?;
    mesh.face_sets.insert(
        "top".into(),
        soil_fem_3d::boundary::surface_owners(&mesh)
            .into_iter()
            .filter(|(face, owners)| {
                owners.len() == 1 && face.iter().all(|&n| mesh.nodes[n][2] == 1.)
            })
            .map(|(face, _)| face)
            .collect(),
    );
    let ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.schema = "soil-fem-mixed-mcc-stages/1".into();
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let geo = run.geostatic.as_ref().unwrap();
    assert!(
        geo.complete
            && geo
                .attempts
                .iter()
                .any(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    assert!(
        geo.last_committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|p| p.abs() > 1e-6))
    );
    assert!(
        geo.last_committed
            .mixed_mcc
            .as_ref()
            .unwrap()
            .relative_volume_equation_residual
            < 1e-9
    );
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let resumed = Checkpoint::seal(&mesh, partial)?.resume(&mesh, ground, case, None)?;
    assert!(
        serde_json::to_vec(&run.last_committed)? == serde_json::to_vec(&resumed.last_committed)?,
        "geostatic restart changed full mixed state"
    );
    Ok(())
}

#[test]
fn unsupported_water_and_underintegrated_cells_are_explicitly_rejected() -> Result<()> {
    let (mesh, mut ground, case) = fixture()?;
    ground.initialization.water_level_m = Some(0.);
    assert!(stages::solve(&mesh, ground, case).is_err());
    let (mut mesh, ground, case) = fixture()?;
    mesh.quadrature_order = Some(2);
    assert!(stages::solve(&mesh, ground, case).is_err());
    Ok(())
}
