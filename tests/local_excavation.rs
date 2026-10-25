use soil_fem_3d::{
    Result,
    ground::{self, Ground},
    mesh::{Mesh, Region},
    stages::{self, Stage, StageCase},
};
use std::collections::BTreeMap;

fn specimen() -> Result<(Mesh, Ground, StageCase)> {
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
        region_ids: vec![],
        regions: vec![],
        node_sets: BTreeMap::new(),
        excavation_sets: BTreeMap::new(),
    };
    let xy = [0., 3., 5., 7., 10.];
    let zs = [0., 4., 6., 8.];
    for z in zs {
        for y in xy {
            for x in xy {
                mesh.nodes.push([x, y, z]);
            }
        }
    }
    let mut cut = vec![];
    for z in 0usize..3 {
        for y in 0usize..4 {
            for x in 0usize..4 {
                let n = |a, b, c| x + a + 5 * (y + b) + 25 * (z + c);
                let c = [
                    n(0, 0, 0),
                    n(1, 0, 0),
                    n(1, 1, 0),
                    n(0, 1, 0),
                    n(0, 0, 1),
                    n(1, 0, 1),
                    n(1, 1, 1),
                    n(0, 1, 1),
                ];
                for t in [
                    [0, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
                ] {
                    if z == 2 && (1..3).contains(&x) && (1..3).contains(&y) {
                        cut.push(mesh.tetrahedra.len());
                    }
                    mesh.tetrahedra.push(t.map(|i| c[i]));
                    mesh.material_ids.push(if z == 0 { 1 } else { 2 });
                    mesh.region_ids.push(if z == 0 { 1 } else { 2 });
                }
            }
        }
    }
    mesh.excavation_sets.insert("Pit".into(), cut);
    for (i, name) in ["LowerClay", "UpperClay"].iter().enumerate() {
        mesh.regions.push(Region {
            id: i as u32 + 1,
            body: (*name).into(),
            material_id: i as u32 + 1,
            bounds_m: [[0., 0., 4. * i as f64], [10., 10., 4. * (i + 1) as f64]],
            volume_m3: 400.,
            source_sha256: "a".repeat(64),
        });
    }
    for (name, z) in [("Pit:floor", 6.), ("Pit:surrounding_surface", 8.)] {
        mesh.node_sets.insert(
            name.into(),
            mesh.nodes
                .iter()
                .enumerate()
                .filter(|(_, p)| {
                    p[2] == z
                        && (z == 8. || ((3. ..=7.).contains(&p[0]) && (3. ..=7.).contains(&p[1])))
                })
                .map(|(i, _)| i)
                .collect(),
        );
    }
    mesh.node_sets
        .insert("Pit:wall".into(), mesh.node_sets["Pit:floor"].clone());
    let ground = serde_json::from_str(include_str!("../examples/ground-local.json"))?;
    let case = serde_json::from_str(include_str!("../examples/stages-local.json"))?;
    Ok((mesh, ground, case))
}

#[test]
fn checkpoint_roundtrip_reproduces_plastic_and_excavation_history_after_every_increment()
-> Result<()> {
    use soil_fem_3d::{
        checkpoint::Checkpoint,
        linalg::{GmresOptions, LinearSolver},
    };
    let (mesh, ground, mut case) = specimen()?;
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions::default(),
    };
    case.stages.insert(
        0,
        Stage::Pressure {
            name: "preload".into(),
            target_pressure_kpa: 90.,
            initial_increment: Some(0.25),
        },
    );
    case.stages.insert(
        1,
        Stage::Pressure {
            name: "unload".into(),
            target_pressure_kpa: 30.,
            initial_increment: Some(0.25),
        },
    );
    case.verification_reject_once_stages = vec!["preload".into(), "excavate_pit".into()];
    let reference = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(reference.complete, "{:?}", reference.failure_reason);
    let mut interrupted = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let mut checkpoints = 0;
    let mut partial_excavation = false;
    while !interrupted.complete {
        partial_excavation |= interrupted.last_committed.stage_index == 3
            && interrupted.last_committed.stage_progress < 1.;
        let sealed = Checkpoint::seal(&mesh, interrupted)?;
        let bytes = serde_json::to_vec_pretty(&sealed)?;
        let restored: Checkpoint = serde_json::from_slice(&bytes)?;
        interrupted = restored.resume(&mesh, ground.clone(), case.clone(), Some(1))?;
        checkpoints += 1;
        assert!(checkpoints < 80);
    }
    assert!(partial_excavation && checkpoints > 10);
    assert_eq!(
        serde_json::to_vec(&interrupted.last_committed)?,
        serde_json::to_vec(&reference.last_committed)?
    );
    assert_eq!(interrupted.frames.len(), reference.frames.len());
    assert_eq!(interrupted.attempts.len(), reference.attempts.len());
    for (a, b) in interrupted.attempts.iter().zip(&reference.attempts) {
        assert_eq!(
            (a.accepted, a.progress_from, a.progress_to),
            (b.accepted, b.progress_from, b.progress_to)
        );
    }
    Ok(())
}

#[test]
fn checkpoints_reject_changed_inputs_binary_and_corrupted_load_history() -> Result<()> {
    use soil_fem_3d::checkpoint::Checkpoint;
    let (mesh, ground, case) = specimen()?;
    let sealed = Checkpoint::seal(
        &mesh,
        stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?,
    )?;
    let original = serde_json::to_vec(&sealed)?;
    let mut changed_mesh = mesh.clone();
    changed_mesh.crs = "DIFFERENT".into();
    assert!(sealed.verify(&changed_mesh, &ground, &case).is_err());
    let mut changed_ground = ground.clone();
    changed_ground.initialization.surface_surcharge_kpa += 1.;
    assert!(sealed.verify(&mesh, &changed_ground, &case).is_err());
    let mut changed_case = case.clone();
    changed_case.settings.relative_force_tolerance = 1e-7;
    assert!(sealed.verify(&mesh, &ground, &changed_case).is_err());
    let mut bad: Checkpoint = serde_json::from_slice(&original)?;
    bad.identity.solver_binary_sha256 = "0".repeat(64);
    assert!(bad.verify(&mesh, &ground, &case).is_err());
    let mut bad: Checkpoint = serde_json::from_slice(&original)?;
    bad.run.last_committed.excavation_release_load_kn[0] += 1.;
    assert!(bad.verify(&mesh, &ground, &case).is_err());
    // Even a recomputed checksum cannot bypass the load-path consistency checks.
    bad.run_sha256 = soil_fem_3d::checkpoint::hash_json(&bad.run)?;
    assert!(
        bad.resume(&mesh, ground.clone(), case.clone(), None)
            .is_err()
    );
    assert!(stages::solve_until(&mesh, ground.clone(), case.clone(), 0).is_err());
    assert!(Checkpoint::seal(&mesh, stages::solve(&mesh, ground, case)?).is_err());
    Ok(())
}

#[test]
fn sparse_newton_matches_dense_for_plastic_local_excavation_and_reports_true_residual() -> Result<()>
{
    use soil_fem_3d::linalg::{GmresOptions, LinearSolver};
    let (mesh, ground, mut case) = specimen()?;
    case.stages.insert(
        0,
        Stage::Pressure {
            name: "preload".into(),
            target_pressure_kpa: 90.,
            initial_increment: Some(0.25),
        },
    );
    case.stages.insert(
        1,
        Stage::Pressure {
            name: "unload".into(),
            target_pressure_kpa: 30.,
            initial_increment: Some(0.25),
        },
    );
    let dense = stages::solve(&mesh, ground.clone(), case.clone())?;
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions {
            restart: 20,
            ..Default::default()
        },
    };
    case.settings.relative_linear_tolerance = 1e-11;
    let sparse = stages::solve(&mesh, ground, case)?;
    assert!(
        dense.complete && sparse.complete,
        "{:?} {:?}",
        dense.failure_reason,
        sparse.failure_reason
    );
    for (a, b) in dense.frames.iter().zip(&sparse.frames) {
        let norm: f64 = a.state.displacement_m.iter().flatten().map(|v| v * v).sum();
        let error: f64 = a
            .state
            .displacement_m
            .iter()
            .flatten()
            .zip(b.state.displacement_m.iter().flatten())
            .map(|(x, y)| (x - y).powi(2))
            .sum();
        assert!((error / norm.max(1e-24)).sqrt() < 1e-6);
        assert_eq!(a.state.active_elements, b.state.active_elements);
        for (sa, sb) in a.state.material_states.iter().zip(&b.state.material_states) {
            assert!(
                (sa.as_mcc().unwrap().preconsolidation_pressure_kpa
                    / sb.as_mcc().unwrap().preconsolidation_pressure_kpa
                    - 1.)
                    .abs()
                    < 1e-7
            );
            assert!((sa.void_ratio() - sb.void_ratio()).abs() < 1e-8);
        }
    }
    let reports: Vec<_> = sparse
        .attempts
        .iter()
        .flat_map(|a| &a.linear_solves)
        .collect();
    assert!(!reports.is_empty());
    assert!(
        reports
            .iter()
            .all(|r| r.method == "csr_gmres_ilu0" && r.true_relative_residual <= 1e-11)
    );
    assert!(
        reports
            .iter()
            .any(|r| r.stored_nonzeros < r.dimension * r.dimension / 2)
    );
    Ok(())
}

#[test]
fn sparse_local_cutback_retries_identically_and_linear_failure_is_not_committed() -> Result<()> {
    use soil_fem_3d::linalg::{GmresOptions, LinearSolver};
    let (mesh, ground, mut case) = specimen()?;
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions::default(),
    };
    case.verification_reject_once_stages = vec!["excavate_pit".into()];
    let injected = stages::solve(&mesh, ground.clone(), case.clone())?;
    case.verification_reject_once_stages.clear();
    if let Stage::Excavate {
        initial_increment, ..
    } = &mut case.stages[0]
    {
        *initial_increment = Some(0.125);
    }
    let clean = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(
        injected.complete && clean.complete,
        "{:?} {:?}",
        injected.failure_reason,
        clean.failure_reason
    );
    assert_eq!(
        serde_json::to_vec(&injected.last_committed)?,
        serde_json::to_vec(&clean.last_committed)?
    );
    assert!(
        injected
            .attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions {
            restart: 1,
            maximum_iterations: 1,
            ..Default::default()
        },
    };
    case.settings.minimum_increment = 0.25;
    if let Stage::Excavate {
        initial_increment, ..
    } = &mut case.stages[0]
    {
        *initial_increment = Some(0.25);
    }
    let failed = stages::solve(&mesh, ground, case)?;
    assert!(!failed.complete);
    assert!(
        failed
            .failure_reason
            .as_deref()
            .unwrap_or("")
            .contains("GMRES iteration limit")
    );
    assert_eq!(
        serde_json::to_vec(&failed.last_committed)?,
        serde_json::to_vec(&failed.frames[0].state)?
    );
    Ok(())
}

#[test]
fn local_cut_keeps_surrounding_soil_and_removes_floor_surcharge() -> Result<()> {
    let (mesh, ground, case) = specimen()?;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let end = &run.frames[1];
    let audit = &run.excavation_audits[0];
    assert_eq!(audit.removed_element_count, 24);
    assert!((audit.removed_volume_m3 - 32.).abs() < 1e-10);
    assert!(audit.zero_release_relative_force_residual < 1e-12);
    assert_eq!(
        end.state.active_elements.iter().filter(|a| **a).count(),
        264
    );
    assert!((end.applied_force_kn[2] + 14624. + 84. * 30.).abs() < 1e-8);
    assert!(
        end.monitors["excavation_floor"]
            .mean_displacement_z_m
            .unwrap_or(0.)
            > 0.
    );
    assert!(
        end.monitors["surrounding_surface"]
            .minimum_displacement_z_m
            .is_some()
    );
    assert!(
        end.state
            .displacement_m
            .iter()
            .any(|u| u[0].abs() + u[1].abs() > 1e-5)
    );
    assert!(
        end.state
            .excavation_release_load_kn
            .iter()
            .all(|x| *x == 0.)
    );
    for &i in &mesh.excavation_sets["Pit"] {
        assert_eq!(
            serde_json::to_string(&end.state.material_states[i])?,
            serde_json::to_string(&run.frames[0].state.material_states[i])?
        );
    }
    Ok(())
}

#[test]
fn local_cutback_restores_plastic_history_loads_and_active_elements() -> Result<()> {
    let (mesh, ground, mut case) = specimen()?;
    case.stages.insert(
        0,
        Stage::Pressure {
            name: "preload".into(),
            target_pressure_kpa: 90.,
            initial_increment: Some(0.25),
        },
    );
    case.stages.insert(
        1,
        Stage::Pressure {
            name: "unload".into(),
            target_pressure_kpa: 30.,
            initial_increment: Some(0.25),
        },
    );
    case.verification_reject_once_stages = vec!["preload".into(), "excavate_pit".into()];
    let mut clean = case.clone();
    clean.verification_reject_once_stages.clear();
    if let Stage::Pressure {
        initial_increment, ..
    } = &mut clean.stages[0]
    {
        *initial_increment = Some(0.125);
    }
    if let Stage::Excavate {
        initial_increment, ..
    } = &mut clean.stages[2]
    {
        *initial_increment = Some(0.125);
    }
    let a = stages::solve(&mesh, ground.clone(), case)?;
    let b = stages::solve(&mesh, ground, clean)?;
    assert!(
        a.complete && b.complete,
        "{:?} {:?}",
        a.failure_reason,
        b.failure_reason
    );
    assert!(
        a.frames[1]
            .state
            .material_states
            .iter()
            .any(|s| { s.plastic_strain().iter().any(|x| x.abs() > 1e-6) })
    );
    let rejects: Vec<_> = a.attempts.iter().filter(|t| !t.accepted).collect();
    assert_eq!(rejects.len(), 2);
    assert!(
        rejects
            .iter()
            .all(|t| t.committed_state_unchanged_on_rejection
                && t.rejected_trial_material_changed
                && t.rejected_trial_load_changed)
    );
    assert!(rejects.iter().any(|t| t.rejected_trial_active_set_changed));
    assert_eq!(
        serde_json::to_string(&a.last_committed)?,
        serde_json::to_string(&b.last_committed)?
    );
    Ok(())
}

#[test]
fn submerged_cut_has_hydrostatic_wall_and_floor_tractions() -> Result<()> {
    let (mesh, mut ground, case) = specimen()?;
    ground.initialization.water_level_m = Some(8.);
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let end = &run.frames[1];
    let water_force = 16. * 2. * 9.80665;
    assert!((end.applied_force_kn[2] + 14624. + 84. * 30. + water_force).abs() < 1e-8);
    let lateral: f64 = end
        .state
        .applied_nodal_load_kn
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 3 != 2)
        .map(|(_, f)| f.abs())
        .sum();
    assert!(lateral > 10.);
    assert!(end.force_balance_error_kn.iter().all(|f| f.abs() < 1e-5));
    Ok(())
}

#[test]
fn buried_cavity_duplicate_cut_and_ambiguous_pressure_are_rejected() -> Result<()> {
    let (mut mesh, ground, mut case) = specimen()?;
    case.pressure_on_original_surface_only = false;
    assert!(stages::solve(&mesh, ground.clone(), case.clone()).is_err());
    case.pressure_on_original_surface_only = true;
    case.stages.push(case.stages[0].clone());
    if let Stage::Excavate { name, .. } = &mut case.stages[1] {
        *name = "repeat".into();
    }
    assert!(stages::solve(&mesh, ground.clone(), case.clone()).is_err());
    case.stages.pop();
    // Interior slab at z=4..6 cannot be removed while the overlying soil remains.
    mesh.excavation_sets.insert(
        "Pit".into(),
        (96..192)
            .filter(|&i| {
                let p = mesh.points(&mesh.tetrahedra[i]);
                let x = p.iter().map(|p| p[0]).sum::<f64>() / 4.;
                let y = p.iter().map(|p| p[1]).sum::<f64>() / 4.;
                x > 3. && x < 7. && y > 3. && y < 7.
            })
            .collect(),
    );
    assert!(stages::solve(&mesh, ground, case).is_err());
    Ok(())
}

#[test]
fn local_failure_keeps_the_last_accepted_state() -> Result<()> {
    let (mesh, ground, mut case) = specimen()?;
    case.settings.maximum_newton_iterations = 1;
    case.settings.minimum_increment = 0.25;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(!run.complete);
    assert_eq!(run.frames.len(), 1);
    assert_eq!(
        serde_json::to_string(&run.last_committed)?,
        serde_json::to_string(&run.frames[0].state)?
    );
    Ok(())
}

#[test]
fn malformed_excavation_sets_and_ground_mapping_are_rejected() -> Result<()> {
    let (mesh, mut ground, _) = specimen()?;
    ground.excavation_bodies.clear();
    assert!(ground::initialize(&mesh, ground).is_err());
    for indices in [vec![], vec![0, 0], vec![mesh.tetrahedra.len()]] {
        let mut bad = mesh.clone();
        bad.excavation_sets.insert("Pit".into(), indices);
        assert!(bad.validate().is_err());
    }
    Ok(())
}

#[test]
fn sequential_local_cuts_can_advance_from_the_existing_floor() -> Result<()> {
    let (mut mesh, mut ground, mut case) = specimen()?;
    let deep: Vec<_> = (96..192)
        .filter(|&i| {
            let p = mesh.points(&mesh.tetrahedra[i]);
            let x = p.iter().map(|p| p[0]).sum::<f64>() / 4.;
            let y = p.iter().map(|p| p[1]).sum::<f64>() / 4.;
            x > 3. && x < 7. && y > 3. && y < 7.
        })
        .collect();
    mesh.excavation_sets.insert("DeepPit".into(), deep);
    ground.excavation_bodies.push("DeepPit".into());
    case.stages.push(Stage::Excavate {
        name: "deepen_pit".into(),
        bodies: vec!["DeepPit".into()],
        remaining_surface_pressure_kpa: 30.,
        initial_increment: None,
    });
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.frames.len(), 3);
    assert_eq!(run.excavation_audits.len(), 2);
    for audit in &run.excavation_audits {
        assert!((audit.removed_volume_m3 - 32.).abs() < 1e-9);
        assert!(audit.zero_release_relative_force_residual < 1e-8);
    }
    assert_eq!(
        run.last_committed
            .active_elements
            .iter()
            .filter(|a| **a)
            .count(),
        240
    );
    Ok(())
}
