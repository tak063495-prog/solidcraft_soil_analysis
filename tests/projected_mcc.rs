use nalgebra::DVector;
use soil_fem_3d::{
    Result,
    boundary::{Conditions, Constraint, NodalLoad},
    checkpoint::Checkpoint,
    ground::{self, Ground},
    linalg::{GmresOptions, LinearSolver},
    mesh::Mesh,
    stages::{self, Kinematics, Stage, StageCase},
};
use std::collections::BTreeMap;

#[test]
fn pressure_unload_then_display_reset_preserves_exact_load_and_mixed_history() -> Result<()> {
    let (mesh, mut ground, mut case) = fixture(false)?;
    ground.initialization.surface_surcharge_kpa = 30.;
    for (i, material) in ground.materials.iter_mut().enumerate() {
        material.unit_weight_kn_m3 = 18. + 2. * i as f64;
    }
    for layer in &mut ground.layers {
        layer.consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 2. };
    }
    case.boundary.as_mut().unwrap().constraints.pop();
    case.settings.initial_increment = 0.25;
    case.stages = vec![
        Stage::Pressure {
            name: "preload".into(),
            target_pressure_kpa: 120.,
            initial_increment: None,
        },
        Stage::Pressure {
            name: "unload".into(),
            target_pressure_kpa: 30.,
            initial_increment: None,
        },
        Stage::ResetDisplay {
            name: "zero_display".into(),
        },
        Stage::Pressure {
            name: "reload".into(),
            target_pressure_kpa: 31.,
            initial_increment: None,
        },
    ];
    case.verification_reject_once_stages = vec!["zero_display".into()];
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    let before = &full.frames[2].state;
    let after = &full.frames[3].state;
    assert_eq!(
        serde_json::to_vec(&before.applied_nodal_load_kn)?,
        serde_json::to_vec(&after.applied_nodal_load_kn)?
    );
    assert_eq!(
        serde_json::to_vec(&before.material_states)?,
        serde_json::to_vec(&after.material_states)?
    );
    assert_eq!(
        serde_json::to_vec(&before.mixed_mcc)?,
        serde_json::to_vec(&after.mixed_mcc)?
    );
    assert_eq!(before.displacement_m, after.displacement_m);
    assert_eq!(
        after.display_reference_displacement_m,
        before.displacement_m
    );
    assert!(full.attempts.iter().any(|a| a.stage == "zero_display"
        && !a.accepted
        && a.committed_state_unchanged_on_rejection));
    let count = full
        .attempts
        .iter()
        .take_while(|a| a.stage != "zero_display")
        .filter(|a| a.accepted)
        .count()
        + 1;
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
    let sealed = Checkpoint::seal(&mesh, partial)?;
    let serialized: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&sealed)?)?;
    let resumed = serialized.resume(&mesh, ground, case, None)?;
    assert_eq!(
        serde_json::to_vec(&full.frames)?,
        serde_json::to_vec(&resumed.frames)?
    );
    Ok(())
}

fn fixture(local: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh =
        serde_json::from_str(include_str!("../examples/hardening-column-tet10-mesh.json"))?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-mixed-mcc.json"))?;
    for material in &mut ground.materials {
        material.unit_weight_kn_m3 = 1e-8;
        if let soil_fem_3d::soil::Model::Mcc(p) = &mut material.model {
            *p = soil_fem_3d::mcc::Parameters {
                critical_state_m: 1.2,
                lambda: 0.2,
                kappa: 0.04,
                poisson: 0.3,
                reference_pressure_kpa: 1.,
                specific_volume_ncl_at_reference: 3.5,
                integration: soil_fem_3d::mcc::Integration::BackwardEuler,
            };
        }
    }
    for layer in &mut ground.layers {
        layer.k0_x = 1.;
        layer.k0_y = 1.;
        layer.consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 1. };
    }
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-mixed-mcc.json"))?;
    case.schema = "soil-fem-projected-mcc-stages/1".into();
    case.settings.kinematics = Kinematics::Explicit;
    case.settings.linear_solver = LinearSolver::DenseLu;
    case.settings.maximum_free_dofs = 2000;
    case.settings.initial_increment = 0.5;
    case.settings.minimum_increment = 1. / 256.;
    case.settings.relative_force_tolerance = 1e-9;
    case.settings.relative_linear_tolerance = 1e-11;
    case.settings.material_maximum_substep = 2e-4;
    case.boundary = Some(Conditions {
        constraints: [
            ("all", 0, 0.),
            ("all", 1, 0.),
            ("bottom", 2, 0.),
            ("top", 2, -0.005),
        ]
        .into_iter()
        .map(|(name, axis, value)| Constraint {
            node_set: name.into(),
            direction: std::array::from_fn(|i| if i == axis { 1. } else { 0. }),
            value_m: value,
            gradient: [0.; 3],
            reference_point_m: [0.; 3],
        })
        .collect(),
        ..Conditions::default()
    });
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
        Stage::Excavate {
            name: "cut".into(),
            bodies: vec![if local { "LocalCut" } else { "UpperClay" }.into()],
            remaining_surface_pressure_kpa: if local { 0. } else { 100. },
            initial_increment: None,
        },
    ];
    if local {
        case.pressure_on_original_surface_only = true;
        let exposed = mesh
            .tetrahedra
            .iter()
            .position(|tet| tet.iter().filter(|&&n| mesh.nodes[n][2] == 1.).count() == 3)
            .unwrap();
        mesh.excavation_sets
            .insert("LocalCut".into(), vec![exposed]);
        ground.excavation_bodies = vec!["LocalCut".into()];
        let point = mesh
            .nodes
            .iter()
            .position(|p| {
                (p[2] - 0.5).abs() < 1e-10 && p[0] > 0.1 && p[0] < 0.9 && p[1] > 0.1 && p[1] < 0.9
            })
            .unwrap();
        mesh.node_sets.insert("point".into(), vec![point]);
        case.boundary.as_mut().unwrap().nodal_loads.push(NodalLoad {
            node_set: "point".into(),
            total_force_kn: [0., 0., -0.3],
        });
    }
    Ok((mesh, ground, case))
}

#[test]
fn affine_plastic_cycle_dense_sparse_and_removed_history() -> Result<()> {
    let (mesh, ground, case) = fixture(false)?;
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.schema, "soil-fem-projected-mcc-stage-result/1");
    for (u, x) in run.frames[1].state.displacement_m.iter().zip(&mesh.nodes) {
        assert!((u[2] + 0.005 * x[2]).abs() < 1e-10);
    }
    assert!(
        run.frames[1]
            .state
            .material_states
            .iter()
            .any(|s| s.plastic_strain()[2] > 1e-4)
    );
    for frame in &run.frames {
        assert!(frame.relative_force_residual < 1e-9);
        assert!(
            frame
                .state
                .mixed_mcc
                .as_ref()
                .unwrap()
                .projected_volume
                .is_some()
        );
    }
    let before = &run.frames[2].state;
    let after = &run.last_committed;
    for (i, point) in run
        .integration_points
        .iter()
        .enumerate()
        .filter(|(_, q)| !after.active_elements[q.element_index])
    {
        assert!(point.element_index >= 6);
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
    let mut sparse = case;
    sparse.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions::default(),
    };
    let other = stages::solve(&mesh, ground, sparse)?;
    assert!(other.complete, "{:?}", other.failure_reason);
    for (a, b) in run
        .last_committed
        .displacement_m
        .iter()
        .flatten()
        .zip(other.last_committed.displacement_m.iter().flatten())
    {
        assert!((a - b).abs() < 1e-10);
    }
    Ok(())
}

#[test]
fn local_death_freezes_nonzero_compatibility_offset_and_restarts() -> Result<()> {
    let (mesh, ground, mut case) = fixture(true)?;
    case.verification_reject_once_stages = vec!["load".into(), "unload".into(), "cut".into()];
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let before = &run.frames[2].state;
    let after = &run.last_committed;
    let projected = after
        .mixed_mcc
        .as_ref()
        .unwrap()
        .projected_volume
        .as_ref()
        .unwrap();
    let mut expected: BTreeMap<_, f64> = projected
        .compatibility_offsets
        .iter()
        .map(|d| ((d.node_index, d.material_id), 0.))
        .collect();
    let mut cursor = 0;
    for (e, tet) in mesh.tetrahedra.iter().enumerate() {
        let q = mesh.integration(e)?;
        if after.active_elements[e] {
            let u = DVector::from_iterator(
                30,
                mesh.cell_nodes(e)
                    .iter()
                    .flat_map(|&n| before.displacement_m[n]),
            );
            for (j, gp) in q.iter().enumerate() {
                let div = (&gp.b * &u).iter().take(3).sum::<f64>();
                let theta = before
                    .mixed_mcc
                    .as_ref()
                    .unwrap()
                    .cumulative_volume_strain_compression_positive[cursor + j];
                let n = [
                    1. - gp.natural.iter().sum::<f64>(),
                    gp.natural[0],
                    gp.natural[1],
                    gp.natural[2],
                ];
                for a in 0..4 {
                    *expected.get_mut(&(tet[a], mesh.material_ids[e])).unwrap() +=
                        gp.weight_m3 * n[a] * (div + theta);
                }
            }
        } else {
            for j in 0..q.len() {
                assert_eq!(
                    serde_json::to_vec(&before.material_states[cursor + j])?,
                    serde_json::to_vec(&after.material_states[cursor + j])?
                );
            }
        }
        cursor += q.len();
    }
    assert!(
        projected
            .compatibility_offsets
            .iter()
            .any(|d| d.integral_m3.abs() > 1e-12)
    );
    for row in &projected.compatibility_offsets {
        assert_eq!(
            row.integral_m3,
            expected[&(row.node_index, row.material_id)]
        );
    }
    for attempt in run.attempts.iter().filter(|a| !a.accepted) {
        assert!(attempt.committed_state_unchanged_on_rejection);
    }
    assert!(
        run.attempts
            .iter()
            .any(|a| !a.accepted && a.rejected_trial_plastic_changed)
    );
    let mut control = case.clone();
    control.verification_reject_once_stages.clear();
    control.settings.initial_increment = 0.25;
    let control = stages::solve(&mesh, ground.clone(), control)?;
    assert!(control.complete, "{:?}", control.failure_reason);
    assert_eq!(
        serde_json::to_vec(&run.last_committed)?,
        serde_json::to_vec(&control.last_committed)?
    );
    let cut_start = run
        .attempts
        .iter()
        .filter(|a| a.accepted && a.stage != "cut")
        .count();
    for stop in [1, cut_start + 1] {
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), stop)?;
        let resumed =
            Checkpoint::seal(&mesh, partial)?.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert_eq!(
            serde_json::to_vec(&run.last_committed)?,
            serde_json::to_vec(&resumed.last_committed)?
        );
    }
    Ok(())
}

#[test]
fn altered_volume_offset_and_formulation_are_rejected_after_resealing() -> Result<()> {
    let (mesh, ground, case) = fixture(true)?;
    let baseline = stages::solve_until(&mesh, ground.clone(), case.clone(), 5)?;
    assert!(!baseline.complete);
    let original = serde_json::to_value(baseline)?;
    for kind in 0..4 {
        let mut value = original.clone();
        let mixed = &mut value["last_committed"]["mixed_mcc"];
        match kind {
            0 => {
                mixed["projected_volume"]["volume_dofs"][0]["compression_positive"] =
                    serde_json::json!(0.03)
            }
            1 => {
                mixed["projected_volume"]["compatibility_offsets"][0]["integral_m3"] =
                    serde_json::json!(0.03)
            }
            2 => mixed["projected_volume"]["formulation"] = serde_json::json!("unknown"),
            _ => mixed["projected_volume"] = serde_json::Value::Null,
        }
        let changed = Checkpoint::seal(&mesh, serde_json::from_value(value)?)?;
        assert!(
            changed
                .resume(&mesh, ground.clone(), case.clone(), None)
                .is_err(),
            "tamper {kind}"
        );
    }
    Ok(())
}
