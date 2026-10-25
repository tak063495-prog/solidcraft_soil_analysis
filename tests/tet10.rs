use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::{self, Ground},
    mesh::{Mesh, Region},
    stages::{self, Stage, StageCase},
};
use std::collections::BTreeMap;

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
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
        material_ids: vec![],
        region_ids: vec![],
        regions: vec![],
        node_sets: BTreeMap::new(),
        face_sets: BTreeMap::new(),
        excavation_sets: BTreeMap::new(),
    };
    let xy = [0., 3., 7., 10.];
    let z = [0., 4., 6., 8.];
    for k in z {
        for y in xy {
            for x in xy {
                mesh.nodes.push([x, y, k]);
            }
        }
    }
    let mut pit = Vec::new();
    for k in 0usize..3 {
        for j in 0usize..3 {
            for i in 0usize..3 {
                let node = |x, y, z| i + x + 4 * (j + y + 4 * (k + z));
                let c = [
                    node(0, 0, 0),
                    node(1, 0, 0),
                    node(1, 1, 0),
                    node(0, 1, 0),
                    node(0, 0, 1),
                    node(1, 0, 1),
                    node(1, 1, 1),
                    node(0, 1, 1),
                ];
                for t in [
                    [0, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
                ] {
                    if i == 1 && j == 1 && k == 2 {
                        pit.push(mesh.tetrahedra.len());
                    }
                    mesh.tetrahedra.push(t.map(|a| c[a]));
                    mesh.material_ids.push(if k == 0 { 1 } else { 2 });
                    mesh.region_ids.push(if k == 0 { 1 } else { 2 });
                }
            }
        }
    }
    mesh.excavation_sets.insert("Pit".into(), pit);
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
    mesh.node_sets.insert(
        "Pit:floor".into(),
        mesh.nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[2] == 6. && (3. ..=7.).contains(&p[0]) && (3. ..=7.).contains(&p[1]))
            .map(|(i, _)| i)
            .collect(),
    );
    mesh.node_sets
        .insert("Pit:wall".into(), mesh.node_sets["Pit:floor"].clone());
    mesh.node_sets.insert(
        "Pit:surrounding_surface".into(),
        mesh.nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[2] == 8.)
            .map(|(i, _)| i)
            .collect(),
    );
    let ground = serde_json::from_str(include_str!("../examples/ground-local.json"))?;
    let case = serde_json::from_str(include_str!("../examples/stages-local-sparse.json"))?;
    for (name, a, v) in [
        ("xmin", 0, 0.),
        ("xmax", 0, 10.),
        ("ymin", 1, 0.),
        ("ymax", 1, 10.),
        ("bottom", 2, 0.),
        ("top", 2, 8.),
    ] {
        mesh.node_sets.insert(
            name.into(),
            mesh.nodes
                .iter()
                .enumerate()
                .filter(|(_, p)| p[a] == v)
                .map(|(i, _)| i)
                .collect(),
        );
        mesh.face_sets.insert(
            name.into(),
            soil_fem_3d::boundary::surface_owners(&mesh)
                .iter()
                .filter(|(f, o)| o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][a] == v))
                .map(|(&f, _)| f)
                .collect(),
        );
    }
    Ok((mesh.elevate_straight()?, ground, case))
}

fn boundary_constraint(
    set: &str,
    direction: [f64; 3],
    value: f64,
) -> soil_fem_3d::boundary::Constraint {
    soil_fem_3d::boundary::Constraint {
        node_set: set.into(),
        direction,
        value_m: value,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    }
}

#[test]
fn general_surface_load_connects_to_mcc_newton_and_matches_pressure_stage() -> Result<()> {
    use soil_fem_3d::boundary::{Conditions, SurfaceLoad};
    let (mesh, ground, mut case) = fixture()?;
    case.stages = vec![Stage::Pressure {
        name: "pressure".into(),
        target_pressure_kpa: 31.,
        initial_increment: None,
    }];
    let reference = stages::solve(&mesh, ground.clone(), case.clone())?;
    case.boundary = Some(Conditions {
        constraints: vec![],
        nodal_loads: vec![],
        surface_loads: vec![SurfaceLoad::Pressure {
            face_set: "top".into(),
            pressure_kpa: 1.,
            gradient_kpa_per_m: [0.; 3],
            reference_point_m: [0.; 3],
        }],
    });
    case.stages = vec![Stage::Boundary {
        name: "face_pressure".into(),
        target_factor: 1.,
        initial_increment: None,
    }];
    let run = stages::solve(&mesh, ground, case)?;
    assert!(
        reference.complete && run.complete,
        "{:?} {:?}",
        reference.failure_reason,
        run.failure_reason
    );
    for (a, b) in reference
        .last_committed
        .displacement_m
        .iter()
        .flatten()
        .zip(run.last_committed.displacement_m.iter().flatten())
    {
        assert!((a - b).abs() < 1e-9);
    }
    assert_eq!(run.last_committed.boundary_factor, 1.);
    Ok(())
}

#[test]
fn prescribed_displacement_mcc_stages_restore_offset_and_plastic_history_on_restart() -> Result<()>
{
    use soil_fem_3d::boundary::Conditions;
    let (mesh, ground, mut case) = fixture()?;
    case.boundary = Some(Conditions {
        constraints: vec![boundary_constraint("top", [0., 0., 1.], -0.02)],
        ..Default::default()
    });
    case.stages = vec![
        Stage::Boundary {
            name: "compression".into(),
            target_factor: 1.,
            initial_increment: None,
        },
        Stage::Boundary {
            name: "reverse".into(),
            target_factor: 0.,
            initial_increment: None,
        },
    ];
    case.verification_reject_once_stages = vec!["compression".into()];
    let reference = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(reference.complete, "{:?}", reference.failure_reason);
    assert!(
        reference.frames[1]
            .state
            .material_states
            .iter()
            .any(|s| { s.plastic_strain().iter().any(|v| v.abs() > 1e-7) })
    );
    for &node in &mesh.node_sets["top"] {
        assert!((reference.frames[1].state.displacement_m[node][2] + 0.02).abs() < 1e-11);
    }
    let saved = stages::solve_until(&mesh, ground.clone(), case.clone(), 2)?;
    let sealed = Checkpoint::seal(&mesh, saved)?;
    let restored: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&sealed)?)?;
    let resumed = restored.resume(&mesh, ground, case, None)?;
    assert!(resumed.complete, "{:?}", resumed.failure_reason);
    assert_eq!(
        serde_json::to_vec(&resumed.last_committed)?,
        serde_json::to_vec(&reference.last_committed)?
    );
    assert!(
        resumed
            .frames
            .iter()
            .all(|f| f.maximum_constraint_error_m < 1e-10)
    );
    Ok(())
}

#[test]
fn nonaxis_displacement_constraints_and_explicit_supports_work_in_mcc() -> Result<()> {
    use soil_fem_3d::boundary::Conditions;
    let (mesh, ground, mut case) = fixture()?;
    case.boundary = Some(Conditions {
        constraints: vec![boundary_constraint("top", [0.3, 0.2, 1.], -0.005)],
        ..Default::default()
    });
    case.stages = vec![Stage::Boundary {
        name: "inclined".into(),
        target_factor: 1.,
        initial_increment: None,
    }];
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for &n in &mesh.node_sets["top"] {
        let u = run.last_committed.displacement_m[n];
        assert!((0.3 * u[0] + 0.2 * u[1] + u[2] + 0.005).abs() < 1e-10);
    }
    case.settings.kinematics = stages::Kinematics::Explicit;
    let mut constraints = vec![];
    for (set, a) in [
        ("xmin", 0),
        ("xmax", 0),
        ("ymin", 1),
        ("ymax", 1),
        ("bottom", 2),
    ] {
        let mut d = [0.; 3];
        d[a] = 1.;
        constraints.push(boundary_constraint(set, d, 0.));
    }
    case.boundary = Some(Conditions {
        constraints,
        ..Default::default()
    });
    case.stages = vec![Stage::Pressure {
        name: "pressure".into(),
        target_pressure_kpa: 31.,
        initial_increment: None,
    }];
    let explicit = stages::solve(&mesh, ground.clone(), case.clone())?;
    case.settings.kinematics = stages::Kinematics::SideRollers;
    case.boundary = None;
    let builtin = stages::solve(&mesh, ground, case)?;
    assert!(
        explicit.complete && builtin.complete,
        "{:?} {:?}",
        explicit.failure_reason,
        builtin.failure_reason
    );
    for (a, b) in explicit
        .last_committed
        .displacement_m
        .iter()
        .flatten()
        .zip(builtin.last_committed.displacement_m.iter().flatten())
    {
        assert!((a - b).abs() < 1e-9);
    }
    Ok(())
}

#[test]
fn quadratic_k0_state_is_in_nodal_equilibrium_and_stores_every_material_point() -> Result<()> {
    let (mesh, mut ground, _) = fixture()?;
    ground.initialization.water_level_m = Some(8.);
    let initial = ground::initialize(&mesh, ground)?;
    assert_eq!(initial.point_states.len(), 27 * mesh.tetrahedra.len());
    assert!(initial.true_relative_free_force_residual < 1e-11);
    assert!((initial.body_volumes_m3["LowerClay"] - 400.).abs() < 1e-9);
    assert!(initial.equilibrium_error_kn.iter().all(|f| f.abs() < 1e-8));
    for (i, p) in initial.point_states.iter().enumerate() {
        assert_eq!(p.element_index, i / 27);
        assert_eq!(p.local_index, i % 27);
        assert!(p.weight_m3 > 0.);
    }
    Ok(())
}

#[test]
fn quadratic_mcc_excavation_cutback_and_checkpoint_restore_all_quadrature_history() -> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    case.verification_reject_once_stages = vec!["excavate_pit".into()];
    let reference = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(reference.complete, "{:?}", reference.failure_reason);
    assert_eq!(reference.schema, "soil-fem-stage-result/2");
    assert_eq!(
        reference.last_committed.material_states.len(),
        27 * mesh.tetrahedra.len()
    );
    let saved = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let checkpoint = Checkpoint::seal(&mesh, saved)?;
    let bytes = serde_json::to_vec(&checkpoint)?;
    let restored: Checkpoint = serde_json::from_slice(&bytes)?;
    let resumed = restored.resume(&mesh, ground, case, None)?;
    assert!(resumed.complete, "{:?}", resumed.failure_reason);
    assert_eq!(
        serde_json::to_vec(&resumed.last_committed)?,
        serde_json::to_vec(&reference.last_committed)?
    );
    assert!(
        resumed.frames[1].monitors["excavation_floor"]
            .mean_displacement_z_m
            .unwrap_or(0.)
            > 0.
    );
    let removed = &mesh.excavation_sets["Pit"];
    for &ei in removed {
        for qi in 0..27 {
            assert_eq!(
                serde_json::to_vec(&resumed.frames[0].state.material_states[27 * ei + qi])?,
                serde_json::to_vec(&resumed.last_committed.material_states[27 * ei + qi])?
            );
        }
    }
    assert!(resumed.excavation_audits[0].zero_release_relative_force_residual < 1e-11);
    assert!(
        resumed
            .attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_state_unchanged_on_rejection)
    );
    Ok(())
}

#[test]
fn quadratic_uniform_constrained_compression_matches_theory() -> Result<()> {
    let (mesh, mut ground, mut case) = fixture()?;
    ground.initialization.water_level_m = None;
    ground.initialization.surface_surcharge_kpa = 100.;
    for material in &mut ground.materials {
        material.unit_weight_kn_m3 = 1e-8;
    }
    case.settings.kinematics = stages::Kinematics::Oedometer;
    case.stages = vec![Stage::Pressure {
        name: "tiny_pressure".into(),
        target_pressure_kpa: 100.0001,
        initial_increment: Some(1.),
    }];
    let initial = ground::initialize(&mesh, ground.clone())?;
    let p = &initial.point_states[0];
    let k = p.bulk_modulus_kpa;
    let g = p.shear_modulus_kpa;
    let expected = -8. * 0.0001 / (k + 4. * g / 3.);
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let value = run.frames[1].monitors["surrounding_surface"]
        .mean_displacement_z_m
        .unwrap_or(0.);
    assert!(
        (value / expected - 1.).abs() < 1e-5,
        "computed={value:e},expected={expected:e}"
    );
    Ok(())
}

#[test]
fn quadratic_mesh_rejects_inconsistent_edges_and_hidden_inversion() -> Result<()> {
    let (mesh, _, _) = fixture()?;
    let mut bad = mesh.clone();
    let duplicate = bad.nodes.len();
    bad.nodes.push(bad.nodes[bad.tetrahedra10[0][4]]);
    bad.tetrahedra10[0][4] = duplicate;
    assert!(bad.validate().is_err());
    let mut bad = mesh.clone();
    bad.quadrature_order = Some(2);
    assert!(bad.validate().is_err());
    let mut bad = mesh.clone();
    let corner = bad.tetrahedra10[0][0];
    bad.tetrahedra10[0][4] = corner;
    assert!(bad.validate().is_err());
    Ok(())
}
