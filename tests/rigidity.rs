use nalgebra::{Rotation3, Vector3};
use soil_fem_3d::{Result, joints, mesh::Mesh, rigidity, stages::StageCase};

fn rows(mesh: &Mesh, fixed: impl Fn(usize) -> bool) -> Vec<Vec<(usize, f64)>> {
    let mut next = 0;
    (0..mesh.nodes.len() * 3)
        .map(|g| {
            if fixed(g / 3) {
                vec![]
            } else {
                let r = vec![(next, 1.)];
                next += 1;
                r
            }
        })
        .collect()
}

#[test]
fn contact_support_counts_open_normal_stick_and_one_sided_limit_modes() -> Result<()> {
    for quadratic in [false, true] {
        for rotated in [false, true] {
            let mut mesh: Mesh =
                serde_json::from_str(include_str!("../examples/joint-column-mesh.json"))?;
            if quadratic {
                mesh = mesh.elevate_straight()?;
            }
            let rotation = if rotated {
                Rotation3::from_euler_angles(0.4, -0.3, 0.7)
            } else {
                Rotation3::identity()
            };
            let lower: std::collections::BTreeSet<_> = mesh
                .region_ids
                .iter()
                .enumerate()
                .filter(|(_, r)| **r == 1)
                .flat_map(|(e, _)| mesh.cell_nodes(e))
                .collect();
            for x in &mut mesh.nodes {
                *x = (rotation * Vector3::from(*x)).into();
            }
            let active = vec![true; mesh.tetrahedra.len()];
            let transform = rows(&mesh, |n| lower.contains(&n));
            let geometry = rigidity::Geometry::new(&mesh, &active)?;
            let case: StageCase =
                serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
            for (mu, jump, missing) in [
                (0.3, [0., 0., 0.], 6),
                (0.3, [0., 0., 0.001], 6),
                (0., [0., 0., -0.0001], 3),
                (0.3, [0., 0., -0.0001], 0),
                (0.3, [0.00015, 0., -0.0001], 1),
                (0.3, [0.001, 0., -0.0001], 1),
            ] {
                let mut configuration = case.joint_system.clone().unwrap();
                configuration.joints[0].material.friction_coefficient = mu;
                let model = joints::Model::new(&mesh, configuration)?;
                let jump: [f64; 3] = (rotation * Vector3::from(jump)).into();
                let u: Vec<_> = (0..mesh.nodes.len())
                    .map(|n| if lower.contains(&n) { [0.; 3] } else { jump })
                    .collect();
                let eval = model.evaluate(&model.initial_state()?, &u, &active, None, false)?;
                let mut space = geometry.constrained(&transform)?;
                model.add_rigid_support(&eval.state, &active, 1e-8, &mut space)?;
                assert_eq!(space.report().continuum_components, 2);
                assert_eq!(
                    space.report().unrestrained_modes,
                    missing,
                    "mu={mu}, jump={jump:?}, rotated={rotated}, quadratic={quadratic}"
                );
                assert_eq!(space.null_modes().len(), missing);
                assert_eq!(space.require_supported().is_ok(), missing == 0);
            }
        }
    }
    Ok(())
}

#[test]
fn shared_edge_and_vertex_retain_real_hinge_modes() -> Result<()> {
    for (tetrahedra, nodes, missing) in [
        (
            vec![[0, 1, 2, 3], [0, 1, 4, 5]],
            vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [0., -1., 0.],
                [0., 0., -1.],
            ],
            7,
        ),
        (
            vec![[0, 1, 2, 3], [0, 5, 4, 6]],
            vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [-1., 0., 0.],
                [0., -1., 0.],
                [0., 0., -1.],
            ],
            9,
        ),
        (
            vec![[0, 1, 2, 3], [0, 2, 1, 4]],
            vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [0., 0., -1.],
            ],
            6,
        ),
    ] {
        let mesh = Mesh {
            schema: "soil-fem-mesh/1".into(),
            units: "m-kN-kPa".into(),
            axis: "z_up".into(),
            crs: "local".into(),
            origin_m: [0.; 3],
            nodes,
            material_ids: vec![1; tetrahedra.len()],
            tetrahedra,
            tetrahedra10: vec![],
            quadrature_order: None,
            node_sets: Default::default(),
            face_sets: Default::default(),
            regions: vec![],
            region_ids: vec![],
            excavation_sets: Default::default(),
        };
        let geometry = rigidity::Geometry::new(&mesh, &[true, true])?;
        let space = geometry.constrained(&rows(&mesh, |_| false))?;
        assert_eq!(space.report().unrestrained_modes, missing);
    }
    Ok(())
}

#[test]
fn contact_alone_never_restrains_global_motion_after_large_coordinate_translation() -> Result<()> {
    for offset in [Vector3::zeros(), Vector3::new(1e6, -2e6, 3e6)] {
        let mut mesh: Mesh =
            serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
        let rotation = Rotation3::from_euler_angles(0.2, -0.3, 0.4);
        for x in &mut mesh.nodes {
            *x = (rotation * Vector3::from(*x) + offset).into();
        }
        let case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
        let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
        let upper: std::collections::BTreeSet<_> = mesh
            .region_ids
            .iter()
            .enumerate()
            .filter(|(_, r)| **r == 2)
            .flat_map(|(e, _)| mesh.cell_nodes(e))
            .collect();
        let u: Vec<_> = (0..mesh.nodes.len())
            .map(|n| {
                if upper.contains(&n) {
                    [0., 0., -1e-4]
                } else {
                    [0.; 3]
                }
            })
            .collect();
        let active = vec![true; mesh.tetrahedra.len()];
        let geometry = rigidity::Geometry::new(&mesh, &active)?;
        let mut space = geometry.constrained(&rows(&mesh, |_| false))?;
        assert_eq!(
            space.report().continuum_components,
            2,
            "unexpected component split"
        );
        assert_eq!(
            space.report().unrestrained_modes,
            12,
            "unexpected shared-node ties"
        );
        let state = model
            .evaluate(&model.initial_state()?, &u, &active, None, false)?
            .state;
        model.add_rigid_support(&state, &active, 1e-8, &mut space)?;
        assert_eq!(space.report().unrestrained_modes, 6, "offset={offset:?}");
    }
    Ok(())
}

#[test]
fn weak_early_rigid_rows_are_repivoted_without_inventing_extra_constraints() -> Result<()> {
    let mesh: Mesh =
        serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
    let active: Vec<_> = mesh.region_ids.iter().map(|&r| r == 1).collect();
    let geometry = rigidity::Geometry::new(&mesh, &active)?;
    let mut space = geometry.constrained(&rows(&mesh, |_| false))?;
    let mut first = nalgebra::DVector::zeros(6);
    first[0] = 1.;
    space.add(first.clone())?;
    let mut weak = first;
    weak[1] = 1e-9;
    weak[2] = 1e-16;
    space.add(weak)?;
    assert_eq!(
        space.report().restrained_modes,
        2,
        "a resolvable weak constraint was dropped"
    );
    let mut stronger = nalgebra::DVector::zeros(6);
    stronger[1] = 1.;
    space.add(stronger)?;
    assert_eq!(
        space.report().restrained_modes,
        2,
        "weak-pivot roundoff invented a third constraint"
    );
    assert_eq!(space.null_modes().len(), 4);
    Ok(())
}

#[test]
fn inactive_contact_histories_are_validated_before_reporting() -> Result<()> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/joint-column-mesh.json"))?;
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
    let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
    let u = vec![[0.; 3]; mesh.nodes.len()];
    let inactive = vec![false; mesh.tetrahedra.len()];
    let initial = model.initial_state()?;
    let mut corrupt = initial.clone();
    corrupt.points[0].traction_kpa[0] = f64::NAN;
    assert!(
        model
            .evaluate(&corrupt, &u, &inactive, None, false)
            .is_err()
    );
    let mut corrupt = initial;
    corrupt.points[0].normal = [1., 0., 0.];
    assert!(
        model
            .evaluate(&corrupt, &u, &inactive, None, false)
            .is_err()
    );
    Ok(())
}

#[test]
fn rigid_audit_rejects_malformed_geometry_and_nonlocal_or_nonorthogonal_transforms() -> Result<()> {
    let mesh: Mesh =
        serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
    let active = vec![true; mesh.tetrahedra.len()];
    let geometry = rigidity::Geometry::new(&mesh, &active)?;
    let transform = rows(&mesh, |_| false);
    assert!(geometry.constrained(&transform[..3]).is_err());
    let mut bad = transform.clone();
    bad[0][0].1 = f64::NAN;
    assert!(geometry.constrained(&bad).is_err());
    let mut bad = transform.clone();
    bad[3][0].0 = bad[0][0].0;
    assert!(geometry.constrained(&bad).is_err());
    let mut bad = transform;
    bad[0][0].1 = 0.5;
    assert!(geometry.constrained(&bad).is_err());
    let mut malformed = mesh;
    malformed.tetrahedra10.pop();
    assert!(rigidity::Geometry::new(&malformed, &active).is_err());
    Ok(())
}

#[test]
fn explicit_contact_supported_column_releases_upper_side_without_artificial_clamps() -> Result<()> {
    let mesh: Mesh =
        serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
    let ground = serde_json::from_str(include_str!("../examples/ground-joint-mc.json"))?;
    let case: StageCase = serde_json::from_str(include_str!(
        "../examples/stages-joint-explicit-normal-column.json"
    ))?;
    let run = soil_fem_3d::stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for f in &run.frames {
        assert_eq!(f.joint_support.as_ref().unwrap().unrestrained_modes, 0);
        assert!(f.maximum_constraint_error_m < 1e-12);
    }
    let e = 1e5;
    let nu = 0.3;
    let c = e * (1. - nu) / ((1. + nu) * (1. - 2. * nu));
    let initial = &run.frames[0];
    for &n in &mesh.node_sets["top"] {
        assert!((initial.state.displacement_m[n][2] - (-10. / c - 110. / 1e6)).abs() < 1e-10);
    }
    let final_frame = run.frames.last().unwrap();
    assert_eq!(
        final_frame
            .joint_support
            .as_ref()
            .unwrap()
            .continuum_components,
        1
    );
    assert!(
        final_frame
            .joint_results
            .as_ref()
            .unwrap()
            .points
            .iter()
            .all(|p| !p.active)
    );
    assert!(run.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
    assert!(
        run.attempts
            .iter()
            .any(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    Ok(())
}
