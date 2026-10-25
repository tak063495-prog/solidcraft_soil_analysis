use nalgebra::{Rotation3, Vector3};
use soil_fem_3d::{
    Result,
    boundary::{self, Conditions, Constraint, NodalLoad, SurfaceLoad},
    linalg::{GmresOptions, LinearSolver},
    mechanics::{self, Case, Material},
    mesh::{Mesh, Region},
};
use std::collections::BTreeMap;

fn mesh(length: f64) -> Result<Mesh> {
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
        node_sets: BTreeMap::new(),
        face_sets: BTreeMap::new(),
        regions: vec![Region {
            id: 1,
            body: "Block".into(),
            material_id: 1,
            bounds_m: [[0., -0.5, -0.5], [length, 0.5, 0.5]],
            volume_m3: length,
            source_sha256: "a".repeat(64),
        }],
        region_ids: vec![],
        excavation_sets: BTreeMap::new(),
    };
    for z in [-0.5, 0.5] {
        for y in [-0.5, 0.5] {
            for x in [0., length / 2., length] {
                mesh.nodes.push([x, y, z]);
            }
        }
    }
    for x in 0usize..2 {
        let n = |a: usize, b: usize, c: usize| x + a + 3 * (b + 2 * c);
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
            mesh.tetrahedra.push(t.map(|i| c[i]));
            mesh.material_ids.push(1);
            mesh.region_ids.push(1);
        }
    }
    for (name, axis, value) in [
        ("left", 0, 0.),
        ("right", 0, length),
        ("ymin", 1, -0.5),
        ("ymax", 1, 0.5),
        ("bottom", 2, -0.5),
        ("top", 2, 0.5),
    ] {
        mesh.node_sets.insert(
            name.into(),
            mesh.nodes
                .iter()
                .enumerate()
                .filter(|(_, p)| p[axis] == value)
                .map(|(i, _)| i)
                .collect(),
        );
        mesh.face_sets.insert(
            name.into(),
            boundary::surface_owners(&mesh)
                .iter()
                .filter(|(f, owners)| {
                    owners.len() == 1 && f.iter().all(|&n| mesh.nodes[n][axis] == value)
                })
                .map(|(&f, _)| f)
                .collect(),
        );
    }
    mesh.elevate_straight()
}
fn constraint(name: &str, direction: [f64; 3], value: f64) -> Constraint {
    Constraint {
        node_set: name.into(),
        direction,
        value_m: value,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    }
}
fn case() -> Case {
    Case {
        schema: "soil-fem-mechanics/1".into(),
        units: "m-kN-kPa".into(),
        materials: vec![Material {
            id: 1,
            young_kpa: 1000.,
            poisson: 0.3,
            body_force_kn_m3: [0.; 3],
        }],
        boundary: Conditions::default(),
        relative_linear_tolerance: 1e-11,
        maximum_free_dofs: 10000,
        linear_solver: LinearSolver::SparseGmresIlu0 {
            options: GmresOptions::default(),
        },
    }
}
fn compression_constraints() -> Vec<Constraint> {
    vec![
        constraint("left", [1., 0., 0.], 0.),
        constraint("ymin", [0., 1., 0.], 0.),
        constraint("bottom", [0., 0., 1.], 0.),
    ]
}
#[test]
fn tet10_cantilever_pure_bending_matches_exact_quadratic_solution() -> Result<()> {
    let mesh = mesh(10.)?;
    let mut case = case();
    case.materials[0].poisson = 0.;
    case.boundary.constraints = (0..3)
        .map(|a| {
            let mut d = [0.; 3];
            d[a] = 1.;
            constraint("left", d, 0.)
        })
        .collect();
    let curvature = 1e-4;
    let mut gradient = [[0.; 3]; 3];
    gradient[0][1] = -case.materials[0].young_kpa * curvature;
    case.boundary.surface_loads = vec![SurfaceLoad::Traction {
        face_set: "right".into(),
        traction_kpa: [0.; 3],
        gradient_kpa_per_m: gradient,
        reference_point_m: [0.; 3],
    }];
    let solution = mechanics::solve(&mesh, case)?;
    for (p, u) in mesh.nodes.iter().zip(&solution.displacement_m) {
        let expected = [-curvature * p[0] * p[1], curvature * p[0] * p[0] / 2., 0.];
        for a in 0..3 {
            assert!(
                (u[a] - expected[a]).abs() < 1e-9,
                "u={u:?}, expected={expected:?}"
            );
        }
    }
    assert!(solution.linear.dimension > 0);
    assert!(
        solution
            .force_balance_error_kn
            .iter()
            .all(|f| f.abs() < 1e-9)
    );
    assert!(
        solution
            .moment_balance_error_kn_m
            .iter()
            .all(|f| f.abs() < 1e-8)
    );
    Ok(())
}
#[test]
fn signed_surface_pressure_and_rotated_direction_constraints_are_objective() -> Result<()> {
    let mut mesh = mesh(1.)?;
    let mut case = case();
    case.boundary.constraints = compression_constraints();
    case.boundary.surface_loads = vec![SurfaceLoad::Pressure {
        face_set: "top".into(),
        pressure_kpa: 1.,
        gradient_kpa_per_m: [0.; 3],
        reference_point_m: [0.; 3],
    }];
    let reference = mechanics::solve(&mesh, case.clone())?;
    for (p, u) in mesh.nodes.iter().zip(&reference.displacement_m) {
        assert!((u[0] - 0.0003 * p[0]).abs() < 1e-10);
        assert!((u[1] - 0.0003 * (p[1] + 0.5)).abs() < 1e-10);
        assert!((u[2] + 0.001 * (p[2] + 0.5)).abs() < 1e-10);
    }
    let rotation = Rotation3::from_euler_angles(0.4, -0.3, 0.7);
    for p in &mut mesh.nodes {
        *p = (rotation * Vector3::from_row_slice(p)).into();
    }
    for side in 0..2 {
        for a in 0..3 {
            mesh.regions[0].bounds_m[side][a] = mesh.nodes.iter().map(|p| p[a]).fold(
                if side == 0 {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                },
                |x, y| if side == 0 { x.min(y) } else { x.max(y) },
            );
        }
    }
    for c in &mut case.boundary.constraints {
        c.direction = (rotation * Vector3::from_row_slice(&c.direction)).into();
    }
    let solution = mechanics::solve(&mesh, case)?;
    for (a, b) in reference
        .displacement_m
        .iter()
        .zip(&solution.displacement_m)
    {
        assert!(
            (rotation * Vector3::from_row_slice(a) - Vector3::from_row_slice(b)).norm() < 1e-10
        );
    }
    assert!(solution.maximum_constraint_error_m < 1e-12);
    Ok(())
}
#[test]
fn nonzero_displacement_lifting_recovers_reaction_and_energy() -> Result<()> {
    let mesh = mesh(1.)?;
    let mut case = case();
    case.boundary.constraints = compression_constraints();
    case.boundary
        .constraints
        .push(constraint("top", [0., 0., 1.], -0.001));
    // Redundant scaled equations are accepted only when physically consistent.
    case.boundary
        .constraints
        .push(constraint("top", [0., 0., 2.], -0.002));
    let solution = mechanics::solve(&mesh, case.clone())?;
    let top_force: f64 = mesh.node_sets["top"]
        .iter()
        .map(|&n| solution.reaction_kn[n][2])
        .sum();
    assert!((top_force + 1.).abs() < 1e-9);
    assert!((solution.strain_energy_kn_m - 0.0005).abs() < 1e-11);
    assert!(solution.maximum_constraint_error_m < 1e-12);
    case.boundary
        .constraints
        .push(constraint("top", [0., 0., 1.], -0.002));
    assert!(mechanics::solve(&mesh, case).is_err());
    Ok(())
}
#[test]
fn nodal_force_totals_and_underconstrained_components_are_checked() -> Result<()> {
    let mesh = mesh(1.)?;
    let mut conditions = Conditions::default();
    conditions.nodal_loads.push(NodalLoad {
        node_set: "top".into(),
        total_force_kn: [1., 2., -3.],
    });
    let f = conditions.loads(&mesh, &vec![true; mesh.tetrahedra.len()])?;
    for a in 0..3 {
        assert!(
            ((a..f.len()).step_by(3).map(|i| f[i]).sum::<f64>() - [1., 2., -3.][a]).abs() < 1e-12
        );
    }
    let mut c = case();
    assert!(mechanics::solve(&mesh, c.clone()).is_err());
    c.boundary.constraints = compression_constraints();
    let mut detached = mesh.clone();
    let start = detached.nodes.len();
    let elements = detached.tetrahedra.len();
    detached
        .nodes
        .extend(mesh.nodes.iter().map(|p| [p[0] + 3., p[1], p[2]]));
    detached
        .tetrahedra
        .extend(mesh.tetrahedra.iter().map(|t| t.map(|n| n + start)));
    detached
        .tetrahedra10
        .extend(mesh.tetrahedra10.iter().map(|t| t.map(|n| n + start)));
    detached.material_ids.extend(vec![1; elements]);
    detached.region_ids.extend(vec![2; elements]);
    detached.regions.push(Region {
        id: 2,
        body: "Detached".into(),
        material_id: 1,
        bounds_m: [[3., -0.5, -0.5], [4., 0.5, 0.5]],
        volume_m3: 1.,
        source_sha256: "b".repeat(64),
    });
    assert!(mechanics::solve(&detached, c).is_err());
    Ok(())
}

#[test]
fn two_material_series_bar_has_correct_interface_force_and_compliance() -> Result<()> {
    let mut mesh = mesh(1.)?;
    mesh.material_ids = (0..mesh.tetrahedra.len())
        .map(|i| if i < 6 { 1 } else { 2 })
        .collect();
    mesh.region_ids = mesh.material_ids.clone();
    mesh.regions[0].bounds_m[1][0] = 0.5;
    mesh.regions[0].volume_m3 = 0.5;
    mesh.regions.push(Region {
        id: 2,
        body: "Right".into(),
        material_id: 2,
        bounds_m: [[0.5, -0.5, -0.5], [1., 0.5, 0.5]],
        volume_m3: 0.5,
        source_sha256: "b".repeat(64),
    });
    let mut c = case();
    c.materials[0].poisson = 0.;
    c.materials.push(Material {
        id: 2,
        young_kpa: 2000.,
        poisson: 0.,
        body_force_kn_m3: [0.; 3],
    });
    c.boundary.constraints = compression_constraints();
    c.boundary.surface_loads.push(SurfaceLoad::Pressure {
        face_set: "right".into(),
        pressure_kpa: 1.,
        gradient_kpa_per_m: [0.; 3],
        reference_point_m: [0.; 3],
    });
    let result = mechanics::solve(&mesh, c)?;
    for (p, u) in mesh.nodes.iter().zip(&result.displacement_m) {
        let expected = -p[0].min(0.5) / 1000. - (p[0] - 0.5).max(0.) / 2000.;
        assert!((u[0] - expected).abs() < 1e-10);
    }
    assert!(
        result
            .stress_tension_positive_kpa
            .iter()
            .all(|s| (s[0] + 1.).abs() < 1e-8)
    );
    Ok(())
}
