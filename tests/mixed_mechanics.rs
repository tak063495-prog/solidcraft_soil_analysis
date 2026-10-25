use nalgebra::Vector3;
use soil_fem_3d::{
    Result,
    boundary::{self, Conditions, Constraint, SurfaceLoad},
    elastic,
    linalg::{GmresOptions, LinearSolver},
    mechanics::{Case, Material},
    mesh::{Mesh, Region},
    mixed_mechanics,
};
use std::collections::BTreeMap;

fn cube(n: usize, layered: bool) -> Result<Mesh> {
    let mut m = Mesh {
        schema: "soil-fem-mesh/2".into(),
        units: "m-kN-kPa".into(),
        axis: "z_up".into(),
        crs: "LOCAL".into(),
        origin_m: [0.; 3],
        nodes: Vec::new(),
        tetrahedra: Vec::new(),
        tetrahedra10: Vec::new(),
        quadrature_order: None,
        material_ids: Vec::new(),
        node_sets: BTreeMap::new(),
        face_sets: BTreeMap::new(),
        regions: vec![Region {
            id: 1,
            body: "lower".into(),
            material_id: 1,
            bounds_m: [[0.; 3], [1., 1., if layered { 0.5 } else { 1. }]],
            volume_m3: if layered { 0.5 } else { 1. },
            source_sha256: "a".repeat(64),
        }],
        region_ids: Vec::new(),
        excavation_sets: BTreeMap::new(),
    };
    if layered {
        m.regions.push(Region {
            id: 2,
            body: "upper".into(),
            material_id: 2,
            bounds_m: [[0., 0., 0.5], [1.; 3]],
            volume_m3: 0.5,
            source_sha256: "b".repeat(64),
        });
    }
    for z in 0..=n {
        for y in 0..=n {
            for x in 0..=n {
                m.nodes.push([
                    x as f64 / n as f64,
                    y as f64 / n as f64,
                    z as f64 / n as f64,
                ]);
            }
        }
    }
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let node = |a, b, c| x + a + (n + 1) * (y + b + (n + 1) * (z + c));
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
                    m.tetrahedra.push(t.map(|i| c[i]));
                    let id = if layered && 2 * z >= n { 2 } else { 1 };
                    m.material_ids.push(id);
                    m.region_ids.push(id);
                }
            }
        }
    }
    for (name, axis, value) in [
        ("left", 0, 0.),
        ("front", 1, 0.),
        ("bottom", 2, 0.),
        ("top", 2, 1.),
    ] {
        m.node_sets.insert(
            name.into(),
            m.nodes
                .iter()
                .enumerate()
                .filter(|(_, p)| p[axis] == value)
                .map(|(i, _)| i)
                .collect(),
        );
        m.face_sets.insert(
            name.into(),
            boundary::surface_owners(&m)
                .iter()
                .filter(|(face, owners)| {
                    owners.len() == 1 && face.iter().all(|&i| m.nodes[i][axis] == value)
                })
                .map(|(&f, _)| f)
                .collect(),
        );
    }
    m.node_sets.insert(
        "boundary".into(),
        m.nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p.iter().any(|&v| v == 0. || v == 1.))
            .map(|(i, _)| i)
            .collect(),
    );
    m.node_sets
        .insert("all".into(), (0..m.nodes.len()).collect());
    let mut elevated = m.elevate_straight()?;
    elevated.node_sets.insert(
        "boundary".into(),
        elevated
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p.iter().any(|&v| v == 0. || v == 1.))
            .map(|(i, _)| i)
            .collect(),
    );
    Ok(elevated)
}
fn constraint(name: &str, axis: usize) -> Constraint {
    let mut direction = [0.; 3];
    direction[axis] = 1.;
    Constraint {
        node_set: name.into(),
        direction,
        value_m: 0.,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    }
}
fn case(poisson: f64) -> Case {
    Case {
        schema: "soil-fem-mechanics/1".into(),
        units: "m-kN-kPa".into(),
        materials: vec![Material {
            id: 1,
            young_kpa: 1000.,
            poisson,
            body_force_kn_m3: [0.; 3],
        }],
        boundary: Conditions::default(),
        relative_linear_tolerance: 1e-11,
        maximum_free_dofs: 2000,
        linear_solver: LinearSolver::DenseLu,
    }
}

#[test]
fn affine_strain_pressure_and_rigid_motion_are_exact_on_straight_and_curved_cells() -> Result<()> {
    for curved in [false, true] {
        for rigid in [false, true] {
            let mut m = cube(1, false)?;
            if curved {
                let i = m
                    .nodes
                    .iter()
                    .position(|p| *p == [0.5; 3])
                    .expect("interior edge midpoint");
                m.nodes[i] = [0.53, 0.48, 0.52];
            }
            let mut c = case(0.4999);
            let gradient = if rigid {
                [
                    [0., -0.003, 0.004],
                    [0.003, 0., -0.002],
                    [-0.004, 0.002, 0.],
                ]
            } else {
                [
                    [-0.001, 0.0003, 0.0002],
                    [0.0001, 0.002, -0.0004],
                    [0.0006, 0.0005, 0.0004],
                ]
            };
            for (a, row) in gradient.iter().enumerate() {
                let mut d = constraint("boundary", a);
                d.gradient = *row;
                d.value_m = 0.01 * (a as f64 + 1.);
                c.boundary.constraints.push(d);
            }
            let result = mixed_mechanics::solve(&m, c.clone())?;
            let strain = [
                gradient[0][0],
                gradient[1][1],
                gradient[2][2],
                gradient[0][1] + gradient[1][0],
                gradient[1][2] + gradient[2][1],
                gradient[0][2] + gradient[2][0],
            ];
            let stress =
                elastic::constitutive(1000., 0.4999)? * nalgebra::SVector::<f64, 6>::from(strain);
            let pressure = -stress.iter().take(3).sum::<f64>() / 3.;
            for (xyz, u) in m.nodes.iter().zip(&result.mechanics.displacement_m) {
                for a in 0..3 {
                    let expected = 0.01 * (a as f64 + 1.)
                        + Vector3::from(gradient[a]).dot(&Vector3::from(*xyz));
                    assert!(
                        (u[a] - expected).abs() < 2e-10,
                        "curved={curved}, rigid={rigid}, u={u:?}"
                    );
                }
            }
            for sig in &result.mechanics.stress_tension_positive_kpa {
                for a in 0..6 {
                    assert!(
                        (sig[a] - stress[a]).abs() < 2e-7,
                        "sig={sig:?}, exact={stress:?}"
                    );
                }
            }
            assert!(
                result
                    .pressure_dofs
                    .iter()
                    .all(|p| (p.compression_positive_kpa - pressure).abs() < 2e-7)
            );
            assert!(result.pressure_equation_relative_residual < 1e-10);
        }
    }
    Ok(())
}

#[test]
fn near_incompressible_uniaxial_stress_sparse_and_dense_match_theory() -> Result<()> {
    let m = cube(2, false)?;
    let mut c = case(0.4999);
    c.boundary.constraints = vec![
        constraint("left", 0),
        constraint("front", 1),
        constraint("bottom", 2),
    ];
    c.boundary.surface_loads = vec![SurfaceLoad::Pressure {
        face_set: "top".into(),
        pressure_kpa: 2.,
        gradient_kpa_per_m: [0.; 3],
        reference_point_m: [0.; 3],
    }];
    let dense = mixed_mechanics::solve(&m, c.clone())?;
    c.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions {
            restart: 100,
            maximum_iterations: 5000,
            ..Default::default()
        },
    };
    let sparse = mixed_mechanics::solve(&m, c)?;
    for (i, xyz) in m.nodes.iter().enumerate() {
        let expected = [
            2. * 0.4999 / 1000. * xyz[0],
            2. * 0.4999 / 1000. * xyz[1],
            -2. / 1000. * xyz[2],
        ];
        for (a, &value) in expected.iter().enumerate() {
            assert!((dense.mechanics.displacement_m[i][a] - value).abs() < 1e-10);
            assert!((sparse.mechanics.displacement_m[i][a] - value).abs() < 1e-10);
        }
    }
    assert!(
        sparse
            .pressure_dofs
            .iter()
            .all(|p| (p.compression_positive_kpa - 2. / 3.).abs() < 1e-7)
    );
    Ok(())
}

#[test]
fn two_material_column_preserves_pressure_jump_and_series_compliance() -> Result<()> {
    let m = cube(2, true)?;
    let mut c = case(0.2);
    c.materials.push(Material {
        id: 2,
        young_kpa: 1000.,
        poisson: 0.4999,
        body_force_kn_m3: [0.; 3],
    });
    c.boundary.constraints = vec![
        constraint("all", 0),
        constraint("all", 1),
        constraint("bottom", 2),
    ];
    c.boundary.surface_loads = vec![SurfaceLoad::Pressure {
        face_set: "top".into(),
        pressure_kpa: 1.,
        gradient_kpa_per_m: [0.; 3],
        reference_point_m: [0.; 3],
    }];
    let solution = mixed_mechanics::solve(&m, c)?;
    let compliance = |nu: f64| (1. + nu) * (1. - 2. * nu) / (1000. * (1. - nu));
    for (xyz, u) in m.nodes.iter().zip(&solution.mechanics.displacement_m) {
        let exact =
            -xyz[2].min(0.5) * compliance(0.2) - (xyz[2] - 0.5).max(0.) * compliance(0.4999);
        assert!((u[2] - exact).abs() < 1e-10);
    }
    for p in &solution.pressure_dofs {
        let nu = if p.material_id == 1 { 0.2 } else { 0.4999 };
        let exact = (1. + nu) / (3. * (1. - nu));
        assert!(
            (p.compression_positive_kpa - exact).abs() < 1e-8,
            "p={}, expected={exact}",
            p.compression_positive_kpa
        );
    }
    let interface = solution
        .pressure_dofs
        .iter()
        .filter(|p| m.nodes[p.node_index][2] == 0.5)
        .count();
    assert_eq!(interface, 18);
    Ok(())
}

#[test]
fn unsupported_linear_cells_underintegration_and_exact_incompressibility_are_rejected() -> Result<()>
{
    let mut m = cube(1, false)?;
    let c = case(0.3);
    m.quadrature_order = Some(1);
    assert!(mixed_mechanics::solve(&m, c.clone()).is_err());
    m.quadrature_order = Some(3);
    let mut bad = c.clone();
    bad.materials[0].poisson = 0.5;
    assert!(mixed_mechanics::solve(&m, bad).is_err());
    m.tetrahedra10.clear();
    m.schema = "soil-fem-mesh/2".into();
    assert!(mixed_mechanics::solve(&m, c).is_err());
    Ok(())
}
