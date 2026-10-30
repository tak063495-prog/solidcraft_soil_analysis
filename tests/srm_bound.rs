use nalgebra::{Matrix3, Rotation3, Unit, Vector3};
use soil_fem_3d::{
    Result,
    boundary::{self, Conditions, Constraint, SurfaceLoad},
    mc, mechanics,
    mesh::{Mesh, Region},
    srm_bound::{self, Assessment, Certificate, Problem},
};
use std::collections::BTreeMap;

fn fixture(stress: Matrix3<f64>) -> Result<(Mesh, Problem)> {
    let nodes: Vec<_> = (0..8)
        .map(|i| [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64])
        .collect();
    let mut tetrahedra = Vec::new();
    for [a, b, c] in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let mut cell = [
            0,
            1 << a,
            (1 << a) | (1 << b),
            (1 << a) | (1 << b) | (1 << c),
        ];
        let p: Vec<_> = cell.iter().map(|&i| Vector3::from(nodes[i])).collect();
        if Matrix3::from_columns(&[p[1] - p[0], p[2] - p[0], p[3] - p[0]]).determinant() < 0. {
            cell.swap(1, 2);
        }
        tetrahedra.push(cell);
    }
    let source_sha256 = soil_fem_3d::checkpoint::hash_json(&(&nodes, &tetrahedra))?;
    let mut mesh = Mesh {
        schema: "soil-fem-mesh/2".into(),
        units: "m-kN-kPa".into(),
        axis: "z_up".into(),
        crs: "LOCAL".into(),
        origin_m: [0.; 3],
        nodes,
        tetrahedra,
        material_ids: vec![1; 6],
        tetrahedra10: vec![],
        quadrature_order: None,
        node_sets: BTreeMap::from([
            ("origin".into(), vec![0]),
            ("x".into(), vec![1]),
            ("y".into(), vec![2]),
            ("all".into(), (0..8).collect()),
        ]),
        face_sets: BTreeMap::new(),
        regions: vec![Region {
            id: 1,
            body: "TestCube".into(),
            material_id: 1,
            bounds_m: [[0.; 3], [1.; 3]],
            volume_m3: 1.,
            source_sha256,
        }],
        region_ids: vec![1; 6],
        excavation_sets: BTreeMap::new(),
    };
    mesh.validate()?;
    let faces = boundary::surface_owners(&mesh);
    let mut loads = Vec::new();
    for axis in 0..3 {
        for side in 0..2 {
            let name = format!("face-{axis}-{side}");
            let selected = faces
                .iter()
                .filter(|(f, o)| {
                    o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][axis] == side as f64)
                })
                .map(|(f, _)| *f)
                .collect();
            mesh.face_sets.insert(name.clone(), selected);
            let mut normal = Vector3::zeros();
            normal[axis] = if side == 0 { -1. } else { 1. };
            let traction = -stress * normal;
            loads.push(SurfaceLoad::Traction {
                face_set: name,
                traction_kpa: traction.into(),
                gradient_kpa_per_m: [[0.; 3]; 3],
                reference_point_m: [0.; 3],
            });
        }
    }
    let mut constraints = Vec::new();
    for (name, axes) in [("origin", vec![0, 1, 2]), ("x", vec![1, 2]), ("y", vec![2])] {
        for axis in axes {
            let mut direction = [0.; 3];
            direction[axis] = 1.;
            constraints.push(Constraint {
                node_set: name.into(),
                direction,
                value_m: 0.,
                gradient: [0.; 3],
                reference_point_m: [0.; 3],
            });
        }
    }
    let problem = Problem {
        schema: "soil-fem-mc-mean-stress-bound/1".into(),
        material_id: 1,
        water_policy: soil_fem_3d::srm_bound::WaterPolicy::Dry,
        material: mc::Parameters {
            young_kpa: 100000.,
            poisson: 0.3,
            cohesion_kpa: 10.,
            friction_deg: 0.,
            dilation_deg: 0.,
            tension_cutoff_kpa: None,
            flow: mc::Flow::NonAssociated,
        },
        boundary: Conditions {
            constraints,
            nodal_loads: vec![],
            surface_loads: loads,
        },
        body_force_kn_m3: [0.; 3],
        active_elements: None,
    };
    Ok((mesh, problem))
}
fn evaluated(mesh: &Mesh, p: &Problem) -> Result<Certificate> {
    match srm_bound::analyze(mesh, p)? {
        Assessment::Evaluated { certificate } => Ok(*certificate),
        Assessment::NotApplicable { reason } => panic!("unexpected unavailable witness: {reason}"),
    }
}
fn matrix(c: &Certificate) -> Matrix3<f64> {
    let s = c.mean_effective_stress_compression_positive_kpa;
    Matrix3::new(s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2])
}

#[test]
fn uniaxial_closed_form_excludes_only_after_the_strength_limit_without_a_solver_failure()
-> Result<()> {
    let (mesh, mut p) = fixture(Matrix3::from_diagonal(&Vector3::new(0., 0., 15.)))?;
    let below = evaluated(&mesh, &p)?;
    assert!(!below.mean_stress_excludes_admissible_discrete_equilibrium);
    assert!((below.mean_mc_yield_margin_kpa + 5.).abs() < 1e-10);
    p.material.cohesion_kpa = 10. / 1.5;
    let above = evaluated(&mesh, &p)?;
    assert!(above.mean_stress_excludes_admissible_discrete_equilibrium);
    assert!((above.mean_mc_yield_margin_kpa - (15. - 20. / 1.5)).abs() < 1e-10);
    assert!(!above.lower_bound_or_stability_proven && !above.continuum_collapse_factor_certified);
    assert!(
        above
            .gauge_reactions
            .iter()
            .all(|r| r.magnitude_kn.abs() < 1e-10)
    );
    Ok(())
}

#[test]
fn mean_stress_agrees_with_an_independent_elastic_global_fem_and_tet10() -> Result<()> {
    let expected = Matrix3::new(6., 1., 3., 1., 8., 2., 3., 2., 15.);
    let (mesh, p) = fixture(expected)?;
    for m in [mesh.clone(), mesh.elevate_straight()?] {
        let certificate = evaluated(&m, &p)?;
        assert!((matrix(&certificate) - expected).norm() < 1e-10);
        let case = mechanics::Case {
            schema: "soil-fem-mechanics/1".into(),
            units: "m-kN-kPa".into(),
            materials: vec![mechanics::Material {
                id: 1,
                young_kpa: 100000.,
                poisson: 0.3,
                body_force_kn_m3: [0.; 3],
            }],
            boundary: p.boundary.clone(),
            relative_linear_tolerance: 1e-12,
            maximum_free_dofs: 1000,
            linear_solver: Default::default(),
        };
        let solution = mechanics::solve(&m, case)?;
        let mut mean = [0.; 6];
        let mut volume = 0.;
        for (q, s) in solution
            .integration_points
            .iter()
            .zip(&solution.stress_tension_positive_kpa)
        {
            for axis in 0..6 {
                mean[axis] -= q.weight_m3 * s[axis];
            }
            volume += q.weight_m3;
        }
        for (a, b) in mean
            .iter()
            .zip(certificate.mean_effective_stress_compression_positive_kpa)
        {
            assert!((a / volume - b).abs() < 1e-8);
        }
        assert!(solution.true_relative_free_force_residual < 1e-11);
    }
    Ok(())
}

#[test]
fn certificate_is_rotation_translation_and_nonassociated_flow_invariant() -> Result<()> {
    let stress = Matrix3::new(6., 1., 3., 1., 8., 2., 3., 2., 15.);
    let (mut mesh, mut p) = fixture(stress)?;
    p.material.friction_deg = 20.;
    let baseline = evaluated(&mesh, &p)?;
    let q = Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(1., 2., -1.)), 0.73)
        .into_inner();
    let shift = Vector3::new(1000., -2000., 3000.);
    for x in &mut mesh.nodes {
        *x = (q * Vector3::from(*x) + shift).into();
    }
    mesh.regions[0].bounds_m = std::array::from_fn(|side| {
        std::array::from_fn(|axis| {
            mesh.nodes.iter().map(|x| x[axis]).fold(
                if side == 0 {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                },
                if side == 0 { f64::min } else { f64::max },
            )
        })
    });
    mesh.regions[0].source_sha256 =
        soil_fem_3d::checkpoint::hash_json(&(&mesh.nodes, &mesh.tetrahedra))?;
    for c in &mut p.boundary.constraints {
        c.direction = (q * Vector3::from(c.direction)).into();
    }
    for load in &mut p.boundary.surface_loads {
        if let SurfaceLoad::Traction { traction_kpa, .. } = load {
            *traction_kpa = (q * Vector3::from(*traction_kpa)).into();
        }
    }
    let rotated = evaluated(&mesh, &p)?;
    assert!((matrix(&rotated) - q * stress * q.transpose()).norm() < 1e-8);
    assert!((rotated.mean_mc_yield_margin_kpa - baseline.mean_mc_yield_margin_kpa).abs() < 1e-8);
    assert_eq!(
        rotated.mean_stress_excludes_admissible_discrete_equilibrium,
        baseline.mean_stress_excludes_admissible_discrete_equilibrium
    );
    Ok(())
}

#[test]
fn body_weight_and_all_six_reactions_enter_the_first_moment() -> Result<()> {
    let (mesh, mut p) = fixture(Matrix3::zeros())?;
    p.body_force_kn_m3 = [0., 0., -20.];
    let c = evaluated(&mesh, &p)?;
    assert!((c.volume_m3 - 1.).abs() < 1e-12);
    assert!((matrix(&c) - Matrix3::from_diagonal(&Vector3::new(0., 0., 10.))).norm() < 1e-10);
    let total = c.gauge_reactions.iter().fold(Vector3::zeros(), |a, r| {
        a + Vector3::from(r.unit_direction) * r.magnitude_kn
    });
    assert!((total - Vector3::new(0., 0., 20.)).norm() < 1e-10);
    assert!(c.relative_resultant_force_and_moment_residual < 1e-12);
    Ok(())
}

#[test]
fn distributed_or_singular_supports_cannot_be_misreported_as_physical_uppers() -> Result<()> {
    let (mesh, mut p) = fixture(Matrix3::from_diagonal(&Vector3::new(0., 0., 100.)))?;
    p.boundary.constraints[0].node_set = "all".into();
    assert!(matches!(
        srm_bound::analyze(&mesh, &p)?,
        Assessment::NotApplicable { .. }
    ));
    let (_, mut p) = fixture(Matrix3::from_diagonal(&Vector3::new(0., 0., 100.)))?;
    p.boundary.constraints[3] = p.boundary.constraints[0].clone();
    assert!(matches!(
        srm_bound::analyze(&mesh, &p)?,
        Assessment::NotApplicable { .. }
    ));
    Ok(())
}

#[test]
fn rankine_exclusion_is_distinct_from_mc_and_not_an_admissibility_sufficiency_claim() -> Result<()>
{
    let (mesh, mut p) = fixture(Matrix3::from_diagonal(&Vector3::new(0., 0., -2.)))?;
    let pure = evaluated(&mesh, &p)?;
    assert!(!pure.mean_stress_excludes_admissible_discrete_equilibrium);
    p.material.tension_cutoff_kpa = Some(1.);
    let capped = evaluated(&mesh, &p)?;
    assert!(capped.mean_stress_excludes_admissible_discrete_equilibrium);
    assert!(capped.mean_mc_yield_margin_kpa < 0.);
    assert!((capped.mean_rankine_yield_margin_kpa.unwrap() - 1.).abs() < 1e-10);
    assert!(!pure.lower_bound_or_stability_proven);
    Ok(())
}

#[test]
fn invalid_or_heterogeneous_inputs_never_become_an_upper_candidate() -> Result<()> {
    let (mut mesh, mut p) = fixture(Matrix3::zeros())?;
    p.active_elements = Some(vec![true]);
    assert!(srm_bound::analyze(&mesh, &p).is_err());
    p.active_elements = None;
    p.body_force_kn_m3[0] = f64::NAN;
    assert!(srm_bound::analyze(&mesh, &p).is_err());
    p.body_force_kn_m3 = [0.; 3];
    mesh.material_ids[0] = 2;
    mesh.region_ids[0] = 2;
    mesh.regions[0].volume_m3 = 5. / 6.;
    mesh.regions.push(Region {
        id: 2,
        body: "OtherMaterial".into(),
        material_id: 2,
        volume_m3: 1. / 6.,
        ..mesh.regions[0].clone()
    });
    mesh.validate()?;
    assert!(matches!(
        srm_bound::analyze(&mesh, &p)?,
        Assessment::NotApplicable { .. }
    ));
    Ok(())
}

#[test]
fn davis_uses_its_effective_strength_instead_of_the_untransformed_cohesion() -> Result<()> {
    let (mesh, mut p) = fixture(Matrix3::from_diagonal(&Vector3::new(0., 0., 30.)))?;
    p.material.friction_deg = 30.;
    let original = evaluated(&mesh, &p)?;
    assert!(!original.mean_stress_excludes_admissible_discrete_equilibrium);
    p.material.flow = mc::Flow::Davis;
    let changed = evaluated(&mesh, &p)?;
    let beta = 30_f64.to_radians().cos();
    let phi = (beta * 30_f64.to_radians().tan()).atan();
    let expected = 30. * (1. - phi.sin()) - 2. * 10. * beta * phi.cos();
    assert!((changed.mean_mc_yield_margin_kpa - expected).abs() < 1e-10);
    assert!(changed.mean_stress_excludes_admissible_discrete_equilibrium);
    Ok(())
}

#[test]
fn unknown_water_policy_and_overflow_do_not_create_physical_exclusion() -> Result<()> {
    let (mesh, mut p) = fixture(Matrix3::zeros())?;
    let mut raw = serde_json::to_value(&p)?;
    raw["water_policy"] = serde_json::json!("fixed_pore_pressure");
    assert!(serde_json::from_value::<Problem>(raw).is_err());
    for load in &mut p.boundary.surface_loads {
        if let SurfaceLoad::Traction { traction_kpa, .. } = load {
            *traction_kpa = [f64::MAX; 3];
        }
    }
    assert!(
        !matches!(srm_bound::analyze(&mesh,&p),Ok(Assessment::Evaluated{certificate}) if certificate.mean_stress_excludes_admissible_discrete_equilibrium)
    );
    Ok(())
}
