use nalgebra::{DVector, Matrix3, Rotation3, Unit, Vector3};
use serde_json::json;
use soil_fem_3d::{
    Result,
    boundary::{Conditions, SurfaceLoad},
    ground::Ground,
    mc, mechanics,
    mesh::Mesh,
    soil, srm_bound, srm_search,
    srm_virtual_work::{self as work, Assessment, Certificate, Problem},
    stages::{self, StageCase},
};

fn boundary(mesh: &Mesh, p: &Problem) -> Conditions {
    let mut b = Conditions {
        constraints: p.constraints.clone(),
        nodal_loads: vec![],
        surface_loads: vec![],
    };
    assert!(mesh.face_sets.contains_key("face-2-1"));
    b.surface_loads.push(SurfaceLoad::Traction {
        face_set: "face-2-1".into(),
        traction_kpa: [0., 0., -30.],
        gradient_kpa_per_m: [[0.; 3]; 3],
        reference_point_m: [0.; 3],
    });
    b
}
fn fixture(quadratic: bool) -> Result<(Mesh, Problem)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-virtual-work-mesh.json"))?;
    let mesh = if quadratic {
        mesh.elevate_straight()?
    } else {
        mesh
    };
    let mut p: Problem =
        serde_json::from_str(include_str!("../examples/mc-virtual-work-problem.json"))?;
    p.physical_nodal_force_kn = boundary(&mesh, &p).loads(&mesh, &p.active_elements)?;
    let phi = p.materials[0].parameters.friction_deg.to_radians();
    let ratio = (1. + phi.sin()) / (1. - phi.sin());
    p.virtual_displacement_m = mesh
        .nodes
        .iter()
        .map(|x| [(ratio + 0.2) * x[0], 0., -x[2]])
        .collect();
    Ok((mesh, p))
}
fn factor(mesh: &Mesh, p: &mut Problem, f: f64, slack: f64) {
    let params = &mut p.materials[0].parameters;
    params.cohesion_kpa = 10. / f;
    params.friction_deg = (30_f64.to_radians().tan() / f).atan().to_degrees();
    let phi = params.friction_deg.to_radians();
    let ratio = (1. + phi.sin()) / (1. - phi.sin());
    p.virtual_displacement_m = mesh
        .nodes
        .iter()
        .map(|x| [(ratio + slack) * x[0], 0., -x[2]])
        .collect();
}
fn cert(mesh: &Mesh, p: &Problem) -> Result<Certificate> {
    match work::analyze(mesh, p)? {
        Assessment::Evaluated { certificate } => Ok(*certificate),
        Assessment::NotApplicable { reason } => panic!("unexpected unavailable witness: {reason}"),
    }
}
fn tensor(e: &[f64]) -> Matrix3<f64> {
    Matrix3::new(
        e[0],
        e[3] / 2.,
        e[5] / 2.,
        e[3] / 2.,
        e[1],
        e[4] / 2.,
        e[5] / 2.,
        e[4] / 2.,
        e[2],
    )
}

#[test]
fn distributed_support_bounds_match_the_closed_form_for_tet4_and_tet10() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, mut p) = fixture(quadratic)?;
        for f in [1., 1.5] {
            factor(&mesh, &mut p, f, 0.2);
            let audit = work::audit(&mesh, &p)?;
            let Assessment::Evaluated { certificate: c } = audit.assessment else {
                panic!("evaluated")
            };
            let phi = p.materials[0].parameters.friction_deg.to_radians();
            let ratio = (1. + phi.sin()) / (1. - phi.sin());
            let expected = 10. / 30_f64.to_radians().tan() * (ratio - 1. + 0.2);
            assert!((c.external_virtual_work_kn_m - 30.).abs() < 1e-11);
            assert!((c.maximum_admissible_internal_virtual_work_kn_m - expected).abs() < 1e-10);
            assert_eq!(c.excludes_admissible_discrete_equilibrium, f == 1.5);
            assert!(
                audit.complete_point_audit
                    && audit.point_checks.iter().all(|q| q.accepted_finite_support)
            );
            assert_eq!(c.integration_points, if quadratic { 162 } else { 6 });
            assert_eq!(audit.point_checks.len(), c.integration_points);
            assert!(!c.continuum_collapse_or_stability_certified);
        }
        let mean = srm_bound::NodalProblem {
            schema: "soil-fem-mc-nodal-mean-stress-bound/1".into(),
            material_id: 1,
            water_policy: p.water_policy,
            material: p.materials[0].parameters.clone(),
            constraints: p.constraints.clone(),
            physical_nodal_force_kn: p.physical_nodal_force_kn.clone(),
            active_elements: p.active_elements.clone(),
        };
        assert!(matches!(
            srm_bound::analyze_nodal(&mesh, &mean)?,
            srm_bound::Assessment::NotApplicable { .. }
        ));
    }
    Ok(())
}

#[test]
fn actual_global_elastic_fem_and_independent_affine_solution_have_identical_virtual_work()
-> Result<()> {
    for quadratic in [false, true] {
        let (mesh, mut p) = fixture(quadratic)?;
        factor(&mesh, &mut p, 1., 0.2);
        let case = mechanics::Case {
            schema: "soil-fem-mechanics/1".into(),
            units: "m-kN-kPa".into(),
            materials: vec![mechanics::Material {
                id: 1,
                young_kpa: 100000.,
                poisson: 0.3,
                body_force_kn_m3: [0.; 3],
            }],
            boundary: boundary(&mesh, &p),
            relative_linear_tolerance: 1e-12,
            maximum_free_dofs: 1000,
            linear_solver: Default::default(),
        };
        let solution = mechanics::solve(&mesh, case)?;
        for (x, u) in mesh.nodes.iter().zip(&solution.displacement_m) {
            assert!((u[0] - 0.3 * 1.3 * 30. / 100000. * x[0]).abs() < 1e-13);
            assert!(u[1].abs() < 1e-13);
            assert!((u[2] + (1. - 0.3 * 0.3) * 30. / 100000. * x[2]).abs() < 1e-13);
        }
        let mut internal = 0.;
        let mut global = 0;
        for element in 0..mesh.tetrahedra.len() {
            let local = DVector::from_iterator(
                3 * mesh.cell_nodes(element).len(),
                mesh.cell_nodes(element)
                    .iter()
                    .flat_map(|&n| p.virtual_displacement_m[n]),
            );
            for q in mesh.integration(element)? {
                let e = &q.b * &local;
                let s = &solution.stress_tension_positive_kpa[global];
                assert!(
                    (s[0]).abs() < 1e-8 && (s[1] + 9.).abs() < 1e-8 && (s[2] + 30.).abs() < 1e-8
                );
                internal += q.weight_m3 * (0..6).map(|i| s[i] * e[i]).sum::<f64>();
                global += 1;
            }
        }
        let c = cert(&mesh, &p)?;
        assert!((internal - c.external_virtual_work_kn_m).abs() < 1e-9);
        assert!(!c.excludes_admissible_discrete_equilibrium);
    }
    Ok(())
}

#[test]
fn full_3d_tensor_rotation_preserves_work_capacity_and_eigenvalue_order() -> Result<()> {
    let (mesh, mut p) = fixture(true)?;
    factor(&mesh, &mut p, 1., 0.2);
    let original: srm_bound::Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    p.constraints = original.boundary.constraints.clone();
    let mut values = Vec::new();
    for rotation in [
        Matrix3::identity(),
        Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(1., 2., -1.)), 0.73)
            .into_inner(),
    ] {
        let rate =
            rotation * Matrix3::from_diagonal(&Vector3::new(3.2, 0., -1.)) * rotation.transpose();
        let stress =
            rotation * Matrix3::from_diagonal(&Vector3::new(0., 9., 30.)) * rotation.transpose();
        let skew = Matrix3::new(
            0.,
            rate[(1, 0)],
            rate[(2, 0)],
            -rate[(1, 0)],
            0.,
            rate[(2, 1)],
            -rate[(2, 0)],
            -rate[(2, 1)],
            0.,
        );
        let gradient = rate + skew;
        p.virtual_displacement_m = mesh
            .nodes
            .iter()
            .map(|x| (gradient * Vector3::from(*x)).into())
            .collect();
        let mut b = original.boundary.clone();
        for load in &mut b.surface_loads {
            if let SurfaceLoad::Traction {
                face_set,
                traction_kpa,
                ..
            } = load
            {
                let parts: Vec<_> = face_set.split('-').collect();
                let axis: usize = parts[1].parse().unwrap();
                let mut normal = Vector3::zeros();
                normal[axis] = if parts[2] == "0" { -1. } else { 1. };
                *traction_kpa = (-stress * normal).into();
            }
        }
        p.physical_nodal_force_kn = b.loads(&mesh, &p.active_elements)?;
        values.push(cert(&mesh, &p)?);
    }
    assert!(
        (values[0].external_virtual_work_kn_m - values[1].external_virtual_work_kn_m).abs() < 1e-10
    );
    assert!(
        (values[0].maximum_admissible_internal_virtual_work_kn_m
            - values[1].maximum_admissible_internal_virtual_work_kn_m)
            .abs()
            < 1e-10
    );
    assert_eq!(
        values[0].excludes_admissible_discrete_equilibrium,
        values[1].excludes_admissible_discrete_equilibrium
    );
    Ok(())
}

#[test]
fn positive_scaling_and_streamed_point_observer_preserve_the_exclusion() -> Result<()> {
    let (mesh, p) = fixture(true)?;
    let base = cert(&mesh, &p)?;
    for scale in [1e-6, 1., 1e6] {
        let mut scaled = p.clone();
        for u in &mut scaled.virtual_displacement_m {
            for x in u {
                *x *= scale;
            }
        }
        let c = cert(&mesh, &scaled)?;
        assert!(c.excludes_admissible_discrete_equilibrium);
        assert!(
            (c.external_virtual_work_kn_m / scale - base.external_virtual_work_kn_m).abs() < 1e-10
        );
        assert!(
            (c.maximum_admissible_internal_virtual_work_kn_m / scale
                - base.maximum_admissible_internal_virtual_work_kn_m)
                .abs()
                < 1e-10
        );
    }
    let audit = work::audit(&mesh, &p)?;
    let mut count = 0;
    let observed = work::analyze_with_observer(&mesh, &p, &mut |_| count += 1)?;
    let Assessment::Evaluated { certificate } = observed else {
        panic!("evaluated")
    };
    assert_eq!(count, audit.point_checks.len());
    assert_eq!(certificate.point_checks_sha256, base.point_checks_sha256);
    Ok(())
}

#[test]
fn even_tiny_supported_displacement_and_unsupported_oblique_reactions_cannot_make_an_upper()
-> Result<()> {
    let (mesh, p) = fixture(false)?;
    let mut bad = p.clone();
    bad.virtual_displacement_m[0][0] = 1e-15;
    assert!(matches!(
        work::analyze(&mesh, &bad)?,
        Assessment::NotApplicable { .. }
    ));
    bad = p.clone();
    bad.constraints[1].direction = [1., 1., 0.];
    assert!(matches!(
        work::analyze(&mesh, &bad)?,
        Assessment::NotApplicable { .. }
    ));
    bad = p;
    bad.constraints[0].value_m = 1e-12;
    assert!(matches!(
        work::analyze(&mesh, &bad)?,
        Assessment::NotApplicable { .. }
    ));
    Ok(())
}

#[test]
fn one_bad_gp_is_rejected_even_when_the_mean_strain_is_strictly_dual_feasible() -> Result<()> {
    let (mesh, mut p) = fixture(true)?;
    let center = mesh
        .nodes
        .iter()
        .position(|x| *x == [0.5, 0.5, 0.5])
        .unwrap();
    p.virtual_displacement_m[center][0] += 1.;
    let mut mean = Matrix3::zeros();
    let mut volume = 0.;
    for element in 0..mesh.tetrahedra.len() {
        let local = DVector::from_iterator(
            3 * mesh.cell_nodes(element).len(),
            mesh.cell_nodes(element)
                .iter()
                .flat_map(|&n| p.virtual_displacement_m[n]),
        );
        for q in mesh.integration(element)? {
            mean += tensor((&q.b * &local).as_slice()) * q.weight_m3;
            volume += q.weight_m3;
        }
    }
    mean /= volume;
    let phi = p.materials[0].parameters.friction_deg.to_radians();
    let ratio = (1. + phi.sin()) / (1. - phi.sin());
    let eigen = mean.symmetric_eigen().eigenvalues;
    assert!(mean.trace() + (ratio - 1.) * eigen.iter().map(|x| x.min(0.)).sum::<f64>() > 0.1);
    let audit = work::audit(&mesh, &p)?;
    assert!(matches!(audit.assessment, Assessment::NotApplicable { .. }));
    assert!(!audit.complete_point_audit);
    assert!(
        audit
            .point_checks
            .iter()
            .any(|q| !q.accepted_finite_support)
    );
    Ok(())
}

#[test]
fn heterogeneous_material_capacities_are_integrated_without_an_average_yield_domain() -> Result<()>
{
    let (mut mesh, mut p) = fixture(true)?;
    mesh.regions[0].volume_m3 = 0.5;
    mesh.regions[0].body = "FirstWedge".into();
    let mut second = mesh.regions[0].clone();
    second.id = 2;
    second.material_id = 2;
    second.body = "SecondWedge".into();
    mesh.regions.push(second);
    for i in 3..6 {
        mesh.material_ids[i] = 2;
        mesh.region_ids[i] = 2;
    }
    mesh.validate()?;
    let mut other = p.materials[0].clone();
    other.id = 2;
    other.parameters.cohesion_kpa = 2.;
    other.parameters.friction_deg = 10.;
    p.materials.push(other);
    let audit = work::audit(&mesh, &p)?;
    let Assessment::Evaluated { certificate: c } = audit.assessment else {
        panic!("evaluated")
    };
    let trace = p.virtual_displacement_m[1][0] - 1.;
    let expected = 0.5 * trace * (10. / 30_f64.to_radians().tan() + 2. / 10_f64.to_radians().tan());
    assert!((c.maximum_admissible_internal_virtual_work_kn_m - expected).abs() < 1e-10);
    assert!(c.excludes_admissible_discrete_equilibrium);
    assert_eq!(
        audit
            .point_checks
            .iter()
            .filter(|q| q.material_id == 1)
            .count(),
        81
    );
    assert_eq!(
        audit
            .point_checks
            .iter()
            .filter(|q| q.material_id == 2)
            .count(),
        81
    );
    Ok(())
}

#[test]
fn davis_effective_strength_and_rankine_superset_are_explicit_and_conservative() -> Result<()> {
    let (mesh, mut p) = fixture(false)?;
    factor(&mesh, &mut p, 1., 0.2);
    p.materials[0].parameters.flow = mc::Flow::Davis;
    let beta = 30_f64.to_radians().cos();
    let phi = (beta * 30_f64.to_radians().tan()).atan();
    let ratio = (1. + phi.sin()) / (1. - phi.sin());
    p.virtual_displacement_m = mesh
        .nodes
        .iter()
        .map(|x| [(ratio + 0.02) * x[0], 0., -x[2]])
        .collect();
    let c = cert(&mesh, &p)?;
    assert!(c.excludes_admissible_discrete_equilibrium);
    assert!(
        (c.maximum_admissible_internal_virtual_work_kn_m
            - 10. / 30_f64.to_radians().tan() * (ratio - 1. + 0.02))
            .abs()
            < 1e-10
    );
    let mut original = p.clone();
    original.materials[0].parameters.flow = mc::Flow::NonAssociated;
    assert!(matches!(
        work::analyze(&mesh, &original)?,
        Assessment::NotApplicable { .. }
    ));
    p.materials[0].parameters.tension_cutoff_kpa = Some(1.);
    let capped = cert(&mesh, &p)?;
    assert_eq!(
        capped.maximum_admissible_internal_virtual_work_kn_m,
        c.maximum_admissible_internal_virtual_work_kn_m
    );
    assert!(capped.pure_mc_superset_used_when_rankine_is_present);
    Ok(())
}

fn column(cohesion: f64, phi: f64, gamma: f64) -> Result<(Mesh, Ground, StageCase)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut g: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    for m in &mut g.materials {
        m.unit_weight_kn_m3 = gamma;
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            panic!("MC")
        };
        p.cohesion_kpa = cohesion;
        p.friction_deg = phi;
        p.dilation_deg = 0.;
        p.tension_cutoff_kpa = None;
        p.flow = mc::Flow::NonAssociated;
    }
    let mut c: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    c.stages = serde_json::from_value(
        json!([{"kind":"pressure","name":"load","target_pressure_kpa":100.0}]),
    )?;
    c.verification_reject_once_stages = vec![];
    Ok((mesh, g, c))
}

#[test]
fn actual_frictional_plastic_fem_internal_work_matches_external_work_and_is_not_excluded()
-> Result<()> {
    let (mesh, g, c) = column(1., 10., 20.)?;
    let run = stages::solve(&mesh, g.clone(), c.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.last_committed
            .material_states
            .iter()
            .any(|s| matches!(s,soil::State::MohrCoulomb(p) if p.equivalent_plastic_strain>0.))
    );
    let virtual_u: Vec<_> = mesh.nodes.iter().map(|x| [0., 0., x[2]]).collect();
    let problem = Problem {
        schema: "soil-fem-mc-virtual-work-bound/1".into(),
        water_policy: srm_bound::WaterPolicy::Dry,
        materials: g
            .materials
            .iter()
            .map(|m| {
                let soil::Model::MohrCoulomb(p) = &m.model else {
                    panic!("MC")
                };
                work::Material {
                    id: m.id,
                    parameters: p.clone(),
                }
            })
            .collect(),
        constraints: c.boundary.as_ref().unwrap().constraints.clone(),
        physical_nodal_force_kn: run.last_committed.applied_nodal_load_kn.clone(),
        active_elements: run.last_committed.active_elements.clone(),
        virtual_displacement_m: virtual_u,
    };
    let evaluated = cert(&mesh, &problem)?;
    assert!(!evaluated.excludes_admissible_discrete_equilibrium);
    let mut internal = 0.;
    let mut global = 0;
    for element in 0..mesh.tetrahedra.len() {
        let local = DVector::from_iterator(
            3 * mesh.cell_nodes(element).len(),
            mesh.cell_nodes(element)
                .iter()
                .flat_map(|&n| problem.virtual_displacement_m[n]),
        );
        for q in mesh.integration(element)? {
            let soil::State::MohrCoulomb(state) = &run.last_committed.material_states[global]
            else {
                panic!("MC")
            };
            let strain = &q.b * &local;
            internal -= q.weight_m3
                * (0..6)
                    .map(|i| state.effective_stress_compression_positive_kpa[i] * strain[i])
                    .sum::<f64>();
            global += 1;
        }
    }
    let frame = run.frames.last().unwrap();
    assert!(frame.relative_force_residual <= c.settings.relative_force_tolerance);
    let f = DVector::from_vec(problem.physical_nodal_force_kn.clone());
    let v = DVector::from_iterator(
        3 * mesh.nodes.len(),
        problem.virtual_displacement_m.iter().flatten().copied(),
    );
    let limit = frame.relative_force_residual * f.norm().max(1.) * v.norm()
        + 1e-10 * internal.abs().max(1.);
    assert!((internal - evaluated.external_virtual_work_kn_m).abs() <= limit);
    assert!(frame.maximum_normalized_yield <= 1e-8);
    Ok(())
}

#[test]
fn invalid_data_zero_strain_and_tresca_nonzero_trace_never_become_an_upper() -> Result<()> {
    let (mesh, p) = fixture(false)?;
    let mut bad = p.clone();
    bad.physical_nodal_force_kn.pop();
    assert!(work::analyze(&mesh, &bad).is_err());
    bad = p.clone();
    bad.materials.clear();
    assert!(work::analyze(&mesh, &bad).is_err());
    bad = p.clone();
    bad.materials.push(bad.materials[0].clone());
    assert!(work::analyze(&mesh, &bad).is_err());
    bad = p.clone();
    bad.virtual_displacement_m[1][0] = f64::NAN;
    assert!(work::analyze(&mesh, &bad).is_err());
    bad = p.clone();
    bad.physical_nodal_force_kn[15] = f64::MAX;
    assert!(work::analyze(&mesh, &bad).is_err());
    bad = p.clone();
    bad.active_elements = vec![true, false, false, false, false, false];
    assert!(work::analyze(&mesh, &bad).is_err());
    for phi in [0., 1e-10] {
        bad = p.clone();
        bad.materials[0].parameters.friction_deg = phi;
        assert!(matches!(
            work::analyze(&mesh, &bad)?,
            Assessment::NotApplicable { .. }
        ));
    }
    bad = p.clone();
    bad.virtual_displacement_m.fill([0.; 3]);
    assert!(matches!(
        work::analyze(&mesh, &bad)?,
        Assessment::NotApplicable { .. }
    ));
    let mut raw = serde_json::to_value(p)?;
    raw["water_policy"] = json!("unknown");
    assert!(serde_json::from_value::<Problem>(raw).is_err());
    Ok(())
}

#[test]
fn fresh_search_uses_a_verified_distributed_support_upper_and_retains_honest_numerical_unknowns()
-> Result<()> {
    let (mesh, mut p) = fixture(true)?;
    factor(&mesh, &mut p, 1., 0.2);
    let mut g: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = g["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    g["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
    g["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    g["materials"] = json!([{"id":1,"name":"Frictional unconfined cube","unit_weight_kn_m3":1e-10,"mohr_coulomb":p.materials[0].parameters}]);
    g["bodies"] = json!([{"body":"TestCube","material_id":1}]);
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(boundary(&mesh, &p));
    case.monitors.clear();
    case.monitor_sets.clear();
    case.stages =
        serde_json::from_value(json!([{"kind":"boundary","name":"load","target_factor":1.0}]))?;
    case.verification_reject_once_stages = vec![];
    let mut spec: srm_search::Specification =
        serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    spec.schema = "soil-fem-srm-search/2".into();
    spec.maximum_factor = 1.5;
    spec.required_factor = 1.05;
    factor(&mesh, &mut p, 1.105, 1e-4);
    spec.virtual_displacements_m = vec![p.virtual_displacement_m.clone()];
    let report = srm_search::run(&mesh, serde_json::from_value(g)?, case, spec)?;
    assert_eq!(
        report.requirement,
        srm_search::Requirement::EquilibriumRequirementMet
    );
    let b = report
        .interval
        .as_ref()
        .expect("a physically bound interval");
    let critical = 2. * (10. * (30. * 30_f64.to_radians().tan() + 10.)).sqrt() / 30.;
    assert!(b.attained_equilibrium_factor < critical && b.independently_excluded_factor > critical);
    for t in &report.trials {
        if t.outcome == srm_search::Outcome::PhysicallyExcluded {
            assert!(matches!(
                t.physical_exclusion,
                srm_bound::Assessment::NotApplicable { .. }
            ));
            assert!(t.virtual_work.iter().any(|v|matches!(&v.assessment,Assessment::Evaluated{certificate} if certificate.excludes_admissible_discrete_equilibrium)));
        }
        if t.outcome == srm_search::Outcome::Indeterminate {
            assert_ne!(t.factor, b.independently_excluded_factor);
        }
    }
    assert!(
        report.unreduced_source_history_unchanged
            && !report.continuum_stability_or_design_safety_certified
    );
    Ok(())
}

#[test]
fn real_nonlinear_failure_with_unavailable_virtual_fields_stays_indeterminate_without_an_upper()
-> Result<()> {
    let (mesh, mut g, mut c) = column(40., 0., 1e-10)?;
    g.initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    c.settings.maximum_newton_iterations = 3;
    c.settings.initial_increment = 1.;
    c.settings.minimum_increment = 1.;
    let mut spec: srm_search::Specification =
        serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    spec.schema = "soil-fem-srm-search/2".into();
    spec.material_ids = g.materials.iter().map(|m| m.id).collect();
    spec.maximum_factor = 4.;
    spec.required_factor = 4.;
    spec.virtual_displacements_m = vec![vec![[0.; 3]; mesh.nodes.len()]];
    let report = srm_search::run(&mesh, g, c, spec)?;
    assert!(report.interval.is_none());
    assert_eq!(report.requirement, srm_search::Requirement::Indeterminate);
    assert!(
        report
            .trials
            .iter()
            .all(|t| t.outcome != srm_search::Outcome::PhysicallyExcluded)
    );
    let failed = report.trials.iter().find(|t| t.factor == 4.).unwrap();
    assert_eq!(failed.outcome, srm_search::Outcome::Indeterminate);
    assert!(failed.reason.is_some());
    assert!(
        failed
            .virtual_work
            .iter()
            .all(|v| matches!(v.assessment, Assessment::NotApplicable { .. }))
    );
    Ok(())
}
