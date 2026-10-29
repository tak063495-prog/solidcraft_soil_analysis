use nalgebra::{Matrix3, Rotation3, SVector, Vector3};
use soil_fem_3d::{
    Result,
    mc::{Flow, Parameters, Return, State},
};

fn parameters(c: f64, phi: f64, psi: f64) -> Parameters {
    Parameters {
        young_kpa: 30_000.,
        poisson: 0.25,
        cohesion_kpa: c,
        friction_deg: phi,
        dilation_deg: psi,
        tension_cutoff_kpa: None,
        flow: Flow::NonAssociated,
    }
}
fn state(stress: [f64; 6]) -> State {
    State {
        effective_stress_compression_positive_kpa: stress,
        plastic_strain_compression_positive: [0.003, -0.001, -0.002, 0., 0., 0.],
        void_ratio: 1.,
        equivalent_plastic_strain: 0.02,
        plastic_work_kpa: 1.,
        last_return: Return::Elastic,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:e} != {b:e}");
}
fn normal(p: &Parameters) -> Matrix3<f64> {
    let g = p.young_kpa / (2. * (1. + p.poisson));
    let l = p.young_kpa * p.poisson / ((1. + p.poisson) * (1. - 2. * p.poisson));
    Matrix3::from_fn(|i, j| l + if i == j { 2. * g } else { 0. })
}
fn tensor(s: [f64; 6]) -> Matrix3<f64> {
    Matrix3::new(s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2])
}
fn stress(m: Matrix3<f64>) -> [f64; 6] {
    [
        m[(0, 0)],
        m[(1, 1)],
        m[(2, 2)],
        m[(0, 1)],
        m[(1, 2)],
        m[(2, 0)],
    ]
}

#[test]
fn zero_strain_tresca_weakening_matches_closed_form_and_preserves_prior_history() -> Result<()> {
    let from = parameters(20., 0., 0.);
    let to = parameters(10., 0., 0.);
    let old = state([120., 100., 80., 0., 0., 0.]);
    let before = serde_json::to_vec(&old)?;
    from.validate_state(&old)?;
    assert!(to.integrate_partition(&old, [0.; 6], 1).is_err());
    let u = to.integrate_transition_partition(&from, &old, [0.; 6], 1)?;
    let a = 20. / (4. * 12_000.);
    for (actual, expected) in u
        .state
        .effective_stress_compression_positive_kpa
        .iter()
        .zip([110., 100., 90., 0., 0., 0.])
    {
        near(*actual, expected, 1e-10);
    }
    for i in 0..6 {
        near(
            u.state.plastic_strain_compression_positive[i],
            old.plastic_strain_compression_positive[i] + [a, 0., -a, 0., 0., 0.][i],
            1e-13,
        );
    }
    near(
        u.state.equivalent_plastic_strain,
        old.equivalent_plastic_strain + 2. * a / 3_f64.sqrt(),
        1e-13,
    );
    near(
        u.state.plastic_work_kpa,
        old.plastic_work_kpa + 20. * a,
        1e-12,
    );
    near(u.state.void_ratio, old.void_ratio, 1e-15);
    assert_eq!(u.plastic_substeps, 1);
    assert_eq!(serde_json::to_vec(&old)?, before);
    to.validate_state(&u.state)?;
    Ok(())
}

#[test]
fn rankine_tension_parameter_jump_matches_uniaxial_multiplier_without_clearing_history()
-> Result<()> {
    let mut from = parameters(100., 0., 0.);
    from.tension_cutoff_kpa = Some(20.);
    let mut to = from.clone();
    to.tension_cutoff_kpa = Some(5.);
    let old = state([0., 0., -10., 0., 0., 0.]);
    from.validate_state(&old)?;
    let u = to.integrate_transition_partition(&from, &old, [0.; 6], 1)?;
    let a = 5. / 36_000.;
    for (x, y) in u
        .state
        .effective_stress_compression_positive_kpa
        .iter()
        .zip([12_000. * a, 12_000. * a, -5., 0., 0., 0.])
    {
        near(*x, y, 1e-10);
    }
    near(
        u.state.plastic_strain_compression_positive[2],
        old.plastic_strain_compression_positive[2] - a,
        1e-13,
    );
    near(
        u.state.plastic_work_kpa,
        old.plastic_work_kpa + 5. * a,
        1e-12,
    );
    near(
        u.state.equivalent_plastic_strain,
        old.equivalent_plastic_strain + 2. * a / 3.,
        1e-13,
    );
    assert_eq!(u.state.last_return, Return::TensionFace);
    to.validate_state(&u.state)?;
    Ok(())
}

#[test]
fn nonassociated_face_rotation_and_algorithmic_tangent_match_independent_formulas() -> Result<()> {
    let from = parameters(20., 30., 10.);
    let to = parameters(10., 20., 10.);
    let old = state([180., 100., 70., 0., 0., 0.]);
    from.validate_state(&old)?;
    let n = Vector3::new(
        1. - 20_f64.to_radians().sin(),
        0.,
        -1. - 20_f64.to_radians().sin(),
    );
    let m = Vector3::new(
        1. - 10_f64.to_radians().sin(),
        0.,
        -1. - 10_f64.to_radians().sin(),
    );
    let s = Vector3::new(180., 100., 70.);
    let d = normal(&to);
    let multiplier = (n.dot(&s) - 20. * 20_f64.to_radians().cos()) / n.dot(&(d * m));
    assert!(multiplier > 0.);
    let expected = s - d * m * multiplier;
    let u = to.integrate_transition_partition(&from, &old, [0.; 6], 1)?;
    for i in 0..3 {
        near(
            u.state.effective_stress_compression_positive_kpa[i],
            expected[i],
            1e-9,
        );
    }
    let r = Rotation3::from_euler_angles(0.31, -0.27, 0.41).into_inner();
    let mut rotated = old.clone();
    rotated.effective_stress_compression_positive_kpa =
        stress(r * tensor(old.effective_stress_compression_positive_kpa) * r.transpose());
    let ep = tensor([
        old.plastic_strain_compression_positive[0],
        old.plastic_strain_compression_positive[1],
        old.plastic_strain_compression_positive[2],
        0.,
        0.,
        0.,
    ]);
    let er = r * ep * r.transpose();
    let ev = stress(er);
    rotated.plastic_strain_compression_positive =
        std::array::from_fn(|i| ev[i] * if i < 3 { 1. } else { 2. });
    let v = to.integrate_transition_partition(&from, &rotated, [0.; 6], 1)?;
    assert!(
        (tensor(v.state.effective_stress_compression_positive_kpa)
            - r * tensor(u.state.effective_stress_compression_positive_kpa) * r.transpose())
        .norm()
            < 1e-8
    );
    let (_, tangent) = to.stress_tangent_transition_partition(&from, &rotated, [0.; 6], 1, 1e-7)?;
    let direction = [0.27, -0.13, 0.08, 0.21, -0.19, 0.09];
    let h = 4e-8;
    let a = to
        .integrate_transition_partition(&from, &rotated, direction.map(|x| h * x), 1)?
        .state;
    let b = to
        .integrate_transition_partition(&from, &rotated, direction.map(|x| -h * x), 1)?
        .state;
    let direct = SVector::<f64, 6>::from_fn(|i, _| {
        (a.effective_stress_compression_positive_kpa[i]
            - b.effective_stress_compression_positive_kpa[i])
            / (2. * h)
    });
    assert!((tangent * SVector::from_row_slice(&direction) - direct).norm() / direct.norm() < 3e-5);
    Ok(())
}

#[test]
fn unchanged_parameters_match_ordinary_updates_and_davis_matches_explicit_effective_model()
-> Result<()> {
    let p = parameters(20., 30., 10.);
    let old = state([180., 100., 70., 0., 0., 0.]);
    let de = [1e-5, -2e-6, 0., 3e-6, 0., 0.];
    let a = p.integrate_partition(&old, de, 4)?;
    let b = p.integrate_transition_partition(&p, &old, de, 4)?;
    assert_eq!(serde_json::to_vec(&a)?, serde_json::to_vec(&b)?);
    let plastic_de = [0.02, -0.004, -0.006, 0.003, 0., 0.];
    let a = p.integrate_partition(&old, plastic_de, 32)?;
    let b = p.integrate_transition_partition(&p, &old, plastic_de, 32)?;
    assert!(a.plastic_substeps > 0);
    assert_eq!(serde_json::to_vec(&a)?, serde_json::to_vec(&b)?);
    let mut from = p.clone();
    from.flow = Flow::Davis;
    let mut to = parameters(10., 20., 10.);
    to.flow = Flow::Davis;
    let effective = |p: &Parameters| {
        let phi = p.friction_deg.to_radians();
        let psi = p.dilation_deg.to_radians();
        let beta = phi.cos() * psi.cos() / (1. - phi.sin() * psi.sin());
        let angle = (beta * phi.tan()).atan().to_degrees();
        Parameters {
            cohesion_kpa: beta * p.cohesion_kpa,
            friction_deg: angle,
            dilation_deg: angle,
            flow: Flow::Associated,
            ..p.clone()
        }
    };
    let a = to.integrate_transition_partition(&from, &old, [0.; 6], 1)?;
    let b = effective(&to).integrate_transition_partition(&effective(&from), &old, [0.; 6], 1)?;
    for i in 0..6 {
        near(
            a.state.effective_stress_compression_positive_kpa[i],
            b.state.effective_stress_compression_positive_kpa[i],
            1e-9,
        );
        near(
            a.state.plastic_strain_compression_positive[i],
            b.state.plastic_strain_compression_positive[i],
            1e-12,
        );
    }
    Ok(())
}

#[test]
fn invalid_origins_parameter_switches_and_genuine_no_solution_leave_inputs_unchanged() -> Result<()>
{
    let from = parameters(20., 30., 0.);
    let old = state([-10., -10., -10., 0., 0., 0.]);
    from.validate_state(&old)?;
    let saved = serde_json::to_vec(&old)?;
    let mut to = from.clone();
    to.cohesion_kpa = 0.;
    assert!(
        to.integrate_transition_partition(&from, &old, [0.; 6], 1)
            .is_err()
    );
    assert_eq!(serde_json::to_vec(&old)?, saved);
    for steps in [0, 100_001] {
        assert!(
            from.integrate_transition_partition(&from, &old, [0.; 6], steps)
                .is_err()
        );
    }
    assert!(
        from.integrate_transition_partition(&from, &old, [f64::NAN; 6], 1)
            .is_err()
    );
    for option in 0..3 {
        let mut changed = from.clone();
        match option {
            0 => changed.young_kpa *= 2.,
            1 => changed.poisson = 0.2,
            _ => {
                changed.flow = Flow::Associated;
                changed.dilation_deg = changed.friction_deg;
            }
        };
        changed.validate()?;
        assert!(
            changed
                .integrate_transition_partition(&from, &old, [0.; 6], 1)
                .is_err()
        );
    }
    let invalid = state([1000., 0., 0., 0., 0., 0.]);
    assert!(
        from.integrate_transition_partition(&from, &invalid, [0.; 6], 1)
            .is_err()
    );
    assert!(
        from.stress_tangent_transition_partition(&from, &old, [0.; 6], 1, 0.)
            .is_err()
    );
    assert_eq!(serde_json::to_vec(&old)?, saved);
    Ok(())
}
