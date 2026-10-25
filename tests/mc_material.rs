use nalgebra::{Matrix3, SVector, Vector3};
use soil_fem_3d::{
    Result, elastic,
    mc::{self, Flow, Parameters, Return, State},
};

fn parameters(psi: f64) -> Parameters {
    Parameters {
        young_kpa: 10000.,
        poisson: 0.25,
        cohesion_kpa: 10.,
        friction_deg: 30.,
        dilation_deg: psi,
        tension_cutoff_kpa: None,
        flow: Flow::NonAssociated,
    }
}
fn initial() -> State {
    State {
        effective_stress_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
        plastic_strain_compression_positive: [0.; 6],
        void_ratio: 0.7,
        equivalent_plastic_strain: 0.,
        plastic_work_kpa: 0.,
        last_return: Return::Elastic,
    }
}
fn increment(p: &Parameters, s: &State, trial: [f64; 6]) -> [f64; 6] {
    let d = elastic::constitutive(p.young_kpa, p.poisson).unwrap();
    let x = d
        .lu()
        .solve(&SVector::from_fn(|i, _| {
            trial[i] - s.effective_stress_compression_positive_kpa[i]
        }))
        .unwrap();
    std::array::from_fn(|i| x[i])
}

#[test]
fn mc_face_edges_apex_and_zero_dilation_solvability_are_explicit() -> Result<()> {
    let p = parameters(10.);
    let s = initial();
    for (trial, mode) in [
        ([300., 100., 50., 0., 0., 0.], Return::ShearFace),
        ([300., 50., 50., 0., 0., 0.], Return::CompressionEdge),
        ([200., 200., 20., 0., 0., 0.], Return::ExtensionEdge),
        ([-100., -110., -120., 0., 0., 0.], Return::ShearApex),
    ] {
        let u = p.integrate_partition(&s, increment(&p, &s, trial), 1)?;
        assert_eq!(u.state.last_return, mode);
        assert!(p.normalized_yield(&u.state).abs() < 1e-10);
        assert!(u.state.plastic_work_kpa >= 0.);
        if mode == Return::ShearApex {
            let expected = -p.cohesion_kpa / p.friction_deg.to_radians().tan();
            for v in &u.state.effective_stress_compression_positive_kpa[..3] {
                assert!((v - expected).abs() < 1e-9);
            }
        }
    }
    let p = parameters(0.);
    let before = serde_json::to_vec(&s)?;
    assert!(
        p.integrate_partition(&s, increment(&p, &s, [-100., -110., -120., 0., 0., 0.]), 1)
            .is_err()
    );
    assert_eq!(before, serde_json::to_vec(&s)?);
    let u = p.integrate_partition(&s, increment(&p, &s, [300., 50., 50., 0., 0., 0.]), 1)?;
    assert!(
        u.state.plastic_strain_compression_positive[..3]
            .iter()
            .sum::<f64>()
            .abs()
            < 1e-12
    );
    Ok(())
}

#[test]
fn rankine_tension_surfaces_return_without_silent_apex_clipping() -> Result<()> {
    let mut p = parameters(0.);
    p.tension_cutoff_kpa = Some(0.);
    let s = initial();
    let u = p.integrate_partition(&s, increment(&p, &s, [-100., -110., -120., 0., 0., 0.]), 1)?;
    assert_eq!(u.state.last_return, Return::TensionApex);
    assert!(
        u.state
            .effective_stress_compression_positive_kpa
            .iter()
            .all(|v| v.abs() < 1e-9)
    );
    p.cohesion_kpa = 50.;
    let u = p.integrate_partition(&s, increment(&p, &s, [100., 20., -10., 0., 0., 0.]), 1)?;
    assert_eq!(u.state.last_return, Return::TensionFace);
    p.tension_cutoff_kpa = Some(1000.);
    assert!(p.validate().is_err());
    Ok(())
}

fn rotate_stress(q: &Matrix3<f64>, s: &[f64; 6]) -> [f64; 6] {
    mc::voigt(&(q * mc::tensor(s) * q.transpose()))
}
fn rotate_strain(q: &Matrix3<f64>, s: &[f64; 6]) -> [f64; 6] {
    let mut t = *s;
    for v in &mut t[3..] {
        *v /= 2.;
    }
    let mut t = rotate_stress(q, &t);
    for v in &mut t[3..] {
        *v *= 2.;
    }
    t
}
#[test]
fn mc_rotation_repeated_eigenvalues_and_plastic_history_are_objective() -> Result<()> {
    let p = parameters(7.);
    let s = initial();
    let q = nalgebra::Rotation3::new(Vector3::new(0.37, -0.52, 0.29)).into_inner();
    for trial in [
        [300., 100., 50., 0., 0., 0.],
        [300., 50., 50., 0., 0., 0.],
        [200., 200., 20., 0., 0., 0.],
        [-100., -110., -120., 0., 0., 0.],
        [101., 101., 101., 0., 0., 0.],
    ] {
        let de = increment(&p, &s, trial);
        let a = p.integrate_partition(&s, de, 1)?.state;
        let b = p.integrate_partition(&s, rotate_strain(&q, &de), 1)?.state;
        let expected = rotate_stress(&q, &a.effective_stress_compression_positive_kpa);
        assert!(
            expected
                .iter()
                .zip(b.effective_stress_compression_positive_kpa)
                .all(|(x, y)| (x - y).abs() < 1e-8)
        );
        let ep = rotate_strain(&q, &a.plastic_strain_compression_positive);
        assert!(
            ep.iter()
                .zip(b.plastic_strain_compression_positive)
                .all(|(x, y)| (x - y).abs() < 1e-10)
        );
        assert!((a.plastic_work_kpa - b.plastic_work_kpa).abs() < 1e-9);
    }
    Ok(())
}
#[test]
fn full_mc_tangent_is_nonsymmetric_and_matches_independent_direction() -> Result<()> {
    let p = parameters(5.);
    let s = initial();
    let de = increment(&p, &s, [300., 100., 50., 0., 0., 0.]);
    let (_, d) = p.stress_tangent_partition(&s, de, 1, 1e-7)?;
    assert!((d - d.transpose()).norm() > 100.);
    let direction = SVector::<f64, 6>::new(0.7, -0.2, 0.3, 0.4, -0.1, 0.25);
    let h = 2e-8;
    let a = p
        .integrate_partition(&s, std::array::from_fn(|i| de[i] + h * direction[i]), 1)?
        .state;
    let b = p
        .integrate_partition(&s, std::array::from_fn(|i| de[i] - h * direction[i]), 1)?
        .state;
    let diff = SVector::<f64, 6>::from_fn(|i, _| {
        (a.effective_stress_compression_positive_kpa[i]
            - b.effective_stress_compression_positive_kpa[i])
            / (2. * h)
    });
    assert!((d * direction - diff).norm() / diff.norm() < 1e-7);
    let mut associated = p.clone();
    associated.flow = Flow::Associated;
    associated.dilation_deg = associated.friction_deg;
    let (_, a) = associated.stress_tangent_partition(&s, de, 1, 1e-7)?;
    assert!((a - a.transpose()).norm() / a.norm() < 1e-7);
    Ok(())
}
#[test]
fn davis_is_an_explicit_weaker_associated_comparison_and_tresca_is_supported() -> Result<()> {
    let s = initial();
    let mut p = parameters(0.);
    p.flow = Flow::Davis;
    let (c, phi, psi) = p.effective_strength();
    assert!((c - 10. * 30f64.to_radians().cos()).abs() < 1e-12);
    assert_eq!(phi, psi);
    let de = increment(&p, &s, [300., 100., 50., 0., 0., 0.]);
    let a = p.integrate_partition(&s, de, 1)?.state;
    let mut q = p.clone();
    q.flow = Flow::Associated;
    q.cohesion_kpa = c;
    q.friction_deg = phi.to_degrees();
    q.dilation_deg = psi.to_degrees();
    let b = q.integrate_partition(&s, de, 1)?.state;
    for i in 0..6 {
        assert!(
            (a.effective_stress_compression_positive_kpa[i]
                - b.effective_stress_compression_positive_kpa[i])
                .abs()
                < 1e-9
        );
    }
    let mut p = parameters(0.);
    p.friction_deg = 0.;
    let a = p
        .integrate_partition(&s, increment(&p, &s, [200., 100., 80., 0., 0., 0.]), 1)?
        .state;
    let (v, _) = mc::principal(&a.effective_stress_compression_positive_kpa);
    assert!((v[0] - v[2] - 2. * p.cohesion_kpa).abs() < 1e-9);
    Ok(())
}
