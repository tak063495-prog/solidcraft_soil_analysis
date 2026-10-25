use nalgebra::{Matrix3, SVector, Vector3};
use soil_fem_3d::{
    Result,
    hardening::{Parameters, Return},
    mc, mcc,
};

fn parameters() -> Parameters {
    serde_json::from_str(include_str!("../examples/hardening-parameters.json")).unwrap()
}
fn hydro(p: f64) -> [f64; 6] {
    [p, p, p, 0., 0., 0.]
}
fn rotate(q: &Matrix3<f64>, a: [f64; 6], strain: bool) -> [f64; 6] {
    let mut a = a;
    if strain {
        for x in &mut a[3..] {
            *x *= 0.5;
        }
    }
    let mut b = mc::voigt(&(q * mc::tensor(&a) * q.transpose()));
    if strain {
        for x in &mut b[3..] {
            *x *= 2.;
        }
    }
    b
}

#[test]
fn pressure_dependent_unloading_and_isotropic_hardening_follow_independent_integrals() -> Result<()>
{
    let p = parameters();
    let seed = p.initial_state(hydro(100.), 0.7, 100.)?;
    let k = p.moduli(&seed)?.0;
    let low = p.initial_state(hydro(25.), 0.7, 25.)?;
    assert!((p.moduli(&low)?.0 / k - 0.5).abs() < 1e-14);
    let overconsolidated = p.initial_state(hydro(100.), 0.7, 400.)?;
    assert!((p.moduli(&overconsolidated)?.0 / k - 2.).abs() < 1e-14);
    let volume = 0.02;
    let effective = k * p.cap_hardening_modulus_ref_kpa / (k + p.cap_hardening_modulus_ref_kpa);
    let expected = 100. * (1. + effective * volume / 200.).powi(2);
    let coarse = p
        .integrate_partition(
            &seed,
            [volume / 3., volume / 3., volume / 3., 0., 0., 0.],
            500,
        )?
        .state;
    let fine = p
        .integrate_partition(
            &seed,
            [volume / 3., volume / 3., volume / 3., 0., 0., 0.],
            1000,
        )?
        .state;
    let pc = mcc::invariants(&coarse.effective_stress_compression_positive_kpa).0;
    let pf = mcc::invariants(&fine.effective_stress_compression_positive_kpa).0;
    eprintln!("hydrostatic IDH p expected={expected}, coarse={pc}, fine={pf}");
    assert!((pf - expected).abs() < (pc - expected).abs() * 0.55);
    assert!((pf - expected).abs() / expected < 2e-4);
    assert!((fine.cap_pressure_kpa - pf).abs() < 1e-7);
    assert!(fine.cap_plastic_volumetric_strain > 0.01 && fine.shear_plastic_strain == 0.);
    let unloaded = p
        .integrate_partition(&fine, [-0.001, -0.001, -0.001, 0., 0., 0.], 1000)?
        .state;
    assert_eq!(unloaded.cap_pressure_kpa, fine.cap_pressure_kpa);
    assert_eq!(
        unloaded.cap_plastic_volumetric_strain,
        fine.cap_plastic_volumetric_strain
    );
    let reloaded = p
        .integrate_partition(&unloaded, [0.0005, 0.0005, 0.0005, 0., 0., 0.], 500)?
        .state;
    assert_eq!(reloaded.cap_pressure_kpa, fine.cap_pressure_kpa);
    assert!((reloaded.void_ratio + 1. - 1.7 * (-(volume - 0.003 + 0.0015)).exp()).abs() < 1e-11);
    Ok(())
}

#[test]
fn shear_hardening_constant_pressure_has_a_closed_form_and_elastic_reversal() -> Result<()> {
    let mut p = parameters();
    p.dilation_ratio = 0.;
    let seed = p.initial_state(hydro(100.), 0.7, 10000.)?;
    let history = 0.002;
    let target = 100. * p.shear_ratio(history);
    let g = p.moduli(&seed)?.1;
    let eq = history + target / (3. * g);
    let update = p
        .integrate_partition(&seed, [eq, -eq / 2., -eq / 2., 0., 0., 0.], 1)?
        .state;
    assert_eq!(update.last_return, Return::Shear);
    assert!(
        (mcc::invariants(&update.effective_stress_compression_positive_kpa).1 - target).abs()
            < 1e-7
    );
    assert!((update.shear_plastic_strain - history).abs() < 1e-10);
    assert_eq!(update.cap_plastic_volumetric_strain, 0.);
    let relax = target / (12. * g);
    let reverse = p
        .integrate_partition(&update, [-relax, relax / 2., relax / 2., 0., 0., 0.], 1)?
        .state;
    assert_eq!(reverse.last_return, Return::Elastic);
    assert_eq!(reverse.shear_plastic_strain, update.shear_plastic_strain);
    let back = p
        .integrate_partition(&reverse, [relax, -relax / 2., -relax / 2., 0., 0., 0.], 1)?
        .state;
    for (a, b) in back
        .effective_stress_compression_positive_kpa
        .iter()
        .zip(update.effective_stress_compression_positive_kpa)
    {
        assert!((a - b).abs() < 1e-8);
    }
    Ok(())
}

#[test]
fn shear_cap_corner_rotation_and_complete_increment_tangent_are_consistent() -> Result<()> {
    let p = parameters();
    let seed = p.initial_state(hydro(100.), 0.7, 120.)?;
    let rotation = nalgebra::Rotation3::new(Vector3::new(0.31, -0.42, 0.27)).into_inner();
    let mut seen = std::collections::BTreeSet::new();
    for strain in [
        [0.001, 0.001, 0.001, 0., 0., 0.],
        [0.004, -0.002, -0.002, 0., 0., 0.],
        [0.002, -0.001, -0.001, 0., 0., 0.],
        [0.01, 0., 0., 0.002, 0.001, -0.0007],
        [0.002, 0.002, 0.004, 0.004, 0., 0.],
    ] {
        let out = p.integrate_partition(&seed, strain, 1)?.state;
        seen.insert(format!("{:?}", out.last_return));
        let rotated = p
            .integrate_partition(&seed, rotate(&rotation, strain, true), 1)?
            .state;
        let expected = rotate(
            &rotation,
            out.effective_stress_compression_positive_kpa,
            false,
        );
        assert!(
            expected
                .iter()
                .zip(rotated.effective_stress_compression_positive_kpa)
                .all(|(a, b)| (a - b).abs() < 1e-7)
        );
        let ep = rotate(&rotation, out.plastic_strain_compression_positive, true);
        assert!(
            ep.iter()
                .zip(rotated.plastic_strain_compression_positive)
                .all(|(a, b)| (a - b).abs() < 1e-9)
        );
        assert!(out.plastic_work_kpa >= 0. && p.normalized_yield(&out) < 2e-8);
    }
    eprintln!("IDH active sets observed: {seen:?}");
    assert!(seen.contains("ShearCap") && seen.contains("Cap") && seen.contains("Shear"));
    let de = [0.004, 0.001, -0.0003, 0.002, 0.001, 0.0004];
    let (_, a) = p.stress_tangent_partition(&seed, de, 8, 1e-7)?;
    let (_, b) = p.stress_tangent_partition(&seed, de, 8, 5e-8)?;
    assert!((a - b).norm() / b.norm() < 1e-4);
    let direction = SVector::<f64, 6>::from_row_slice(&[0.2, -0.1, 0.15, 0.3, -0.12, 0.1]);
    let h = 2e-7;
    let plus = std::array::from_fn(|i| de[i] + h * direction[i]);
    let minus = std::array::from_fn(|i| de[i] - h * direction[i]);
    let up = p.integrate_partition(&seed, plus, 8)?.state;
    let down = p.integrate_partition(&seed, minus, 8)?.state;
    let derivative = SVector::<f64, 6>::from_fn(|i, _| {
        (up.effective_stress_compression_positive_kpa[i]
            - down.effective_stress_compression_positive_kpa[i])
            / (2. * h)
    });
    assert!((a * direction - derivative).norm() / derivative.norm() < 1e-4);
    assert!((a - a.transpose()).norm() / a.norm() > 1e-3);
    Ok(())
}

#[test]
fn invalid_pressure_history_and_failed_trials_do_not_modify_committed_idh_state() -> Result<()> {
    let p = parameters();
    let seed = p.initial_state(hydro(100.), 0.7, 110.)?;
    let bytes = serde_json::to_vec(&seed)?;
    assert!(
        p.initial_state([300., 0., 0., 0., 0., 0.], 0.7, 1e12)
            .is_err(),
        "a distant cap must not hide shear inadmissibility"
    );
    assert!(
        p.integrate_partition(&seed, [-0.1, -0.1, -0.1, 0., 0., 0.], 1)
            .is_err()
    );
    assert_eq!(bytes, serde_json::to_vec(&seed)?);
    let mut bad = seed.clone();
    bad.cap_pressure_kpa += 1.;
    assert!(p.validate_state(&bad).is_err());
    let mut bad = seed.clone();
    bad.plastic_strain_compression_positive[0] = 0.01;
    assert!(p.validate_state(&bad).is_err());
    let mut bad = p.clone();
    bad.poisson = -1.;
    assert!(bad.validate().is_err());
    let mut bad = p.clone();
    bad.shear_ratio_ultimate = 3.;
    assert!(bad.validate().is_err());
    Ok(())
}

#[test]
fn closed_nonproportional_elastic_cycles_leave_no_residual_stress_or_work() -> Result<()> {
    for exponent in [0., 0.5, 1.] {
        let mut p = parameters();
        p.stiffness_exponent = exponent;
        let seed = p.initial_state(hydro(100.), 0.7, 300.)?;
        let mut state = seed.clone();
        let mut work = 0.;
        for de in [
            [0.001 / 3., 0.001 / 3., 0.001 / 3., 0., 0., 0.],
            [0., 0., 0., 0.0002, 0., 0.],
            [-0.001 / 3., -0.001 / 3., -0.001 / 3., 0., 0., 0.],
            [0., 0., 0., -0.0002, 0., 0.],
        ] {
            let update = p.integrate_partition(&state, de, 100)?;
            assert_eq!(update.plastic_substeps, 0);
            work += de
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    0.5 * (state.effective_stress_compression_positive_kpa[i]
                        + update.state.effective_stress_compression_positive_kpa[i])
                        * e
                })
                .sum::<f64>();
            state = update.state;
        }
        assert!(
            state
                .effective_stress_compression_positive_kpa
                .iter()
                .zip(seed.effective_stress_compression_positive_kpa)
                .all(|(a, b)| (a - b).abs() < 1e-9)
        );
        assert!((state.void_ratio - seed.void_ratio).abs() < 1e-10);
        assert_eq!(state.cap_pressure_kpa, seed.cap_pressure_kpa);
        assert!(work.abs() < 1e-10);
    }
    Ok(())
}

#[test]
fn a_dilative_return_can_recover_a_negative_elastic_predictor_without_pressure_clipping()
-> Result<()> {
    let p = parameters();
    let seed = p.initial_state(hydro(100.), 0.7, 300.)?;
    let (k, g) = p.moduli(&seed)?;
    let volume = -101. / k;
    let deviator = 1000. / (3. * g);
    let de = [
        volume / 3. + deviator,
        volume / 3. - deviator / 2.,
        volume / 3. - deviator / 2.,
        0.,
        0.,
        0.,
    ];
    let result = p.integrate_partition(&seed, de, 1)?.state;
    assert_eq!(result.last_return, Return::Shear);
    let (mean, _) = mcc::invariants(&result.effective_stress_compression_positive_kpa);
    assert!(mean > p.minimum_mean_effective_pressure_kpa);
    assert!(p.normalized_yield(&result) < 1e-8);
    assert_eq!(seed.effective_stress_compression_positive_kpa, hydro(100.));
    Ok(())
}
