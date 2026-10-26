use nalgebra::{Matrix3, SVector};
use soil_fem_3d::{Result, mcc, mcc_split};

fn parameters() -> mcc::Parameters {
    mcc::Parameters {
        critical_state_m: 1.2,
        lambda: 0.2,
        kappa: 0.04,
        poisson: 0.3,
        reference_pressure_kpa: 1.,
        specific_volume_ncl_at_reference: 3.5,
        integration: soil_fem_3d::mcc::Integration::BackwardEuler,
    }
}

fn state(pc: f64) -> Result<mcc::State> {
    Ok(mcc::State {
        effective_stress_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
        void_ratio: parameters().reference_void_ratio(100., pc)?,
        preconsolidation_pressure_kpa: pc,
        plastic_strain_compression_positive: [0.; 6],
    })
}

#[test]
fn elastic_split_has_exact_bulk_and_deviatoric_moduli() -> Result<()> {
    let material = parameters();
    let old = state(400.)?;
    let r = mcc_split::evaluate(
        &material,
        &old,
        [2e-5, -1e-5, -1e-5, 3e-5, 0., 0.],
        1e-5,
        1,
        1e-7,
    )?;
    let (bulk, shear) = material.moduli(&old)?;
    assert_eq!(r.update.plastic_substeps, 0);
    assert!((r.mean_pressure_compression_positive_kpa - (100. + bulk * 1e-5)).abs() < 1e-12);
    assert!((r.tangent.mean_by_volume - bulk).abs() / bulk < 1e-9);
    assert!(
        r.tangent
            .deviatoric_by_volume
            .iter()
            .all(|x| x.abs() < 1e-6)
    );
    assert!(r.tangent.mean_by_deviatoric.iter().all(|x| x.abs() < 1e-6));
    assert!((r.tangent.deviatoric_by_deviatoric[(0, 0)] - 4. * shear / 3.).abs() / shear < 1e-9);
    assert!((r.tangent.deviatoric_by_deviatoric[(3, 3)] - shear).abs() / shear < 1e-9);
    Ok(())
}

#[test]
fn plastic_split_reconstructs_all_history_and_unsymmetrized_full_tangent() -> Result<()> {
    let material = parameters();
    let normal = SVector::<f64, 6>::from_row_slice(&[1., 1., 1., 0., 0., 0.]);
    for pc in [100., 250.] {
        let old = state(pc)?;
        for (dev, volume) in [
            ([0.006, -0.003, -0.003, 0.03, 0.004, -0.002], 0.006),
            ([0.003, -0.001, -0.002, 0.04, 0.002, 0.004], -0.002),
        ] {
            let r = mcc_split::evaluate(&material, &old, dev, volume, 160, 1e-7)?;
            let mut full_strain = dev;
            for x in &mut full_strain[..3] {
                *x += volume / 3.;
            }
            let (direct, tangent) =
                material.stress_tangent_partition(&old, full_strain, 160, 1e-7)?;
            assert_eq!(
                serde_json::to_vec(&r.update.state)?,
                serde_json::to_vec(&direct.state)?
            );
            assert!(r.update.plastic_substeps > 0);
            let dev_stress =
                SVector::<f64, 6>::from_row_slice(&r.deviatoric_stress_compression_positive_kpa);
            let reconstructed = dev_stress + normal * r.mean_pressure_compression_positive_kpa;
            let expected = SVector::<f64, 6>::from_row_slice(
                &direct.state.effective_stress_compression_positive_kpa,
            );
            assert!((reconstructed - expected).norm() < 1e-12);
            let t = &r.tangent;
            let reconstructed_tangent = t.deviatoric_by_deviatoric
                + SVector::<f64, 6>::from_row_slice(&t.deviatoric_by_volume) * normal.transpose()
                + normal * SVector::<f64, 6>::from_row_slice(&t.mean_by_deviatoric).transpose()
                + normal * normal.transpose() * t.mean_by_volume;
            assert!((reconstructed_tangent - tangent).norm() / tangent.norm() < 1e-12);
        }
    }
    Ok(())
}

#[test]
fn split_blocks_predict_independent_directional_derivatives() -> Result<()> {
    let material = parameters();
    let old = state(100.)?;
    let dev = [0.006, -0.003, -0.003, 0.03, 0.004, -0.002];
    let volume = 0.006;
    let r = mcc_split::evaluate(&material, &old, dev, volume, 160, 1e-7)?;
    for (direction, dv) in [([0.5, -0.2, -0.3, 0.4, -0.2, 0.1], 0.7), ([0.; 6], 1.)] {
        let h = 5e-8;
        let plus = std::array::from_fn(|i| dev[i] + h * direction[i]);
        let minus = std::array::from_fn(|i| dev[i] - h * direction[i]);
        let a = mcc_split::evaluate(&material, &old, plus, volume + h * dv, 160, 1e-7)?;
        let b = mcc_split::evaluate(&material, &old, minus, volume - h * dv, 160, 1e-7)?;
        let actual_dev =
            (SVector::<f64, 6>::from_row_slice(&a.deviatoric_stress_compression_positive_kpa)
                - SVector::<f64, 6>::from_row_slice(&b.deviatoric_stress_compression_positive_kpa))
                / (2. * h);
        let expected_dev = r.tangent.deviatoric_by_deviatoric
            * SVector::<f64, 6>::from_row_slice(&direction)
            + SVector::<f64, 6>::from_row_slice(&r.tangent.deviatoric_by_volume) * dv;
        assert!((actual_dev - expected_dev).norm() / actual_dev.norm().max(1.) < 1e-5);
        let actual_p = (a.mean_pressure_compression_positive_kpa
            - b.mean_pressure_compression_positive_kpa)
            / (2. * h);
        let expected_p = SVector::<f64, 6>::from_row_slice(&r.tangent.mean_by_deviatoric)
            .dot(&SVector::<f64, 6>::from_row_slice(&direction))
            + r.tangent.mean_by_volume * dv;
        assert!((actual_p - expected_p).abs() / actual_p.abs().max(1.) < 1e-5);
    }
    Ok(())
}

fn tensor(s: [f64; 6], strain: bool) -> Matrix3<f64> {
    let a = if strain { 0.5 } else { 1. };
    Matrix3::new(
        s[0],
        a * s[3],
        a * s[5],
        a * s[3],
        s[1],
        a * s[4],
        a * s[5],
        a * s[4],
        s[2],
    )
}

fn voigt(s: Matrix3<f64>, strain: bool) -> [f64; 6] {
    let a = if strain { 2. } else { 1. };
    [
        s[(0, 0)],
        s[(1, 1)],
        s[(2, 2)],
        a * s[(0, 1)],
        a * s[(1, 2)],
        a * s[(2, 0)],
    ]
}

#[test]
fn rotated_plastic_response_and_failed_trials_preserve_committed_state() -> Result<()> {
    let material = parameters();
    let old = state(100.)?;
    let snapshot = serde_json::to_vec(&old)?;
    let dev = [0.006, -0.003, -0.003, 0.03, 0.004, -0.002];
    let q = nalgebra::Rotation3::from_euler_angles(0.23, -0.31, 0.47).into_inner();
    let rotated = voigt(q * tensor(dev, true) * q.transpose(), true);
    let a = mcc_split::evaluate(&material, &old, dev, 0.006, 160, 1e-7)?;
    let b = mcc_split::evaluate(&material, &old, rotated, 0.006, 160, 1e-7)?;
    let expected = q * tensor(a.deviatoric_stress_compression_positive_kpa, false) * q.transpose();
    assert!((expected - tensor(b.deviatoric_stress_compression_positive_kpa, false)).norm() < 1e-8);
    assert!(
        (a.mean_pressure_compression_positive_kpa - b.mean_pressure_compression_positive_kpa).abs()
            < 1e-8
    );
    assert!(
        (a.update.state.preconsolidation_pressure_kpa
            - b.update.state.preconsolidation_pressure_kpa)
            .abs()
            < 1e-8
    );
    assert!(mcc_split::evaluate(&material, &old, [0.; 6], -100., 1, 1e-7).is_err());
    assert!(mcc_split::evaluate(&material, &old, [0.01; 6], 0., 1, 1e-7).is_err());
    assert!(mcc_split::evaluate(&material, &old, [0.; 6], f64::NAN, 1, 1e-7).is_err());
    assert_eq!(snapshot, serde_json::to_vec(&old)?);
    Ok(())
}
