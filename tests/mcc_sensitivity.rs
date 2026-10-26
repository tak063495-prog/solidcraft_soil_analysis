use nalgebra::{Matrix3, SVector};
use soil_fem_3d::{Result, mcc};

fn model() -> mcc::Parameters {
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
fn isotropic(pc: f64) -> Result<mcc::State> {
    Ok(mcc::State {
        effective_stress_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
        void_ratio: model().reference_void_ratio(100., pc)?,
        preconsolidation_pressure_kpa: pc,
        plastic_strain_compression_positive: [0.; 6],
    })
}

#[test]
fn primal_histories_and_substeps_are_exact_on_nc_oc_loading_unloading_and_shear() -> Result<()> {
    let m = model();
    for pc in [100., 200., 400.] {
        let old = isotropic(pc)?;
        for de in [
            [0.003, 0.004, 0.005, 0.012, -0.002, 0.001],
            [0.003, 0.003, 0.003, 0., 0., 0.],
            [0., 0., 0., 0.02, 0.01, -0.005],
            [-0.0003, -0.0003, -0.0003, 0.0001, 0., 0.],
            [0.; 6],
        ] {
            for steps in [1, 37, 200] {
                let expected = m.integrate_partition(&old, de, steps)?;
                let (actual, _) = m.implicit_stress_tangent_partition(&old, de, steps)?;
                assert!(
                    serde_json::to_vec(&expected.state)? == serde_json::to_vec(&actual.state)?,
                    "primal changed pc={pc}, steps={steps}"
                );
                assert_eq!(expected.plastic_substeps, actual.plastic_substeps);
                assert_eq!(expected.substeps, actual.substeps);
            }
        }
    }
    Ok(())
}

#[test]
fn complete_tangent_matches_two_finite_difference_steps_and_hydrostatic_limit() -> Result<()> {
    let m = model();
    for (pc, de, steps) in [
        (400., [1e-5, -1e-5, 2e-5, 1e-5, 0., 0.], 37),
        (100., [0.003, 0.003, 0.003, 0., 0., 0.], 200),
        (100., [0.006, -0.003, -0.003, 0.03, 0.004, -0.002], 200),
        (250., [0.004, -0.001, -0.002, 0.08, 0.002, 0.004], 400),
    ] {
        let old = isotropic(pc)?;
        let (_, analytic) = m.implicit_stress_tangent_partition(&old, de, steps)?;
        for h in [2e-7, 1e-7] {
            let (_, difference) = m.stress_tangent_partition(&old, de, steps, h)?;
            let error = (analytic - difference).norm() / difference.norm();
            println!("pc={pc} steps={steps} h={h:e} tangent_relative_error={error:e}");
            assert!(error < 1e-4, "complete tangent error={error:e}");
        }
    }
    Ok(())
}

#[test]
fn critical_state_zero_volumetric_flow_and_elastic_reversal_have_finite_consistent_tangents()
-> Result<()> {
    let m = model();
    let mut old = isotropic(200.)?;
    let q = m.critical_state_m * 100.;
    old.effective_stress_compression_positive_kpa =
        [100. + 2. * q / 3., 100. - q / 3., 100. - q / 3., 0., 0., 0.];
    m.validate_state(&old)?;
    let de = [0.003, -0.0015, -0.0015, 0.002, 0., 0.];
    let (_, analytic) = m.implicit_stress_tangent_partition(&old, de, 100)?;
    let (_, difference) = m.stress_tangent_partition(&old, de, 100, 1e-7)?;
    assert!((analytic - difference).norm() / difference.norm() < 1e-4);
    let load = m
        .integrate_partition(&isotropic(100.)?, [0.005, 0.005, 0.005, 0., 0., 0.], 100)?
        .state;
    let de = [-0.001, -0.001, -0.001, 0., 0., 0.];
    let (unload, analytic) = m.implicit_stress_tangent_partition(&load, de, 100)?;
    let (_, difference) = m.stress_tangent_partition(&load, de, 100, 1e-7)?;
    assert_eq!(unload.plastic_substeps, 0);
    assert_eq!(
        unload.state.preconsolidation_pressure_kpa,
        load.preconsolidation_pressure_kpa
    );
    assert!((analytic - difference).norm() / difference.norm() < 1e-5);
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
fn nonisotropic_history_and_rotated_direction_preserve_objectivity() -> Result<()> {
    let m = model();
    let mut old = m
        .integrate_partition(
            &isotropic(200.)?,
            [0.004, -0.002, -0.002, 0.025, 0.003, -0.005],
            200,
        )?
        .state;
    let de = [0.003, -0.001, -0.002, 0.035, -0.002, 0.004];
    let direction = [0.1, -0.2, 0.3, 0.4, 0.1, -0.2];
    let (_, a) = m.implicit_stress_tangent_partition(&old, de, 200)?;
    let q = nalgebra::Rotation3::from_euler_angles(0.2, -0.4, 0.3).into_inner();
    old.effective_stress_compression_positive_kpa = voigt(
        q * tensor(old.effective_stress_compression_positive_kpa, false) * q.transpose(),
        false,
    );
    old.plastic_strain_compression_positive = voigt(
        q * tensor(old.plastic_strain_compression_positive, true) * q.transpose(),
        true,
    );
    let rd = voigt(q * tensor(de, true) * q.transpose(), true);
    let rotated_direction = voigt(q * tensor(direction, true) * q.transpose(), true);
    let (_, b) = m.implicit_stress_tangent_partition(&old, rd, 200)?;
    let sa = a * SVector::<f64, 6>::from_row_slice(&direction);
    let sb = b * SVector::<f64, 6>::from_row_slice(&rotated_direction);
    let expected = q * tensor(std::array::from_fn(|i| sa[i]), false) * q.transpose();
    let actual = tensor(std::array::from_fn(|i| sb[i]), false);
    assert!((actual - expected).norm() / expected.norm() < 1e-7);
    Ok(())
}

#[test]
fn invalid_trials_and_zero_partitions_never_mutate_committed_history() -> Result<()> {
    let m = model();
    let old = isotropic(200.)?;
    let bytes = serde_json::to_vec(&old)?;
    assert!(
        m.implicit_stress_tangent_partition(&old, [0.; 6], 0)
            .is_err()
    );
    assert!(
        m.implicit_stress_tangent_partition(&old, [f64::NAN; 6], 10)
            .is_err()
    );
    assert!(
        m.implicit_stress_tangent_partition(&old, [-100., -100., -100., 0., 0., 0.], 1)
            .is_err()
    );
    assert!(bytes == serde_json::to_vec(&old)?);
    Ok(())
}

#[test]
fn report_material_tangent_timing_without_changing_acceptance_or_primal_values() -> Result<()> {
    let m = model();
    let old = isotropic(100.)?;
    let de = [0.006, -0.003, -0.003, 0.03, 0.004, -0.002];
    let mut measurements = Vec::new();
    for _ in 0..5 {
        let start = std::time::Instant::now();
        let expected = m.stress_tangent_partition(&old, de, 1000, 1e-7)?;
        let finite = start.elapsed().as_secs_f64();
        let start = std::time::Instant::now();
        let actual = m.implicit_stress_tangent_partition(&old, de, 1000)?;
        let implicit = start.elapsed().as_secs_f64();
        assert!(serde_json::to_vec(&expected.0.state)? == serde_json::to_vec(&actual.0.state)?);
        assert!((expected.1 - actual.1).norm() / expected.1.norm() < 1e-4);
        measurements.push((finite, implicit));
    }
    measurements.sort_by(|a, b| a.0.total_cmp(&b.0));
    let finite = measurements[2].0;
    measurements.sort_by(|a, b| a.1.total_cmp(&b.1));
    let implicit = measurements[2].1;
    println!(
        "MCC tangent 1000 substeps: finite_difference_median_seconds={finite:e}; implicit_median_seconds={implicit:e}; speedup={:e}",
        finite / implicit
    );
    Ok(())
}
