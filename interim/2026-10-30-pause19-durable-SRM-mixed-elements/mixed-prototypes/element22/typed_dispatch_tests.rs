use super::*;
use soil_fem_3d::{mcc, mcc_split, soil};
fn mcc_input() -> Result<(mcc::Parameters, mcc::State)> {
    let p = mcc::Parameters {
        critical_state_m: 1.2,
        lambda: 0.2,
        kappa: 0.04,
        poisson: 0.3,
        reference_pressure_kpa: 1.,
        specific_volume_ncl_at_reference: 3.5,
        integration: mcc::Integration::BackwardEuler,
    };
    let s = mcc::State {
        effective_stress_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
        void_ratio: p.reference_void_ratio(100., 250.)?,
        preconsolidation_pressure_kpa: 250.,
        plastic_strain_compression_positive: [0.; 6],
    };
    Ok((p, s))
}
fn mcc_bits(s: &mcc::State) -> Vec<u64> {
    s.effective_stress_compression_positive_kpa
        .iter()
        .chain(&s.plastic_strain_compression_positive)
        .chain([&s.void_ratio, &s.preconsolidation_pressure_kpa])
        .map(|v| v.to_bits())
        .collect()
}
#[test]
fn typed_mcc_dispatch_keeps_all_state_and_original_implicit_tangent_bits() -> Result<()> {
    let (p, s) = mcc_input()?;
    let before = mcc_bits(&s);
    let model = soil::Model::Mcc(p.clone());
    let old = soil::State::Mcc(s.clone());
    let dev = [0.003, -0.001, -0.002, 0.01, 0.002, 0.004];
    let v = -0.002;
    let direct = mcc_split::evaluate_implicit(&p, &s, dev, v, 16)?;
    let r = soil_split::evaluate(&model, &old, dev, v, 16, soil_split::Method::ImplicitMcc)?;
    let soil::State::Mcc(state) = r.values.update.state else {
        panic!("MCC type changed")
    };
    assert_eq!(mcc_bits(&state), mcc_bits(&direct.update.state));
    assert_eq!(
        (r.values.update.substeps, r.values.update.plastic_substeps),
        (direct.update.substeps, direct.update.plastic_substeps)
    );
    assert_eq!(
        r.values.mean_pressure_compression_positive_kpa.to_bits(),
        direct.mean_pressure_compression_positive_kpa.to_bits()
    );
    assert_eq!(
        r.values
            .deviatoric_stress_compression_positive_kpa
            .map(f64::to_bits),
        direct
            .deviatoric_stress_compression_positive_kpa
            .map(f64::to_bits)
    );
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(
                r.tangent.deviatoric_by_deviatoric[i][j].to_bits(),
                direct.tangent.deviatoric_by_deviatoric[(i, j)].to_bits()
            );
        }
    }
    assert_eq!(
        r.tangent.deviatoric_by_volume.map(f64::to_bits),
        direct.tangent.deviatoric_by_volume.map(f64::to_bits)
    );
    assert_eq!(
        r.tangent.mean_by_deviatoric.map(f64::to_bits),
        direct.tangent.mean_by_deviatoric.map(f64::to_bits)
    );
    assert_eq!(
        r.tangent.mean_by_volume.to_bits(),
        direct.tangent.mean_by_volume.to_bits()
    );
    let integrated = soil_split::integrate(&model, &old, dev, v, 16)?;
    let unchanged = mcc_split::integrate(&p, &s, dev, v, 16)?;
    let soil::State::Mcc(state) = integrated.update.state else {
        panic!("MCC type changed")
    };
    assert_eq!(mcc_bits(&state), mcc_bits(&unchanged.update.state));
    assert_eq!(mcc_bits(&s), before);
    Ok(())
}
#[test]
fn typed_mc_dispatch_preserves_mc_state_and_all_split_block_bits() -> Result<()> {
    let p = parameters(5.);
    let s = initial();
    let before = bits(&s);
    let (dev, v) = split(increment(&p, &s, [300., 100., 50., 0., 0., 0.]));
    let direct = mc_split::evaluate(&p, &s, dev, v, 1, 1e-7)?;
    let r = soil_split::evaluate(
        &soil::Model::MohrCoulomb(p),
        &soil::State::MohrCoulomb(s.clone()),
        dev,
        v,
        1,
        soil_split::Method::FiniteDifference(1e-7),
    )?;
    let soil::State::MohrCoulomb(state) = r.values.update.state else {
        panic!("MC converted to another state")
    };
    assert_eq!(bits(&state), bits(&direct.values.update.state));
    assert_eq!(state.last_return, direct.values.update.state.last_return);
    assert_eq!(
        (r.values.update.substeps, r.values.update.plastic_substeps),
        (
            direct.values.update.substeps,
            direct.values.update.plastic_substeps
        )
    );
    assert_eq!(
        r.values.mean_pressure_compression_positive_kpa.to_bits(),
        direct
            .values
            .mean_pressure_compression_positive_kpa
            .to_bits()
    );
    assert_eq!(
        r.values
            .deviatoric_stress_compression_positive_kpa
            .map(f64::to_bits),
        direct
            .values
            .deviatoric_stress_compression_positive_kpa
            .map(f64::to_bits)
    );
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(
                r.tangent.deviatoric_by_deviatoric[i][j].to_bits(),
                direct.tangent.deviatoric_by_deviatoric[i][j].to_bits()
            );
        }
    }
    assert_eq!(
        r.tangent.deviatoric_by_volume.map(f64::to_bits),
        direct.tangent.deviatoric_by_volume.map(f64::to_bits)
    );
    assert_eq!(
        r.tangent.mean_by_deviatoric.map(f64::to_bits),
        direct.tangent.mean_by_deviatoric.map(f64::to_bits)
    );
    assert_eq!(
        r.tangent.mean_by_volume.to_bits(),
        direct.tangent.mean_by_volume.to_bits()
    );
    assert_eq!(bits(&s), before);
    Ok(())
}
#[test]
fn typed_mismatches_and_unavailable_mc_implicit_method_are_rejected_without_conversions()
-> Result<()> {
    let (p, s) = mcc_input()?;
    let mc = parameters(0.);
    let ms = initial();
    let before = bits(&ms);
    let mcc_before = mcc_bits(&s);
    assert!(
        soil_split::integrate(
            &soil::Model::Mcc(p.clone()),
            &soil::State::MohrCoulomb(ms.clone()),
            [0.; 6],
            0.,
            1
        )
        .is_err()
    );
    assert!(
        soil_split::integrate(
            &soil::Model::MohrCoulomb(mc.clone()),
            &soil::State::Mcc(s.clone()),
            [0.; 6],
            0.,
            1
        )
        .is_err()
    );
    assert!(
        soil_split::evaluate(
            &soil::Model::MohrCoulomb(mc.clone()),
            &soil::State::MohrCoulomb(ms.clone()),
            [0.; 6],
            0.,
            1,
            soil_split::Method::ImplicitMcc
        )
        .is_err()
    );
    assert!(
        soil_split::integrate_transition(
            &soil::Model::MohrCoulomb(mc),
            &soil::Model::Mcc(p),
            &soil::State::Mcc(s.clone()),
            [0.; 6],
            0.,
            1
        )
        .is_err()
    );
    assert_eq!(bits(&ms), before);
    assert_eq!(mcc_bits(&s), mcc_before);
    Ok(())
}
#[test]
fn typed_mc_transition_retains_complete_committed_history_and_matches_explicit_split() -> Result<()>
{
    let from = transition_parameters(20., 0., 0.);
    let to = transition_parameters(10., 0., 0.);
    let old = transition_state([120., 100., 80., 0., 0., 0.]);
    let before = bits(&old);
    let direct = mc_split::integrate_transition(&to, &from, &old, [0.; 6], 0., 1)?;
    let r = soil_split::integrate_transition(
        &soil::Model::MohrCoulomb(to),
        &soil::Model::MohrCoulomb(from),
        &soil::State::MohrCoulomb(old.clone()),
        [0.; 6],
        0.,
        1,
    )?;
    let soil::State::MohrCoulomb(state) = r.update.state else {
        panic!("typed transition relabelled history")
    };
    assert_eq!(bits(&state), bits(&direct.update.state));
    assert_eq!(state.last_return, direct.update.state.last_return);
    assert_eq!(
        (r.update.substeps, r.update.plastic_substeps),
        (direct.update.substeps, direct.update.plastic_substeps)
    );
    assert_eq!(bits(&old), before);
    Ok(())
}

#[test]
fn typed_transition_jacobian_keeps_every_nonsymmetric_block_and_rejects_foreign_method()
-> Result<()> {
    let from = transition_parameters(20., 30., 10.);
    let to = transition_parameters(10., 20., 10.);
    let old = transition_state([180., 100., 70., 0., 0., 0.]);
    let before = bits(&old);
    let direct = mc_split::evaluate_transition(&to, &from, &old, [0.; 6], 0., 1, 1e-7)?;
    let a = soil::Model::MohrCoulomb(to);
    let b = soil::Model::MohrCoulomb(from);
    let s = soil::State::MohrCoulomb(old.clone());
    let r = soil_split::evaluate_transition(
        &a,
        &b,
        &s,
        [0.; 6],
        0.,
        1,
        soil_split::Method::FiniteDifference(1e-7),
    )?;
    let soil::State::MohrCoulomb(state) = r.values.update.state else {
        panic!("typed Jacobian relabelled history")
    };
    assert_eq!(bits(&state), bits(&direct.values.update.state));
    assert_eq!(state.last_return, direct.values.update.state.last_return);
    assert_eq!(
        (r.values.update.substeps, r.values.update.plastic_substeps),
        (
            direct.values.update.substeps,
            direct.values.update.plastic_substeps
        )
    );
    assert_eq!(
        r.values.mean_pressure_compression_positive_kpa.to_bits(),
        direct
            .values
            .mean_pressure_compression_positive_kpa
            .to_bits()
    );
    assert_eq!(
        r.values
            .deviatoric_stress_compression_positive_kpa
            .map(f64::to_bits),
        direct
            .values
            .deviatoric_stress_compression_positive_kpa
            .map(f64::to_bits)
    );
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(
                r.tangent.deviatoric_by_deviatoric[i][j].to_bits(),
                direct.tangent.deviatoric_by_deviatoric[i][j].to_bits()
            );
        }
    }
    assert_eq!(
        r.tangent.deviatoric_by_volume.map(f64::to_bits),
        direct.tangent.deviatoric_by_volume.map(f64::to_bits)
    );
    assert_eq!(
        r.tangent.mean_by_deviatoric.map(f64::to_bits),
        direct.tangent.mean_by_deviatoric.map(f64::to_bits)
    );
    assert_eq!(
        r.tangent.mean_by_volume.to_bits(),
        direct.tangent.mean_by_volume.to_bits()
    );
    assert!(
        soil_split::evaluate_transition(
            &a,
            &b,
            &s,
            [0.; 6],
            0.,
            1,
            soil_split::Method::ImplicitMcc
        )
        .is_err()
    );
    let (m, ms) = mcc_input()?;
    assert!(
        soil_split::evaluate_transition(
            &soil::Model::Mcc(m),
            &b,
            &soil::State::Mcc(ms),
            [0.; 6],
            0.,
            1,
            soil_split::Method::FiniteDifference(1e-7)
        )
        .is_err()
    );
    assert_eq!(bits(&old), before);
    Ok(())
}
