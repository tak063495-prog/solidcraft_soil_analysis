mod mc_split;
use soil_fem_3d::{
    Result,
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
fn increment(p: &Parameters, old: &State, trial: [f64; 6]) -> [f64; 6] {
    let ds: [f64; 6] =
        std::array::from_fn(|i| trial[i] - old.effective_stress_compression_positive_kpa[i]);
    std::array::from_fn(|i| {
        if i < 3 {
            (ds[i] - p.poisson * (ds[(i + 1) % 3] + ds[(i + 2) % 3])) / p.young_kpa
        } else {
            2. * (1. + p.poisson) * ds[i] / p.young_kpa
        }
    })
}
fn split(full: [f64; 6]) -> ([f64; 6], f64) {
    let volume = full[..3].iter().sum::<f64>();
    let mut dev: [f64; 6] = std::array::from_fn(|i| full[i] - if i < 3 { volume / 3. } else { 0. });
    // Independent coordinates: the third normal component is constrained.
    // This avoids cancellation creating an invalid nonzero trace near hydrostatic strain.
    dev[2] = -(dev[0] + dev[1]);
    (dev, volume)
}
fn full(dev: [f64; 6], v: f64) -> [f64; 6] {
    std::array::from_fn(|i| dev[i] + if i < 3 { v / 3. } else { 0. })
}
fn bits(s: &State) -> Vec<u64> {
    s.effective_stress_compression_positive_kpa
        .iter()
        .chain(&s.plastic_strain_compression_positive)
        .chain([
            &s.void_ratio,
            &s.equivalent_plastic_strain,
            &s.plastic_work_kpa,
        ])
        .map(|x| x.to_bits())
        .collect()
}
fn same(a: &mc::Update, b: &mc::Update) {
    assert_eq!(bits(&a.state), bits(&b.state));
    assert_eq!(a.state.last_return, b.state.last_return);
    assert_eq!(
        (a.substeps, a.plastic_substeps),
        (b.substeps, b.plastic_substeps)
    );
}
#[test]
fn elastic_bulk_shear_and_near_incompressibility_match_closed_form() -> Result<()> {
    for nu in [0.25, 0.49, 0.4999] {
        let mut p = parameters(0.);
        p.poisson = nu;
        p.cohesion_kpa = 1e6;
        let old = initial();
        let dev = [2e-5, -1e-5, -1e-5, 3e-5, 0., 0.];
        let volume = 1e-6;
        let r = mc_split::evaluate(&p, &old, dev, volume, 1, 1e-7)?;
        let bulk = p.young_kpa / (3. * (1. - 2. * p.poisson));
        let shear = p.young_kpa / (2. * (1. + p.poisson));
        assert_eq!(r.values.update.plastic_substeps, 0);
        assert!(
            (r.values.mean_pressure_compression_positive_kpa - (100. + bulk * volume)).abs() < 1e-7
        );
        assert!((r.tangent.mean_by_volume - bulk).abs() / bulk < 1e-7);
        assert!((r.tangent.deviatoric_by_deviatoric[3][3] - shear).abs() / shear < 1e-7);
        assert!((r.tangent.deviatoric_by_deviatoric[0][0] - 4. * shear / 3.).abs() / shear < 1e-7);
        for i in 0..6 {
            let expected = dev[i] * if i < 3 { 2. * shear } else { shear };
            assert!(
                (r.values.deviatoric_stress_compression_positive_kpa[i] - expected).abs() < 1e-7
            );
        }
    }
    Ok(())
}
#[test]
fn face_edges_and_apex_retain_typed_mc_return_and_all_history_bits() -> Result<()> {
    let p = parameters(10.);
    let old = initial();
    let before = bits(&old);
    for (target, kind) in [
        ([300., 100., 50., 0., 0., 0.], Return::ShearFace),
        ([300., 50., 50., 0., 0., 0.], Return::CompressionEdge),
        ([200., 200., 20., 0., 0., 0.], Return::ExtensionEdge),
        ([-100., -110., -120., 0., 0., 0.], Return::ShearApex),
    ] {
        let (dev, v) = split(increment(&p, &old, target));
        let r = mc_split::integrate(&p, &old, dev, v, 1)?;
        let direct = p.integrate_partition(&old, full(dev, v), 1)?;
        same(&r.update, &direct);
        assert_eq!(r.update.state.last_return, kind);
        assert!(p.normalized_yield(&r.update.state).abs() < 1e-10);
    }
    assert_eq!(bits(&old), before);
    Ok(())
}
#[test]
fn rankine_face_and_apex_are_not_pressure_clipped() -> Result<()> {
    let mut p = parameters(0.);
    p.tension_cutoff_kpa = Some(0.);
    let old = initial();
    for (cohesion, target, kind) in [
        (10., [-100., -110., -120., 0., 0., 0.], Return::TensionApex),
        (50., [100., 20., -10., 0., 0., 0.], Return::TensionFace),
    ] {
        p.cohesion_kpa = cohesion;
        let (dev, v) = split(increment(&p, &old, target));
        let r = mc_split::integrate(&p, &old, dev, v, 1)?;
        same(&r.update, &p.integrate_partition(&old, full(dev, v), 1)?);
        assert_eq!(r.update.state.last_return, kind);
    }
    Ok(())
}
#[test]
fn nonsymmetric_split_reconstructs_full_tangent_and_predicts_independent_directions() -> Result<()>
{
    let p = parameters(5.);
    let old = initial();
    let (dev, v) = split(increment(&p, &old, [300., 100., 50., 0., 0., 0.]));
    let r = mc_split::evaluate(&p, &old, dev, v, 1, 1e-7)?;
    let (_, d) = p.stress_tangent_partition(&old, full(dev, v), 1, 1e-7)?;
    let mut asym = 0.;
    for i in 0..6 {
        for j in 0..6 {
            let ni = if i < 3 { 1. } else { 0. };
            let nj = if j < 3 { 1. } else { 0. };
            let reconstructed = r.tangent.deviatoric_by_deviatoric[i][j]
                + r.tangent.deviatoric_by_volume[i] * nj
                + ni * r.tangent.mean_by_deviatoric[j]
                + ni * nj * r.tangent.mean_by_volume;
            assert!((reconstructed - d[(i, j)]).abs() < 1e-8);
            asym += (d[(i, j)] - d[(j, i)]).powi(2);
        }
    }
    assert!(asym.sqrt() > 100.);
    for (direction, dv) in [([0.5, -0.2, -0.3, 0.4, -0.2, 0.1], 0.7), ([0.; 6], 1.)] {
        let h = 2e-8;
        let a = mc_split::integrate(
            &p,
            &old,
            std::array::from_fn(|i| dev[i] + h * direction[i]),
            v + h * dv,
            1,
        )?;
        let b = mc_split::integrate(
            &p,
            &old,
            std::array::from_fn(|i| dev[i] - h * direction[i]),
            v - h * dv,
            1,
        )?;
        let measured_mean = (a.mean_pressure_compression_positive_kpa
            - b.mean_pressure_compression_positive_kpa)
            / (2. * h);
        let predicted_mean = (0..6)
            .map(|j| r.tangent.mean_by_deviatoric[j] * direction[j])
            .sum::<f64>()
            + r.tangent.mean_by_volume * dv;
        assert!((measured_mean - predicted_mean).abs() < 1e-3);
        for i in 0..6 {
            let measured = (a.deviatoric_stress_compression_positive_kpa[i]
                - b.deviatoric_stress_compression_positive_kpa[i])
                / (2. * h);
            let predicted = (0..6)
                .map(|j| r.tangent.deviatoric_by_deviatoric[i][j] * direction[j])
                .sum::<f64>()
                + r.tangent.deviatoric_by_volume[i] * dv;
            assert!((measured - predicted).abs() < 1e-3);
        }
    }
    Ok(())
}
#[test]
fn shear_work_and_partitioned_plastic_history_use_engineering_strain_convention() -> Result<()> {
    let p = parameters(7.);
    let mut old = initial();
    for de in [
        [0.01, -0.005, -0.005, 0.012, 0.001, -0.002],
        [-0.001, 0., 0., -0.0005, 0., 0.],
    ] {
        let (dev, v) = split(de);
        let r = mc_split::integrate(&p, &old, dev, v, 8)?;
        let direct = p.integrate_partition(&old, full(dev, v), 8)?;
        same(&r.update, &direct);
        let total = (0..6)
            .map(|i| r.update.state.effective_stress_compression_positive_kpa[i] * full(dev, v)[i])
            .sum::<f64>();
        let separated = (0..6)
            .map(|i| r.deviatoric_stress_compression_positive_kpa[i] * dev[i])
            .sum::<f64>()
            + r.mean_pressure_compression_positive_kpa * v;
        assert!((total - separated).abs() < 1e-12);
        old = r.update.state;
    }
    assert!(old.equivalent_plastic_strain > 0.);
    Ok(())
}
#[test]
fn invalid_split_and_genuine_zero_dilation_no_solution_leave_committed_state_unchanged() {
    let p = parameters(0.);
    let old = initial();
    let before = bits(&old);
    for (dev, v, steps) in [
        ([1e-4; 6], 0., 1),
        ([f64::NAN; 6], 0., 1),
        ([0.; 6], f64::INFINITY, 1),
        ([0.; 6], 0., 0),
    ] {
        assert!(mc_split::integrate(&p, &old, dev, v, steps).is_err());
    }
    let (dev, v) = split(increment(&p, &old, [-100., -110., -120., 0., 0., 0.]));
    assert!(mc_split::integrate(&p, &old, dev, v, 1).is_err());
    assert!(mc_split::evaluate(&p, &old, dev, v, 1, 1e-7).is_err());
    assert_eq!(bits(&old), before);
    assert_eq!(old.last_return, Return::Elastic);
}

#[test]
fn arbitrary_rotation_preserves_stress_and_plastic_history_with_engineering_shear() -> Result<()> {
    let axis: [f64; 3] = [0.37, -0.52, 0.29];
    let angle = axis.iter().map(|v| v * v).sum::<f64>().sqrt();
    let n = axis.map(|v| v / angle);
    let (s, c) = angle.sin_cos();
    let skew = [[0., -n[2], n[1]], [n[2], 0., -n[0]], [-n[1], n[0], 0.]];
    let q: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            c * if i == j { 1. } else { 0. } + (1. - c) * n[i] * n[j] + s * skew[i][j]
        })
    });
    let rotate = |v: [f64; 6], engineering: bool| {
        let factor = if engineering { 0.5 } else { 1. };
        let a = [
            [v[0], v[3] * factor, v[5] * factor],
            [v[3] * factor, v[1], v[4] * factor],
            [v[5] * factor, v[4] * factor, v[2]],
        ];
        let b: [[f64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let mut x = 0.;
                for k in 0..3 {
                    for l in 0..3 {
                        x += q[i][k] * a[k][l] * q[j][l];
                    }
                }
                x
            })
        });
        [
            b[0][0],
            b[1][1],
            b[2][2],
            b[0][1] / factor,
            b[1][2] / factor,
            b[2][0] / factor,
        ]
    };
    let p = parameters(7.);
    let mut a = initial();
    let mut b = initial();
    for de in [
        [0.01, -0.005, -0.005, 0.012, 0.001, -0.002],
        [-0.001, 0., 0., -0.0005, 0., 0.],
        [1e-5, 1e-5, 1e-5, 0., 0., 0.],
    ] {
        let (dev, v) = split(de);
        let (rd, rv) = split(rotate(de, true));
        a = mc_split::integrate(&p, &a, dev, v, 8)?.update.state;
        b = mc_split::integrate(&p, &b, rd, rv, 8)?.update.state;
        let expected = rotate(a.effective_stress_compression_positive_kpa, false);
        let plastic = rotate(a.plastic_strain_compression_positive, true);
        assert!(
            (0..6).all(
                |i| (expected[i] - b.effective_stress_compression_positive_kpa[i]).abs() < 1e-8
            )
        );
        assert!(
            (0..6).all(|i| (plastic[i] - b.plastic_strain_compression_positive[i]).abs() < 1e-10)
        );
        assert!((a.void_ratio - b.void_ratio).abs() < 1e-12);
        assert!((a.equivalent_plastic_strain - b.equivalent_plastic_strain).abs() < 1e-10);
        assert!((a.plastic_work_kpa - b.plastic_work_kpa).abs() < 1e-9);
    }
    Ok(())
}

#[test]
fn test_input_coordinates_are_exactly_trace_free_without_relaxing_material_guard() -> Result<()> {
    let p = parameters(0.);
    let old = initial();
    let raw = [
        1e-5 - 1.6940658945086007e-21,
        1e-5,
        1e-5,
        0.,
        8.470329472543003e-22,
        0.,
    ];
    let (dev, v) = split(raw);
    assert_eq!(dev[..3].iter().sum::<f64>(), 0.);
    let restored = full(dev, v);
    assert!((0..6).all(|i| (restored[i] - raw[i]).abs() <= 8. * f64::EPSILON * 1e-5));
    mc_split::integrate(&p, &old, dev, v, 1)?;
    assert!(mc_split::integrate(&p, &old, [1e-22, 0., 0., 0., 0., 0.], 0., 1).is_err());
    Ok(())
}

#[test]
fn remaining_rankine_edge_and_combined_shear_tension_returns_are_explicit() -> Result<()> {
    let old = initial();
    let mut p = parameters(10.);
    p.tension_cutoff_kpa = Some(0.);
    for (c, target, kind) in [
        (50., [100., -10., -20., 0., 0., 0.], Return::TensionEdge),
        (10., [300., 100., -100., 0., 0., 0.], Return::Mixed),
    ] {
        p.cohesion_kpa = c;
        let (dev, v) = split(increment(&p, &old, target));
        let r = mc_split::integrate(&p, &old, dev, v, 1)?;
        assert_eq!(r.update.state.last_return, kind);
        assert!(p.normalized_yield(&r.update.state).abs() < 1e-10);
        same(&r.update, &p.integrate_partition(&old, full(dev, v), 1)?);
    }
    Ok(())
}

#[test]
fn associated_and_explicit_davis_flows_keep_their_own_strength_and_history() -> Result<()> {
    let old = initial();
    let mut p = parameters(30.);
    p.flow = Flow::Associated;
    let (dev, v) = split(increment(&p, &old, [300., 100., 50., 0., 0., 0.]));
    let r = mc_split::evaluate(&p, &old, dev, v, 1, 1e-7)?;
    same(
        &r.values.update,
        &p.integrate_partition(&old, full(dev, v), 1)?,
    );
    let (_, d) = p.stress_tangent_partition(&old, full(dev, v), 1, 1e-7)?;
    let asym = (0..6)
        .flat_map(|i| (0..6).map(move |j| (i, j)))
        .map(|(i, j)| (d[(i, j)] - d[(j, i)]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(asym < 1e-3);
    p.flow = Flow::Davis;
    p.dilation_deg = 0.;
    let r = mc_split::integrate(&p, &old, dev, v, 1)?;
    let mut explicit = p.clone();
    let beta = p.friction_deg.to_radians().cos();
    let angle = (beta * p.friction_deg.to_radians().tan())
        .atan()
        .to_degrees();
    explicit.cohesion_kpa = beta * p.cohesion_kpa;
    explicit.friction_deg = angle;
    explicit.dilation_deg = angle;
    explicit.flow = Flow::Associated;
    let b = mc_split::integrate(&explicit, &old, dev, v, 1)?;
    for i in 0..6 {
        assert!(
            (r.update.state.effective_stress_compression_positive_kpa[i]
                - b.update.state.effective_stress_compression_positive_kpa[i])
                .abs()
                < 1e-9
        );
        assert!(
            (r.update.state.plastic_strain_compression_positive[i]
                - b.update.state.plastic_strain_compression_positive[i])
                .abs()
                < 1e-12
        );
    }
    Ok(())
}

fn transition_parameters(c: f64, phi: f64, psi: f64) -> Parameters {
    Parameters {
        young_kpa: 30000.,
        poisson: 0.25,
        cohesion_kpa: c,
        friction_deg: phi,
        dilation_deg: psi,
        tension_cutoff_kpa: None,
        flow: Flow::NonAssociated,
    }
}
fn transition_state(stress: [f64; 6]) -> State {
    State {
        effective_stress_compression_positive_kpa: stress,
        plastic_strain_compression_positive: [0.003, -0.001, -0.002, 0., 0., 0.],
        void_ratio: 1.,
        equivalent_plastic_strain: 0.02,
        plastic_work_kpa: 1.,
        last_return: Return::Elastic,
    }
}
#[test]
fn zero_volume_tresca_weakening_preserves_prior_history_and_matches_closed_form() -> Result<()> {
    let from = transition_parameters(20., 0., 0.);
    let to = transition_parameters(10., 0., 0.);
    let old = transition_state([120., 100., 80., 0., 0., 0.]);
    let before = bits(&old);
    assert!(mc_split::integrate(&to, &old, [0.; 6], 0., 1).is_err());
    let r = mc_split::integrate_transition(&to, &from, &old, [0.; 6], 0., 1)?;
    let shear = to.young_kpa / (2. * (1. + to.poisson));
    let lambda = 20. / (4. * shear);
    for (i, expected) in [110., 100., 90., 0., 0., 0.].iter().enumerate() {
        assert!(
            (r.update.state.effective_stress_compression_positive_kpa[i] - expected).abs() < 1e-10
        );
        let ep = old.plastic_strain_compression_positive[i] + [lambda, 0., -lambda, 0., 0., 0.][i];
        assert!((r.update.state.plastic_strain_compression_positive[i] - ep).abs() < 1e-13);
    }
    assert!(
        (r.update.state.equivalent_plastic_strain
            - old.equivalent_plastic_strain
            - 2. * lambda / 3f64.sqrt())
        .abs()
            < 1e-13
    );
    assert!((r.update.state.plastic_work_kpa - old.plastic_work_kpa - 20. * lambda).abs() < 1e-12);
    assert_eq!(
        r.update.state.void_ratio.to_bits(),
        old.void_ratio.to_bits()
    );
    assert_eq!(bits(&old), before);
    same(
        &r.update,
        &to.integrate_transition_partition(&from, &old, [0.; 6], 1)?,
    );
    Ok(())
}

#[test]
fn nonassociated_parameter_transition_split_jacobian_matches_independent_trial_directions()
-> Result<()> {
    let from = transition_parameters(20., 30., 10.);
    let to = transition_parameters(10., 20., 10.);
    let old = transition_state([180., 100., 70., 0., 0., 0.]);
    let before = bits(&old);
    let r = mc_split::evaluate_transition(&to, &from, &old, [0.; 6], 0., 1, 1e-7)?;
    same(
        &r.values.update,
        &to.integrate_transition_partition(&from, &old, [0.; 6], 1)?,
    );
    let direction = [0.27, -0.13, 0.08, 0.21, -0.19, 0.09];
    let (dd, dv) = split(direction);
    let h = 4e-8;
    let a = to.integrate_transition_partition(&from, &old, direction.map(|v| v * h), 1)?;
    let b = to.integrate_transition_partition(&from, &old, direction.map(|v| -v * h), 1)?;
    let mean = (0..6)
        .map(|j| r.tangent.mean_by_deviatoric[j] * dd[j])
        .sum::<f64>()
        + r.tangent.mean_by_volume * dv;
    let mut err2 = 0.;
    let mut norm2 = 0.;
    for i in 0..6 {
        let predicted = (0..6)
            .map(|j| r.tangent.deviatoric_by_deviatoric[i][j] * dd[j])
            .sum::<f64>()
            + r.tangent.deviatoric_by_volume[i] * dv
            + if i < 3 { mean } else { 0. };
        let actual = (a.state.effective_stress_compression_positive_kpa[i]
            - b.state.effective_stress_compression_positive_kpa[i])
            / (2. * h);
        err2 += (actual - predicted).powi(2);
        norm2 += actual * actual;
    }
    assert!(norm2 > 0. && (err2 / norm2).sqrt() < 3e-5);
    assert_eq!(bits(&old), before);
    Ok(())
}

#[test]
fn unchanged_parameter_transitions_equal_ordinary_map_and_invalid_transitions_never_mutate_old_state()
-> Result<()> {
    let p = parameters(7.);
    let old = initial();
    let (dev, v) = split([0.01, -0.005, -0.005, 0.012, 0.001, -0.002]);
    let r = mc_split::integrate_transition(&p, &p, &old, dev, v, 8)?;
    same(&r.update, &mc_split::integrate(&p, &old, dev, v, 8)?.update);
    let before = bits(&old);
    for mut to in [p.clone(), p.clone(), p.clone()] {
        if to.young_kpa == p.young_kpa {
            to.young_kpa *= 2.;
        }
        assert!(mc_split::integrate_transition(&to, &p, &old, dev, v, 1).is_err());
    }
    let mut to = p.clone();
    to.flow = Flow::Davis;
    assert!(mc_split::integrate_transition(&to, &p, &old, dev, v, 1).is_err());
    assert!(mc_split::integrate_transition(&p, &p, &old, dev, v, 0).is_err());
    assert!(mc_split::evaluate_transition(&p, &p, &old, dev, v, 1, f64::NAN).is_err());
    assert_eq!(bits(&old), before);
    Ok(())
}
