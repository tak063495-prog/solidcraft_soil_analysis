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
    (
        std::array::from_fn(|i| full[i] - if i < 3 { volume / 3. } else { 0. }),
        volume,
    )
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
        let (bulk, shear) = p.moduli();
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
