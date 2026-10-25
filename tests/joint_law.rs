use nalgebra::{Matrix3, Rotation3, Unit, Vector3};
use soil_fem_3d::{
    Result,
    joint_law::{Material, PointState, Status},
};

fn material() -> Material {
    Material {
        normal_stiffness_kpa_m: 1e6,
        tangential_stiffness_kpa_m: 2e5,
        friction_coefficient: 0.3,
    }
}

#[test]
fn normal_compression_stick_and_opening_obey_unilateral_friction_and_work() -> Result<()> {
    let m = material();
    let s = PointState::initial([0., 0., 1.])?;
    let closed = m.integrate(&s, [0.001, 0., -0.001], 1e-5)?;
    assert!((closed.traction_kpa[2] + 1000.).abs() < 1e-10);
    assert!((closed.traction_kpa[0] - 200.).abs() < 1e-9);
    assert_eq!(closed.plastic_slip_m, [0.; 3]);
    assert_eq!(closed.status, Status::Stick);
    let open = m.integrate(&closed, [0.001, 0., 0.002], 1e-5)?;
    assert_eq!(open.status, Status::Open);
    assert_eq!(open.traction_kpa, [0.; 3]);
    assert!(open.plastic_dissipation_kpa_m > 0. && open.numerical_dissipation_kpa_m > 0.);
    assert!(
        (open.endpoint_work_kpa_m
            - open.plastic_dissipation_kpa_m
            - open.numerical_dissipation_kpa_m)
            .abs()
            < 1e-8
    );
    Ok(())
}

#[test]
fn two_direction_reversal_preserves_residual_slip_and_positive_dissipation() -> Result<()> {
    let m = material();
    let mut s = PointState::initial([0., 0., 1.])?;
    for target in [
        [0., 0., -0.001],
        [0.01, 0.006, -0.001],
        [-0.01, -0.003, -0.001],
        [0., 0., -0.001],
    ] {
        let old = s.clone();
        s = m.integrate(&s, target, 1e-5)?;
        assert!(s.accumulated_plastic_slip_m >= old.accumulated_plastic_slip_m);
        assert!(s.plastic_dissipation_kpa_m >= old.plastic_dissipation_kpa_m);
    }
    assert_eq!(s.status, Status::Slip);
    let tau = Vector3::new(s.traction_kpa[0], s.traction_kpa[1], 0.);
    assert!((tau.norm() - 300.).abs() < 1e-8);
    assert!(s.accumulated_plastic_slip_m > Vector3::from(s.plastic_slip_m).norm());
    assert!((Vector3::from(s.plastic_slip_m) + tau / m.tangential_stiffness_kpa_m).norm() < 1e-10);
    Ok(())
}

#[test]
fn open_motion_is_not_plastic_sliding_and_recontact_starts_unloaded() -> Result<()> {
    let m = material();
    let initial = PointState::initial([0., 0., 1.])?;
    let open = m.integrate(&initial, [0., 0., 0.002], 1e-4)?;
    let moved = m.integrate(&open, [0.1, -0.2, 0.003], 1e-4)?;
    assert_eq!(moved.plastic_slip_m, [0.; 3]);
    assert_eq!(moved.accumulated_plastic_slip_m, 0.);
    assert_eq!(moved.endpoint_work_kpa_m, 0.);
    let closed = m.integrate(&moved, [0.1, -0.2, -0.001], 1e-5)?;
    assert!(Vector3::new(closed.traction_kpa[0], closed.traction_kpa[1], 0.).norm() < 1e-8);
    assert_eq!(closed.plastic_slip_m, [0.; 3]);
    let sheared = m.integrate(&closed, [0.1001, -0.1998, -0.001], 1e-5)?;
    assert!(
        (sheared.traction_kpa[0] - 20.).abs() < 1e-7
            && (sheared.traction_kpa[1] - 40.).abs() < 1e-7
    );
    Ok(())
}

#[test]
fn global_rotation_and_restart_preserve_all_interface_history() -> Result<()> {
    let m = material();
    let r = Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(1., 2., -3.)), 0.8)
        .into_inner();
    let mut a = PointState::initial([0., 0., 1.])?;
    let mut b = a.rotated(&r)?;
    for target in [
        [0.004, 0.002, -0.001],
        [-0.003, 0.004, -0.0005],
        [0.01, -0.02, 0.001],
        [0.011, -0.021, -0.001],
    ] {
        a = m.integrate(&a, target, 1e-4)?;
        b = m.integrate(&b, (r * Vector3::from(target)).into(), 1e-4)?;
        let expected = a.rotated(&r)?;
        assert!(
            (Vector3::from(b.traction_kpa) - Vector3::from(expected.traction_kpa)).norm() < 1e-7
        );
        assert!(
            (Vector3::from(b.plastic_slip_m) - Vector3::from(expected.plastic_slip_m)).norm()
                < 1e-10
        );
        assert!((b.plastic_dissipation_kpa_m - a.plastic_dissipation_kpa_m).abs() < 1e-8);
        b = serde_json::from_slice(&serde_json::to_vec(&b)?)?;
    }
    Ok(())
}

#[test]
fn algorithmic_tangent_is_nonsymmetric_and_matches_half_step() -> Result<()> {
    let m = material();
    let initial = PointState::initial([0., 0., 1.])?;
    let s = m.integrate(&initial, [0., 0., -0.001], 1e-4)?;
    let target = [0.003, 0.002, -0.0012];
    let count = m.partition(&s, target, 1e-4)?;
    let d = m.tangent(&s, target, count, 1e-8)?;
    let half = m.tangent(&s, target, count, 5e-9)?;
    assert!((d - d.transpose()).norm() > 1e4);
    assert!((d - half).norm() / d.norm() < 1e-6);
    let direction = Vector3::new(0.4, -0.3, 0.2);
    let h: f64 = 2e-8;
    let plus = m.integrate_partition(&s, (Vector3::from(target) + h * direction).into(), count)?;
    let minus = m.integrate_partition(&s, (Vector3::from(target) - h * direction).into(), count)?;
    let fd = (Vector3::from(plus.traction_kpa) - Vector3::from(minus.traction_kpa)) / (2. * h);
    assert!((fd - d * direction).norm() / fd.norm() < 1e-6);
    assert!((m.tangent(&s, [0.004, 0.002, 0.01], count, 1e-8)? - Matrix3::zeros()).norm() < 1e-8);
    Ok(())
}

#[test]
fn numerical_dissipation_refines_and_failed_history_does_not_mutate_state() -> Result<()> {
    let m = material();
    let initial = PointState::initial([0., 0., 1.])?;
    let a = m.integrate(&initial, [0.003, 0.001, -0.001], 1e-4)?;
    let b = m.integrate(&initial, [0.003, 0.001, -0.001], 5e-5)?;
    assert!(b.numerical_dissipation_kpa_m < a.numerical_dissipation_kpa_m);
    let before = initial.clone();
    assert!(m.integrate(&initial, [1.; 3], 1e-8).is_err());
    assert_eq!(initial, before);
    let mut corrupt = a.clone();
    corrupt.accumulated_plastic_slip_m = 0.;
    assert!(m.validate_state(&corrupt).is_err());
    let mut bad = m;
    bad.normal_stiffness_kpa_m = 0.;
    assert!(bad.integrate(&initial, [0.; 3], 1e-4).is_err());
    Ok(())
}
