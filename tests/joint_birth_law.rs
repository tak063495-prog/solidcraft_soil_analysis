use nalgebra::Vector3;
use soil_fem_3d::{Result, joint_law::Material};
#[test]
fn inherited_energy_work_and_stress_free_reference_obey_independent_increments() -> Result<()> {
    let m = Material {
        normal_stiffness_kpa_m: 1e6,
        tangential_stiffness_kpa_m: 2e5,
        friction_coefficient: 0.3,
    };
    let reference = [0.02, -0.01, 0.03];
    let initial = m.birth([0., 0., 1.], reference, [12., 5., -100.])?;
    assert_eq!(initial.plastic_slip_m, [0.; 3]);
    assert_eq!(initial.accumulated_plastic_slip_m, 0.);
    assert_eq!(initial.plastic_dissipation_kpa_m, 0.);
    assert_eq!(initial.numerical_dissipation_kpa_m, 0.);
    assert_eq!(initial.endpoint_work_kpa_m, 0.);
    let mut s = initial.clone();
    for increment in [
        [1e-5, -2e-5, -1e-5],
        [-1e-5, 2e-5, 0.],
        [0.003, 0.001, -1e-4],
        [0.003, 0.001, 0.001],
    ] {
        s = m.integrate(
            &s,
            std::array::from_fn(|a| reference[a] + increment[a]),
            1e-5,
        )?;
        let balance = initial.birth.as_ref().unwrap().initial_energy_kpa_m + s.endpoint_work_kpa_m
            - m.energy(&s)
            - s.plastic_dissipation_kpa_m
            - s.numerical_dissipation_kpa_m;
        assert!(balance.abs() < 1e-10);
    }
    assert_eq!(s.traction_kpa, [0.; 3]);
    let b = m.birth([0., 0., 1.], reference, [0.; 3])?;
    let compressed = m.integrate(
        &b,
        [reference[0], reference[1], reference[2] - 50. / 1e6],
        1e-5,
    )?;
    assert!((compressed.traction_kpa[2] + 50.).abs() < 1e-8);
    assert!(Vector3::from(compressed.plastic_slip_m).norm() == 0.);
    assert!(m.birth([0., 0., 1.], reference, [0., 0., 1.]).is_err());
    assert!(m.birth([0., 0., 1.], reference, [31., 0., -100.]).is_err());
    Ok(())
}
