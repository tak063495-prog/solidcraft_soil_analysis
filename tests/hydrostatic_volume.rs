use soil_fem_3d::{Result, elements, hydrostatic_volume::affine_load};

fn points() -> Vec<[f64; 3]> {
    let mut p = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    for [i, j] in elements::EDGES {
        p.push(std::array::from_fn(|a| (p[i][a] + p[j][a]) / 2.));
    }
    p
}

#[test]
fn partial_water_virtual_work_matches_analytic_pressure_integral() -> Result<()> {
    for count in [4, 10] {
        let p = points();
        for h in [0., 1e-5, 0.25, 0.75, 0.99999, 1.] {
            let f = affine_load(&p[..count], h, 10.)?.unwrap();
            let integral = 10. * (h * h / 4. - h.powi(3) / 6. + h.powi(4) / 24.);
            for a in 0..3 {
                assert!(f.iter().map(|v| v[a]).sum::<f64>().abs() < 1e-13);
                for b in [0, 1, 2] {
                    let work: f64 = (0..count).map(|i| p[i][b] * f[i][a]).sum();
                    assert!((work + if a == b { integral } else { 0. }).abs() < 1e-13);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn fully_wet_quadratic_force_matches_polynomial_material_quadrature() -> Result<()> {
    let p = points();
    let f = affine_load(&p, 1.5, 10.)?.unwrap();
    let mut expected = [[0.; 3]; 10];
    for q in elements::volume_points(&p, 3)? {
        for (i, value) in expected.iter_mut().enumerate() {
            for (a, component) in value.iter_mut().enumerate() {
                *component -= q.b[(a, 3 * i + a)] * 10. * (1.5 - q.position[2]) * q.weight_m3;
            }
        }
    }
    assert!(
        f.iter()
            .flatten()
            .zip(expected.iter().flatten())
            .all(|(a, b)| (a - b).abs() < 1e-13)
    );
    Ok(())
}

#[test]
fn curved_mapping_and_invalid_data_do_not_take_affine_shortcut() -> Result<()> {
    let mut p = points();
    p[4][2] += 1e-9;
    assert!(affine_load(&p, 0.75, 10.)?.is_none());
    assert!(affine_load(&p[..3], 0.75, 10.).is_err());
    assert!(affine_load(&p, f64::NAN, 10.).is_err());
    assert!(affine_load(&p, 0.75, -1.).is_err());
    Ok(())
}
