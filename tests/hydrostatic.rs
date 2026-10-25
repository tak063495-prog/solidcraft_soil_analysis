use nalgebra::Vector3;
use soil_fem_3d::{
    Result, elements,
    hydrostatic::{self, Integration},
};

fn nodes(f: impl Fn(f64, f64) -> [f64; 3]) -> Vec<[f64; 3]> {
    [
        [0., 0.],
        [1., 0.],
        [0., 1.],
        [0.5, 0.],
        [0.5, 0.5],
        [0., 0.5],
    ]
    .map(|r| f(r[0], r[1]))
    .to_vec()
}
fn multiply(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            c[i + j] += x * y;
        }
    }
    c
}
fn integral(a: &[f64], x: f64) -> f64 {
    a.iter()
        .enumerate()
        .map(|(i, a)| a * x.powi(i as i32 + 1) / (i + 1) as f64)
        .sum()
}

#[test]
fn planar_partial_water_force_and_moment_match_closed_polynomial_integrals() -> Result<()> {
    for count in [3, 6] {
        let p = nodes(|r, s| [0., r, s]);
        let p = &p[..count];
        let h: f64 = 0.37;
        let gamma = 10.;
        let q = hydrostatic::face_load(p, h, gamma, Integration::default())?;
        let mut force = Vector3::zeros();
        let mut moment = Vector3::zeros();
        for (p, f) in p.iter().zip(&q.nodal_force_kn) {
            force += Vector3::from(*f);
            moment += Vector3::from(*p).cross(&Vector3::from(*f));
        }
        let expected = Vector3::new(-gamma * (h * h / 2. - h.powi(3) / 6.), 0., 0.);
        let expected_m = Vector3::new(
            0.,
            -gamma * (h.powi(3) / 6. - h.powi(4) / 12.),
            gamma / 2. * (h * h / 2. - h.powi(3) / 3. + h.powi(4) / 12.),
        );
        assert!((force - expected).norm() < 1e-13);
        assert!((moment - expected_m).norm() < 1e-13);
        assert!(q.domain_error_bound_kn < 1e-13);
    }
    Ok(())
}

#[test]
fn curved_partial_face_matches_independent_root_and_all_six_nodal_integrals() -> Result<()> {
    let e: f64 = 0.2;
    let a = 1. - e;
    let h: f64 = 0.37;
    let gamma = 10.;
    let p = nodes(|r, s| [0., r, e * s + a * s * s]);
    let upper = (-e + (e * e + 4. * a * h).sqrt()) / (2. * a);
    let pressure_j = multiply(&[h, -e, -a], &[e, 2. * a]);
    // Exact integrals of each original Tri6 shape over r in [0,1-s].
    let shape_integrals = [
        vec![1. / 6., -1., 1.5, -2. / 3.],
        vec![1. / 6., -1., 1.5, -2. / 3.],
        vec![0., -1., 3., -2.],
        vec![2. / 3., -2., 2., -2. / 3.],
        vec![0., 2., -4., 2.],
        vec![0., 2., -4., 2.],
    ];
    let expected: Vec<_> = shape_integrals
        .iter()
        .map(|n| [-gamma * integral(&multiply(n, &pressure_j), upper), 0., 0.])
        .collect();
    let mut errors = Vec::new();
    for tol in [1e-7, 1e-10] {
        let q = hydrostatic::face_load(
            &p,
            h,
            gamma,
            Integration {
                relative_tolerance: tol,
                absolute_tolerance_kn: tol / 10.,
                ..Integration::default()
            },
        )?;
        let error = q
            .nodal_force_kn
            .iter()
            .flatten()
            .zip(expected.iter().flatten())
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(
            error <= q.domain_error_bound_kn + 2e-13,
            "{error:e} > {}",
            q.domain_error_bound_kn
        );
        assert!(q.domain_error_bound_kn <= q.permitted_domain_error_kn);
        errors.push(error);
    }
    assert!(errors[1] < errors[0] / 10., "{errors:?}");
    Ok(())
}

#[test]
fn wet_interior_island_is_detected_when_all_six_surface_nodes_are_dry() -> Result<()> {
    let radius2 = 0.01;
    let gamma = 10.;
    let p = nodes(|r, s| {
        [
            r,
            s,
            (r - 1. / 3.).powi(2) + (s - 1. / 3.).powi(2) - radius2,
        ]
    });
    assert!(p.iter().all(|p| p[2] > 0.));
    let q = hydrostatic::face_load(&p, 0., gamma, Integration::default())?;
    let force = q
        .nodal_force_kn
        .iter()
        .fold(Vector3::zeros(), |a, f| a + Vector3::from(*f));
    let expected = Vector3::new(0., 0., -gamma * std::f64::consts::PI * radius2.powi(2) / 2.);
    assert!((force - expected).norm() < 6_f64.sqrt() * q.domain_error_bound_kn + 1e-12);
    assert!(force[2] < -1e-3 && q.visited_cells > 1);
    Ok(())
}

#[test]
fn wet_dry_limits_reversed_normal_and_world_translation_are_consistent() -> Result<()> {
    let p = nodes(|r, s| [0., r, 0.2 * s + 0.8 * s * s]);
    for h in [-1., 2.] {
        let q = hydrostatic::face_load(&p, h, 10., Integration::default())?;
        let mut expected = [[0.; 3]; 6];
        for x in elements::face_points(&p)? {
            for (i, n) in x.shape.iter().enumerate() {
                for (a, value) in expected[i].iter_mut().enumerate() {
                    *value -= 10. * (h - x.position[2]).max(0.) * n * x.area_vector[a];
                }
            }
        }
        assert!(
            q.nodal_force_kn
                .iter()
                .flatten()
                .zip(expected.iter().flatten())
                .all(|(a, b)| (a - b).abs() < 1e-13)
        );
    }
    let q = hydrostatic::face_load(&p, 0.37, 10., Integration::default())?;
    let shift = [1e6, -2e6, 10.];
    let translated: Vec<_> = p
        .iter()
        .map(|p| std::array::from_fn(|a| p[a] + shift[a]))
        .collect();
    let t = hydrostatic::face_load(&translated, 10.37, 10., Integration::default())?;
    assert!(
        q.nodal_force_kn
            .iter()
            .flatten()
            .zip(t.nodal_force_kn.iter().flatten())
            .all(|(a, b)| (a - b).abs() < 1e-8)
    );
    let reversed = [p[0], p[2], p[1], p[5], p[4], p[3]];
    let r = hydrostatic::face_load(&reversed, 0.37, 10., Integration::default())?;
    for (i, j) in [0, 2, 1, 5, 4, 3].into_iter().enumerate() {
        for a in 0..3 {
            assert!((r.nodal_force_kn[i][a] + q.nodal_force_kn[j][a]).abs() < 1e-9);
        }
    }
    Ok(())
}

#[test]
fn invalid_geometry_controls_and_unresolved_domain_are_rejected() -> Result<()> {
    let p = nodes(|r, s| [0., r, 0.2 * s + 0.8 * s * s]);
    assert!(
        hydrostatic::face_load(
            &p,
            0.37,
            10.,
            Integration {
                maximum_cells: 1,
                ..Integration::default()
            }
        )
        .is_err()
    );
    assert!(hydrostatic::face_load(&p, f64::NAN, 10., Integration::default()).is_err());
    assert!(hydrostatic::face_load(&p, 0.3, 0., Integration::default()).is_err());
    assert!(hydrostatic::face_load(&[[0.; 3]; 6], 0.3, 10., Integration::default()).is_err());
    // A wholly dry face must not serialize a NaN error budget as JSON null
    // when finite coordinates overflow the surface-area norm.
    assert!(
        hydrostatic::face_load(
            &[[0., 0., 0.], [1e100, 0., 0.], [0., 1e100, 0.]],
            -1.,
            10.,
            Integration::default()
        )
        .is_err()
    );
    Ok(())
}
