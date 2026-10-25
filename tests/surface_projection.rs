use soil_fem_3d::{Result, surface_projection::upward_area};
#[test]
fn curved_triangle_rain_projection_crossing_zero_normal_matches_exact_integral() -> Result<()> {
    let refs = [
        [0., 0.],
        [1., 0.],
        [0., 1.],
        [0.5, 0.],
        [0.5, 0.5],
        [0., 0.5],
    ];
    let mut points = refs.map(|[r, s]| [r + 0.5 * r * r, (r - 0.5) * s, s]);
    let up = upward_area(&points)?;
    assert!((up.area_m2 - 7. / 192.).abs() <= up.domain_error_bound_m2 + 1e-14);
    assert!(up.domain_error_bound_m2 < 2e-11 && up.visited_cells > 1);
    points.swap(1, 2);
    points.swap(3, 5);
    let reverse = upward_area(&points)?;
    assert!((reverse.area_m2 - 23. / 192.).abs() <= reverse.domain_error_bound_m2 + 1e-14);
    Ok(())
}
#[test]
fn planar_projection_handles_vertical_downward_and_scaled_shifted_geometry() -> Result<()> {
    assert_eq!(
        upward_area(&[[0., 0., 0.], [0., 1., 0.], [0., 0., 1.]])?.area_m2,
        0.
    );
    assert_eq!(
        upward_area(&[[0., 0., 0.], [0., 1., 0.], [1., 0., 0.]])?.area_m2,
        0.
    );
    let q = upward_area(&[[1e8, 1e8, 1e8], [1e8 + 3., 1e8, 1e8], [1e8, 1e8 + 3., 1e8]])?;
    assert!((q.area_m2 - 4.5).abs() < 1e-14);
    assert_eq!(q.domain_error_bound_m2, 0.);
    Ok(())
}
