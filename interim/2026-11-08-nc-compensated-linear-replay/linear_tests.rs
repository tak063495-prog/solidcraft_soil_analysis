#[path = "../src/nc_linear_refinement.rs"]
mod refinement;
use nalgebra::{DMatrix, DVector};

#[test]
fn original_and_compensated_gates_do_not_accept_a_cancelled_wrong_solution() {
    // In binary64, the ordinary first-row dot can erase the unit term.
    // A solver must reject instead of accepting the apparent zero residual.
    let a = DMatrix::from_row_slice(3, 3, &[1e16, 1., -1e16, 1., 0., 0., 0., 0., 1.]);
    let b = DVector::from_column_slice(&[1., 1., 1.]);
    assert!(refinement::solve(&a, &b, 1e-10).is_none());
}

#[test]
fn unsymmetric_equations_and_zero_rhs_keep_the_original_system() {
    let a = DMatrix::from_row_slice(2, 2, &[4., 2., 1., 3.]);
    let b = DVector::from_column_slice(&[1., 2.]);
    let (x, _, r, c) = refinement::solve(&a, &b, 1e-10).unwrap();
    assert!((x[0] + 0.1).abs() < 1e-15 && (x[1] - 0.7).abs() < 1e-15);
    assert!(r <= 1e-10 && c <= 1e-10);
    let (zero, _, r, c) = refinement::solve(&a, &DVector::zeros(2), 1e-10).unwrap();
    assert_eq!(zero, DVector::zeros(2));
    assert_eq!((r, c), (0., 0.));
}

#[test]
fn invalid_inputs_and_singular_systems_are_rejected() {
    let a = DMatrix::identity(2, 2);
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(refinement::solve(&a, &DVector::zeros(2), tolerance).is_none());
    }
    assert!(refinement::solve(&a, &DVector::zeros(3), 1e-10).is_none());
    assert!(refinement::solve(&a, &DVector::from_element(2, f64::NAN), 1e-10).is_none());
    assert!(refinement::solve(&DMatrix::zeros(2, 2), &DVector::from_element(2, 1.), 1e-10).is_none());
}
