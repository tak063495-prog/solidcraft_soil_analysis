//! Bounded LU residual correction. Both original and compensated gates must pass.
use nalgebra::{DMatrix, DVector};

// Ogita/Rump/Oishi-style product error with Neumaier summation of all terms.
fn compensated_residual(a: &DMatrix<f64>, x: &DVector<f64>, b: &DVector<f64>) -> DVector<f64> {
    DVector::from_fn(b.len(), |i, _| {
        let mut sum = b[i];
        let mut error = 0.;
        for j in 0..x.len() {
            let product = -a[(i, j)] * x[j];
            let product_error = (-a[(i, j)]).mul_add(x[j], -product);
            let next = sum + product;
            error += if sum.abs() >= product.abs() {
                (sum - next) + product
            } else {
                (product - next) + sum
            } + product_error;
            sum = next;
        }
        sum + error
    })
}

pub fn solve(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    tolerance: f64,
) -> Option<(DVector<f64>, usize, f64, f64)> {
    if a.nrows() != a.ncols()
        || a.nrows() != b.len()
        || !tolerance.is_finite()
        || tolerance <= 0.
        || !a.iter().chain(b.iter()).all(|v| v.is_finite())
    {
        return None;
    }
    let lu = a.clone().lu();
    let mut x = lu.solve(b)?;
    let scale = b.norm().max(1e-30);
    for corrections in 0..=8 {
        let original = (a * &x - b).norm() / scale;
        let residual = compensated_residual(a, &x, b);
        let accurate = residual.norm() / scale;
        if x.iter().all(|v| v.is_finite())
            && original.is_finite()
            && accurate.is_finite()
            && original <= tolerance
            && accurate <= tolerance
        {
            return Some((x, corrections, original, accurate));
        }
        if corrections < 8 {
            x += lu.solve(&residual)?;
        }
    }
    None
}
