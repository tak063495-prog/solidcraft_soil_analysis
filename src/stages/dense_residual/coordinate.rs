//! Cyclic coordinate minimization of ||b - A x||₂ on a rejected dense solve.
use super::{DVector, System, extended_residual};

pub(super) fn polish(system: &System<'_>, x: &mut DVector<f64>) -> (usize, f64, f64) {
    // Scaling each search direction cancels algebraically in its projection and
    // avoids overflow or underflow when forming squared column norms.
    let columns: Vec<_> = system
        .a
        .column_iter()
        .map(|column| {
            let scale = column.amax();
            let norm_squared = if scale.is_finite() && scale > 0. {
                column.iter().map(|v| (v / scale).powi(2)).sum()
            } else {
                0.
            };
            (scale, norm_squared)
        })
        .collect();
    const MAX_SWEEPS: usize = 64;
    let mut sweep = 0;
    loop {
        sweep += 1;
        let mut r = extended_residual(system.a, x, system.b);
        let mut changed = false;
        for (j, column) in system.a.column_iter().enumerate() {
            let (scale, norm_squared) = columns[j];
            if norm_squared == 0. {
                continue;
            }
            let dot: f64 = column.iter().zip(&r).map(|(a, ri)| (a / scale) * ri).sum();
            let step = (dot / norm_squared) / scale;
            let candidate = x[j] + step;
            let applied = candidate - x[j];
            if applied != 0. && applied.is_finite() && candidate.is_finite() {
                x[j] = candidate;
                r.axpy(-applied, &column, 1.);
                changed = true;
            }
        }
        let (absolute, residual) = system.residuals(x);
        if !changed || !residual.is_finite() || residual <= system.tolerance || sweep == MAX_SWEEPS
        {
            return (sweep, absolute, residual);
        }
    }
}
