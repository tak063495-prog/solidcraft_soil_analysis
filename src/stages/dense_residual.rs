//! Opt-in recovery of a failed original-coordinate dense LU residual check.
use super::{DMatrix, DVector, Error, Instant, LinearReport, Result, dense_refinement};
use nalgebra::{Dyn, LU};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DenseResidualMethod {
    RepeatedRefinement,
    ExtendedResidual,
    Ruiz,
    GmresPolish,
    CoordinatePolish,
}

impl DenseResidualMethod {
    fn name(self) -> &'static str {
        match self {
            Self::RepeatedRefinement => "repeated_refinement",
            Self::ExtendedResidual => "extended_residual",
            Self::Ruiz => "ruiz",
            Self::GmresPolish => "gmres_polish",
            Self::CoordinatePolish => "coordinate_polish",
        }
    }
}

fn singular() -> Error {
    Error::Solve("singular non-symmetric tangent; inspect active supports".into())
}

// Knuth TwoSum, including when |a| < |b|.
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    (s, (a - (s - v)) + (b - v))
}

// Distill both parts of every product through two error-free sums.
// Only the third (lowest) accumulator is summed without compensation.
fn extended_residual(a: &DMatrix<f64>, x: &DVector<f64>, b: &DVector<f64>) -> DVector<f64> {
    DVector::from_fn(b.len(), |i, _| {
        let (mut hi, mut mid, mut lo) = (b[i], 0., 0.);
        for j in 0..x.len() {
            let coefficient = -a[(i, j)];
            let product = coefficient * x[j];
            for term in [product, coefficient.mul_add(x[j], -product)] {
                let (s, t) = two_sum(hi, term);
                let (m, u) = two_sum(mid, t);
                hi = s;
                mid = m;
                lo += u;
            }
        }
        let (s, t) = two_sum(hi, mid);
        s + (t + lo)
    })
}

struct System<'a> {
    a: &'a DMatrix<f64>,
    b: &'a DVector<f64>,
    scale: f64,
    tolerance: f64,
}

impl System<'_> {
    fn residuals(&self, x: &DVector<f64>) -> (f64, f64) {
        let absolute = (self.a * x - self.b).norm();
        (absolute, absolute / self.scale)
    }
}

// Simultaneous infinity-norm row/column equilibration (Ruiz, Alg. 2.1).
fn equilibrate(system: &System<'_>) -> Result<(DVector<f64>, usize)> {
    let mut a = system.a.clone();
    let n = a.nrows();
    let mut rhs = system.b.clone();
    let mut col_scale = DVector::from_element(n, 1.);
    let mut sweeps = 0;
    let mut previous = f64::INFINITY;
    for _ in 0..8 {
        let mut rows = DVector::<f64>::zeros(n);
        let mut cols = DVector::<f64>::zeros(n);
        for (j, col) in a.column_iter().enumerate() {
            for (i, &v) in col.iter().enumerate() {
                rows[i] = rows[i].max(v.abs());
                cols[j] = cols[j].max(v.abs());
            }
        }
        if rows
            .iter()
            .chain(cols.iter())
            .any(|v| !v.is_finite() || *v == 0.)
        {
            return Err(singular());
        }
        let deviation = rows
            .iter()
            .chain(cols.iter())
            .map(|v| (v - 1.).abs())
            .fold(0., f64::max);
        if deviation <= 0.05 || deviation >= previous {
            break;
        }
        previous = deviation;
        rows.apply(|v| *v = v.sqrt().recip());
        cols.apply(|v| *v = v.sqrt().recip());
        for (j, mut col) in a.column_iter_mut().enumerate() {
            for (i, v) in col.iter_mut().enumerate() {
                *v = (*v * rows[i]) * cols[j];
            }
        }
        rhs.component_mul_assign(&rows);
        col_scale.component_mul_assign(&cols);
        sweeps += 1;
    }
    let mut y = a.lu().solve(&rhs).ok_or_else(singular)?;
    y.component_mul_assign(&col_scale);
    Ok((y, sweeps))
}

fn correction(
    system: &System<'_>,
    x: &DVector<f64>,
    lu: &LU<f64, Dyn, Dyn>,
    method: DenseResidualMethod,
) -> Result<DVector<f64>> {
    let a = system.a;
    let b = system.b;
    let rhs = match method {
        DenseResidualMethod::RepeatedRefinement | DenseResidualMethod::CoordinatePolish => {
            dense_refinement::compensated_residual(a, x, b)
        }
        DenseResidualMethod::ExtendedResidual => extended_residual(a, x, b),
        DenseResidualMethod::GmresPolish => b - a * x,
        DenseResidualMethod::Ruiz => unreachable!("equilibration uses its own factors"),
    };
    let mut direction = lu.solve(&rhs).ok_or_else(singular)?;
    if method == DenseResidualMethod::GmresPolish {
        // One right-preconditioned GMRES(1) cycle: minimize ||r - alpha A z||.
        let image = a * &direction;
        let norm = image.norm();
        if norm == 0. || !norm.is_finite() {
            return Err(Error::Solve("dense GMRES polish breakdown".into()));
        }
        let alpha = (image / norm).dot(&rhs) / norm;
        direction *= alpha;
    }
    Ok(direction)
}

pub(super) fn solve(
    a: &DMatrix<f64>,
    rhs: &DVector<f64>,
    tolerance: f64,
    method: DenseResidualMethod,
) -> Result<(DVector<f64>, LinearReport)> {
    let started = Instant::now();
    let system = System {
        a,
        b: rhs,
        scale: rhs.norm().max(1e-30),
        tolerance,
    };
    let lu = a.clone().lu();
    let mut x = lu.solve(rhs).ok_or_else(singular)?;
    let (mut absolute, mut residual) = system.residuals(&x);
    let initial = residual;
    let mut steps = 0;
    let mut sweeps = 0;
    let mut min_pivot = f64::NAN;
    if initial.is_finite() && initial > system.tolerance {
        min_pivot = (0..a.nrows())
            .map(|i| lu.lu_internal()[(i, i)].abs())
            .fold(f64::INFINITY, f64::min);
        if method == DenseResidualMethod::Ruiz {
            drop(lu);
            (x, sweeps) = equilibrate(&system)?;
            steps = 1;
            (absolute, residual) = system.residuals(&x);
        } else if method == DenseResidualMethod::CoordinatePolish {
            x += correction(&system, &x, &lu, method)?;
            steps = 1;
            (absolute, residual) = system.residuals(&x);
            if residual.is_finite() && residual > tolerance {
                let count;
                (count, absolute, residual) = coordinate::polish(&system, &mut x);
                steps += count;
            }
        } else {
            let limit = if method == DenseResidualMethod::ExtendedResidual {
                1
            } else {
                5
            };
            let mut best = initial;
            let mut stale = 0;
            for _ in 0..limit {
                x += correction(&system, &x, &lu, method)?;
                steps += 1;
                (absolute, residual) = system.residuals(&x);
                if !residual.is_finite() || residual <= system.tolerance {
                    break;
                }
                if residual < best {
                    best = residual;
                    stale = 0;
                } else {
                    stale += 1;
                    if stale == 2 {
                        break;
                    }
                }
            }
        }
    }
    let name = method.name();
    let absolute = (method == DenseResidualMethod::CoordinatePolish).then_some(absolute);
    if !residual.is_finite() || residual > tolerance || !x.iter().all(|v| v.is_finite()) {
        if let Some(absolute) = absolute {
            return Err(Error::Solve(format!(
                "true linear residual {residual:e} ({name}, absolute {absolute:e}, initial {initial:e}, steps {steps}, unscaled min abs U pivot {min_pivot:e})"
            )));
        }
        return Err(Error::Solve(format!(
            "true linear residual {residual:e} ({name}, initial {initial:e}, steps {steps}, scale sweeps {sweeps}, unscaled min abs U pivot {min_pivot:e})"
        )));
    }
    let n = a.nrows();
    let vectors = match method {
        DenseResidualMethod::Ruiz
        | DenseResidualMethod::GmresPolish
        | DenseResidualMethod::CoordinatePolish => 5,
        _ => 3,
    };
    let report = LinearReport {
        method: name.into(),
        dimension: n,
        stored_nonzeros: a.len(),
        iterations: 1 + steps,
        restarts: if method == DenseResidualMethod::GmresPolish {
            steps.saturating_sub(1)
        } else {
            0
        },
        true_absolute_residual: absolute,
        true_relative_residual: residual,
        matrix_storage_bytes: a.len() * size_of::<f64>(),
        workspace_storage_bytes: (a.len() + vectors * n) * size_of::<f64>()
            + n * size_of::<usize>(),
        elapsed_seconds: started.elapsed().as_secs_f64(),
    };
    Ok((x, report))
}

mod coordinate;

#[cfg(test)]
mod tests;
