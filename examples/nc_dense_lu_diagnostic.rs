//! Offline diagnosis of a saved dense system; never accepts a nonlinear step.
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::io;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    matrix: Vec<Vec<f64>>,
    rhs: Vec<f64>,
    refinement_steps: usize,
}
#[derive(Serialize)]
struct Observation {
    refinement: usize,
    solution: Vec<f64>,
    true_relative_residual: f64,
    normwise_backward_error: f64,
    relative_matvec_roundoff_bound: f64,
}
fn observe(a: &DMatrix<f64>, b: &DVector<f64>, x: &DVector<f64>, step: usize) -> Observation {
    let norm = (a * x - b).norm();
    let eps_n = a.ncols() as f64 * f64::EPSILON;
    let gamma = eps_n / (1. - eps_n);
    let scale = a.map(f64::abs) * x.map(f64::abs) + b.map(f64::abs);
    Observation {
        refinement: step,
        solution: x.as_slice().to_vec(),
        true_relative_residual: norm / b.norm().max(1e-30),
        normwise_backward_error: norm / (a.norm() * x.norm() + b.norm()).max(1e-30),
        relative_matvec_roundoff_bound: gamma * scale.norm() / b.norm().max(1e-30),
    }
}
fn diagnose(input: Input) -> Result<Vec<Observation>, Box<dyn std::error::Error>> {
    let n = input.rhs.len();
    if n == 0
        || n > 2000
        || input.refinement_steps > 32
        || input.matrix.len() != n
        || input.matrix.iter().any(|r| r.len() != n)
        || input
            .matrix
            .iter()
            .flatten()
            .chain(&input.rhs)
            .any(|x| !x.is_finite())
    {
        return Err("invalid finite square diagnostic system or refinement count".into());
    }
    let a = DMatrix::from_fn(n, n, |i, j| input.matrix[i][j]);
    let b = DVector::from_vec(input.rhs);
    let lu = a.clone().lu();
    let mut x = lu.solve(&b).ok_or("singular diagnostic matrix")?;
    let mut result = vec![observe(&a, &b, &x, 0)];
    for step in 1..=input.refinement_steps {
        // Reuse exactly the original LU and the original, unmodified system.
        let correction = lu.solve(&(&b - &a * &x)).ok_or("singular correction")?;
        x += correction;
        if !x.iter().all(|v| v.is_finite()) {
            return Err("nonfinite diagnostic correction".into());
        }
        result.push(observe(&a, &b, &x, step));
    }
    Ok(result)
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input: Input = serde_json::from_reader(io::stdin().lock())?;
    serde_json::to_writer(io::stdout().lock(), &diagnose(input)?)?;
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn original_nonsymmetric_system_is_retained() -> Result<(), Box<dyn std::error::Error>> {
        let rows = diagnose(Input {
            matrix: vec![vec![2., 3.], vec![0., 4.]],
            rhs: vec![8., 8.],
            refinement_steps: 3,
        })?;
        assert_eq!(rows.len(), 4);
        for row in rows {
            assert_eq!(row.solution, vec![1., 2.]);
            assert_eq!(row.true_relative_residual, 0.);
        }
        Ok(())
    }
    #[test]
    fn malformed_and_singular_systems_are_rejected() {
        for matrix in [
            vec![vec![1.]],
            vec![vec![1., 2.], vec![2., 4.]],
            vec![vec![f64::NAN, 0.], vec![0., 1.]],
        ] {
            assert!(
                diagnose(Input {
                    matrix,
                    rhs: vec![1., 1.],
                    refinement_steps: 2
                })
                .is_err()
            );
        }
    }
}
