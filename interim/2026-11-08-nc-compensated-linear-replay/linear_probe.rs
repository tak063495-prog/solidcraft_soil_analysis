//! Saved-matrix probe, no constitutive or continuation solve.
#[path = "../src/nc_linear_refinement.rs"]
mod refinement;
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::{fs, path::Path};
#[derive(Deserialize)]
struct Request {
    matrix_columns: Vec<Vec<f64>>,
    reduced_residual: Vec<f64>,
    label: String,
}
#[derive(Serialize)]
struct Response {
    label: String,
    baseline_x: Option<Vec<f64>>,
    baseline_relative_residual: Option<f64>,
    candidate_x: Option<Vec<f64>>,
    corrections: Option<usize>,
    original_relative_residual: Option<f64>,
    compensated_relative_residual: Option<f64>,
    accepted: bool,
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("linear_probe REQUEST NEW_OUTPUT".into());
    }
    if Path::new(&args[2]).exists() {
        return Err("output exists".into());
    }
    let requests: Vec<Request> = serde_json::from_slice(&fs::read(&args[1])?)?;
    let mut responses = vec![];
    for r in requests {
        let n = r.reduced_residual.len();
        if r.matrix_columns.len() != n || !r.matrix_columns.iter().all(|c| c.len() == n) {
            return Err("invalid matrix size".into());
        }
        let a = DMatrix::from_fn(n, n, |i, j| r.matrix_columns[j][i]);
        let b = DVector::from_vec(r.reduced_residual);
        let baseline = a.clone().lu().solve(&b);
        let baseline_relative_residual = baseline
            .as_ref()
            .map(|x| (&a * x - &b).norm() / b.norm().max(1e-30));
        let candidate = refinement::solve(&a, &b, 1e-10);
        let (candidate_x, corrections, original, accurate) = match candidate {
            Some((x, c, o, h)) => (Some(x.as_slice().to_vec()), Some(c), Some(o), Some(h)),
            None => (None, None, None, None),
        };
        responses.push(Response {
            label: r.label,
            baseline_x: baseline.map(|x| x.as_slice().to_vec()),
            baseline_relative_residual,
            accepted: candidate_x.is_some(),
            candidate_x,
            corrections,
            original_relative_residual: original,
            compensated_relative_residual: accurate,
        });
    }
    fs::write(&args[2], serde_json::to_vec_pretty(&responses)?)?;
    Ok(())
}
