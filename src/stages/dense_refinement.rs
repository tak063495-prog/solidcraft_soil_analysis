//! One residual correction with the existing partial-pivot LU factors.
use super::{DMatrix, DVector, Error, Instant, LinearReport, Result};

// FMA recovers the product roundoff. Neumaier accumulation retains low parts
// in the subtraction b_i - sum_j A_ij x_j (also for nonsymmetric tangents).
pub(super) fn compensated_residual(
    a: &DMatrix<f64>,
    x: &DVector<f64>,
    b: &DVector<f64>,
) -> DVector<f64> {
    DVector::from_fn(b.len(), |i, _| {
        let mut sum = b[i];
        let mut error = 0.;
        for j in 0..x.len() {
            let coefficient = -a[(i, j)];
            let product = coefficient * x[j];
            let next = sum + product;
            let addition_error = if sum.abs() >= product.abs() {
                (sum - next) + product
            } else {
                (product - next) + sum
            };
            error += addition_error + coefficient.mul_add(x[j], -product);
            sum = next;
        }
        sum + error
    })
}

pub(super) fn solve(
    a: &DMatrix<f64>,
    rhs: &DVector<f64>,
    tolerance: f64,
) -> Result<(DVector<f64>, LinearReport)> {
    let started = Instant::now();
    let lu = a.clone().lu();
    let singular =
        || Error::Solve("singular non-symmetric tangent; inspect active supports".into());
    let mut x = lu.solve(rhs).ok_or_else(singular)?;
    let scale = rhs.norm().max(1e-30);
    let original_residual = |x: &DVector<f64>| (a * x - rhs).norm() / scale;
    let initial = original_residual(&x);
    let correct = initial.is_finite() && initial > tolerance;
    let residual = if correct {
        let mut correction = compensated_residual(a, &x, rhs);
        if !lu.solve_mut(&mut correction) {
            return Err(singular());
        }
        x += correction;
        original_residual(&x)
    } else {
        initial
    };
    if !residual.is_finite() || residual > tolerance || !x.iter().all(|v| v.is_finite()) {
        return Err(Error::Solve(format!(
            "true linear residual {residual:e} (dense_lu_refinement, initial {initial:e})"
        )));
    }
    let report = LinearReport {
        method: "dense_lu_refinement".into(),
        dimension: a.nrows(),
        stored_nonzeros: a.len(),
        iterations: 1 + usize::from(correct),
        restarts: 0,
        true_relative_residual: residual,
        matrix_storage_bytes: a.len() * size_of::<f64>(),
        workspace_storage_bytes: (a.len() + 3 * a.nrows()) * size_of::<f64>(),
        elapsed_seconds: started.elapsed().as_secs_f64(),
    };
    Ok((x, report))
}

#[cfg(test)]
mod tests {
    use super::super::{Settings, StageCase, Tangent};
    use super::*;
    use serde::Deserialize;
    use std::{fs, path::Path};

    fn case() -> Result<StageCase> {
        Ok(serde_json::from_str(include_str!(
            "../../examples/stages-free-contact.json"
        ))?)
    }

    fn settings() -> Result<Settings> {
        let mut s = case()?.settings;
        s.linear_solver = Default::default();
        Ok(s)
    }

    #[derive(Deserialize)]
    struct SavedSystem {
        matrix: Vec<Vec<f64>>,
        rhs: Vec<f64>,
    }

    #[test]
    fn saved_nc_tangents_pass_original_gate_after_one_correction() -> Result<()> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/linear-refinement");
        for (name, count) in [("nc12.json", 1), ("nc96.jsonl", 7)] {
            let data = fs::read(root.join(name))?;
            let systems = serde_json::Deserializer::from_slice(&data).into_iter::<SavedSystem>();
            let mut seen = 0;
            for saved in systems {
                let saved = saved?;
                seen += 1;
                let rhs = DVector::from_vec(saved.rhs);
                let a = DMatrix::from_fn(rhs.len(), rhs.len(), |i, j| saved.matrix[i][j]);
                let matrix = Tangent::Dense(a);
                let mut s = settings()?;
                assert_eq!(s.relative_linear_tolerance, 1e-10);
                assert!(
                    serde_json::to_value(&s)?
                        .get("dense_lu_refinement")
                        .is_none()
                );
                let Err(Error::Solve(reason)) = matrix.solve(&rhs, &s) else {
                    panic!("{name}: expected the original residual rejection")
                };
                assert!(reason.starts_with("true linear residual "));
                s.dense_lu_refinement = true;
                let s: Settings = serde_json::from_value(serde_json::to_value(&s)?)?;
                let (x, report) = matrix.solve(&rhs, &s)?;
                let actual = (matrix.action(&x)? - &rhs).norm() / rhs.norm();
                println!("{name}: {reason}; corrected {actual:e}");
                assert_eq!(report.true_relative_residual, actual);
                assert!(actual <= s.relative_linear_tolerance);
                assert_eq!(report.iterations, 2);
            }
            assert_eq!(seen, count);
        }
        Ok(())
    }

    #[test]
    fn exact_dense_solve_and_singular_rejection_are_preserved() -> Result<()> {
        let mut s = settings()?;
        let matrix = Tangent::Dense(DMatrix::from_row_slice(2, 2, &[2., 3., 0., 4.]));
        let b = DVector::from_vec(vec![8., 8.]);
        let (baseline, _) = matrix.solve(&b, &s)?;
        s.dense_lu_refinement = true;
        let (x, report) = matrix.solve(&b, &s)?;
        assert_eq!(x, baseline);
        assert_eq!(report.iterations, 1);
        assert!(Tangent::Dense(DMatrix::zeros(2, 2)).solve(&b, &s).is_err());
        Ok(())
    }

    #[test]
    fn sparse_solver_rejects_dense_only_option() -> Result<()> {
        let mut case = case()?;
        case.validate()?;
        case.settings.dense_lu_refinement = true;
        assert!(
            matches!(case.validate(), Err(Error::Input(reason)) if reason.contains("requires dense_lu"))
        );
        Ok(())
    }
}
