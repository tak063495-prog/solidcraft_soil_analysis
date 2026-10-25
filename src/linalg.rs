//! Deterministic CSR assembly and row-equilibrated, right-preconditioned
//! restarted GMRES with ILU(0). No symmetry, diagonal shift or dense fallback.
use serde::{Deserialize, Serialize};
use std::{sync::Arc, time::Instant};

#[derive(Debug, thiserror::Error, Clone, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum LinearError {
    #[error("invalid sparse system: {message}")]
    Invalid { message: String },
    #[error("ILU(0) pivot failure at row {row}: scaled pivot={pivot:e}")]
    Pivot { row: usize, pivot: f64 },
    #[error("GMRES Arnoldi breakdown after {iterations} iterations; true residual={residual:e}")]
    Breakdown { iterations: usize, residual: f64 },
    #[error("GMRES iteration limit {iterations}; true residual={residual:e}")]
    IterationLimit { iterations: usize, residual: f64 },
    #[error("linear solve cancelled after {iterations} iterations")]
    Cancelled { iterations: usize },
}
type LinearResult<T> = std::result::Result<T, LinearError>;
fn invalid(message: &str) -> LinearError {
    LinearError::Invalid {
        message: message.into(),
    }
}

#[derive(Debug, Clone)]
pub struct CsrPattern {
    n: usize,
    row_ptr: Vec<usize>,
    columns: Vec<usize>,
    diagonal: Vec<usize>,
}
impl CsrPattern {
    /// Full element connectivity, including constrained/tied DOF reduction.
    /// Duplicate DOFs are allowed here; numerical contributions are still added.
    pub fn from_connectivity(
        n: usize,
        elements: impl IntoIterator<Item = Vec<usize>>,
    ) -> LinearResult<Arc<Self>> {
        if n == 0 || n > 1_000_000 {
            return Err(invalid("dimension must be in 1..=1000000"));
        }
        let mut rows: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
        let mut contributions = 0usize;
        for mut indices in elements {
            contributions = contributions
                .checked_add(indices.len().saturating_mul(indices.len()))
                .ok_or_else(|| invalid("symbolic assembly budget overflow"))?;
            if indices.len() > 256 || contributions > 200_000_000 {
                return Err(invalid("symbolic connectivity/assembly budget exceeded"));
            }
            if indices.iter().any(|&i| i >= n) {
                return Err(invalid("connectivity references an unknown DOF"));
            }
            indices.sort_unstable();
            indices.dedup();
            for &i in &indices {
                rows[i].extend(&indices);
            }
        }
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut diagonal = Vec::with_capacity(n);
        for (i, row) in rows.iter_mut().enumerate() {
            row.sort_unstable();
            row.dedup();
            let local = row
                .binary_search(&i)
                .map_err(|_| invalid("missing diagonal"))?;
            diagonal.push(columns.len() + local);
            columns.extend(row.iter().copied());
            if columns.len() > 100_000_000 {
                return Err(invalid("CSR nonzero budget exceeded"));
            }
            row_ptr.push(columns.len());
        }
        Ok(Arc::new(Self {
            n,
            row_ptr,
            columns,
            diagonal,
        }))
    }
    pub fn dimension(&self) -> usize {
        self.n
    }
    pub fn nonzeros(&self) -> usize {
        self.columns.len()
    }
    pub fn storage_bytes(&self) -> usize {
        (self.row_ptr.len() + self.columns.len() + self.diagonal.len()) * size_of::<usize>()
    }
}

#[derive(Debug, Clone)]
pub struct CsrMatrix {
    pattern: Arc<CsrPattern>,
    values: Vec<f64>,
}
impl CsrMatrix {
    pub fn zeros(pattern: Arc<CsrPattern>) -> Self {
        let values = vec![0.; pattern.nonzeros()];
        Self { pattern, values }
    }
    pub fn pattern(&self) -> &CsrPattern {
        &self.pattern
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn row_ptr(&self) -> &[usize] {
        &self.pattern.row_ptr
    }
    pub fn columns(&self) -> &[usize] {
        &self.pattern.columns
    }
    pub fn add(&mut self, row: usize, col: usize, value: f64) -> LinearResult<()> {
        if row >= self.pattern.n || col >= self.pattern.n || !value.is_finite() {
            return Err(invalid("out-of-range/non-finite assembly contribution"));
        }
        let start = self.pattern.row_ptr[row];
        let end = self.pattern.row_ptr[row + 1];
        let index = self.pattern.columns[start..end]
            .binary_search(&col)
            .map_err(|_| invalid("assembly contribution outside symbolic pattern"))?
            + start;
        let next = self.values[index] + value;
        if !next.is_finite() {
            return Err(invalid("non-finite assembled value"));
        }
        self.values[index] = next;
        Ok(())
    }
    pub fn matvec(&self, x: &[f64]) -> LinearResult<Vec<f64>> {
        if x.len() != self.pattern.n || !x.iter().all(|x| x.is_finite()) {
            return Err(invalid("invalid matvec input"));
        }
        let mut y = vec![0.; x.len()];
        self.multiply_into(x, &mut y);
        if !y.iter().all(|x| x.is_finite()) {
            return Err(invalid("matvec overflow"));
        }
        Ok(y)
    }
    fn multiply_into(&self, x: &[f64], out: &mut [f64]) {
        for (i, y) in out.iter_mut().enumerate() {
            *y = (self.pattern.row_ptr[i]..self.pattern.row_ptr[i + 1])
                .map(|p| self.values[p] * x[self.pattern.columns[p]])
                .sum();
        }
    }
    pub fn storage_bytes(&self) -> usize {
        self.pattern.storage_bytes() + self.values.len() * size_of::<f64>()
    }
    /// Temporary zero-increment coordinates, preserving the sparse pattern.
    /// The caller must still check the original force equations.
    pub(crate) fn zero_increment_gauges(
        &self,
        indices: &[usize],
        diagonal: f64,
    ) -> LinearResult<Self> {
        if !diagonal.is_finite() || diagonal <= 0. || indices.iter().any(|&i| i >= self.pattern.n) {
            return Err(invalid("invalid zero-increment gauge"));
        }
        let mut pinned = vec![false; self.pattern.n];
        for &i in indices {
            if pinned[i] {
                return Err(invalid("duplicate zero-increment gauge"));
            }
            pinned[i] = true;
        }
        let mut result = self.clone();
        for i in 0..self.pattern.n {
            for p in self.pattern.row_ptr[i]..self.pattern.row_ptr[i + 1] {
                if pinned[i] || pinned[self.pattern.columns[p]] {
                    result.values[p] = 0.;
                }
            }
            if pinned[i] {
                result.values[self.pattern.diagonal[i]] = diagonal;
            }
        }
        Ok(result)
    }
    pub fn true_residual(&self, x: &[f64], b: &[f64]) -> LinearResult<f64> {
        if b.len() != self.pattern.n || !b.iter().all(|v| v.is_finite()) {
            return Err(invalid("invalid right-hand side"));
        }
        let ax = self.matvec(x)?;
        let residual: Vec<_> = b.iter().zip(ax).map(|(b, a)| b - a).collect();
        Ok(norm(&residual) / norm(b).max(f64::MIN_POSITIVE))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GmresOptions {
    pub restart: usize,
    pub maximum_iterations: usize,
    pub pivot_tolerance: f64,
}
impl Default for GmresOptions {
    fn default() -> Self {
        Self {
            restart: 60,
            maximum_iterations: 3000,
            pivot_tolerance: 1e-14,
        }
    }
}
impl GmresOptions {
    pub fn validate(&self) -> LinearResult<()> {
        if self.restart == 0
            || self.restart > 500
            || self.maximum_iterations == 0
            || self.maximum_iterations > 1_000_000
            || !self.pivot_tolerance.is_finite()
            || !(1e-16..=1e-6).contains(&self.pivot_tolerance)
        {
            return Err(invalid(
                "invalid restart, iteration limit or ILU pivot threshold",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(tag = "method", rename_all = "snake_case", deny_unknown_fields)]
pub enum LinearSolver {
    #[default]
    DenseLu,
    SparseGmresIlu0 {
        options: GmresOptions,
    },
}
impl LinearSolver {
    pub fn is_dense(&self) -> bool {
        matches!(self, Self::DenseLu)
    }
    pub fn validate(&self) -> LinearResult<()> {
        if let Self::SparseGmresIlu0 { options } = self {
            options.validate()?;
        }
        Ok(())
    }
    pub fn dimension_limit(&self) -> usize {
        if self.is_dense() { 2000 } else { 300_000 }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LinearReport {
    pub method: String,
    pub dimension: usize,
    pub stored_nonzeros: usize,
    pub iterations: usize,
    pub restarts: usize,
    pub true_relative_residual: f64,
    pub matrix_storage_bytes: usize,
    /// Explicitly allocated solve arrays, excluding matrix and allocator overhead.
    pub workspace_storage_bytes: usize,
    pub elapsed_seconds: f64,
}

struct Ilu0<'a> {
    pattern: &'a CsrPattern,
    factors: Vec<f64>,
}
impl<'a> Ilu0<'a> {
    fn factor(matrix: &'a CsrMatrix, threshold: f64) -> LinearResult<Self> {
        let p = &matrix.pattern;
        let mut values = matrix.values.clone();
        for i in 0..p.n {
            for ij in p.row_ptr[i]..p.diagonal[i] {
                let j = p.columns[ij];
                let pivot = values[p.diagonal[j]];
                if !pivot.is_finite() || pivot.abs() < threshold {
                    return Err(LinearError::Pivot { row: j, pivot });
                }
                let multiplier = values[ij] / pivot;
                values[ij] = multiplier;
                for jk in p.diagonal[j] + 1..p.row_ptr[j + 1] {
                    let col = p.columns[jk];
                    if let Ok(index) = p.columns[ij + 1..p.row_ptr[i + 1]].binary_search(&col) {
                        let ik = ij + 1 + index;
                        values[ik] -= multiplier * values[jk];
                    }
                }
            }
            let pivot = values[p.diagonal[i]];
            if !pivot.is_finite() || pivot.abs() < threshold {
                return Err(LinearError::Pivot { row: i, pivot });
            }
        }
        if !values.iter().all(|v| v.is_finite()) {
            return Err(invalid("ILU factor overflow"));
        }
        Ok(Self {
            pattern: p,
            factors: values,
        })
    }
    fn apply(&self, rhs: &[f64]) -> LinearResult<Vec<f64>> {
        let p = self.pattern;
        let mut x = rhs.to_vec();
        for i in 0..p.n {
            for ij in p.row_ptr[i]..p.diagonal[i] {
                x[i] -= self.factors[ij] * x[p.columns[ij]];
            }
        }
        for i in (0..p.n).rev() {
            for ij in p.diagonal[i] + 1..p.row_ptr[i + 1] {
                x[i] -= self.factors[ij] * x[p.columns[ij]];
            }
            x[i] /= self.factors[p.diagonal[i]];
        }
        if !x.iter().all(|v| v.is_finite()) {
            return Err(invalid("ILU triangular solve overflow"));
        }
        Ok(x)
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(a: &[f64]) -> f64 {
    let scale = a.iter().map(|x| x.abs()).fold(0., f64::max);
    if scale == 0. {
        0.
    } else {
        scale * a.iter().map(|x| (x / scale).powi(2)).sum::<f64>().sqrt()
    }
}
fn candidate(
    base: &[f64],
    z: &[Vec<f64>],
    h: &[Vec<f64>],
    g: &[f64],
    m: usize,
) -> LinearResult<Vec<f64>> {
    let mut y = g[..m].to_vec();
    for i in (0..m).rev() {
        for j in i + 1..m {
            y[i] -= h[i][j] * y[j];
        }
        if h[i][i].abs() < f64::MIN_POSITIVE {
            return Err(invalid("singular GMRES Hessenberg factor"));
        }
        y[i] /= h[i][i];
    }
    let mut x = base.to_vec();
    for j in 0..m {
        for (i, value) in x.iter_mut().enumerate() {
            *value += z[j][i] * y[j];
        }
    }
    Ok(x)
}

pub fn gmres_ilu0(
    matrix: &CsrMatrix,
    rhs: &[f64],
    tolerance: f64,
    options: &GmresOptions,
) -> LinearResult<(Vec<f64>, LinearReport)> {
    gmres_ilu0_with_cancel(matrix, rhs, tolerance, options, || false)
}

pub fn gmres_ilu0_with_cancel(
    matrix: &CsrMatrix,
    rhs: &[f64],
    tolerance: f64,
    options: &GmresOptions,
    cancelled: impl Fn() -> bool,
) -> LinearResult<(Vec<f64>, LinearReport)> {
    let started = Instant::now();
    options.validate()?;
    let n = matrix.pattern.n;
    if rhs.len() != n
        || !rhs.iter().chain(&matrix.values).all(|v| v.is_finite())
        || !tolerance.is_finite()
        || !(1e-15..=1e-2).contains(&tolerance)
    {
        return Err(invalid(
            "non-finite/dimension mismatch/invalid relative tolerance",
        ));
    }
    if cancelled() {
        return Err(LinearError::Cancelled { iterations: 0 });
    }
    let mut x = vec![0.; n];
    let bnorm = norm(rhs);
    if !bnorm.is_finite() {
        return Err(invalid("right-hand side norm overflow"));
    }
    let restart = options.restart.min(n).min(options.maximum_iterations);
    let report = |iterations, restarts, residual| LinearReport {
        method: "csr_gmres_ilu0".into(),
        dimension: n,
        stored_nonzeros: matrix.values.len(),
        iterations,
        restarts,
        true_relative_residual: residual,
        matrix_storage_bytes: matrix.storage_bytes(),
        workspace_storage_bytes: (2 * matrix.values.len()
            + (2 * restart + 12) * n
            + (restart + 1) * restart
            + 4 * restart
            + 1)
            * size_of::<f64>(),
        elapsed_seconds: started.elapsed().as_secs_f64(),
    };
    if bnorm == 0. {
        return Ok((x, report(0, 0, 0.)));
    }
    // Row equilibration changes only the preconditioned system. Acceptance is
    // always checked using the original A and b, with no diagonal regularization.
    let mut scaled = matrix.clone();
    let mut b = rhs.to_vec();
    for (i, bi) in b.iter_mut().enumerate() {
        let range = matrix.pattern.row_ptr[i]..matrix.pattern.row_ptr[i + 1];
        let maximum = matrix.values[range.clone()]
            .iter()
            .map(|v| v.abs())
            .fold(0., f64::max);
        if maximum == 0. {
            return Err(LinearError::Pivot { row: i, pivot: 0. });
        }
        for j in range {
            scaled.values[j] /= maximum;
        }
        *bi /= maximum;
    }
    if !b.iter().all(|x| x.is_finite()) {
        return Err(invalid("row scaling overflow"));
    }
    let ilu = Ilu0::factor(&scaled, options.pivot_tolerance)?;
    let mut iterations = 0;
    let mut cycles = 0;
    loop {
        if cancelled() {
            return Err(LinearError::Cancelled { iterations });
        }
        let mut r = vec![0.; n];
        scaled.multiply_into(&x, &mut r);
        for i in 0..n {
            r[i] = b[i] - r[i];
        }
        let beta = norm(&r);
        let residual = matrix.true_residual(&x, rhs)?;
        if residual <= tolerance {
            return Ok((x, report(iterations, cycles, residual)));
        }
        if iterations >= options.maximum_iterations {
            return Err(LinearError::IterationLimit {
                iterations,
                residual,
            });
        }
        if beta == 0. || !beta.is_finite() {
            return Err(LinearError::Breakdown {
                iterations,
                residual,
            });
        }
        let mut v = vec![r.iter().map(|r| r / beta).collect::<Vec<_>>()];
        let mut z = Vec::with_capacity(restart);
        let mut h = vec![vec![0.; restart]; restart + 1];
        let mut cs = vec![0.; restart];
        let mut sn = vec![0.; restart];
        let mut g = vec![0.; restart + 1];
        g[0] = beta;
        for j in 0..restart {
            if cancelled() {
                return Err(LinearError::Cancelled { iterations });
            }
            z.push(ilu.apply(&v[j])?);
            let mut w = vec![0.; n];
            scaled.multiply_into(&z[j], &mut w);
            let wnorm = norm(&w);
            // Twice modified Gram-Schmidt avoids accepting an optimistic
            // recursive residual when Arnoldi vectors lose orthogonality.
            for _ in 0..2 {
                for i in 0..=j {
                    let projection = dot(&v[i], &w);
                    h[i][j] += projection;
                    for (wi, vi) in w.iter_mut().zip(&v[i]) {
                        *wi -= projection * vi;
                    }
                }
            }
            let next_norm = norm(&w);
            h[j + 1][j] = next_norm;
            let breakdown = next_norm <= 1e-14 * wnorm.max(f64::MIN_POSITIVE);
            if !breakdown {
                v.push(w.iter().map(|w| w / next_norm).collect());
            }
            for i in 0..j {
                let a = h[i][j];
                let b = h[i + 1][j];
                h[i][j] = cs[i] * a + sn[i] * b;
                h[i + 1][j] = -sn[i] * a + cs[i] * b;
            }
            let rho = h[j][j].hypot(h[j + 1][j]);
            if rho == 0. || !rho.is_finite() {
                return Err(LinearError::Breakdown {
                    iterations,
                    residual,
                });
            }
            cs[j] = h[j][j] / rho;
            sn[j] = h[j + 1][j] / rho;
            h[j][j] = rho;
            h[j + 1][j] = 0.;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            iterations += 1;
            let end_cycle = j + 1 == restart || iterations == options.maximum_iterations;
            // The Hessenberg estimate is only an early-check trigger.
            if g[j + 1].abs() <= tolerance * norm(&b)
                || end_cycle
                || breakdown
                || iterations % 10 == 0
            {
                let next = candidate(&x, &z, &h, &g, j + 1)?;
                let true_residual = matrix.true_residual(&next, rhs)?;
                if true_residual <= tolerance {
                    return Ok((next, report(iterations, cycles, true_residual)));
                }
                if breakdown {
                    return Err(LinearError::Breakdown {
                        iterations,
                        residual: true_residual,
                    });
                }
                if end_cycle {
                    x = next;
                    break;
                }
            }
        }
        cycles += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    fn from_dense(a: &DMatrix<f64>) -> LinearResult<CsrMatrix> {
        let n = a.nrows();
        let connections = (0..n)
            .flat_map(|i| (0..n).filter_map(move |j| (a[(i, j)] != 0.).then_some(vec![i, j])));
        let mut matrix = CsrMatrix::zeros(CsrPattern::from_connectivity(n, connections)?);
        for i in 0..n {
            for j in 0..n {
                if a[(i, j)] != 0. {
                    matrix.add(i, j, a[(i, j)])?;
                }
            }
        }
        Ok(matrix)
    }

    #[test]
    fn duplicate_assembly_preserves_unsymmetric_entries_and_matches_dense() -> LinearResult<()> {
        let pattern = CsrPattern::from_connectivity(3, [vec![0, 1], vec![1, 2], vec![0, 1, 1]])?;
        let mut a = CsrMatrix::zeros(pattern);
        for (i, j, v) in [
            (0, 0, 4.),
            (0, 1, 3.),
            (1, 0, -2.),
            (1, 1, 3.),
            (1, 2, 1.),
            (2, 1, -1.),
            (2, 2, 2.),
            (0, 0, 1.),
        ] {
            a.add(i, j, v)?;
        }
        assert_eq!(a.matvec(&[1., 2., 3.])?, vec![11., 7., 4.]);
        let expected = vec![1., -2., 0.5];
        let b = a.matvec(&expected)?;
        let (x, report) = gmres_ilu0(&a, &b, 1e-12, &GmresOptions::default())?;
        assert!(x.iter().zip(expected).all(|(x, y)| (x - y).abs() < 1e-12));
        assert!(report.true_relative_residual < 1e-12);
        assert!(a.add(0, 2, 1.).is_err());
        Ok(())
    }

    #[test]
    fn restarted_nonnormal_indefinite_system_matches_independent_pivoted_lu() -> LinearResult<()> {
        let n = 80;
        let mut dense = DMatrix::zeros(n, n);
        for i in 0..n {
            dense[(i, i)] = if i % 2 == 0 { 4. } else { -4. };
            if i + 1 < n {
                dense[(i, i + 1)] = 1.7;
                dense[(i + 1, i)] = -0.9;
            }
            if i + 9 < n {
                dense[(i, i + 9)] = 0.6;
                dense[(i + 9, i)] = -0.4;
            }
        }
        let a = from_dense(&dense)?;
        let b = DVector::from_fn(n, |i, _| ((i * 17) as f64).sin());
        let exact = dense
            .lu()
            .solve(&b)
            .ok_or_else(|| invalid("dense reference singular"))?;
        let options = GmresOptions {
            restart: 3,
            maximum_iterations: 500,
            ..Default::default()
        };
        let (x, report) = gmres_ilu0(&a, b.as_slice(), 1e-11, &options)?;
        assert!(report.restarts > 0);
        assert!((DVector::from_vec(x) - exact).norm() / b.norm() < 1e-10);
        assert!(report.true_relative_residual < 1e-11);
        Ok(())
    }

    #[test]
    fn row_scaling_does_not_replace_original_residual_acceptance() -> LinearResult<()> {
        let mut dense = DMatrix::zeros(60, 60);
        for i in 0..60 {
            let scale = 10f64.powi((i % 17) as i32 - 8);
            dense[(i, i)] = 5. * scale;
            for offset in [1, 7] {
                if i + offset < 60 {
                    dense[(i, i + offset)] = -scale;
                }
            }
            for offset in [1, 7] {
                if i >= offset {
                    dense[(i, i - offset)] = 0.8 * scale;
                }
            }
        }
        let a = from_dense(&dense)?;
        let expected: Vec<_> = (0..60).map(|i| (i as f64).cos()).collect();
        let b = a.matvec(&expected)?;
        let (x, report) = gmres_ilu0(&a, &b, 1e-12, &GmresOptions::default())?;
        assert!(a.true_residual(&x, &b)? < 1e-12);
        assert!(report.true_relative_residual < 1e-12);
        assert!(x.iter().zip(expected).all(|(a, b)| (a - b).abs() < 1e-9));
        Ok(())
    }

    #[test]
    fn invalid_pivot_limit_and_cancellation_have_distinct_failures() -> LinearResult<()> {
        assert!(CsrPattern::from_connectivity(0, []).is_err());
        assert!(CsrPattern::from_connectivity(2, [vec![2]]).is_err());
        let singular = from_dense(&DMatrix::from_row_slice(2, 2, &[1., 1., 1., 1.]))?;
        assert!(matches!(
            gmres_ilu0(&singular, &[1., 2.], 1e-10, &GmresOptions::default()),
            Err(LinearError::Pivot { .. })
        ));
        let dense = DMatrix::from_fn(20, 20, |i, j| {
            if i == j {
                4.
            } else if (i + j) % 7 == 0 {
                0.5
            } else {
                0.
            }
        });
        let a = from_dense(&dense)?;
        let b: Vec<_> = (0..20).map(|i| i as f64).collect();
        let options = GmresOptions {
            restart: 1,
            maximum_iterations: 1,
            ..Default::default()
        };
        assert!(matches!(
            gmres_ilu0(&a, &b, 1e-14, &options),
            Err(LinearError::IterationLimit { .. })
        ));
        assert!(matches!(
            gmres_ilu0_with_cancel(&a, &b, 1e-10, &options, || true),
            Err(LinearError::Cancelled { .. })
        ));
        assert!(a.matvec(&[1.]).is_err());
        let (x, r) = gmres_ilu0(&a, &[0.; 20], 1e-12, &options)?;
        assert_eq!(x, vec![0.; 20]);
        assert_eq!(r.iterations, 0);
        Ok(())
    }
}
