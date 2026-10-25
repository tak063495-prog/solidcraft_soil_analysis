//! Keep independent GP volume first, then use deterministic RCM on the u/p
//! graph. This is a symmetric permutation, not a change to the equations or
//! a regularization. The original unpermuted residual remains mandatory.
use super::*;
use std::collections::VecDeque;

pub(super) struct Ordering {
    new_to_old: Vec<usize>,
    old_to_new: Vec<usize>,
    pattern: Arc<CsrPattern>,
}

impl Ordering {
    pub(super) fn new(pattern: Arc<CsrPattern>, volume_count: usize) -> Result<Self> {
        let original = CsrMatrix::zeros(pattern);
        let n = original.pattern().dimension();
        if volume_count > n {
            return Err(Error::Input(
                "mixed ordering volume range exceeds dimension".into(),
            ));
        }
        let mut graph = vec![Vec::new(); n - volume_count];
        for (i, neighbors) in graph.iter_mut().enumerate() {
            let row = i + volume_count;
            for &column in &original.columns()[original.row_ptr()[row]..original.row_ptr()[row + 1]]
            {
                if column >= volume_count && column != row {
                    neighbors.push(column - volume_count);
                }
            }
        }
        let degree: Vec<_> = graph.iter().map(Vec::len).collect();
        let mut visited = vec![false; graph.len()];
        let mut new_to_old: Vec<_> = (0..volume_count).collect();
        while let Some(root) = (0..graph.len())
            .filter(|&i| !visited[i])
            .min_by_key(|&i| (degree[i], i))
        {
            let mut queue = VecDeque::from([root]);
            let mut component = Vec::new();
            visited[root] = true;
            while let Some(node) = queue.pop_front() {
                component.push(node + volume_count);
                let mut neighbors: Vec<_> = graph[node]
                    .iter()
                    .copied()
                    .filter(|&i| !visited[i])
                    .collect();
                neighbors.sort_unstable_by_key(|&i| (degree[i], i));
                for next in neighbors {
                    if !visited[next] {
                        visited[next] = true;
                        queue.push_back(next);
                    }
                }
            }
            new_to_old.extend(component.into_iter().rev());
        }
        let mut old_to_new = vec![0; n];
        for (new, &old) in new_to_old.iter().enumerate() {
            old_to_new[old] = new;
        }
        // The FEM pattern is symmetric even when its numerical entries are
        // unsymmetric. Pair connectivity preserves precisely this pattern.
        let pairs = (0..n).flat_map(|i| {
            let mapping = &old_to_new;
            original.columns()[original.row_ptr()[i]..original.row_ptr()[i + 1]]
                .iter()
                .filter(move |&&j| j >= i)
                .map(move |&j| vec![mapping[i], mapping[j]])
        });
        let pattern = CsrPattern::from_connectivity(n, pairs)?;
        if pattern.nonzeros() != original.pattern().nonzeros() {
            return Err(Error::Input(
                "mixed ordering requires a symmetric symbolic pattern".into(),
            ));
        }
        Ok(Self {
            new_to_old,
            old_to_new,
            pattern,
        })
    }

    pub(super) fn solve(
        &self,
        matrix: &Tangent,
        rhs: &DVector<f64>,
        settings: &Settings,
    ) -> Result<(DVector<f64>, LinearReport)> {
        let Tangent::Sparse(original) = matrix else {
            return matrix.solve(rhs, settings);
        };
        let LinearSolver::SparseGmresIlu0 { options } = &settings.linear_solver else {
            return Err(Error::Input("mixed sparse ordering/solver mismatch".into()));
        };
        let mut reordered = CsrMatrix::zeros(self.pattern.clone());
        for i in 0..rhs.len() {
            for k in original.row_ptr()[i]..original.row_ptr()[i + 1] {
                reordered.add(
                    self.old_to_new[i],
                    self.old_to_new[original.columns()[k]],
                    original.values()[k],
                )?;
            }
        }
        let b: Vec<_> = self.new_to_old.iter().map(|&i| rhs[i]).collect();
        let (ordered, mut report) =
            linalg::gmres_ilu0(&reordered, &b, settings.relative_linear_tolerance, options)?;
        let mut x = DVector::zeros(rhs.len());
        for (new, &old) in self.new_to_old.iter().enumerate() {
            x[old] = ordered[new];
        }
        let residual = original.true_residual(x.as_slice(), rhs.as_slice())?;
        if residual > settings.relative_linear_tolerance || !residual.is_finite() {
            return Err(Error::Solve(format!(
                "mixed unpermuted true linear residual {residual:e}"
            )));
        }
        report.method = "csr_gmres_ilu0_volume_first_rcm_up".into();
        report.true_relative_residual = residual;
        // Include both matrices and the two permutation vectors in this
        // linear-solve accounting; the full FEM/material storage is separate.
        report.matrix_storage_bytes += original.storage_bytes();
        report.workspace_storage_bytes += 2 * rhs.len() * size_of::<usize>();
        Ok((x, report))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn volume_first_rcm_preserves_unsymmetric_saddle_equations_and_matches_dense() -> Result<()> {
        let nt = 5;
        let n = 17;
        let mut rows = Vec::new();
        for i in 0..6 {
            rows.push(vec![nt + i, nt + 6 + i]);
            if i + 1 < 6 {
                rows.push(vec![nt + i, nt + i + 1, nt + 6 + i, nt + 6 + i + 1]);
            }
        }
        for i in 0..nt {
            rows.push(vec![i, nt + i, nt + 6 + i]);
        }
        let pattern = CsrPattern::from_connectivity(n, rows)?;
        let ordering = Ordering::new(pattern.clone(), nt)?;
        assert_eq!(&ordering.new_to_old[..nt], &[0, 1, 2, 3, 4]);
        let mut sparse = CsrMatrix::zeros(pattern.clone());
        let mut dense = DMatrix::zeros(n, n);
        for i in 0..n {
            let cols = sparse.columns()[sparse.row_ptr()[i]..sparse.row_ptr()[i + 1]].to_vec();
            for j in cols {
                let value = if i == j {
                    if i < nt {
                        1200.
                    } else if i < nt + 6 {
                        8.
                    } else {
                        -0.2
                    }
                } else {
                    0.1 * (1. + (i + 2 * j) as f64).sin()
                };
                sparse.add(i, j, value)?;
                dense[(i, j)] = value;
            }
        }
        let rhs = DVector::from_fn(n, |i, _| (i as f64).cos());
        let case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-mixed-mcc.json"))?;
        let mut settings = case.settings;
        settings.linear_solver = LinearSolver::SparseGmresIlu0 {
            options: Default::default(),
        };
        let expected = dense.lu().solve(&rhs).unwrap();
        let (x, report) = ordering.solve(&Tangent::Sparse(sparse), &rhs, &settings)?;
        assert!((&x - expected).norm() < 1e-8);
        assert!(report.true_relative_residual <= settings.relative_linear_tolerance);
        assert_eq!(report.method, "csr_gmres_ilu0_volume_first_rcm_up");
        let again = Ordering::new(pattern, nt)?;
        assert_eq!(again.new_to_old, ordering.new_to_old);
        Ok(())
    }
}
