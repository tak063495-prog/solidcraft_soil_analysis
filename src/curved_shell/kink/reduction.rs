//! Eliminate the compatibility equations in physical coordinates. Keeping the
//! prescribed columns in each row also recovers the constraint forces there.
use crate::{Error, Result};
use nalgebra::{DMatrix, DVector};
use std::collections::{BTreeMap, BTreeSet};

pub(super) type Row = BTreeMap<usize, f64>;
pub(crate) struct Reduction {
    pivots: Vec<(usize, Row)>,
    columns: Vec<Vec<(usize, f64)>>,
}

impl Reduction {
    pub(super) fn new(rows: Vec<(String, Row)>, free: &[usize]) -> Result<Self> {
        let mut remaining: BTreeSet<_> = free.iter().copied().collect();
        let mut pivots: Vec<(usize, Row)> = Vec::new();
        let total = rows.len();
        let mut dependent = BTreeSet::new();
        for (owner, mut row) in rows {
            for (p, previous) in &pivots {
                let coeff = row.get(p).copied().unwrap_or(0.);
                subtract(&mut row, previous, coeff);
                row.remove(p);
            }
            let pivot = row
                .iter()
                .filter(|(j, _)| remaining.contains(j))
                .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()));
            let Some((&j, &a)) = pivot.filter(|(_, a)| a.abs() > 128. * f64::EPSILON) else {
                dependent.insert(owner);
                continue;
            };
            for v in row.values_mut() {
                *v /= a;
            }
            for (_, previous) in &mut pivots {
                let coeff = previous.get(&j).copied().unwrap_or(0.);
                subtract(previous, &row, coeff);
                previous.remove(&j);
            }
            remaining.remove(&j);
            pivots.push((j, row));
        }
        if !dependent.is_empty() {
            return Err(Error::Input(format!(
                "joint {dependent:?}: constraint rank {}/{} with prescribed or shared DOFs",
                pivots.len(),
                total
            )));
        }
        let columns = remaining
            .into_iter()
            .map(|j| {
                let mut col = vec![(j, 1.)];
                for (p, row) in &pivots {
                    if let Some(a) = row.get(&j) {
                        col.push((*p, -a));
                    }
                }
                col
            })
            .collect();
        Ok(Self { pivots, columns })
    }
    pub(crate) fn enforce(&self, q: &mut DVector<f64>, initial: &[f64]) {
        for (p, row) in &self.pivots {
            q[*p] = initial[*p]
                - row
                    .iter()
                    .filter(|(j, _)| *j != p)
                    .map(|(j, a)| a * (q[*j] - initial[*j]))
                    .sum::<f64>();
        }
    }
    pub(crate) fn solve(
        &self,
        k: &DMatrix<f64>,
        rhs: &DVector<f64>,
        q: &mut DVector<f64>,
    ) -> Result<()> {
        let n = self.columns.len();
        if n == 0 {
            return Ok(());
        }
        let reduced = DMatrix::from_fn(n, n, |i, j| {
            self.columns[i]
                .iter()
                .map(|(a, v)| {
                    v * self.columns[j]
                        .iter()
                        .map(|(b, w)| w * k[(*a, *b)])
                        .sum::<f64>()
                })
                .sum()
        });
        let load = DVector::from_iterator(
            n,
            self.columns
                .iter()
                .map(|col| col.iter().map(|(j, w)| rhs[*j] * w).sum()),
        );
        let x = super::super::path::solve_supported(reduced, &load)?;
        for (col, u) in self.columns.iter().zip(x.iter()) {
            for (j, w) in col {
                q[*j] += w * u;
            }
        }
        Ok(())
    }
    pub(crate) fn internal_transfer(&self, residual: &DVector<f64>) -> DVector<f64> {
        let mut f = DVector::zeros(residual.len());
        for (p, row) in &self.pivots {
            for (j, w) in row {
                f[*j] += w * residual[*p];
            }
        }
        f
    }
}

fn subtract(row: &mut Row, source: &Row, coeff: f64) {
    if coeff != 0. {
        for (j, a) in source {
            *row.entry(*j).or_default() -= a * coeff;
        }
    }
}
