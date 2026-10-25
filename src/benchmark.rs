//! Actual Tet4 stiffness assembly with a manufactured discrete solution.
//! This measures the linear FEM backend, not full staged plastic analysis.
use crate::{
    Error, Result, elastic,
    linalg::{self, CsrMatrix, CsrPattern, GmresOptions, LinearReport},
};
use serde::{Deserialize, Serialize};
use std::time::Instant;

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GridCase {
    pub subdivisions: usize,
    pub young_kpa: f64,
    pub poisson: f64,
    pub relative_linear_tolerance: f64,
    pub options: GmresOptions,
    #[serde(default)]
    pub export_system: bool,
}
#[derive(Debug, Serialize)]
pub struct ExportedSystem {
    pub row_ptr: Vec<usize>,
    pub columns: Vec<usize>,
    pub values: Vec<f64>,
    pub rhs: Vec<f64>,
    pub expected_solution: Vec<f64>,
    pub computed_solution: Vec<f64>,
}
#[derive(Debug, Serialize)]
pub struct GridReport {
    pub schema: &'static str,
    pub scope: &'static str,
    pub case: GridCase,
    pub nodes: usize,
    pub tetrahedra: usize,
    pub free_dofs: usize,
    pub symbolic_seconds: f64,
    pub assembly_seconds: f64,
    pub total_seconds: f64,
    pub linear: LinearReport,
    pub relative_solution_error: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub system: Option<ExportedSystem>,
}

pub fn tet4_grid(case: GridCase) -> Result<GridReport> {
    let started = Instant::now();
    let n = case.subdivisions;
    if !(1..=40).contains(&n) || (case.export_system && n > 6) {
        return Err(Error::Input(
            "grid subdivisions must be 1..40; exporting the full system is limited to 6".into(),
        ));
    }
    case.options.validate()?;
    let d = elastic::constitutive(case.young_kpa, case.poisson)?;
    let width = n + 1;
    let index = |x, y, z| x + width * (y + width * z);
    let mut points = Vec::new();
    let mut map = Vec::<[Option<usize>; 3]>::new();
    let mut free = 0;
    let mut expected = Vec::new();
    for z in 0..=n {
        for y in 0..=n {
            for x in 0..=n {
                let p = [
                    x as f64 / n as f64,
                    y as f64 / n as f64,
                    z as f64 / n as f64,
                ];
                points.push(p);
                let fixed = [x == 0 || x == n, y == 0 || y == n, z == 0];
                let mut dofs = [None; 3];
                let values = [
                    p[0] * (1. - p[0]) * (1. + p[1]) * p[2],
                    p[1] * (1. - p[1]) * (1. + p[2]) * p[0],
                    p[2] * (1. + p[0] * p[1]),
                ];
                for a in 0..3 {
                    if !fixed[a] {
                        dofs[a] = Some(free);
                        free += 1;
                        expected.push(values[a]);
                    }
                }
                map.push(dofs);
            }
        }
    }
    let mut elements = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let c = [
                    index(x, y, z),
                    index(x + 1, y, z),
                    index(x + 1, y + 1, z),
                    index(x, y + 1, z),
                    index(x, y, z + 1),
                    index(x + 1, y, z + 1),
                    index(x + 1, y + 1, z + 1),
                    index(x, y + 1, z + 1),
                ];
                for t in [
                    [0, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
                ] {
                    elements.push(t.map(|i| c[i]));
                }
            }
        }
    }
    let symbolic = Instant::now();
    let pattern = CsrPattern::from_connectivity(
        free,
        elements
            .iter()
            .map(|t| t.iter().flat_map(|&i| map[i]).flatten().collect()),
    )?;
    let symbolic_seconds = symbolic.elapsed().as_secs_f64();
    let mut matrix = CsrMatrix::zeros(pattern);
    let assembly = Instant::now();
    for tet in &elements {
        let p = tet.map(|i| points[i]);
        let (_, _, ke) = elastic::element(&p, &d)?;
        let dofs: [Option<usize>; 12] = std::array::from_fn(|a| map[tet[a / 3]][a % 3]);
        for i in 0..12 {
            if let Some(gi) = dofs[i] {
                for j in 0..12 {
                    if let Some(gj) = dofs[j] {
                        matrix.add(gi, gj, ke[(i, j)])?;
                    }
                }
            }
        }
    }
    let assembly_seconds = assembly.elapsed().as_secs_f64();
    let rhs = matrix.matvec(&expected)?;
    let (x, linear) =
        linalg::gmres_ilu0(&matrix, &rhs, case.relative_linear_tolerance, &case.options)?;
    let relative_solution_error = (x
        .iter()
        .zip(&expected)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        / expected.iter().map(|v| v * v).sum::<f64>())
    .sqrt();
    if relative_solution_error > 1e-7 {
        return Err(Error::Solve(format!(
            "manufactured grid solution error {relative_solution_error:e}"
        )));
    }
    let system = case.export_system.then(|| ExportedSystem {
        row_ptr: matrix.row_ptr().to_vec(),
        columns: matrix.columns().to_vec(),
        values: matrix.values().to_vec(),
        rhs,
        expected_solution: expected,
        computed_solution: x,
    });
    Ok(GridReport {
        schema: "soil-fem-tet4-linear-benchmark/1",
        scope: "actual Tet4 linear stiffness assembly and solve; excludes plastic integration, staging, meshing and file IO",
        case,
        nodes: points.len(),
        tetrahedra: elements.len(),
        free_dofs: free,
        symbolic_seconds,
        assembly_seconds,
        total_seconds: started.elapsed().as_secs_f64(),
        linear,
        relative_solution_error,
        system,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn real_tet4_grid_assembly_solves_manufactured_field() -> Result<()> {
        let report = tet4_grid(GridCase {
            subdivisions: 3,
            young_kpa: 1000.,
            poisson: 0.3,
            relative_linear_tolerance: 1e-11,
            options: GmresOptions::default(),
            export_system: true,
        })?;
        assert_eq!(report.tetrahedra, 162);
        assert!(report.relative_solution_error < 1e-9);
        assert!(report.linear.stored_nonzeros < report.free_dofs.pow(2));
        assert!(report.system.is_some());
        Ok(())
    }
}
