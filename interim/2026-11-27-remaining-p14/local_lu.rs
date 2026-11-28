// Standalone observation of exactly the SMatrix/SVector LU used by mcc_midpoint.
use nalgebra::{SMatrix, SVector};
use serde::{Deserialize, Serialize};
use std::io;

#[derive(Deserialize)]
struct System {
    matrix: [[f64; 8]; 8],
    b: [f64; 8],
}

#[derive(Serialize)]
struct Observation {
    dx: [f64; 8],
    pivots: [f64; 8],
    permutation: [usize; 8],
}

fn main() {
    let input: System = serde_json::from_reader(io::stdin()).unwrap();
    let matrix = SMatrix::<f64, 8, 8>::from_fn(|i, j| input.matrix[i][j]);
    let lu = matrix.lu();
    let dx = lu.solve(&SVector::from(input.b)).unwrap();
    let mut labels = SVector::<usize, 8>::from_fn(|i, _| i);
    lu.p().permute_rows(&mut labels);
    let observation = Observation {
        dx: dx.into(),
        pivots: lu.u().diagonal().into(),
        permutation: labels.into(),
    };
    serde_json::to_writer(io::stdout(), &observation).unwrap();
}
