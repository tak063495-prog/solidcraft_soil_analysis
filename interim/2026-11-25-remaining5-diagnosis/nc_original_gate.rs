//! Read saved candidate vectors and evaluate the original nalgebra binary64 product.
use nalgebra::{DMatrix, DVector};
use serde::Deserialize;
use std::{collections::BTreeMap, env, error::Error, fs};

#[derive(Deserialize)]
struct Record {
    a_column_major: Vec<f64>,
    b: Vec<f64>,
    relative: f64,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    let record: Record = serde_json::from_slice(&fs::read(&args[0])?)?;
    let vectors: BTreeMap<String, Vec<f64>> = serde_json::from_slice(&fs::read(&args[1])?)?;
    let b = DVector::from_vec(record.b);
    let a = DMatrix::from_column_slice(b.len(), b.len(), &record.a_column_major);
    let mut rows = BTreeMap::new();
    for (name, values) in vectors {
        let x = DVector::from_vec(values);
        let rel = (&a * x - &b).norm() / b.norm().max(1e-30);
        rows.insert(name, rel);
    }
    assert_eq!(rows["saved"], record.relative);
    println!("{}", serde_json::to_string_pretty(&rows)?);
    Ok(())
}
