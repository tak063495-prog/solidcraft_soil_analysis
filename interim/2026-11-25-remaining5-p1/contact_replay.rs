//! Replay a pair of quadrature-point paths with the production Coulomb update.
use serde_json::{Value, json};
use soil_fem_3d::joint_law::{Material, PointState};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let row: Value = serde_json::from_str(&line?)?;
        let law: Material = serde_json::from_value(row["material"].clone())?;
        let parent: PointState = serde_json::from_value(row["parent"].clone())?;
        let jumps: [[f64; 3]; 3] = serde_json::from_value(row["jumps"].clone())?;
        let n: [usize; 3] = serde_json::from_value(row["counts"].clone())?;
        let coarse = law.integrate_partition(&parent, jumps[0], n[0])?;
        let mid = law.integrate_partition(&parent, jumps[1], n[1])?;
        let fine = law.integrate_partition(&mid, jumps[2], n[2])?;
        println!(
            "{}",
            json!([
                coarse.plastic_slip_m[0],
                mid.plastic_slip_m[0],
                fine.plastic_slip_m[0]
            ])
        );
    }
    Ok(())
}
