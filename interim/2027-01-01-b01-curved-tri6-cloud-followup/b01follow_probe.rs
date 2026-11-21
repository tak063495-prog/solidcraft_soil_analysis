// Explicit observation executable, unmodified production kernel and guard.
use nalgebra::DVector;
use serde_json::{Value, json};
use soil_fem_3d::plates::{Section, mindlin::Kernel};
use std::io::{self, Read};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let requests: Vec<Value> = serde_json::from_str(&input)?;
    let mut results = Vec::new();
    for request in requests {
        let points: [[f64; 3]; 6] = serde_json::from_value(request["points"].clone())?;
        let section: Section = serde_json::from_value(request["section"].clone())?;
        match Kernel::new(points, &section) {
            Ok(kernel) => {
                let rows: Vec<Vec<f64>> = (0..36)
                    .map(|i| (0..36).map(|j| kernel.stiffness[(i, j)]).collect())
                    .collect();
                let axes: Vec<Vec<f64>> = (0..3)
                    .map(|i| (0..3).map(|j| kernel.axes[(i, j)]).collect())
                    .collect();
                let mut fields = Vec::new();
                for field in request["fields"].as_array().ok_or("fields required")? {
                    let values: Vec<f64> = serde_json::from_value(field["d"].clone())?;
                    if values.len() != 36 {
                        return Err("36 DOFs required".into());
                    }
                    fields.push(json!({"name":field["name"], "recovery":kernel.recover(
                        request["id"].as_str().ok_or("id required")?.into(),
                        &DVector::from_vec(values))}));
                }
                results.push(json!({"id":request["id"], "accepted":true,
                    "stiffness":rows, "axes":axes, "area_m2":kernel.area_m2, "fields":fields}));
            }
            Err(error) => results.push(json!({"id":request["id"], "accepted":false,
                "error":error.to_string()})),
        }
    }
    println!("{}", serde_json::to_string(&results)?);
    Ok(())
}
