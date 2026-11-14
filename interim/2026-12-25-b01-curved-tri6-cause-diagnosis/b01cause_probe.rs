// Standalone internal probe. The generated module exists only under work/.
use nalgebra::DVector;
use serde_json::{Value, json};
pub use soil_fem_3d::{Error, Result, elastic};
use std::io::{self, Read};
include!("b01cause_flat.rs");

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let requests: Vec<Value> = serde_json::from_str(&input)?;
    let mut results = Vec::new();
    for request in requests {
        let points: [[f64; 3]; 6] = serde_json::from_value(request["points"].clone())?;
        let section: soil_fem_3d::plates::Section =
            serde_json::from_value(request["section"].clone())?;
        let production = soil_fem_3d::plates::mindlin::Kernel::new(points, &section);
        let guard = match &production {
            Ok(_) => json!({"accepted":true,"error":null}),
            Err(e) => json!({"accepted":false,"error":e.to_string()}),
        };
        let isolated_section = serde_json::from_value(request["section"].clone())?;
        match b01cause_flat::mindlin::Kernel::new(points, &isolated_section) {
            Ok(kernel) => {
                let k = &kernel.stiffness;
                let rows: Vec<Vec<f64>> = (0..36)
                    .map(|i| (0..36).map(|j| k[(i, j)]).collect())
                    .collect();
                let production_equal = production.as_ref().ok().map(|p| p.stiffness == *k);
                let axes: Vec<Vec<f64>> = (0..3)
                    .map(|i| (0..3).map(|j| kernel.axes[(i, j)]).collect())
                    .collect();
                let mut fields = Vec::new();
                for displacement in request["fields"].as_array().ok_or("fields required")? {
                    let values: Vec<f64> = serde_json::from_value(displacement["d"].clone())?;
                    if values.len() != 36 {
                        return Err("36 DOFs required".into());
                    }
                    let value = kernel.recover("internal-only".into(), &DVector::from_vec(values));
                    fields.push(json!({"name":displacement["name"],"recovery":value}));
                }
                results.push(json!({"id":request["id"],"production_guard":guard,
                    "production_stiffness_exact_when_accepted":production_equal,
                    "stiffness":rows,"axes":axes,"area_m2":kernel.area_m2,
                    "fields":fields}));
            }
            Err(e) => results.push(json!({"id":request["id"],"production_guard":guard,
                "internal_error":e.to_string()})),
        }
    }
    println!("{}", serde_json::to_string(&results)?);
    Ok(())
}
