// Explicit standalone observation. No guard transplant, no production edits.
use nalgebra::DVector;
use serde_json::{Value, json};
use soil_fem_3d::{beams, mesh::Mesh, plates};
use std::io::{self, Read};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut text = String::new();
    io::stdin().read_to_string(&mut text)?;
    let requests: Vec<Value> = serde_json::from_str(&text)?;
    let mut rows = Vec::new();
    for r in requests {
        if r["kind"] == "model_guard" {
            let mesh: Mesh = serde_json::from_value(r["mesh"].clone())?;
            let config: beams::Configuration = serde_json::from_value(r["configuration"].clone())?;
            let result = match beams::Model::new(&mesh, config) {
                Ok(_) => json!({"accepted": true, "error": null}),
                Err(e) => json!({"accepted": false, "error": e.to_string()}),
            };
            rows.push(json!({"id": r["id"], "production_model_guard": result}));
            continue;
        }
        let points: [[f64; 3]; 6] = serde_json::from_value(r["points"].clone())?;
        let section: plates::Section = serde_json::from_value(r["section"].clone())?;
        let e = r["e"].as_f64().ok_or("finite e required")?;
        match plates::offset::Kernel::new(points, &section, e) {
            Ok(kernel) => {
                let matrix: Vec<Vec<f64>> = (0..36)
                    .map(|i| (0..36).map(|j| kernel.stiffness[(i, j)]).collect())
                    .collect();
                let neutral: Vec<Vec<f64>> = (0..36)
                    .map(|i| (0..36).map(|j| kernel.neutral.stiffness[(i, j)]).collect())
                    .collect();
                let traction: [f64; 3] = serde_json::from_value(r["traction"].clone())?;
                let load = kernel.uniform_load(traction);
                let mut fields = Vec::new();
                for field in r["fields"].as_array().ok_or("fields required")? {
                    let values: Vec<f64> = serde_json::from_value(field["d"].clone())?;
                    if values.len() != 36 {
                        return Err("36 DOFs required".into());
                    }
                    let d = DVector::from_vec(values);
                    fields.push(json!({"name":field["name"],
                        "neutral_d":kernel.neutral_kinematics(&d).as_slice(),
                        "recovery":kernel.recover("b03cause-internal".into(), &d)}));
                }
                rows.push(json!({"id":r["id"], "stiffness":matrix,
                    "neutral_stiffness":neutral, "load":load.as_slice(),
                    "fields":fields, "error":null}));
            }
            Err(e) => rows.push(json!({"id":r["id"], "error":e.to_string()})),
        }
    }
    println!("{}", serde_json::to_string(&rows)?);
    Ok(())
}
