//! Opt-in public-kernel cone/birth probe. Rejection is retained, never repaired.
use soil_fem_3d::{Result, joint_law::Material, joints, mesh::Mesh, stages::StageCase};
use std::io::{self, BufRead, Write};
fn request(v: serde_json::Value) -> Result<serde_json::Value> {
    if v["kind"] == "local" {
        let mat: Material = serde_json::from_value(v["material"].clone())?;
        let normal = serde_json::from_value(v["normal"].clone())?;
        let q = serde_json::from_value(v["traction"].clone())?;
        return Ok(match mat.birth(normal, [0.; 3], q) {
            Ok(old) => {
                let target = serde_json::from_value(v["target"].clone())?;
                match mat.integrate_partition(&old, target, 4) {
                    Ok(state) => serde_json::json!({"accepted":true,"initial":old,"state":state}),
                    Err(e) => {
                        serde_json::json!({"accepted":true,"initial":old,"integration_error":e.to_string()})
                    }
                }
            }
            Err(e) => serde_json::json!({"accepted":false,"error":e.to_string(),"work":null}),
        });
    }
    let mesh: Mesh = serde_json::from_value(v["mesh"].clone())?;
    assert_eq!(mesh.tetrahedra.len(), 24);
    let case: StageCase = serde_json::from_value(v["case"].clone())?;
    let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
    let state = model.initial_state()?;
    let u = vec![[0.; 3]; mesh.nodes.len()];
    let active = vec![true; mesh.tetrahedra.len()];
    let rows: Vec<_> = (0..3 * u.len()).map(|g| vec![(g, 1.)]).collect();
    let force: Vec<f64> = serde_json::from_value(v["residual_kn"].clone())?;
    Ok(
        match model.birth(&state, &u, &active, &["Seam".into()], &force, &rows) {
            Ok(born) => {
                let mut target = u.clone();
                for x in &mut target {
                    x[0] += 1e-6;
                }
                let evaluated = model.evaluate(
                    &born,
                    &target,
                    &active,
                    Some(&vec![4; model.point_count()]),
                    false,
                )?;
                serde_json::json!({"accepted":true,"initial":born,"state":evaluated.state,"layout":model.layout()})
            }
            Err(e) => serde_json::json!({"accepted":false,"error":e.to_string(),"work":null}),
        },
    )
}
fn main() {
    let mut out = io::stdout().lock();
    for line in io::stdin().lock().lines() {
        let answer = match request(serde_json::from_str(&line.unwrap()).unwrap()) {
            Ok(v) => v,
            Err(e) => serde_json::json!({"error":e.to_string(),"probe_error":true}),
        };
        serde_json::to_writer(&mut out, &answer).unwrap();
        writeln!(out).unwrap();
        out.flush().unwrap();
    }
}
