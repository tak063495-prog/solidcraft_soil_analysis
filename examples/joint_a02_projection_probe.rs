//! Diagnostic public-API requests only; no production input or solver changes.
use soil_fem_3d::{Result, joints, mesh::Mesh, stages::StageCase};
use std::fs;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 4, "MESH CASE REQUEST OUTPUT");
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let case: StageCase = serde_json::from_slice(&fs::read(&args[1])?)?;
    let request: serde_json::Value = serde_json::from_slice(&fs::read(&args[2])?)?;
    let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
    let layout = model.layout();
    let u = vec![[0.; 3]; mesh.nodes.len()];
    let mut active = vec![true; mesh.tetrahedra.len()];
    let rows: Vec<_> = (0..3 * u.len()).map(|g| vec![(g, 1.)]).collect();
    let ids = vec!["Seam".into()];
    let mut force: Vec<f64> = serde_json::from_value(request["residual_kn"].clone())?;
    assert_eq!(force.len(), rows.len());
    let initial = model.initial_state()?;
    let mut state = initial.clone();
    let operation = request["operation"].as_str().unwrap();
    if operation == "partial_rebirth" {
        state = model.birth(&initial, &u, &active, &ids, &vec![0.; rows.len()], &rows)?;
        state = model.death(&state, &ids)?;
    }
    if let Some(owner) = request["inactive_owner"].as_u64() {
        active[owner as usize] = false;
    }
    let result = if operation == "body_birth" {
        let mut before = active.clone();
        before[layout[0].owner_b] = false;
        model.prepare_body_birth(&state, &before, &active)
    } else {
        model.birth(&state, &u, &active, &ids, &force, &rows)
    };
    let (accepted, error, output) = match result {
        Ok(born) => {
            let repeat = model.birth(&state, &u, &active, &ids, &force, &rows)?;
            assert_eq!(serde_json::to_vec(&born)?, serde_json::to_vec(&repeat)?);
            let value = model.evaluate(&born, &u, &active, None, false)?;
            let maximum = value
                .internal_kn
                .iter()
                .zip(&force)
                .map(|(a, b)| (a - b).abs())
                .fold(0_f64, f64::max);
            assert!(maximum <= 1e-10);
            (
                true,
                None,
                serde_json::json!({"state":born,"nodal_force_kn":value.internal_kn,
                "maximum_force_error_kn":maximum,"exact_repeat":true}),
            )
        }
        Err(err) => (false, Some(err.to_string()), serde_json::Value::Null),
    };
    // Keep the exact requested vector; never clip, smooth, or redistribute it.
    force.shrink_to_fit();
    let output = serde_json::json!({"pid":std::process::id(),"accepted":accepted,
        "error":error,"result":output,"requested_residual_kn":force,"layout":layout});
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[3])?;
    serde_json::to_writer_pretty(file, &output)?;
    if !accepted {
        std::process::exit(1);
    }
    Ok(())
}
