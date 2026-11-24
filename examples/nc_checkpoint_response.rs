//! Evaluate the original force residual and specified tangent actions at a saved Run.
use soil_fem_3d::{mesh::Mesh, soil, stages::Run};
use std::{env, error::Error, fs, io::BufWriter};

// Only unassembled point stresses cross this interface. Python independently
// differences them, constructs the element operators and scatters to the mesh.
fn material_stencils(run: &Run) -> Result<Vec<serde_json::Value>, Box<dyn Error>> {
    let mut rows = Vec::new();
    let old = &run.last_committed;
    let h = run.case.settings.tangent_perturbation;
    for (i, gp) in run.integration_points.iter().enumerate() {
        let e = gp.element_index;
        if !old.active_elements[e] {
            continue;
        }
        let material = run
            .ground
            .materials
            .iter()
            .find(|m| m.id == old.element_material_ids[e])
            .ok_or("unknown checkpoint material")?;
        let (soil::Model::Mcc(mut law), soil::State::Mcc(state)) =
            (material.model.clone(), &old.material_states[i])
        else {
            return Err("expected MCC checkpoint".into());
        };
        if let Some(integration) = run.case.settings.mcc_integration {
            law.integration = integration;
        }
        let mut stresses = Vec::new();
        let mut plastic = Vec::new();
        for j in 0..6 {
            for width in [-h, h] {
                let mut strain = [0.; 6];
                strain[j] = width;
                let update = law.integrate_partition(state, strain, 1)?;
                stresses.push(update.state.effective_stress_compression_positive_kpa);
                plastic.push(update.plastic_substeps);
            }
        }
        rows.push(serde_json::json!({"point":i,"stress":stresses,"plastic":plastic}));
    }
    Ok(rows)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    if let [mode, run, output] = args.as_slice()
        && mode == "--local-stress-stencil"
    {
        let run: Run = serde_json::from_slice(&fs::read(run)?)?;
        serde_json::to_writer(
            BufWriter::new(fs::File::create_new(output)?),
            &material_stencils(&run)?,
        )?;
        return Ok(());
    }
    let [mesh, run, directions, output] = args.as_slice() else {
        return Err("usage: nc_checkpoint_response MESH RUN DIRECTIONS NEW_OUTPUT".into());
    };
    let mesh: Mesh = serde_json::from_slice(&fs::read(mesh)?)?;
    let run: Run = serde_json::from_slice(&fs::read(run)?)?;
    let directions = serde_json::from_slice::<Vec<Vec<f64>>>(&fs::read(directions)?)?;
    let response = run.checkpoint_response(&mesh, &directions)?;
    serde_json::to_writer(BufWriter::new(fs::File::create_new(output)?), &response)?;
    Ok(())
}
