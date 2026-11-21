//! Read-only public production local kernel, with no global solve or adoption.
use soil_fem_3d::joint_law::{Material, PointState};
use std::{fs, io::Write};
fn main() -> soil_fem_3d::Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let input: serde_json::Value = serde_json::from_slice(&fs::read(&args[0])?)?;
    let mesh: soil_fem_3d::mesh::Mesh = serde_json::from_value(input["mesh"].clone())?;
    let config: soil_fem_3d::joints::Configuration =
        serde_json::from_value(input["configuration"].clone())?;
    let model = soil_fem_3d::joints::Model::new(&mesh, config.clone())?;
    let old: soil_fem_3d::joints::State = serde_json::from_value(input["old_state"].clone())?;
    let u: Vec<[f64; 3]> = serde_json::from_value(input["displacement"].clone())?;
    let counts: Vec<usize> = serde_json::from_value(input["counts"].clone())?;
    let eval = model.evaluate(
        &old,
        &u,
        &vec![true; mesh.tetrahedra.len()],
        Some(&counts),
        true,
    )?;
    let mut values = vec![];
    for (i, info) in model.layout().iter().enumerate() {
        let definition = config
            .joints
            .iter()
            .find(|d| d.id == info.joint_id)
            .unwrap();
        values.push(
            serde_json::json!({"material":definition.material,"old":old.points[i],
            "target":eval.state.points[i].jump_m,"count":counts[i],"h":config.tangent_step_m}),
        );
    }
    let mut records = vec![];
    for value in &values {
        let mat: Material = serde_json::from_value(value["material"].clone())?;
        let old: PointState = serde_json::from_value(value["old"].clone())?;
        let target: [f64; 3] = serde_json::from_value(value["target"].clone())?;
        let count = value["count"].as_u64().unwrap() as usize;
        let h = value["h"].as_f64().unwrap();
        let state = mat.integrate_partition(&old, target, count)?;
        let tangent = mat.tangent(&old, target, count, h)?;
        let mut probes = vec![];
        for a in 0..3 {
            for sign in [-1., 1.] {
                let mut x = target;
                x[a] += sign * h;
                probes.push(mat.integrate_partition(&old, x, count)?);
            }
        }
        records.push(serde_json::json!({"state":state,"tangent_columns":(0..3).map(|i|tangent.column(i).as_slice().to_vec()).collect::<Vec<_>>(),"fd_probes":probes}));
    }
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[1])?;
    serde_json::to_writer(
        &mut f,
        &serde_json::json!({"records":records,"local_inputs":values,
        "contact_entries":eval.tangent,"layout":model.layout(),"jump_coefficients":model.jump_coefficients()}),
    )?;
    writeln!(f)?;
    Ok(())
}
