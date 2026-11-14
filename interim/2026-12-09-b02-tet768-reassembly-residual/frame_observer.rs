// Diagnostic executable: only calls existing public production APIs.
use serde::de::DeserializeOwned;
use soil_fem_3d::{
    Result, beams,
    ground::Ground,
    mesh::Mesh,
    stages::{self, Run, StageCase},
};
use std::{
    fs::{self, File, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};

fn read<T: DeserializeOwned>(p: &str) -> Result<T> {
    Ok(serde_json::from_slice(&fs::read(p)?)?)
}
fn write(p: &Path, v: &impl serde::Serialize) -> Result<()> {
    let f = OpenOptions::new().write(true).create_new(true).open(p)?;
    serde_json::to_writer(BufWriter::new(f), v)?;
    Ok(())
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let m: Mesh = read(&a[2])?;
    let g: Ground = read(&a[3])?;
    let c: StageCase = read(&a[4])?;
    let out = Path::new(&a[5]);
    fs::create_dir(out)?;
    if a[1] == "plain" {
        let r = stages::solve(&m, g, c)?;
        assert!(r.complete, "{:?}", r.failure_reason);
        return write(&out.join("run.json"), &r);
    }
    let model = beams::Model::new(&m, c.beam_system.clone().unwrap())?;
    let quadrature: Vec<_> = (0..m.tetrahedra.len())
        .map(|i| m.integration(i))
        .collect::<Result<_>>()?;
    let ops: Vec<_> = quadrature.iter().enumerate().map(|(ei, points)| {
        serde_json::json!({"nodes": m.cell_nodes(ei), "weights": points.iter().map(|q| q.weight_m3).collect::<Vec<_>>(),
        "b_column_major": points.iter().map(|q| q.b.as_slice().to_vec()).collect::<Vec<_>>()})
    }).collect();
    let cfg = &model.configuration;
    let node_index = |id: &String| cfg.nodes.iter().position(|n| &n.id == id).unwrap();
    let mut elements = Vec::new();
    for (i, p) in cfg.plates.iter().enumerate() {
        let nodes: Vec<_> = p.nodes.iter().map(node_index).collect();
        let positions: Vec<_> = nodes.iter().map(|n| model.positions[*n]).collect();
        let section = cfg
            .plate_sections
            .iter()
            .find(|s| s.id == p.section_id)
            .unwrap();
        let kernel = soil_fem_3d::plates::Element::new_with_offset(
            &positions,
            section,
            p.formulation,
            p.attachment_offset_m,
        )?;
        elements.push(serde_json::json!({"kind":"plate", "index":i, "nodes":nodes, "stiffness_column_major":kernel.stiffness().as_slice()}));
    }
    for (i, b) in cfg.members.iter().enumerate() {
        assert!(!b.rotational_releases.iter().flatten().any(|x| *x));
        assert_eq!(b.initial_axial_force_kn, 0.);
        let nodes: Vec<_> = b.nodes.iter().map(node_index).collect();
        let section = cfg.sections.iter().find(|s| s.id == b.section_id).unwrap();
        let kernel = soil_fem_3d::beams::Kernel::new(
            model.positions[nodes[0]],
            model.positions[nodes[1]],
            b.local_y,
            section,
            b.rotational_releases,
        )?;
        elements.push(serde_json::json!({"kind":"beam", "index":i, "nodes":nodes, "stiffness_column_major":kernel.stiffness.as_slice()}));
    }
    write(
        &out.join("operators.json"),
        &serde_json::json!({"soil":ops, "structures":elements}),
    )?;
    let mut history = BufWriter::new(File::create(out.join("accepted.jsonl"))?);
    let r: Run = read(&a[6])?;
    assert!(r.complete, "{:?}", r.failure_reason);
    let mut n = 0;
    for (frame_index, frame) in r.frames.iter().enumerate() {
        let geo = frame_index == 0;
        let s = &frame.state;
        let progress = 1.;
        n += if geo {
            r.geostatic
                .as_ref()
                .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count())
        } else {
            r.attempts
                .iter()
                .filter(|a| a.accepted && a.stage == frame.name)
                .count()
        };
        let e = model.evaluate(
            s.beams.as_ref().unwrap(),
            &s.displacement_m,
            &vec![true; m.nodes.len()],
            s.boundary_factor,
        )?;
        // Production quadrature + actual saved production stress. This is
        // observation/reassembly, never a new independent solution.
        let mut soil_internal = vec![0.; 3 * m.nodes.len()];
        let mut point = 0;
        for (ei, points) in quadrature.iter().enumerate() {
            let nodes = m.cell_nodes(ei);
            for q in points {
                let stress = s.material_states[point].stress();
                for (i, node) in nodes.iter().enumerate() {
                    for a in 0..3 {
                        soil_internal[3 * node + a] -= q.weight_m3
                            * (0..6).map(|j| q.b[(j, 3 * i + a)] * stress[j]).sum::<f64>();
                    }
                }
                point += 1;
            }
        }
        serde_json::to_writer(
            &mut history,
            &serde_json::json!({
                "accepted_index": n, "stage_index": if geo {0} else {s.stage_index},
                "progress": progress, "state": s, "beam_results": e.report,
                "reevaluated_beams": e.state, "soil_internal_kn": soil_internal, "interface_dofs": e.interface_dofs,
                "interface_internal_kn": e.interface_internal_kn.as_slice()
            }),
        )?;
        history.write_all(b"\n")?;
        history.flush()?;
    }
    write(&out.join("run.json"), &r)?;
    Ok(())
}
