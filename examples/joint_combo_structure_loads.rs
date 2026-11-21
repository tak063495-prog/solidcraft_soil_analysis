//! Observation only: existing frame/plate physical nodal loads at accepted states.
use soil_fem_3d::{Result, beams, mesh::Mesh, plates, stages};
use std::{
    collections::BTreeMap,
    fs::{self, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
};

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 4, "MESH CASE ACCEPTED_JSONL NEW_LOADS_JSONL");
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&args[1])?)?;
    let cfg = case
        .beam_system
        .as_ref()
        .expect("frame/plate input required");
    let model = beams::Model::new(&mesh, cfg.clone())?;
    let ids: BTreeMap<_, _> = cfg
        .nodes
        .iter()
        .enumerate()
        .map(|(i, n)| (&n.id, i))
        .collect();
    let sections: BTreeMap<_, _> = cfg.sections.iter().map(|s| (&s.id, s)).collect();
    let plate_sections: BTreeMap<_, _> = cfg.plate_sections.iter().map(|s| (&s.id, s)).collect();
    let mut output = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&args[3])?,
    );
    for line in BufReader::new(fs::File::open(&args[2])?).lines() {
        let state: stages::Committed = serde_json::from_str(&line?)?;
        let history = state.beams.as_ref().expect("missing accepted frame state");
        let evaluated = model.evaluate(
            history,
            &state.displacement_m,
            &vec![true; mesh.nodes.len()],
            state.boundary_factor,
        )?;
        let positions = &evaluated.report.node_positions_m;
        let mut loads = vec![0.; 6 * cfg.nodes.len()];
        for (i, member) in cfg.members.iter().enumerate() {
            if !history.active_members[i] {
                continue;
            }
            assert_eq!(
                member.rotational_releases, [[false; 3]; 2],
                "bounded observation excludes releases"
            );
            let ends = member.nodes.each_ref().map(|n| ids[n]);
            let section = sections[&member.section_id];
            let kernel = beams::Kernel::new(
                positions[ends[0]],
                positions[ends[1]],
                member.local_y,
                section,
                member.rotational_releases,
            )?;
            let mut line = std::array::from_fn(|a| {
                member.constant_line_load_kn_m[a]
                    + state.boundary_factor * member.boundary_line_load_kn_m[a]
            });
            line[2] -= section.unit_weight_kn_m3 * section.area_m2;
            let force = kernel.transform.transpose() * kernel.uniform_load(line);
            for end in 0..2 {
                for a in 0..6 {
                    loads[6 * ends[end] + a] += force[6 * end + a];
                }
            }
        }
        for (i, plate) in cfg.plates.iter().enumerate() {
            if !history.active_plates[i] {
                continue;
            }
            assert_eq!(
                plate.formulation,
                plates::Formulation::CstDkt,
                "bounded DKT/CST observation"
            );
            assert_eq!(plate.attachment_offset_m, 0.);
            let nodes: Vec<_> = plate.nodes.iter().map(|n| ids[n]).collect();
            let section = plate_sections[&plate.section_id];
            let kernel = plates::Element::new(
                &nodes.iter().map(|&n| positions[n]).collect::<Vec<_>>(),
                section,
                plate.formulation,
            )?;
            let mut traction = std::array::from_fn(|a| {
                plate.constant_traction_kpa[a]
                    + state.boundary_factor * plate.boundary_traction_kpa[a]
            });
            traction[2] -= section.unit_weight_kn_m3 * section.thickness_m;
            let force = kernel.uniform_load(traction);
            for (end, &node) in nodes.iter().enumerate() {
                for a in 0..6 {
                    loads[6 * node + a] += force[6 * end + a];
                }
            }
        }
        for load in &cfg.nodal_loads {
            let n = ids[&load.node_id];
            for a in 0..3 {
                loads[6 * n + a] +=
                    load.constant_force_kn[a] + state.boundary_factor * load.boundary_force_kn[a];
                loads[6 * n + 3 + a] += load.constant_moment_kn_m[a]
                    + state.boundary_factor * load.boundary_moment_kn_m[a];
            }
        }
        let mut attached = vec![0.; 3 * mesh.nodes.len()];
        for (g, value) in model.attached_load(history, state.boundary_factor)? {
            attached[g] += value;
        }
        serde_json::to_writer(
            &mut output,
            &serde_json::json!({"factor":state.boundary_factor,"structure_nodal_force_moment":loads,"attached_nodal_load_kn":attached}),
        )?;
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
