//! A-03 hydrostatic accepted-state observation and public-kernel reassembly.
//! No Newton diagnostics or production diagnostic-budget changes.
use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Kind, Observer, Progress},
    stages,
};
use std::{
    cell::{Cell, RefCell},
    fs::{self, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};
fn write(path: &Path, value: &impl serde::Serialize) -> Result<()> {
    let mut f = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?);
    serde_json::to_writer(&mut f, value)?;
    f.flush()?;
    Ok(())
}
struct Trace {
    stream: RefCell<BufWriter<std::fs::File>>,
    accepted: Cell<usize>,
    diagnostics: RefCell<BufWriter<std::fs::File>>,
    last: RefCell<Option<serde_json::Value>>,
    diagnostic: bool,
}
impl Observer for Trace {
    fn on_progress(&self, p: &Progress) -> Result<()> {
        let mut f = self.stream.borrow_mut();
        serde_json::to_writer(&mut *f, &serde_json::json!({"kind":"progress","event":p}))?;
        writeln!(f)?;
        if p.kind == Kind::Accepted {
            self.accepted.set(self.accepted.get() + 1);
            if let Some(value) = self.last.borrow_mut().take() {
                let mut stream = self.diagnostics.borrow_mut();
                serde_json::to_writer(&mut *stream, &value)?;
                writeln!(stream)?;
                stream.flush()?;
            }
        }
        Ok(())
    }
    fn diagnostic_enabled(&self, _stage: usize, _progress: f64) -> bool {
        self.diagnostic
    }
    fn on_diagnostic(&self, value: &serde_json::Value) -> Result<()> {
        if value["kind"] == "iterate" {
            *self.last.borrow_mut() = Some(serde_json::json!({
                "stage": value["stage"], "progress": value["progress"],
                "partitions": value["partitions"], "strain_increment": value["strain_increment"],
                "relative_force_residual": value["relative_force_residual"],
                "internal": value["internal"], "target": value["target"],
            }));
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(self.accepted.get() >= 1)
    }
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert!(
        args.len() == 4 || args.len() == 6 || (args.len() == 5 && args[4] == "--whole"),
        "MESH GROUND CASE NEW_FOLDER [--resume CHECKPOINT]"
    );
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&args[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&args[2])?)?;
    assert!(mesh.tetrahedra.len() <= 48 || std::env::var_os("A03_USER_PC").is_some());
    let folder = Path::new(&args[3]);
    fs::create_dir(folder)?;
    if args.len() == 5 {
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete);
        write(&folder.join("run.json"), &run)?;
        return Ok(());
    }
    let trace = Trace {
        stream: RefCell::new(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(folder.join("progress.jsonl"))?,
        )),
        accepted: Cell::new(0),
        diagnostics: RefCell::new(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(folder.join("constitutive.jsonl"))?,
        )),
        last: RefCell::new(None),
        diagnostic: false,
    };
    let resumed = args.len() == 6;
    let beams = case
        .beam_system
        .clone()
        .map(|cfg| soil_fem_3d::beams::Model::new(&mesh, cfg))
        .transpose()?;
    let anchors = case
        .anchor_system
        .clone()
        .map(|cfg| soil_fem_3d::anchors::Model::new(&mesh, cfg))
        .transpose()?;
    let mut reports = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(folder.join("structure-reports.jsonl"))?,
    );
    assert!(ground.initialization.water_level_m.is_some());
    assert!(beams.is_none() && anchors.is_none());
    let mut run = if resumed {
        assert_eq!(args[4], "--resume");
        let cp: Checkpoint = serde_json::from_slice(&fs::read(&args[5])?)?;
        cp.resume_observed(&mesh, ground.clone(), case.clone(), &trace)?
    } else {
        stages::solve_observed(&mesh, ground.clone(), case.clone(), &trace)?
    };
    let mut stream = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(folder.join("accepted.jsonl"))?,
    );
    let mut count = 0;
    loop {
        let actual = run.attempts.iter().filter(|a| a.accepted).count()
            + run
                .geostatic
                .as_ref()
                .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count());
        if trace.accepted.get() == 0 {
            assert!(run.complete);
            break;
        }
        assert_eq!(trace.accepted.get(), 1);
        count += 1;
        serde_json::to_writer(&mut stream, &run.last_committed)?;
        writeln!(stream)?;
        stream.flush()?;
        let state = &run.last_committed;
        let diagnostic = observe(&mesh, &ground, &case, state, false)?;
        let mut output = trace.diagnostics.borrow_mut();
        serde_json::to_writer(&mut *output, &diagnostic)?;
        writeln!(output)?;
        output.flush()?;
        drop(output);
        let beam = beams
            .as_ref()
            .zip(state.beams.as_ref())
            .map(|(m, s)| {
                m.evaluate(
                    s,
                    &state.displacement_m,
                    &vec![true; mesh.nodes.len()],
                    state.boundary_factor,
                )
                .map(|v| v.report)
            })
            .transpose()?;
        let anchor = anchors
            .as_ref()
            .zip(state.anchors.as_ref())
            .map(|(m, s)| {
                m.evaluate(s, &state.displacement_m, &state.active_elements, None)
                    .map(|v| v.report)
            })
            .transpose()?;
        serde_json::to_writer(
            &mut reports,
            &serde_json::json!({"beams":beam,"anchors":anchor}),
        )?;
        writeln!(reports)?;
        reports.flush()?;
        if run.complete {
            break;
        }
        assert_eq!(run.failure_reason.as_deref(), Some("checkpoint_requested"));
        let cp = Checkpoint::seal(&mesh, run)?;
        if !resumed && cp.run.last_committed.stage_progress == 0.5 {
            write(&folder.join(format!("checkpoint-{actual:04}.json")), &cp)?;
        }
        trace.accepted.set(0);
        run = cp.resume_observed(&mesh, ground.clone(), case.clone(), &trace)?;
    }
    let mut preparations = Vec::new();
    for audit in &run.transition_audits {
        if case
            .stages
            .iter()
            .any(|s| matches!(s, stages::Stage::Birth { name, .. } if name == &audit.stage))
        {
            preparations.push(serde_json::json!({"stage":audit.stage,
                "diagnostic":observe(&mesh, &ground, &case, &audit.prepared_state, true)?}));
        }
    }
    write(&folder.join("preparations.json"), &preparations)?;
    trace.stream.borrow_mut().flush()?;
    write(&folder.join("run.json"), &run)?;
    write(
        &folder.join("verification.json"),
        &serde_json::json!({"accepted_states_written":count,"complete":run.complete,"resumed_from_file":resumed,"process_id":std::process::id()}),
    )?;
    assert!(run.complete);
    println!("complete; wrote {count} accepted states");
    Ok(())
}

// This is a production-kernel observation of committed state, not the independent solver.
fn observe(
    mesh: &Mesh,
    ground: &Ground,
    case: &stages::StageCase,
    state: &stages::Committed,
    prepared: bool,
) -> Result<serde_json::Value> {
    use soil_fem_3d::{elements, joint_water, joints};
    use std::collections::{BTreeMap, BTreeSet};
    let joints = joints::Model::new(mesh, case.joint_system.clone().unwrap())?;
    let wet = joint_water::Model::new(mesh, &joints, ground)?;
    let level = state
        .hydrostatic_water_level_m
        .or(ground.initialization.water_level_m)
        .unwrap();
    let gamma = ground.initialization.water_unit_weight_kn_m3;
    assert!(state.pore_pressure_kpa.is_none());
    let nd = mesh.nodes.len() * 3;
    let mut effective = vec![0.; nd];
    let mut pore = vec![0.; nd];
    let mut gravity = vec![0.; nd];
    let mut pressure = Vec::new();
    let mut head = Vec::new();
    let mut start = 0;
    for e in 0..mesh.tetrahedra.len() {
        let q = mesh.integration(e)?;
        let nodes = mesh.cell_nodes(e);
        let mid = state
            .element_material_ids
            .get(e)
            .copied()
            .unwrap_or(mesh.material_ids[e]);
        let weight = ground
            .materials
            .iter()
            .find(|m| m.id == mid)
            .unwrap()
            .unit_weight_kn_m3;
        for (i, p) in q.iter().enumerate() {
            let water = gamma * (level - p.position[2]);
            assert!(water >= 0.);
            pressure.push(water);
            head.push(p.position[2] + water / gamma);
            if state.active_elements[e] {
                let stress = nalgebra::SVector::<f64, 6>::from_row_slice(
                    state.material_states[start + i].stress(),
                );
                let ps = nalgebra::SVector::<f64, 6>::new(water, water, water, 0., 0., 0.);
                let fe = -p.weight_m3 * p.b.transpose() * stress;
                let fp = -p.weight_m3 * p.b.transpose() * ps;
                for (j, &n) in nodes.iter().enumerate() {
                    for a in 0..3 {
                        effective[3 * n + a] += fe[3 * j + a];
                        pore[3 * n + a] += fp[3 * j + a];
                    }
                    gravity[3 * n + 2] -= p.weight_m3 * weight * p.shape[j];
                }
            }
        }
        start += q.len();
    }
    let contact = joints.evaluate(
        state.joints.as_ref().unwrap(),
        &state.displacement_m,
        &state.active_elements,
        None,
        false,
    )?;
    let internal: Vec<_> = effective
        .iter()
        .zip(&pore)
        .zip(&contact.internal_kn)
        .map(|((a, b), c)| a + b + c)
        .collect();
    let mut topology = BTreeMap::<[usize; 3], Vec<usize>>::new();
    for (e, t) in mesh.tetrahedra.iter().enumerate() {
        for omitted in 0..4 {
            let mut f: Vec<_> = t
                .iter()
                .enumerate()
                .filter_map(|(j, &n)| (j != omitted).then_some(n))
                .collect();
            f.sort_unstable();
            topology.entry(f.try_into().unwrap()).or_default().push(e);
        }
    }
    let cfg = ground.initialization.geostatic.as_ref().unwrap();
    let mut hydro = BTreeSet::new();
    for name in &cfg.hydrostatic_face_sets {
        for face in &mesh.face_sets[name] {
            let mut f = *face;
            f.sort_unstable();
            hydro.insert(f);
        }
    }
    let top = mesh
        .nodes
        .iter()
        .map(|p| p[2])
        .fold(f64::NEG_INFINITY, f64::max);
    let mut surcharge = vec![0.; nd];
    let mut exterior_water = vec![0.; nd];
    let mut exterior = Vec::new();
    let mut water_faces = Vec::new();
    let mut face_ids = Vec::new();
    for (f, owners) in &topology {
        let live: Vec<_> = owners
            .iter()
            .copied()
            .filter(|&e| state.active_elements[e])
            .collect();
        let opposite = joints.counterpart(f);
        let exposed = live.len() == 1 && !opposite.is_some_and(|e| state.active_elements[e]);
        let wet_face = exposed && (owners.len() == 2 || opposite.is_some() || hydro.contains(f));
        exterior.push(exposed);
        water_faces.push(wet_face);
        face_ids.push(*f);
        if exposed {
            let nodes = mesh.face_nodes(live[0], f)?;
            let points: Vec<_> = nodes.iter().map(|&n| mesh.nodes[n]).collect();
            for q in elements::face_points(&points)? {
                let n = q.area_vector;
                let norm = n.iter().map(|x| x * x).sum::<f64>().sqrt();
                let up = n[2] > 0. && n[0].abs() + n[1].abs() < 1e-8 * norm;
                let loaded = up
                    && (!case.pressure_on_original_surface_only
                        || (q.position[2] - top).abs() < 1e-8);
                for (j, &node) in nodes.iter().enumerate() {
                    for a in 0..3 {
                        if loaded {
                            surcharge[3 * node + a] -=
                                q.shape[j] * n[a] * state.surface_pressure_kpa;
                        }
                        if wet_face {
                            exterior_water[3 * node + a] -=
                                q.shape[j] * n[a] * gamma * (level - q.position[2]);
                        }
                    }
                }
            }
        }
    }
    let water = if prepared {
        wet.regenerate_birth(&state.active_elements, Some(level), None)?
    } else {
        state.joint_water.clone().unwrap()
    };
    let base = cfg.base_loads.on_active(mesh, &state.active_elements)?;
    Ok(
        serde_json::json!({"source":"committed-state public production kernels; no independent solve",
        "stage":state.stage_index,"progress":state.stage_progress,
        "internal":internal,"effective_internal":effective,"pore_internal":pore,"contact_internal":contact.internal_kn,
        "gravity_load":gravity,"base_load":base,"surcharge_load":surcharge,"exterior_water_load":exterior_water,
        "joint_water_load":water.nodal_load_kn,"joint_water_pressure":water.physical_pressure_kpa,
        "joint_face_water_force":water.face_nodal_force_a_kn,"pore_pressure":pressure,"head":head,
        "exterior_faces":exterior,"wet_exterior_faces":water_faces,"face_ids":face_ids}),
    )
}
