use serde_json::{Value, json};
use soil_fem_3d::{Result, ground::Ground, mesh::Mesh, one_way};

fn fixture(policy: Value) -> Result<(Mesh, Ground, one_way::Case)> {
    let mesh = serde_json::from_str(include_str!("../examples/joint-water-column-mesh.json"))?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-joint-water-mc.json"))?;
    ground.initialization.water_level_m = None;
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .hydrostatic_face_sets
        .clear();
    let joint: Value =
        serde_json::from_str(include_str!("../examples/stages-joint-water-cycle.json"))?;
    let mut case: Value = serde_json::from_str(include_str!("../examples/one-way-mc.json"))?;
    for key in [
        "settings",
        "boundary",
        "joint_system",
        "monitors",
        "monitor_sets",
        "pressure_on_original_surface_only",
    ] {
        case["mechanics"][key] = joint[key].clone();
    }
    if policy["kind"] == "soil" {
        case["mechanics"]["joint_system"]["structure_bodies"] = json!(["UpperClay"]);
    }
    case["mechanics"]["joint_system"]["joints"][0]["water"] = policy;
    case["flow"]["initial_head"] = json!({"head_m":2.});
    case["flow"]["analysis"]["periods"] = json!([{"name":"step","end_time_s":0.25,"boundaries":[
        {"face_set":"bottom","condition":{"kind":"head","value":{"head_m":2.1}}}]}]);
    case["flow"]["materials"][0]["volume_source_per_s"] = json!(1e-6);
    case["flow"]["settings"]["initial_time_step_s"] = json!(0.125);
    case["flow"]["settings"]["maximum_time_step_s"] = json!(0.125);
    case["exterior_periods"] = json!([{"name":"step","reservoirs":[]}]);
    Ok((mesh, ground, serde_json::from_value(case)?))
}

#[test]
fn same_endpoint_continuous_rt0_and_opposing_outflows() -> Result<()> {
    let (mesh, ground, case) = fixture(json!({"kind":"continuous"}))?;
    let run = one_way::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for frame in run.frames.iter().filter(|f| f.phase == "transient") {
        let points = &frame.state.fluid.joint_traces;
        assert!(!points.is_empty());
        for point in points {
            let a = point.side_a.as_ref().unwrap();
            let b = point.side_b.as_ref().unwrap();
            assert_ne!(a.owner, b.owner);
            assert_eq!((a.face, a.total_head_m), (b.face, b.total_head_m));
            assert_eq!(point.time_s, frame.state.hydraulic.time_s);
            assert_eq!(
                a.signed_pressure_kpa,
                10. * (a.total_head_m - point.elevation_m)
            );
            assert!((a.outward_flow_m3_s + b.outward_flow_m3_s).abs() < 1e-14);
        }
    }
    Ok(())
}

#[test]
fn soil_requires_declared_side_and_keeps_structure_dry() -> Result<()> {
    let (mesh, ground, case) = fixture(json!({"kind":"soil","soil_side":"a"}))?;
    let run = one_way::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for p in &run.last_committed.fluid.joint_traces {
        assert!(p.side_b.is_none());
        assert!(p.side_a.as_ref().unwrap().outward_flow_m3_s.abs() < 1e-14);
    }
    let (mesh, g, c) = fixture(json!({"kind":"soil","soil_side":"b"}))?;
    assert!(
        one_way::solve(&mesh, g, c)
            .unwrap_err()
            .to_string()
            .contains("declared soil side")
    );
    Ok(())
}

#[test]
fn prescribed_independent_wall_integral_and_duplicate_exterior() -> Result<()> {
    let (mesh, ground, mut case) = fixture(json!({"kind":"prescribed","level_m":2.}))?;
    let run = one_way::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let field = run.last_committed.mechanical.joint_water.as_ref().unwrap();
    let mesh_faces = &mesh.face_sets["Seam:a"];
    for (key, face) in mesh_faces.iter().zip(&field.face_nodal_force_a_kn) {
        let a = nalgebra::Vector3::from(mesh.nodes[key[0]]);
        let b = nalgebra::Vector3::from(mesh.nodes[key[1]]);
        let c = nalgebra::Vector3::from(mesh.nodes[key[2]]);
        let area = (b - a).cross(&(c - a)).norm() / 2.;
        assert_eq!(a.z, b.z);
        for force in face {
            assert!((force[2].abs() - 10. * (2. - a.z) * area / 3.).abs() < 1e-12);
        }
    }
    for p in &run.last_committed.fluid.joint_traces {
        assert!((p.side_a.as_ref().unwrap().total_head_m - 2.).abs() < 1e-10);
        assert!((p.side_b.as_ref().unwrap().total_head_m - 2.).abs() < 1e-10);
    }
    case.initial_exterior_water = vec![one_way::Reservoir {
        face_set: "Seam:b".into(),
        level_m: 2.,
    }];
    assert!(
        one_way::solve(&mesh, ground, case)
            .unwrap_err()
            .to_string()
            .contains("exterior and JOINT")
    );
    Ok(())
}

#[test]
fn rejected_endpoint_restores_contact_and_hydraulics() -> Result<()> {
    let (mesh, ground, mut case) = fixture(json!({"kind":"continuous"}))?;
    case.verification_reject_once_at_times_s = vec![0.125];
    let run = one_way::solve(&mesh, ground, case)?;
    assert!(run.complete);
    let a = run.attempts.iter().find(|a| !a.accepted).unwrap();
    assert_eq!(
        a.committed_pair_sha256_before,
        a.committed_pair_sha256_after
    );
    let times: Vec<_> = run
        .frames
        .iter()
        .filter(|f| f.phase == "transient")
        .map(|f| f.state.hydraulic.time_s)
        .collect();
    assert!(times.windows(2).all(|w| w[0] < w[1]));
    Ok(())
}

#[test]
fn staggered_zero_feedback_reassembles_same_final_pair() -> Result<()> {
    let (mesh, g, mut case) = fixture(json!({"kind":"continuous"}))?;
    let one = one_way::solve(&mesh, g.clone(), case.clone())?;
    case.schema = "soil-fem-staggered-transfer/1".into();
    let repeated = one_way::solve(&mesh, g, case)?;
    assert!(one.complete && repeated.complete);
    assert_eq!(
        serde_json::to_value(&one.last_committed)?,
        serde_json::to_value(&repeated.last_committed)?
    );
    assert!(
        repeated
            .frames
            .iter()
            .filter(|f| f.phase == "transient")
            .all(|f| f.staggered_iterations == Some(2))
    );
    Ok(())
}

#[test]
fn prescribed_history_splits_at_knot_and_loads_both_walls_at_endpoint() -> Result<()> {
    let (mesh, g, mut case) = fixture(json!({"kind":"prescribed"}))?;
    case.joint_reservoir_heads =
        serde_json::from_value(json!({"Seam":[[0.,2.],[0.0625,2.25],[0.25,1.875]]}))?;
    let run = one_way::solve(&mesh, g.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let times: Vec<_> = run
        .frames
        .iter()
        .filter(|f| f.phase == "transient")
        .map(|f| f.state.hydraulic.time_s)
        .collect();
    assert_eq!(times[0], 0.0625);
    assert!(times.iter().any(|t| *t > 0.0625 && *t < 0.25));
    for frame in &run.frames {
        let t = frame.state.hydraulic.time_s;
        let h = if t <= 0.0625 {
            2. + 4. * t
        } else {
            2.25 - 2. * (t - 0.0625)
        };
        for p in &frame.state.fluid.joint_traces {
            assert!((p.applied_pressure_kpa - 10. * (h - p.elevation_m)).abs() < 1e-12);
            for side in [&p.side_a, &p.side_b] {
                assert!((side.as_ref().unwrap().total_head_m - h).abs() < 1e-10);
            }
        }
    }
    for invalid in [
        json!([[0., 2.], [0., 3.]]),
        json!([[0.1, 2.], [0.25, 3.]]),
        json!([[0., 2.], [0.2, 3.]]),
    ] {
        case.joint_reservoir_heads = serde_json::from_value(json!({"Seam":invalid}))?;
        assert!(
            one_way::solve(&mesh, g.clone(), case.clone())
                .unwrap_err()
                .to_string()
                .contains("reservoir head table")
        );
    }
    Ok(())
}

#[test]
fn duplicate_wet_frames_keep_contact_labels_and_record_zero_slip_and_work() -> Result<()> {
    let mesh = serde_json::from_str(include_str!("../examples/joint-water-column-mesh.json"))?;
    let g = serde_json::from_str(include_str!("../examples/ground-joint-water-mc.json"))?;
    let mut c: Value =
        serde_json::from_str(include_str!("../examples/stages-joint-water-cycle.json"))?;
    c["verification_reject_once_stages"] = json!([]);
    let mut stages = vec![
        json!({"kind":"boundary","name":"repeat_initial","target_factor":0.0,"initial_increment":1.0}),
    ];
    for stage in c["stages"].as_array().unwrap() {
        stages.push(stage.clone());
        let mut copy = stage.clone();
        copy["name"] = json!(format!("repeat_{}", stage["name"].as_str().unwrap()));
        stages.push(copy);
    }
    c["stages"] = json!(stages);
    let run = soil_fem_3d::stages::solve(&mesh, g, serde_json::from_value(c)?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let mut observed = std::collections::BTreeSet::new();
    for pair in run
        .frames
        .windows(2)
        .filter(|v| v[1].name.starts_with("repeat_"))
    {
        let previous = pair[0].state.joints.as_ref().unwrap();
        let next = pair[1].state.joints.as_ref().unwrap();
        assert_eq!(next, previous);
        for (p, s) in pair[1]
            .joint_results
            .as_ref()
            .unwrap()
            .points
            .iter()
            .zip(&next.points)
        {
            let increment = p.interval.as_ref().unwrap();
            assert_eq!(increment.status, s.status);
            assert_eq!((increment.slip_m, increment.endpoint_work_kpa_m), (0., 0.));
            observed.insert(serde_json::to_string(&s.status)?);
        }
    }
    assert_eq!(observed.len(), 3, "{observed:?}");
    Ok(())
}

#[test]
fn wet_construction_archives_removed_generation_and_restarts_the_reinstalled_wall() -> Result<()> {
    use soil_fem_3d::one_way_construction as construction;
    let (mesh, ground, base) = fixture(json!({"kind":"continuous"}))?;
    let mut one = serde_json::to_value(base)?;
    one["mechanics"]["joint_system"]["joints"][0]["birth_policy"] = json!("stress_free");
    let period_names = ["hold", "place", "remove", "replace"];
    let periods = period_names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut boundaries = vec![
                json!({"face_set":"bottom","condition":{"kind":"head","value":{"head_m":2.}}}),
            ];
            let mut faces = vec!["xmin", "xmax", "ymin", "ymax"];
            faces.push(if i % 2 == 0 { "Seam:a" } else { "top" });
            for face in faces {
                boundaries.push(
                    json!({"face_set":face,"condition":{"kind":"outward_flux","flux_m_s":0.}}),
                );
            }
            json!({"name":name,"end_time_s":(i+1) as f64 * 0.125,"boundaries":boundaries})
        })
        .collect::<Vec<_>>();
    one["flow"]["analysis"]["periods"] = json!(periods);
    one["exterior_periods"] = json!(
        period_names
            .iter()
            .map(|name| json!({"name":name,"reservoirs":[]}))
            .collect::<Vec<_>>()
    );
    let template: Value = serde_json::from_str(include_str!("../examples/one-way-lifecycle.json"))?;
    let transitions = ["place", "replace"].map(|name| {
        let mut event = template["transitions"][0].clone();
        event["name"] = json!(name);
        event["before_period"] = json!(name);
        event["water_source"] = json!({"reservoir_id":"fill_water"});
        event["action"]["initial_head"] = json!({"head_m":2.});
        event
    });
    let input = json!({"schema":"soil-fem-one-way-construction/1","one_way":one,
        "initially_inactive_bodies":["UpperClay"],"transitions":transitions,
        "excavations":[{"name":"remove","before_period":"remove","bodies":["UpperClay"],
        "water_boundary_policy":"explicit_period_boundaries","cut_traction_policy":"explicit_period_reservoirs",
        "water_destination":{"reservoir_id":"haul_water"}}]});
    let mut missing = input.clone();
    missing["transitions"][0]
        .as_object_mut()
        .unwrap()
        .remove("water_source");
    assert!(
        construction::solve(&mesh, ground.clone(), serde_json::from_value(missing)?)
            .unwrap_err()
            .to_string()
            .contains("named event water source/destination")
    );
    let case: construction::Case = serde_json::from_value(input)?;
    let run = construction::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.transitions.len(), 2);
    for (index, r) in run.transitions.iter().enumerate() {
        let state = r.prepared_mechanical.joints.as_ref().unwrap();
        assert_eq!(state.lifecycle[0].generation, (index + 1) as u64);
        assert_eq!(r.water_source.as_ref().unwrap().reservoir_id, "fill_water");
        assert!(r.water.water_added_m3 > 0.);
        for p in &state.points {
            assert_eq!(p.birth.as_ref().unwrap().initial_traction_kpa, [0.; 3]);
            assert_eq!(p.accumulated_plastic_slip_m, 0.);
        }
    }
    let cut = &run.excavations[0];
    assert_eq!(
        cut.water_destination.as_ref().unwrap().reservoir_id,
        "haul_water"
    );
    assert!(cut.removal.water_removed_m3 > 0.);
    let old = run
        .frames
        .iter()
        .rev()
        .find(|f| f.phase == "transient" && f.state.hydraulic.time_s == cut.time_s)
        .unwrap();
    let archive = &run
        .last_committed
        .mechanical
        .joints
        .as_ref()
        .unwrap()
        .lifecycle[0]
        .archived;
    assert_eq!(archive.len(), 1);
    for (i, p) in &archive[0].points {
        assert_eq!(p, &old.state.mechanical.joints.as_ref().unwrap().points[*i]);
    }
    let partial =
        construction::solve_until(&mesh, ground.clone(), case.clone(), run.frames.len() - 1)?;
    let resumed = construction::Checkpoint::seal(partial)?.resume(&mesh, ground, case, None)?;
    assert_eq!(
        construction::normalized(&run)?,
        construction::normalized(&resumed)?
    );
    Ok(())
}

#[test]
fn continuous_transmission_matches_uncut_column_with_nonzero_crossflow() -> Result<()> {
    use soil_fem_3d::richards;
    let (mesh, ground, mut case) = fixture(json!({"kind":"continuous"}))?;
    let richards::Analysis::Transient { periods } = &mut case.flow.analysis else {
        unreachable!()
    };
    periods[0].boundaries.push(serde_json::from_value(json!({
        "face_set":"top","condition":{"kind":"head","value":{"head_m":2.0}}
    }))?);
    let mut uncut: Mesh =
        serde_json::from_str(include_str!("../examples/seepage-column-tet4-mesh.json"))?;
    uncut.material_ids.clone_from(&mesh.material_ids);
    for (i, cell) in mesh.tetrahedra.iter().enumerate() {
        assert_eq!(
            cell.map(|n| mesh.nodes[n]),
            uncut.tetrahedra[i].map(|n| uncut.nodes[n])
        );
    }
    let reference = richards::solve(
        &uncut,
        case.flow.clone(),
        &vec![true; uncut.tetrahedra.len()],
        None,
    )?;
    let split = one_way::solve(&mesh, ground, case)?;
    assert!(reference.complete && split.complete);
    let steps = split
        .frames
        .iter()
        .filter(|f| f.phase == "transient")
        .collect::<Vec<_>>();
    assert_eq!(steps.len() + 1, reference.frames.len());
    for (actual, expected) in steps.into_iter().zip(&reference.frames[1..]) {
        let h = &actual.state.hydraulic;
        assert_eq!(h.time_s, expected.state.time_s);
        for (i, v) in h.relative_cell_head_m.iter().enumerate() {
            assert!((v.unwrap() - expected.state.relative_cell_head_m[i].unwrap()).abs() < 1e-10);
        }
        let mut crossing = false;
        for p in &actual.state.fluid.joint_traces {
            for s in [&p.side_a, &p.side_b].into_iter().flatten() {
                let flow = expected.state.cell_outflow_m3_s[s.owner].unwrap()[s.local_face];
                assert!((s.outward_flow_m3_s - flow).abs() < 1e-14);
                crossing |= flow.abs() > 1e-12;
            }
        }
        assert!(crossing);
    }
    Ok(())
}
