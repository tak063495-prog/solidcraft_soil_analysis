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
fn joint_policy_is_checked_before_split_displacement_equilibration() -> Result<()> {
    let (mesh, _, mut case) = fixture(json!({"kind": "continuous"}))?;
    let ground: Ground = serde_json::from_str(include_str!("../examples/ground-one-way-mc.json"))?;
    let with_policy = one_way::solve(&mesh, ground.clone(), case.clone())
        .unwrap_err()
        .to_string();
    assert!(with_policy.contains("general geostatic equilibration for split JOINT displacements"));
    assert!(with_policy.contains("construction/structures need their coupled policies"));

    case.mechanics.joint_system.as_mut().unwrap().joints[0].water = None;
    let without_policy = one_way::solve(&mesh, ground, case).unwrap_err().to_string();
    assert!(without_policy.contains("construction/structures need their coupled policies"));
    assert!(!without_policy.contains("general geostatic equilibration"));
    Ok(())
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

fn lifecycle_fixture(policy: Value) -> Result<(Mesh, Ground, Value)> {
    let (mesh, ground, base) = fixture(policy)?;
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
    Ok((mesh, ground, input))
}

#[test]
fn wet_construction_archives_removed_generation_and_restarts_the_reinstalled_wall() -> Result<()> {
    use soil_fem_3d::one_way_construction as construction;
    let (mesh, ground, input) = lifecycle_fixture(json!({"kind":"continuous"}))?;
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
    for (&element, &generation) in &run.transitions[0].placed_generations {
        assert_eq!(generation, 1);
        assert_eq!(cut.removed_generations[&element], 1);
        assert_eq!(run.transitions[1].placed_generations[&element], 2);
    }
    assert!(!cut.removed_hydraulic_points.is_empty());
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
fn split_cut_checks_exterior_level_against_prescribed_head() -> Result<()> {
    use soil_fem_3d::{Error, one_way_construction as construction};
    let (mut mesh, ground, input) = lifecycle_fixture(json!({"kind":"continuous"}))?;
    mesh.face_sets
        .insert("exterior_seam".into(), mesh.face_sets["Seam:a"].clone());
    let head = json!({"kind":"head","value":{"head_m":2.}});
    let pressure = json!({"kind":"pressure","pressure_head_m":1.5});
    let flux = json!({"kind":"outward_flux","flux_m_s":0.});
    for (condition, level, face_set, later, rejected) in [
        (&head, 3., "Seam:a", false, true),
        (&head, 2., "Seam:a", false, false),
        (&head, 3., "exterior_seam", false, true),
        (&head, 3., "exterior_seam", true, true),
        (&pressure, 3., "exterior_seam", false, true),
        (&pressure, 2., "exterior_seam", false, false),
        (&flux, 3., "exterior_seam", false, false),
    ] {
        let mut case = input.clone();
        let periods = case["one_way"]["flow"]["analysis"]["periods"]
            .as_array_mut()
            .unwrap();
        let cut_boundary = periods[2]["boundaries"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|b| b["face_set"] == "Seam:a")
            .unwrap();
        cut_boundary["condition"] = condition.clone();
        if later {
            let mut exposed = periods[2].clone();
            exposed["name"] = json!("still_exposed");
            exposed["end_time_s"] = json!(0.5);
            periods[3]["end_time_s"] = json!(0.625);
            periods.insert(3, exposed);
        }
        let exterior = case["one_way"]["exterior_periods"].as_array_mut().unwrap();
        let period = if later {
            exterior[2]["reservoirs"] = json!([{"face_set":face_set,"level_m":2.}]);
            exterior.insert(
                3,
                json!({
                    "name":"still_exposed",
                    "reservoirs":[{"face_set":face_set,"level_m":level}]
                }),
            );
            "still_exposed"
        } else {
            exterior[2]["reservoirs"] = json!([{"face_set":face_set,"level_m":level}]);
            "remove"
        };
        let result = construction::solve(&mesh, ground.clone(), serde_json::from_value(case)?);
        if rejected {
            let Err(Error::Input(reason)) = result else {
                panic!("expected preflight rejection for {period}, {face_set}: {condition}");
            };
            assert!(
                reason.contains(&format!(
                    "period {period}: exterior reservoir and prescribed hydraulic head disagree on cut face"
                )),
                "{reason}"
            );
        } else {
            let run = result?;
            assert!(run.complete, "{:?}", run.failure_reason);
            assert_eq!(run.excavations.len(), 1);
        }
    }
    Ok(())
}

/// Analytic saturated Tet4 RT0 mass matrix from the uniform tetrahedron covariance.
fn saturated_tet_outflow(
    vertices: [[f64; 3]; 4],
    conductivity: f64,
    head_drop: [f64; 4],
) -> [f64; 4] {
    use nalgebra::{Matrix4, Vector3, Vector4};
    let vertices = vertices.map(Vector3::from);
    let [a, b, c, d] = vertices;
    let volume = (b - a).cross(&(c - a)).dot(&(d - a)) / 6.;
    let centroid = (a + b + c + d) / 4.;
    let covariance_trace = vertices
        .iter()
        .map(|v| (v - centroid).norm_squared())
        .sum::<f64>()
        / 20.;
    let mass = Matrix4::from_fn(|i, j| {
        ((centroid - vertices[i]).dot(&(centroid - vertices[j])) + covariance_trace)
            / (9. * volume * conductivity)
    });
    mass.lu().solve(&Vector4::from(head_drop)).unwrap().into()
}

#[test]
fn independent_unit_tet_matches_birth_rt0_reference() {
    let result = saturated_tet_outflow(
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        1.,
        [-0.25, 0., 0., 0.],
    );
    for (q, expected) in result
        .into_iter()
        .zip([-4. / 3., 1. / 6., 1. / 6., 1. / 6.])
    {
        assert!((q - expected).abs() < 1e-14);
    }
}

fn verify_prescribed_birth(
    policy: Value,
    history: Option<Value>,
    expected_head: f64,
) -> Result<()> {
    use soil_fem_3d::one_way_construction as construction;
    let (mesh, ground, mut input) = lifecycle_fixture(policy)?;
    input["one_way"]["flow"]["materials"][0]["volume_source_per_s"] = json!(0.);
    if let Some(history) = history {
        input["one_way"]["joint_reservoir_heads"] = json!({"Seam":history});
    }
    let case: construction::Case = serde_json::from_value(input)?;
    let conductivity = case.one_way.flow.materials[0].saturated_conductivity_m_s[0][0];
    let theta = case.one_way.flow.materials[1]
        .retention
        .saturated_water_content;
    let run = construction::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let index = run
        .frames
        .iter()
        .position(|f| f.name == "place" && f.phase == "material_transition")
        .unwrap();
    let before = &run.frames[index - 1];
    let birth = &run.frames[index];
    let old = &before.state.hydraulic;
    let new = &birth.state.hydraulic;
    assert_eq!((old.time_s, new.time_s), (0.125, 0.125));
    assert!(birth.hydraulic_balance.is_none());
    assert!(
        birth
            .state
            .fluid
            .joint_traces
            .iter()
            .all(|p| p.time_s == new.time_s)
    );

    let mut owners = std::collections::BTreeSet::new();
    for trace in birth
        .state
        .fluid
        .joint_traces
        .iter()
        .filter_map(|p| p.side_a.as_ref())
    {
        assert_eq!(trace.total_head_m, expected_head);
        assert_eq!(trace.local_face, 0);
        let owner = trace.owner;
        let flow = new.cell_outflow_m3_s[owner].unwrap();
        assert_eq!(trace.outward_flow_m3_s, flow[trace.local_face]);
        if !owners.insert(owner) {
            continue;
        }
        assert_eq!(old.relative_cell_head_m[owner], Some(2.));
        assert_eq!(new.relative_cell_head_m[owner], Some(2.));
        assert_eq!(
            old.cell_storage_offset_m3[owner],
            new.cell_storage_offset_m3[owner]
        );
        assert_eq!(old.cell_outflow_m3_s[owner], Some([0.; 4]));
        let new_face = birth.hydraulic_faces[trace.face].corners;
        let old_face = before
            .hydraulic_faces
            .iter()
            .position(|f| f.corners == new_face)
            .unwrap();
        assert_eq!(old.relative_face_head_m[old_face], 2.);
        let expected = saturated_tet_outflow(
            mesh.tetrahedra[owner].map(|n| mesh.nodes[n]),
            conductivity,
            [2. - expected_head, 0., 0., 0.],
        );
        for (q, expected) in flow.into_iter().zip(expected) {
            assert!((q - expected).abs() < 1e-14, "{owner}: {flow:?}");
        }
        println!("head={expected_head} owner={owner} q={flow:?}");
    }
    assert_eq!(owners, [3, 4].into_iter().collect());

    let water = &run.transitions[0].water;
    assert_eq!(&water.state, new);
    assert!(water.internal_transfers.is_empty());
    assert!(water.changed_cells.iter().all(|c| c.before.is_none()));
    assert_eq!(water.water_removed_m3, 0.);
    let incoming = water.water_added_m3;
    let before_volume = before.state.fluid.water_equivalent_mass_volume_m3;
    let after_volume = birth.state.fluid.water_equivalent_mass_volume_m3;
    let change = after_volume - before_volume;
    let denominator = (incoming + change.abs()) / 2.;
    let error = change - incoming;
    assert!((incoming - 0.5 * theta).abs() <= 1e-12_f64.max(1e-7 * denominator));
    assert_eq!(before.cumulative_net_water_inflow_m3, 0.);
    let cumulative = after_volume - run.initial_water_m3 - incoming;
    for (step_r, cumulative_r, floor) in [(1e-7, 1e-8, 1e-12_f64), (1e-3, 5e-4, 1e-9)] {
        assert!(error.abs() <= floor.max(step_r * denominator));
        // The hold has zero source and flow, so its D and unsigned defect are zero.
        assert!(cumulative.abs().max(error.abs()) <= floor.max(cumulative_r * denominator));
    }
    for frame in run.frames[index..]
        .iter()
        .take_while(|f| f.state.hydraulic.time_s == new.time_s)
    {
        assert_eq!(
            frame.cumulative_net_water_inflow_m3,
            before.cumulative_net_water_inflow_m3
        );
        assert_eq!(frame.cumulative_added_water_m3, incoming);
        assert_eq!(&frame.state.hydraulic, new);
    }
    println!(
        "head={expected_head} birth_water={incoming} D={denominator} e={error} E={cumulative}"
    );

    let partial = construction::solve_until(&mesh, ground.clone(), case.clone(), index + 1)?;
    assert_eq!(&partial.last_committed.hydraulic, new);
    let saved = serde_json::to_vec(&construction::Checkpoint::seal(partial)?)?;
    let checkpoint: construction::Checkpoint = serde_json::from_slice(&saved)?;
    let resumed = checkpoint.resume(&mesh, ground, case, None)?;
    assert_eq!(
        construction::normalized(&run)?,
        construction::normalized(&resumed)?
    );
    Ok(())
}

#[test]
fn fixed_prescribed_birth_recomputes_all_retained_rt0_outflows() -> Result<()> {
    verify_prescribed_birth(json!({"kind":"prescribed","level_m":2.25}), None, 2.25)
}

#[test]
fn scheduled_prescribed_birth_recomputes_all_retained_rt0_outflows() -> Result<()> {
    verify_prescribed_birth(
        json!({"kind":"prescribed"}),
        Some(json!([[0., 2.125], [0.25, 2.375], [0.5, 2.25]])),
        2.25,
    )
}

#[test]
fn unchanged_prescribed_birth_keeps_zero_retained_rt0_outflows() -> Result<()> {
    verify_prescribed_birth(json!({"kind":"prescribed","level_m":2.}), None, 2.)
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
