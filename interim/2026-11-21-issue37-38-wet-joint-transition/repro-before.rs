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
fn issue_37_before_fix() -> Result<()> {
    use soil_fem_3d::one_way_construction as construction;
    let (mesh, ground, mut input) = lifecycle_fixture(json!({"kind":"continuous"}))?;
    let boundaries = input["one_way"]["flow"]["analysis"]["periods"][2]["boundaries"]
        .as_array_mut().unwrap();
    boundaries.iter_mut().find(|b| b["face_set"] == "Seam:a").unwrap()["condition"] =
        json!({"kind":"head","value":{"head_m":2.}});
    input["one_way"]["exterior_periods"][2]["reservoirs"] =
        json!([{"face_set":"Seam:a","level_m":3.}]);
    let run = construction::solve(&mesh, ground, serde_json::from_value(input)?)?;
    println!("complete={} frames={} excavations={} failure={:?}",
        run.complete, run.frames.len(), run.excavations.len(), run.failure_reason);
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.excavations.len(), 1);
    let f = run.frames.iter().find(|f| f.name == "remove" && f.hydraulic_balance.is_some()).unwrap();
    for (index, face) in f.hydraulic_faces.iter().enumerate()
        .filter(|(_, f)| [[4, 6, 7], [4, 5, 7]].contains(&f.corners)) {
        println!("t={} face={:?} head={} exterior_level=3", f.state.hydraulic.time_s,
            face.corners, f.state.hydraulic.relative_face_head_m[index]);
    }
    Ok(())
}

#[test]
fn issue_38_before_fix() -> Result<()> {
    use soil_fem_3d::one_way_construction as construction;
    let (mesh, ground, mut input) =
        lifecycle_fixture(json!({"kind":"prescribed","level_m":2.25}))?;
    input["one_way"]["flow"]["materials"][0]["volume_source_per_s"] = json!(0.);
    let run = construction::solve(&mesh, ground, serde_json::from_value(input)?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let index = run.frames.iter().position(|f| f.name == "place" && f.hydraulic_balance.is_none()).unwrap();
    let before = &run.frames[index - 1];
    let birth = &run.frames[index];
    let old = &before.state.hydraulic;
    let new = &birth.state.hydraulic;
    assert_eq!(old.time_s, new.time_s);
    for (owner, local_face) in [(3, 0), (4, 0)] {
        let before_face = before.hydraulic_faces.iter().position(|f| {
            f.owners.iter().any(|&[cell, j]| before.hydraulic_cell_element_indices[cell] == owner && j == local_face)
        }).unwrap();
        let trace = birth.state.fluid.joint_traces.iter().filter_map(|p| p.side_a.as_ref())
            .find(|t| t.owner == owner && t.local_face == local_face).unwrap();
        let old_head = old.relative_face_head_m[before_face];
        println!("t={} phase={} owner={owner} local={local_face} head_before={old_head} head_after={} q_before={:?} q_after={:?} trace_q={}",
            new.time_s, birth.phase, trace.total_head_m, old.cell_outflow_m3_s[owner],
            new.cell_outflow_m3_s[owner], trace.outward_flow_m3_s);
        assert_eq!(old_head, 2.);
        assert_eq!(trace.total_head_m, 2.25);
        assert_eq!(new.cell_outflow_m3_s[owner], old.cell_outflow_m3_s[owner]);
    }
    println!("water={}", serde_json::to_value(&run.transitions[0].water)?["water_added_m3"]);
    Ok(())
}
