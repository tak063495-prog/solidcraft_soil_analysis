use serde_json::{Value, json};
use soil_fem_3d::{
    Result, ground::Ground, hydrostatic, hydrostatic_volume, mcc, mesh::Mesh,
    soil_birth::Formation, stages,
};

fn fixture(init: Value, level: f64) -> Result<(Mesh, Ground, Value)> {
    let mesh = serde_json::from_str(include_str!(
        "../examples/joint-water-column-tet10-mesh.json"
    ))?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-joint-water-mc.json"))?;
    ground.initialization.water_level_m = Some(level);
    let mut case: Value =
        serde_json::from_str(include_str!("../examples/stages-joint-water-column.json"))?;
    case["joint_system"]["joints"][0]["birth_policy"] = json!("stress_free");
    case["initially_inactive_bodies"] = json!(["UpperClay"]);
    let seed = json!({"effective_stress_compression_positive_kpa":[5.,5.,5.,0.,0.,0.],"plastic_strain_compression_positive":[0.,0.,0.,0.,0.,0.],"void_ratio":0.7,"equivalent_plastic_strain":0.,"plastic_work_kpa":0.,"last_return":"elastic"});
    case["stages"] = json!([{"kind":"birth","name":"fill","bodies":["UpperClay"],"placement":"cad_new_nodes","states":{"UpperClay":seed},"hydraulic_initialization":init}]);
    Ok((mesh, ground, case))
}
fn placement() -> Value {
    json!({"material_id":2,"water_source":{"reservoir_id":"imported_soil"}})
}

#[test]
fn post_darcy_keeps_old_heads_and_initializes_only_new_points() -> Result<()> {
    let init = json!({"kind":"post_darcy","preceding_stage":"flow","initial_head":{"head_m":1.5,"gradient":[0.,0.,0.2],"reference_point_m":[0.,0.,0.5]},"placement_water":placement()});
    let (mesh, ground, mut case) = fixture(init, 1.5)?;
    let flow: Value =
        serde_json::from_str(include_str!("../examples/stages-joint-water-flow.json"))?;
    let mut stage = flow["stages"][0].clone();
    stage["flow"]["heads"] =
        json!([{"face_set":"bottom","head_m":1.6},{"face_set":"Seam:a","head_m":1.5}]);
    stage["exterior_water_level_m"] = json!(1.5);
    case["stages"].as_array_mut().unwrap().insert(0, stage);
    let run = stages::solve(&mesh, ground.clone(), serde_json::from_value(case.clone())?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let old = &run.frames[1].state;
    let placed = &run.transition_audits[0].prepared_state;
    for (i, q) in run.integration_points.iter().enumerate() {
        let pressure = placed.pore_pressure_kpa.as_ref().unwrap()[i];
        if old.active_elements[q.element_index] {
            assert_eq!(pressure, old.pore_pressure_kpa.as_ref().unwrap()[i]);
        } else {
            let z = q.position_m[2];
            assert!((pressure - 10. * (1.5 + 0.2 * (z - 0.5) - z)).abs() < 1e-12);
            assert_eq!(
                placed.material_states[i].stress(),
                &[5., 5., 5., 0., 0., 0.]
            );
        }
    }
    case["stages"][1]["hydraulic_initialization"]["preceding_stage"] = json!("other");
    assert!(stages::solve(&mesh, ground, serde_json::from_value(case)?).is_err());
    Ok(())
}

#[test]
fn two_region_column_reports_separate_placement_inventories_and_weights() -> Result<()> {
    let n = 0.7 / 1.7;
    let init = json!({"kind":"two_region","level_m":0.75,"mixture":{"porosity":n,"upper_saturation":0.25,"gravity_m_s2":10.,"particle_density_kg_m3":2700.,"water_density_kg_m3":1000.},"placement_water":placement()});
    let (mesh, ground, case) = fixture(init, 0.75)?;
    let run = stages::solve(&mesh, ground, serde_json::from_value(case)?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let cells = &run.transition_audits[0].placed_water;
    for (actual, expected) in [
        (cells.iter().map(|c| c.lower_volume_m3).sum::<f64>(), 0.25),
        (cells.iter().map(|c| c.upper_volume_m3).sum(), 0.25),
        (cells.iter().map(|c| c.lower_water_m3).sum(), 0.25 * n),
        (
            cells.iter().map(|c| c.upper_water_m3).sum(),
            0.25 * n * 0.25,
        ),
        (cells.iter().map(|c| c.lower_weight_kn).sum(), 5.),
        (
            cells.iter().map(|c| c.upper_weight_kn).sum(),
            0.25 * (20. - 7.5 * n),
        ),
    ] {
        assert!((actual - expected).abs() < 1e-12, "{actual} {expected}");
    }
    Ok(())
}

#[test]
fn split_tetrahedron_and_face_match_volume_and_pressure_moment() -> Result<()> {
    let points = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    for level in [-0.1_f64, 0., 0.25, 0.75, 1., 1.1] {
        let h = level.clamp(0., 1.);
        let w = hydrostatic_volume::wet_shape_volumes(&points, level)?;
        assert!((w.iter().sum::<f64>() - (1. - (1. - h).powi(3)) / 6.).abs() < 1e-15);
        let face = [points[0], points[1], points[3]];
        let f = hydrostatic::face_load(&face, level, 10., Default::default())?;
        let total: f64 = f.nodal_force_kn.iter().map(|f| f[1]).sum();
        let expected = 10. * (level * (h - h * h / 2.) - h * h / 2. + h.powi(3) / 3.);
        assert!((total - expected).abs() < 1e-12);
        let moment =
            10. * (level * (h * h / 2. - h.powi(3) / 3.) - h.powi(3) / 3. + h.powi(4) / 4.);
        assert!((f.nodal_force_kn[2][1] - moment).abs() < 1e-12);
    }
    Ok(())
}

#[test]
fn mcc_formation_requires_consistent_e_ocr_and_admissible_stress() -> Result<()> {
    let ground: Value =
        serde_json::from_str(include_str!("../examples/ground-joint-water-mcc.json"))?;
    let p: mcc::Parameters = serde_json::from_value(ground["materials"][1]["mcc"].clone())?;
    let seed: mcc::State = serde_json::from_value(
        ground["initialization"]["geostatic"]["seed"]["states"]["UpperClay"].clone(),
    )?;
    let mut f = Formation {
        history: "synthetic drained OC calibration: isotropic 200 to 100 kPa".into(),
        hydraulic_porosity: seed.void_ratio / (1. + seed.void_ratio),
        pc_over_mean_effective: Some(2.),
    };
    f.validate(&p, &seed)?;
    let mut bad = seed.clone();
    bad.void_ratio += 0.1;
    assert!(f.validate(&p, &bad).is_err());
    bad = seed.clone();
    bad.effective_stress_compression_positive_kpa[3] = 200.;
    assert!(f.validate(&p, &bad).is_err());
    f.pc_over_mean_effective = Some(1.);
    assert!(f.validate(&p, &seed).is_err());
    Ok(())
}

#[test]
fn saturated_mcc_places_the_supplied_formation_state_before_release() -> Result<()> {
    let template: Value =
        serde_json::from_str(include_str!("../examples/ground-joint-water-mcc.json"))?;
    let seed = &template["initialization"]["geostatic"]["seed"]["states"]["UpperClay"];
    let ratio = seed["void_ratio"].as_f64().unwrap();
    let init = json!({"kind":"saturated_mcc","level_m":1.5,"placement_water":placement(),
        "formation":{"history":"Synthetic calibration: drained isotropic compression to 200 kPa, unload to 100 kPa",
        "hydraulic_porosity":ratio/(1.+ratio),"pc_over_mean_effective":2.}});
    let (mesh, ground, mut case) = fixture(init, 1.5)?;
    let mut g = serde_json::to_value(ground)?;
    g["materials"][1] = template["materials"][1].clone();
    g["initialization"]["geostatic"]["seed"]["states"]["UpperClay"] = seed.clone();
    case["stages"][0]["states"]["UpperClay"] = seed.clone();
    let run = stages::solve(
        &mesh,
        serde_json::from_value(g)?,
        serde_json::from_value(case)?,
    )?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let placed = &run.transition_audits[0];
    for (i, q) in run.integration_points.iter().enumerate() {
        if placed.selected_elements.contains(&q.element_index) {
            let state = placed.prepared_state.material_states[i].as_mcc().unwrap();
            assert_eq!(state.preconsolidation_pressure_kpa, 200.);
            assert_eq!(state.void_ratio, ratio);
        }
    }
    Ok(())
}
