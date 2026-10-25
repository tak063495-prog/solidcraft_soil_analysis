use serde_json::json;
use soil_fem_3d::{
    Result, checkpoint,
    ground::Ground,
    mesh::Mesh,
    one_way_construction::{self as lifecycle, Case},
    richards,
};

fn fixture(quadratic: bool) -> Result<(Mesh, Ground, Case)> {
    let mesh = serde_json::from_str(if quadratic {
        include_str!("../examples/seepage-column-tet10-mesh.json")
    } else {
        include_str!("../examples/seepage-column-tet4-mesh.json")
    })?;
    let mut ground: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-one-way-mc.json"))?;
    let mut extra = ground["materials"][1].clone();
    extra["id"] = json!(3);
    extra["name"] = json!("improved fill");
    extra["mohr_coulomb"]["young_kpa"] = json!(20000.);
    ground["materials"].as_array_mut().unwrap().push(extra);
    ground["construction_material_ids"] = json!([3]);
    let mut base: serde_json::Value =
        serde_json::from_str(include_str!("../examples/one-way-mc.json"))?;
    let mut material = base["flow"]["materials"][1].clone();
    material["id"] = json!(3);
    material["retention"]["saturated_water_content"] = json!(0.5);
    material["saturated_conductivity_m_s"] =
        json!([[5e-6, 0., 0.], [0., 5e-6, 0.], [0., 0., 5e-6]]);
    base["flow"]["materials"]
        .as_array_mut()
        .unwrap()
        .push(material);
    base["initial_porosity_tolerance"] = json!(0.15);
    let periods=(0..3).map(|i| {
        let mut b=vec![json!({"face_set":"bottom","condition":{"kind":"head","value":{"head_m":0.}}}),json!({"face_set":if i==0 {"floor"}else {"top"},"condition":{"kind":"atmospheric_no_ponding","rainfall_vertical_m_s":0.}})];
        for name in ["xmin","xmax","ymin","ymax"] {b.push(json!({"face_set":name,"condition":{"kind":"outward_flux","flux_m_s":0.}}));}
        json!({"name":(["hold","born","changed"][i]),"end_time_s":(i+1)as f64*100.,"boundaries":b})
    }).collect::<Vec<_>>();
    base["flow"]["analysis"] = json!({"kind":"transient","periods":periods});
    base["exterior_periods"] = json!([{"name":"hold","reservoirs":[]},{"name":"born","reservoirs":[]},{"name":"changed","reservoirs":[]}]);
    let fresh = json!({"effective_stress_compression_positive_kpa":[0.,0.,0.,0.,0.,0.],"plastic_strain_compression_positive":[0.,0.,0.,0.,0.,0.],"void_ratio":0.7,"equivalent_plastic_strain":0.,"plastic_work_kpa":0.,"last_return":"elastic"});
    let case = json!({"schema":"soil-fem-one-way-construction/1","one_way":base,"excavations":[],"initially_inactive_bodies":["UpperClay"],"transitions":[
        {"name":"place","before_period":"born","water_boundary_policy":"explicit_period_boundaries","surface_traction_policy":"explicit_period_reservoirs","action":{"kind":"birth","bodies":["UpperClay"],"states":{"UpperClay":fresh},"placement":"cad_new_nodes","initial_head":{"head_m":0.}}},
        {"name":"improve","before_period":"changed","water_boundary_policy":"explicit_period_boundaries","surface_traction_policy":"explicit_period_reservoirs","action":{"kind":"material_set","bodies":["UpperClay"],"material_id":3,"history":{"mode":"preserve"},"water":{"kind":"preserve_head"}}}
    ]});
    Ok((
        mesh,
        serde_json::from_value(ground)?,
        serde_json::from_value(case)?,
    ))
}
fn saturated(case: &mut Case) {
    case.one_way.flow.initial_head.head_m = 1.;
    let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
        unreachable!()
    };
    for p in periods {
        if let richards::Condition::Head { value } = &mut p.boundaries[0].condition {
            value.head_m = 1.;
        }
        p.boundaries[1].condition = richards::Condition::OutwardFlux { flux_m_s: 0. };
    }
    if let lifecycle::Action::Birth { initial_head, .. } = &mut case.transitions[0].action {
        initial_head.head_m = 1.;
    }
}

#[test]
fn time_zero_cut_keeps_clock_zero_and_nonzero_holding_checkpoint_resumes_exactly() -> Result<()> {
    let mesh: Mesh =
        serde_json::from_str(include_str!("../examples/seepage-column-tet10-mesh.json"))?;
    let ground: Ground = serde_json::from_str(include_str!("../examples/ground-one-way-mcc.json"))?;
    let mut case: Case =
        serde_json::from_str(include_str!("../examples/one-way-construction-mcc.json"))?;
    let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
        unreachable!()
    };
    periods.remove(0);
    case.one_way.exterior_periods.remove(0);
    let direct = lifecycle::solve(&mesh, ground.clone(), case.clone())?;
    assert!(direct.complete, "{:?}", direct.failure_reason);
    assert_eq!(direct.excavations[0].time_s, 0.);
    assert!(
        direct
            .frames
            .iter()
            .filter(|f| f.phase == "excavation_release")
            .all(|f| f.state.hydraulic.time_s == 0. && f.hydraulic_balance.is_none())
    );
    let partial = lifecycle::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    assert_eq!(partial.last_committed.hydraulic.time_s, 0.);
    assert!(
        partial
            .last_committed
            .mechanical
            .excavation_release_load_kn
            .iter()
            .any(|h| *h != 0.)
    );
    let resumed = lifecycle::Checkpoint::seal(partial)?.resume(&mesh, ground, case, None)?;
    assert_eq!(
        lifecycle::normalized(&direct)?,
        lifecycle::normalized(&resumed)?
    );
    Ok(())
}

#[test]
fn wet_stress_free_birth_matches_closed_saturated_oedometer_theory_and_added_water() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, mut case) = fixture(quadratic)?;
        saturated(&mut case);
        case.transitions.truncate(1);
        if let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis {
            periods.truncate(2);
        }
        case.one_way.exterior_periods.truncate(2);
        let run = lifecycle::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        let record = &run.transitions[0];
        assert_eq!(record.time_s, 100.);
        assert!(record.zero_release_relative_force_residual < 1e-12);
        assert!((record.water.water_added_m3 - 0.5 * 0.7 / 1.7).abs() < 1e-13);
        assert!(record.water.water_removed_m3 == 0.);
        let gamma = 16. + 10. * 0.7 / 1.7;
        let expected = -(48.75 + 0.375 * gamma) / 12000.;
        let mean = mesh.node_sets["top"]
            .iter()
            .map(|&n| run.last_committed.mechanical.displacement_m[n][2])
            .sum::<f64>()
            / mesh.node_sets["top"].len() as f64;
        assert!(
            (mean - expected).abs() < 2e-10,
            "Tet10={quadratic} {mean} vs {expected}"
        );
        for frame in run
            .frames
            .iter()
            .filter(|f| f.phase == "material_transition")
        {
            assert_eq!(frame.state.hydraulic.time_s, 100.);
            assert!(frame.hydraulic_balance.is_none());
        }
        assert!(
            run.frames
                .iter()
                .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-10)
        );
    }
    Ok(())
}

#[test]
fn material_set_preserve_head_uses_new_water_law_and_accounts_for_added_water() -> Result<()> {
    let (mesh, ground, case) = fixture(true)?;
    let run = lifecycle::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let record = &run.transitions[1];
    assert!(record.water.water_added_m3 > 0.);
    assert!(record.water.water_removed_m3 == 0.);
    for cell in &record.water.changed_cells {
        assert_eq!(
            cell.before.as_ref().unwrap().relative_head_m,
            cell.after.relative_head_m
        );
        assert_eq!(cell.new_material_id, 3);
        assert!(
            cell.after.water_equivalent_volume_m3
                > cell.before.as_ref().unwrap().water_equivalent_volume_m3
        );
    }
    assert!(
        run.frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-10)
    );
    assert!(
        run.last_committed
            .mechanical
            .element_material_ids
            .iter()
            .skip(6)
            .all(|id| *id == 3)
    );
    Ok(())
}

#[test]
fn material_set_preserve_water_inverts_head_and_retains_unchanged_primary_cells_exactly()
-> Result<()> {
    let (mesh, ground, mut case) = fixture(true)?;
    if let lifecycle::Action::MaterialSet { water, .. } = &mut case.transitions[1].action {
        *water = serde_json::from_value(
            json!({"kind":"preserve_water","reference_head":{"head_m":0.}}),
        )?;
    }
    let run = lifecycle::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let r = &run.transitions[1];
    assert_eq!(r.water.water_added_m3, 0.);
    assert_eq!(r.water.water_removed_m3, 0.);
    assert!(r.water.changed_cells.iter().all(|c| {
        (c.before.as_ref().unwrap().water_equivalent_volume_m3 - c.after.water_equivalent_volume_m3)
            .abs()
            < 1e-13
            && c.before.as_ref().unwrap().relative_head_m != c.after.relative_head_m
    }));
    let old = &run
        .frames
        .iter()
        .rfind(|f| f.phase == "transient" && f.state.hydraulic.time_s == 200.)
        .unwrap()
        .state
        .hydraulic;
    for e in 0..6 {
        assert_eq!(
            r.water.state.relative_cell_head_m[e],
            old.relative_cell_head_m[e]
        );
        assert_eq!(
            r.water.state.cell_storage_offset_m3[e],
            old.cell_storage_offset_m3[e]
        );
        assert_eq!(r.water.state.cell_outflow_m3_s[e], old.cell_outflow_m3_s[e]);
    }
    assert!(
        run.frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-10)
    );
    Ok(())
}

#[test]
fn lifecycle_rejections_restore_ledgers_and_same_clean_release_path_is_exact() -> Result<()> {
    let (mesh, ground, mut case) = fixture(false)?;
    case.verification_reject_once_transitions = vec!["place".into(), "improve".into()];
    let run = lifecycle::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(
        run.attempts
            .iter()
            .filter(|a| !a.accepted && a.candidate_changed)
            .count(),
        2
    );
    assert!(
        run.attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_pair_sha256_before == a.committed_pair_sha256_after)
    );
    case.verification_reject_once_transitions.clear();
    for event in &mut case.transitions {
        event.initial_increment = Some(0.25);
    }
    let clean = lifecycle::solve(&mesh, ground, case)?;
    assert!(clean.complete, "{:?}", clean.failure_reason);
    assert_eq!(
        serde_json::to_value(run.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?,
        serde_json::to_value(clean.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?
    );
    Ok(())
}

#[test]
fn lifecycle_checkpoints_restore_birth_reference_material_ids_porosity_and_mass_sources()
-> Result<()> {
    let (mesh, ground, case) = fixture(false)?;
    let direct = lifecycle::solve(&mesh, ground.clone(), case.clone())?;
    assert!(direct.complete);
    for count in [1, 2, 3, 4, 5, 6] {
        let partial = lifecycle::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        let saved = lifecycle::Checkpoint::seal(partial)?;
        if count == 5 {
            for mode in 0..3 {
                let mut v = serde_json::to_value(&saved)?;
                match mode {
                    0 => v["run"]["runtime"]["porosity_reference"][0] = json!(0.9),
                    1 => {
                        v["run"]["runtime"]["pending_transition"]["record"]["water"]["water_added_m3"] =
                            json!(999.)
                    }
                    _ => {
                        v["run"]["last_committed"]["mechanical"]["element_material_ids"][6] =
                            json!(2)
                    }
                }
                let mut forged: lifecycle::Checkpoint = serde_json::from_value(v)?;
                forged.run_sha256 = checkpoint::hash_json(&forged.run)?;
                assert!(
                    forged
                        .resume(&mesh, ground.clone(), case.clone(), None)
                        .is_err()
                );
            }
        }
        let resumed = saved.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert_eq!(
            lifecycle::normalized(&direct)?,
            lifecycle::normalized(&resumed)?
        );
    }
    Ok(())
}

#[test]
fn impossible_water_preserving_change_keeps_the_previous_accepted_state_and_archive() -> Result<()>
{
    let (mesh, ground, mut case) = fixture(false)?;
    case.one_way.flow.materials[2]
        .retention
        .residual_water_content = 0.45;
    if let lifecycle::Action::MaterialSet { water, .. } = &mut case.transitions[1].action {
        *water = serde_json::from_value(
            json!({"kind":"preserve_water","reference_head":{"head_m":0.}}),
        )?;
    }
    let run = lifecycle::solve(&mesh, ground, case)?;
    assert!(!run.complete);
    assert!(
        run.failure_reason
            .as_ref()
            .unwrap()
            .contains("finite retention head")
    );
    assert_eq!(run.last_committed.hydraulic.time_s, 200.);
    assert_eq!(run.transitions.len(), 1);
    assert_eq!(
        serde_json::to_value(&run.last_committed)?,
        serde_json::to_value(&run.frames.last().unwrap().state)?
    );
    assert_eq!(
        run.attempts.last().unwrap().committed_pair_sha256_before,
        run.attempts.last().unwrap().committed_pair_sha256_after
    );
    Ok(())
}

#[test]
fn zero_storage_saturated_change_requires_a_physical_pressure_reference() -> Result<()> {
    let (mesh, ground, mut case) = fixture(false)?;
    saturated(&mut case);
    case.one_way.flow.materials[2]
        .retention
        .saturated_water_content = 0.7 / 1.7;
    for head in [-1., 1.] {
        let mut input = case.clone();
        if let lifecycle::Action::MaterialSet { water, .. } = &mut input.transitions[1].action {
            *water = serde_json::from_value(
                json!({"kind":"preserve_water","reference_head":{"head_m":head}}),
            )?;
        }
        let run = lifecycle::solve(&mesh, ground.clone(), input)?;
        if head < 0. {
            assert!(!run.complete);
            assert!(run.failure_reason.unwrap().contains("pressure reference"));
        } else {
            assert!(run.complete, "{:?}", run.failure_reason);
            assert!(
                run.transitions[1]
                    .water
                    .changed_cells
                    .iter()
                    .all(|c| c.after.relative_head_m == 1.)
            );
        }
    }
    Ok(())
}

#[test]
fn excavated_mcc_can_be_replaced_while_inactive_and_reborn_as_wet_fresh_mc_without_losing_archive()
-> Result<()> {
    let (mesh, mc_ground, mut case) = fixture(true)?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-one-way-mcc.json"))?;
    for layer in &mut ground.layers {
        layer.consolidation =
            serde_json::from_value(json!({"definition":"mcc_pc_over_mean_effective","ratio":2.}))?;
    }
    ground.materials.push(mc_ground.materials[2].clone());
    ground.construction_material_ids = vec![3];
    case.initially_inactive_bodies.clear();
    case.one_way.flow.materials[0]
        .retention
        .saturated_water_content = 0.6;
    case.one_way.flow.materials[1]
        .retention
        .saturated_water_content = 0.6;
    case.one_way.flow.materials[2]
        .retention
        .saturated_water_content = 0.7 / 1.7;
    let mut v = serde_json::to_value(&case)?;
    v["excavations"] = json!([{"name":"cut","before_period":"hold","bodies":["UpperClay"],"water_boundary_policy":"explicit_period_boundaries","cut_traction_policy":"explicit_period_reservoirs"}]);
    let mut born = v["transitions"][0].clone();
    born["before_period"] = json!("changed");
    let fresh = born["action"]["states"].clone();
    v["transitions"] = json!([
        {"name":"replace_inactive","before_period":"born","water_boundary_policy":"explicit_period_boundaries","surface_traction_policy":"explicit_period_reservoirs","action":{"kind":"material_set","bodies":["UpperClay"],"material_id":3,"history":{"mode":"replace","states":fresh},"water":{"kind":"inactive"}}},born]);
    v["one_way"]["flow"]["analysis"]["periods"][1]["boundaries"][1]["face_set"] = json!("floor");
    let case: Case = serde_json::from_value(v)?;
    let direct = lifecycle::solve(&mesh, ground.clone(), case.clone())?;
    assert!(direct.complete, "{:?}", direct.failure_reason);
    assert_eq!(direct.transitions[0].water.water_added_m3, 0.);
    assert_eq!(direct.transitions[0].water.water_removed_m3, 0.);
    assert!(direct.transitions[0].water.changed_cells.is_empty());
    let original = &direct.initialization.last_committed.mechanical;
    assert!(!direct.excavations[0].removed_material_points.is_empty());
    for p in &direct.excavations[0].removed_material_points {
        assert!(p.state.as_mcc().is_some());
        assert_eq!(p.material_id, 2);
        assert_eq!(
            serde_json::to_value(&p.state)?,
            serde_json::to_value(&original.material_states[p.integration_point_index])?
        );
        assert!(matches!(
            direct.last_committed.mechanical.material_states[p.integration_point_index],
            soil_fem_3d::soil::State::MohrCoulomb(_)
        ));
    }
    assert!(
        direct
            .frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-9)
    );
    let partial = lifecycle::solve_until(&mesh, ground.clone(), case.clone(), 7)?;
    assert!(
        partial
            .runtime
            .as_ref()
            .unwrap()
            .pending_transition
            .is_some()
    );
    let resumed = lifecycle::Checkpoint::seal(partial)?.resume(&mesh, ground, case, None)?;
    assert_eq!(
        lifecycle::normalized(&direct)?,
        lifecycle::normalized(&resumed)?
    );
    Ok(())
}

#[test]
fn rejected_birth_with_genuine_plastic_history_replays_the_same_clean_path_exactly() -> Result<()> {
    let (mesh, mut ground, mut case) = fixture(true)?;
    for material in &mut ground.materials[1..] {
        if let soil_fem_3d::soil::Model::MohrCoulomb(p) = &mut material.model {
            p.cohesion_kpa = 2.;
            p.friction_deg = 10.;
        }
    }
    case.verification_reject_once_transitions = vec!["place".into(), "improve".into()];
    let run = lifecycle::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.attempts
            .iter()
            .any(|a| !a.accepted && a.name == "place" && a.plastic_candidate_changed)
    );
    assert!(
        run.attempts
            .iter()
            .filter(|a| !a.accepted)
            .all(|a| a.committed_pair_sha256_before == a.committed_pair_sha256_after)
    );
    case.verification_reject_once_transitions.clear();
    for event in &mut case.transitions {
        event.initial_increment = Some(0.25);
    }
    let clean = lifecycle::solve(&mesh, ground, case)?;
    assert!(clean.complete, "{:?}", clean.failure_reason);
    assert_eq!(
        serde_json::to_value(run.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?,
        serde_json::to_value(clean.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?
    );
    Ok(())
}

#[test]
fn birth_water_source_includes_declared_positive_pressure_storage() -> Result<()> {
    let (mesh, ground, mut case) = fixture(true)?;
    saturated(&mut case);
    case.transitions.truncate(1);
    for material in &mut case.one_way.flow.materials {
        material.retention.saturated_specific_storage_per_m = 0.1;
    }
    if let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis {
        periods.truncate(2);
    }
    case.one_way.exterior_periods.truncate(2);
    let run = lifecycle::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        (run.transitions[0].water.water_added_m3 - (0.5 * 0.7 / 1.7 + 0.1 * 0.125)).abs() < 1e-13
    );
    assert!(
        run.frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-10)
    );
    Ok(())
}

#[test]
fn new_born_surface_missing_water_law_and_inherited_plastic_history_are_rejected_before_work()
-> Result<()> {
    let (mesh, ground, case) = fixture(false)?;
    let mut missing = case.clone();
    if let richards::Analysis::Transient { periods } = &mut missing.one_way.flow.analysis {
        periods[1].boundaries.pop();
    }
    assert!(lifecycle::solve(&mesh, ground.clone(), missing).is_err());
    let mut missing = case.clone();
    missing.one_way.flow.materials.pop();
    assert!(lifecycle::solve(&mesh, ground.clone(), missing).is_err());
    let mut inherited = serde_json::to_value(&case)?;
    inherited["transitions"][0]["action"]["states"]["UpperClay"]["equivalent_plastic_strain"] =
        json!(0.1);
    assert!(lifecycle::solve(&mesh, ground, serde_json::from_value(inherited)?).is_err());
    Ok(())
}
