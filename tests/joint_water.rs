use serde_json::json;
use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    darcy,
    ground::Ground,
    joint_water, joints,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool, policy: serde_json::Value) -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(if quadratic {
        include_str!("../examples/joint-column-tet10-mesh.json")
    } else {
        include_str!("../examples/joint-column-mesh.json")
    })?;
    let mut g: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-joint-mc.json"))?;
    g["initialization"]["water_unit_weight_kn_m3"] = json!(10.);
    g["initialization"]["water_level_m"] = json!(0.75);
    let mut c: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-normal-column.json"))?;
    c["joint_system"]["joints"][0]["water"] = policy;
    c["verification_reject_once_stages"] = json!([]);
    c["stages"] = json!([{"kind":"hydrostatic_water","name":"rise","target_level_m":1.},{"kind":"pressure","name":"preload","target_pressure_kpa":50.},{"kind":"hydrostatic_water","name":"drawdown","target_level_m":0.5}]);
    Ok((mesh, serde_json::from_value(g)?, serde_json::from_value(c)?))
}
fn models(
    mesh: &Mesh,
    ground: &Ground,
    case: &StageCase,
) -> Result<(joints::Model, joint_water::Model)> {
    let contact = joints::Model::new(mesh, case.joint_system.clone().unwrap())?;
    let water = joint_water::Model::new(mesh, &contact, ground)?;
    Ok((contact, water))
}

#[test]
fn hydrostatic_contact_uses_effective_spring_pressure_and_total_equals_effective_plus_water()
-> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, case) = fixture(quadratic, json!({"kind":"continuous"}))?;
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(run.schema, "soil-fem-stage-result/13");
        for frame in &run.frames {
            let level = frame.state.hydrostatic_water_level_m.unwrap();
            let water = 10. * (level - 0.5).max(0.);
            let total = 110. + frame.state.surface_pressure_kpa;
            let effective = total - water;
            let report = frame.joint_results.as_ref().unwrap();
            let area: f64 = run
                .joint_integration_points
                .iter()
                .map(|p| p.weight_m2)
                .sum();
            let mean: f64 = report
                .points
                .iter()
                .zip(&run.joint_integration_points)
                .map(|(p, q)| p.normal_pressure_kpa * q.weight_m2)
                .sum::<f64>()
                / area;
            assert!(
                (mean - effective).abs() < 2e-5,
                "{} mean {} vs {}",
                frame.name,
                mean,
                effective
            );
            for p in &report.points {
                let w = p.hydraulic.as_ref().unwrap();
                assert!((w.physical_pore_pressure_kpa - water).abs() < 1e-12);
                assert!(
                    (w.total_normal_pressure_kpa - p.normal_pressure_kpa - water).abs() < 1e-12
                );
                assert!((w.friction_limit_kpa - 0.3 * p.normal_pressure_kpa).abs() < 1e-12);
                assert!((p.gap_m + p.normal_pressure_kpa / 1e6).abs() < 2e-11);
            }
            assert!(
                report
                    .water
                    .as_ref()
                    .unwrap()
                    .resultant_force_kn
                    .iter()
                    .all(|v| v.abs() < 1e-12)
            );
            assert!(frame.relative_force_residual < 1e-8);
        }
    }
    Ok(())
}

#[test]
fn pressure_survives_open_contact_and_does_not_enter_friction_twice() -> Result<()> {
    let (mesh, ground, case) = fixture(true, json!({"kind":"continuous"}))?;
    let (contact, water) = models(&mesh, &ground, &case)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let upper: std::collections::BTreeSet<_> = mesh
        .region_ids
        .iter()
        .enumerate()
        .filter(|(_, r)| **r == 2)
        .flat_map(|(e, _)| mesh.cell_nodes(e))
        .collect();
    let u: Vec<_> = (0..mesh.nodes.len())
        .map(|n| {
            if upper.contains(&n) {
                [0.001, -0.0007, 0.0003]
            } else {
                [0.; 3]
            }
        })
        .collect();
    let mut evaluated = contact.evaluate(&contact.initial_state()?, &u, &active, None, false)?;
    let hydraulic = water.hydrostatic(&active, Some(0.75))?;
    water.annotate(&mut evaluated.report, &hydraulic)?;
    for p in &evaluated.report.points {
        let w = p.hydraulic.as_ref().unwrap();
        assert!(p.opening_m > 0.);
        assert_eq!(p.normal_pressure_kpa, 0.);
        assert_eq!(w.physical_pore_pressure_kpa, 2.5);
        assert_eq!(w.total_normal_pressure_kpa, 2.5);
        assert_eq!(w.friction_limit_kpa, 0.);
    }
    let force: f64 = upper
        .iter()
        .map(|&n| hydraulic.nodal_load_kn[3 * n + 2])
        .sum();
    assert!((force - 2.5).abs() < 1e-12);
    assert!(evaluated.internal_kn.iter().all(|v| v.abs() < 1e-12));
    Ok(())
}

#[test]
fn partial_submergence_integrates_actual_sloped_tri3_and_tri6_force_and_moment() -> Result<()> {
    for quadratic in [false, true] {
        let (mut mesh, ground, case) = fixture(false, json!({"kind":"continuous"}))?;
        for p in &mut mesh.nodes {
            if (p[2] - 0.5).abs() < 1e-12 {
                p[2] = 0.2 + 0.6 * p[0];
            }
        }
        if quadratic {
            mesh = mesh.elevate_straight()?;
        }
        let (_, water) = models(&mesh, &ground, &case)?;
        let active = vec![true; mesh.tetrahedra.len()];
        let f = water.hydrostatic(&active, Some(0.5))?;
        let mut lower = std::collections::BTreeSet::new();
        for (e, &r) in mesh.region_ids.iter().enumerate() {
            if r == 1 {
                lower.extend(mesh.cell_nodes(e));
            }
        }
        let sum: [f64; 3] =
            std::array::from_fn(|a| lower.iter().map(|&n| f.nodal_load_kn[3 * n + a]).sum());
        assert!((sum[0] - 0.45).abs() < 2e-10, "{:?}", sum);
        assert!(sum[1].abs() < 1e-12);
        assert!((sum[2] + 0.75).abs() < 2e-10);
        assert!(
            f.resultant_force_kn
                .iter()
                .chain(&f.resultant_moment_kn_m)
                .all(|v| v.abs() < 1e-11)
        );
    }
    Ok(())
}

fn flow_case(mesh: &mut Mesh) -> Result<darcy::Case> {
    for (name, z) in [("bottom", 0.), ("top", 1.)] {
        let faces = soil_fem_3d::boundary::surface_owners(mesh)
            .into_iter()
            .filter(|(f, o)| {
                o.len() == 1 && f.iter().all(|&n| (mesh.nodes[n][2] - z).abs() < 1e-12)
            })
            .map(|(f, _)| f)
            .collect();
        mesh.face_sets.insert(name.into(), faces);
    }
    Ok(serde_json::from_value(
        json!({"schema":"soil-fem-darcy/1","units":"m-s","materials":[{"id":1,"hydraulic_conductivity_m_s":[[1e-4,0.,0.],[0.,1e-4,0.],[0.,0.,1e-4]]},{"id":2,"hydraulic_conductivity_m_s":[[2e-4,0.,0.],[0.,2e-4,0.],[0.,0.,2e-4]]}],"heads":[{"face_set":"bottom","head_m":2.},{"face_set":"top","head_m":3.}]}),
    )?)
}

#[test]
fn continuous_heads_pass_series_flow_through_split_displacement_nodes() -> Result<()> {
    for quadratic in [false, true] {
        let (mut mesh, ground, case) = fixture(quadratic, json!({"kind":"continuous"}))?;
        let flow = flow_case(&mut mesh)?;
        let (contact, water) = models(&mesh, &ground, &case)?;
        let original = mesh.clone();
        let active = vec![true; mesh.tetrahedra.len()];
        let result = water.solve_steady(&flow, &active, None, Some(2.5))?;
        assert!(
            result.solution.complete,
            "{:?}",
            result.solution.failure_reason
        );
        assert!(result.hydraulic_mesh.nodes.len() < mesh.nodes.len());
        assert_eq!(
            serde_json::to_value(&mesh)?,
            serde_json::to_value(original)?
        );
        for p in &result.solution.points {
            assert!((p.darcy_flux_m_s[2] + 1. / 7500.).abs() < 1e-14);
            assert!(p.darcy_flux_m_s[..2].iter().all(|v| v.abs() < 1e-14));
        }
        for p in contact.layout() {
            for (&a, &b) in p.nodes_a.iter().zip(&p.nodes_b) {
                assert_ne!(a, b);
                assert_eq!(
                    result.original_node_to_hydraulic_node[a],
                    result.original_node_to_hydraulic_node[b]
                );
                assert_eq!(
                    result.original_field.head_m[a],
                    result.original_field.head_m[b]
                );
            }
        }
        let pressure = water.from_heads(&active, &result.original_field, Some(2.5))?;
        for p in pressure.physical_pressure_kpa {
            assert!((p - 10. * (2. + 2. / 3. - 0.5)).abs() < 1e-10);
        }
    }
    Ok(())
}

#[test]
fn soil_structure_is_impermeable_and_prescribed_heads_are_not_averaged() -> Result<()> {
    let (mut mesh, ground, mut case) = fixture(false, json!({"kind":"soil","soil_side":"a"}))?;
    let mut flow = flow_case(&mut mesh)?;
    case.joint_system.as_mut().unwrap().structure_bodies = vec!["UpperClay".into()];
    flow.heads.truncate(1);
    flow.materials.truncate(1);
    let (_, water) = models(&mesh, &ground, &case)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let result = water.solve_steady(&flow, &active, None, Some(2.5))?;
    assert!(result.solution.complete);
    assert!(
        result
            .solution
            .points
            .iter()
            .all(|p| p.element_index < 6 && p.darcy_flux_m_s.iter().all(|v| v.abs() < 1e-14))
    );
    let pressure = water.from_heads(&active, &result.original_field, Some(2.5))?;
    assert!(
        pressure
            .physical_pressure_kpa
            .iter()
            .all(|p| (*p - 15.).abs() < 1e-10)
    );
    let (mut mesh, ground, case) = fixture(false, json!({"kind":"prescribed","level_m":2.5}))?;
    let flow = flow_case(&mut mesh)?;
    let (_, water) = models(&mesh, &ground, &case)?;
    let result = water.solve_steady(&flow, &active, None, None)?;
    assert!(result.solution.complete);
    let p = water.from_heads(&active, &result.original_field, None)?;
    assert!(
        p.physical_pressure_kpa
            .iter()
            .all(|p| (*p - 20.).abs() < 1e-10)
    );
    let mut conflict = flow.clone();
    conflict.heads.push(darcy::Head {
        face_set: "Seam:a".into(),
        head_m: 2.6,
        gradient: [0.; 3],
        reference_point_m: [0.; 3],
    });
    assert!(water.solve_steady(&conflict, &active, None, None).is_err());
    assert!(water.hydrostatic(&active, Some(1.)).is_err());
    let mut negative = result.original_field.clone();
    negative.head_m.fill(0.);
    negative.relative_head_m.fill(0.);
    negative.datum_head_m.fill(0.);
    assert!(water.from_heads(&active, &negative, None).is_err());
    let (mesh, mut ground, mut case) = fixture(false, json!({"kind":"soil","soil_side":"a"}))?;
    case.joint_system.as_mut().unwrap().structure_bodies = vec!["UpperClay".into()];
    let mut mcc: Ground = serde_json::from_str(include_str!("../examples/ground-joint-mcc.json"))?;
    let mut replacement = mcc.materials.pop().unwrap();
    replacement.id = 3;
    ground.materials.push(replacement);
    ground.construction_material_ids = vec![3];
    let seed: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-joint-mcc.json"))?;
    case.stages = serde_json::from_value(
        json!([{"kind":"material_set","name":"wrong_structure_model","bodies":["UpperClay"],"material_id":3,"history":{"mode":"replace","states":{"UpperClay":seed["initialization"]["geostatic"]["seed"]["states"]["UpperClay"]}}}]),
    )?;
    let error = stages::solve(&mesh, ground, case)
        .expect_err("soil hardening cannot be a total-stress structure");
    assert!(
        error
            .to_string()
            .contains("impermeable total-stress structure"),
        "{error}"
    );
    Ok(())
}

#[test]
fn wet_stage_cutback_restores_pressure_contact_history_and_checkpoint() -> Result<()> {
    let (mesh, ground, mut case) = fixture(false, json!({"kind":"continuous"}))?;
    case.verification_reject_once_stages = vec!["rise".into()];
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.attempts.iter().any(|a| !a.accepted
        && a.rejected_trial_water_changed
        && a.committed_state_unchanged_on_rejection));
    let mut clean = case.clone();
    clean.verification_reject_once_stages.clear();
    let mut steps = serde_json::to_value(&clean.stages)?;
    steps[0]["initial_increment"] = json!(0.125);
    clean.stages = serde_json::from_value(steps)?;
    let fresh = stages::solve(&mesh, ground.clone(), clean)?;
    assert!(fresh.complete);
    assert_eq!(
        serde_json::to_value(&run.last_committed)?,
        serde_json::to_value(&fresh.last_committed)?
    );
    let geo = run
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), geo + 1)?;
    assert!(partial.last_committed.hydrostatic_water_level_m.unwrap() > 0.75);
    let cp = Checkpoint::seal(&mesh, partial)?;
    let resumed = cp.resume(&mesh, ground.clone(), case.clone(), None)?;
    assert_eq!(
        serde_json::to_value(&run.frames)?,
        serde_json::to_value(&resumed.frames)?
    );
    let mut wrong = case;
    wrong.joint_system.as_mut().unwrap().joints[0].water = None;
    assert!(stages::solve(&mesh, ground, wrong).is_err());
    Ok(())
}

#[test]
fn staged_darcy_interface_policies_preserve_primary_pressure_and_exact_midstage_restart()
-> Result<()> {
    for quadratic in [false, true] {
        for policy in [
            json!({"kind":"continuous"}),
            json!({"kind":"soil","soil_side":"a"}),
            json!({"kind":"prescribed","level_m":2.5}),
        ] {
            let (mut mesh, mut ground, mut case) = fixture(quadratic, policy.clone())?;
            let mut flow = flow_case(&mut mesh)?;
            if policy["kind"] == "soil" {
                case.joint_system.as_mut().unwrap().structure_bodies = vec!["UpperClay".into()];
                flow.heads.truncate(1);
                flow.materials.truncate(1);
            }
            if policy["kind"] == "prescribed" {
                ground.initialization.water_level_m = Some(2.5);
            }
            case.stages = serde_json::from_value(
                json!([{"kind":"steady_seepage","name":"flow","flow":flow,"exterior_water_level_m":0.75},{"kind":"pressure","name":"preload","target_pressure_kpa":20.}]),
            )?;
            case.verification_reject_once_stages = vec!["flow".into()];
            let run = stages::solve(&mesh, ground.clone(), case.clone())?;
            assert!(run.complete, "{:?}", run.failure_reason);
            let expected = if policy["kind"] == "continuous" {
                10. * (2. + 2. / 3. - 0.5)
            } else if policy["kind"] == "soil" {
                15.
            } else {
                20.
            };
            let report = run.frames[1].joint_results.as_ref().unwrap();
            let mean: f64 = report
                .points
                .iter()
                .zip(&run.joint_integration_points)
                .map(|(p, q)| p.normal_pressure_kpa * q.weight_m2)
                .sum();
            assert!(
                (mean - (110. - expected)).abs() < 2e-5,
                "{:?} {}",
                policy,
                mean
            );
            for p in &report.points {
                assert!(
                    (p.hydraulic.as_ref().unwrap().physical_pore_pressure_kpa - expected).abs()
                        < 1e-9
                );
            }
            assert!(run.flow_audits[0].joint_flow.is_some());
            if policy["kind"] == "soil" {
                for (p, &v) in run
                    .integration_points
                    .iter()
                    .zip(run.last_committed.pore_pressure_kpa.as_ref().unwrap())
                {
                    if mesh.region_ids[p.element_index] == 2 {
                        assert_eq!(v, 0.);
                    }
                }
            }
            let geo = run
                .geostatic
                .as_ref()
                .unwrap()
                .attempts
                .iter()
                .filter(|a| a.accepted)
                .count();
            let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), geo + 1)?;
            assert!(partial.last_committed.joint_water.is_some());
            let encoded = serde_json::to_vec(&partial)?;
            let cp = Checkpoint::seal(&mesh, serde_json::from_slice(&encoded)?)?;
            let resumed = cp.resume(&mesh, ground.clone(), case.clone(), None)?;
            assert_eq!(
                serde_json::to_value(&resumed.frames)?,
                serde_json::to_value(&run.frames)?
            );
            let mut bad: stages::Run = serde_json::from_slice(&encoded)?;
            bad.last_committed
                .joint_water
                .as_mut()
                .unwrap()
                .nodal_load_kn[0] += 1.;
            assert!(
                Checkpoint::seal(&mesh, bad)?
                    .resume(&mesh, ground.clone(), case.clone(), None)
                    .is_err()
            );
            let mut bad = partial;
            bad.flow_audits[0]
                .joint_flow
                .as_mut()
                .unwrap()
                .original_node_to_hydraulic_node[0] += 1;
            assert!(
                Checkpoint::seal(&mesh, bad)?
                    .resume(&mesh, ground, case, None)
                    .is_err()
            );
        }
    }
    Ok(())
}

#[test]
fn death_after_wet_seepage_releases_interface_water_and_archives_contact_history() -> Result<()> {
    let (mut mesh, ground, mut case) = fixture(true, json!({"kind":"continuous"}))?;
    let flow = flow_case(&mut mesh)?;
    case.stages = serde_json::from_value(
        json!([{"kind":"steady_seepage","name":"flow","flow":flow,"exterior_water_level_m":0.75},{"kind":"excavate","name":"cut","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.}]),
    )?;
    case.verification_reject_once_stages = vec!["cut".into()];
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
    assert_eq!(run.frames[1].state.joints, run.frames[2].state.joints);
    assert!(
        run.frames[2]
            .state
            .joint_water
            .as_ref()
            .unwrap()
            .nodal_load_kn
            .iter()
            .all(|v| *v == 0.)
    );
    assert!(
        run.frames[2]
            .joint_results
            .as_ref()
            .unwrap()
            .points
            .iter()
            .all(|p| !p.active && p.hydraulic.as_ref().unwrap().physical_pore_pressure_kpa == 0.)
    );
    let geo = run
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let first = run
        .attempts
        .iter()
        .filter(|a| a.accepted && a.stage == "flow")
        .count();
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), geo + first + 1)?;
    assert!(
        partial
            .last_committed
            .excavation_release_load_kn
            .iter()
            .any(|v| v.abs() > 1e-8)
    );
    let resumed = Checkpoint::seal(&mesh, partial)?.resume(&mesh, ground, case, None)?;
    assert_eq!(
        serde_json::to_value(&resumed.frames)?,
        serde_json::to_value(&run.frames)?
    );
    Ok(())
}

#[test]
fn wet_global_sliding_opening_recontact_preserves_water_load_and_rejected_slip_history()
-> Result<()> {
    let (mesh, ground, _) = fixture(false, json!({"kind":"continuous"}))?;
    let mut c: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
    c["joint_system"]["joints"][0]["water"] = json!({"kind":"continuous"});
    c["verification_reject_once_stages"] = json!(["shear", "open", "reverse"]);
    let case: StageCase = serde_json::from_value(c)?;
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.attempts.iter().any(|a| !a.accepted
        && a.rejected_trial_joint_slip_changed
        && a.committed_state_unchanged_on_rejection));
    let open = run.frames.iter().find(|f| f.name == "open").unwrap();
    assert!(
        open.joint_results
            .as_ref()
            .unwrap()
            .points
            .iter()
            .all(|p| p.opening_m > 1e-6
                && p.normal_pressure_kpa == 0.
                && p.hydraulic.as_ref().unwrap().physical_pore_pressure_kpa == 2.5)
    );
    assert!(
        open.state
            .joint_water
            .as_ref()
            .unwrap()
            .nodal_load_kn
            .iter()
            .any(|v| v.abs() > 1e-3)
    );
    let mut clean = case.clone();
    clean.verification_reject_once_stages.clear();
    let mut steps = serde_json::to_value(&clean.stages)?;
    for s in steps.as_array_mut().unwrap() {
        if ["shear", "open", "reverse"].contains(&s["name"].as_str().unwrap()) {
            s["initial_increment"] = json!(
                s.get("initial_increment")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(clean.settings.initial_increment)
                    / 2.
            );
        }
    }
    clean.stages = serde_json::from_value(steps)?;
    let direct = stages::solve(&mesh, ground, clean)?;
    assert!(direct.complete);
    assert_eq!(
        serde_json::to_value(&run.last_committed)?,
        serde_json::to_value(&direct.last_committed)?
    );
    Ok(())
}
