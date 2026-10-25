use nalgebra::{Matrix3, Vector3};
use serde_json::json;
use soil_fem_3d::{
    Result,
    anchors::{Configuration, Model},
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};

fn fixture(active: bool) -> Result<(Mesh, Model)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let head = Vector3::new(0.12, 0.17, 0.2);
    let tip = Vector3::new(0.8, 0.7, 0.9);
    let points: Vec<_> = [0.5, 0.75, 1.]
        .into_iter()
        .map(|t| json!({"position_m":<[f64;3]>::from(head+t*(tip-head))}))
        .collect();
    let config: Configuration = serde_json::from_value(
        json!({"schema":"soil-fem-anchors/1","embedding_tolerance_m":1e-10,
        "anchors":[{"id":"tie","young_kpa":2e8,"area_m2":0.001,"head":{"position_m":<[f64;3]>::from(head)},"bonded_path":points,"initially_active":active}]}),
    )?;
    let model = Model::new(&mesh, config)?;
    Ok((mesh, model))
}

#[test]
fn anchor_affine_tension_rigid_motion_and_distributed_transfer_preserve_work() -> Result<()> {
    let (mesh, model) = fixture(true)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let strain = 0.001;
    let u: Vec<[f64; 3]> = mesh.nodes.iter().map(|p| p.map(|x| strain * x)).collect();
    let e = model.evaluate(&model.initial_state(), &u, &active, None)?;
    let r = &e.report.anchors[0];
    let expected = 2e8 * 0.001 * strain;
    assert!(
        r.segment_forces_kn
            .iter()
            .all(|n| (n - expected).abs() < 1e-9)
    );
    assert!((r.free_length_m - r.bonded_length_m).abs() < 1e-12);
    for force in &r.point_soil_force_kn[1..r.point_soil_force_kn.len() - 1] {
        assert!(Vector3::from(*force).norm() < 1e-9);
    }
    assert!(Vector3::from(e.report.resultant_internal_force_kn).norm() < 1e-9);
    assert!(Vector3::from(e.report.resultant_internal_moment_kn_m).norm() < 1e-9);
    let work: f64 = e
        .internal_kn
        .iter()
        .enumerate()
        .map(|(g, f)| f * u[g / 3][g % 3])
        .sum();
    assert!((work - 2. * r.elastic_energy_kn_m).abs() < 1e-10);
    let rotation = Vector3::new(0.02, -0.03, 0.01);
    let shift = Vector3::new(0.03, 0.02, -0.01);
    let rigid: Vec<[f64; 3]> = mesh
        .nodes
        .iter()
        .map(|p| (shift + rotation.cross(&Vector3::from(*p))).into())
        .collect();
    let e = model.evaluate(&model.initial_state(), &rigid, &active, None)?;
    assert!(e.internal_kn.iter().all(|v| v.abs() < 1e-8));
    Ok(())
}

#[test]
fn jacking_lockoff_slack_and_reloading_follow_free_length_theory() -> Result<()> {
    let (mesh, model) = fixture(true)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let zero = vec![[0.; 3]; mesh.nodes.len()];
    let initial = model.initial_state();
    let names = vec!["tie".into()];
    let control = model.control(&initial, &names, 100., 0.5)?;
    let half = model.evaluate(&initial, &zero, &active, Some(&control))?;
    assert_eq!(half.state.anchors[0].force_control_kn, Some(50.));
    assert_eq!(half.state.anchors[0].segment_forces_kn[0], 50.);
    assert_eq!(half.state.anchors[0].segment_forces_kn[1], 0.);
    let lock = model.control(&initial, &names, 100., 1.)?;
    let locked = model.evaluate(&half.state, &zero, &active, Some(&lock))?;
    assert_eq!(locked.state.anchors[0].force_control_kn, None);
    assert_eq!(locked.state.anchors[0].last_lock_off_force_kn, 100.);
    let reevaluated = model.evaluate(&locked.state, &zero, &active, None)?;
    assert!((reevaluated.state.anchors[0].segment_forces_kn[0] - 100.).abs() < 1e-12);
    for strain in [0.0001, -0.001, 0.0002] {
        let u: Vec<_> = mesh.nodes.iter().map(|p| p.map(|x| strain * x)).collect();
        let actual = model.evaluate(&locked.state, &u, &active, None)?;
        let expected = (100. + 2e8 * 0.001 * strain).max(0.);
        assert!((actual.state.anchors[0].segment_forces_kn[0] - expected).abs() < 1e-9);
        if strain < 0. {
            assert!(!actual.state.anchors[0].segment_taut[0]);
            assert_eq!(actual.report.anchors[0].elastic_energy_kn_m, 0.);
        }
    }
    let control = model.control(&initial, &names, 0., 1.)?;
    let slack_u: Vec<_> = mesh.nodes.iter().map(|p| p.map(|x| -0.001 * x)).collect();
    let zero_control = model.control(&locked.state, &names, 0., 1.)?;
    let slack_lock = model.evaluate(&locked.state, &slack_u, &active, Some(&zero_control))?;
    assert_eq!(
        slack_lock.state.anchors[0].reference_extensions_m[0],
        locked.state.anchors[0].reference_extensions_m[0],
        "a zero-force jack must preserve existing slack instead of shortening the natural length"
    );
    assert!(!slack_lock.state.anchors[0].segment_taut[0]);
    let locked = model.evaluate(&initial, &zero, &active, Some(&control))?;
    assert_eq!(
        locked.state,
        model.evaluate(&locked.state, &zero, &active, None)?.state
    );
    Ok(())
}

#[test]
fn anchor_tangent_and_energy_match_independent_perturbations() -> Result<()> {
    let (mesh, model) = fixture(true)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let u: Vec<_> = mesh.nodes.iter().map(|p| p.map(|x| 0.001 * x)).collect();
    let initial = model.initial_state();
    let e = model.evaluate(&initial, &u, &active, None)?;
    let n = mesh.nodes.len() * 3;
    let matrix = Matrix3::new(0.4, -0.3, 0.2, 0.1, 0.2, -0.3, -0.2, 0.4, 0.3);
    let direction: Vec<[f64; 3]> = mesh
        .nodes
        .iter()
        .map(|p| (matrix * Vector3::from(*p)).into())
        .collect();
    let mut kd = vec![0.; n];
    for &(g, h, k) in &e.tangent {
        kd[g] += k * direction[h / 3][h % 3];
    }
    let h = 1e-7;
    let plus: Vec<_> = u
        .iter()
        .zip(&direction)
        .map(|(p, d)| std::array::from_fn(|a| p[a] + h * d[a]))
        .collect();
    let minus: Vec<_> = u
        .iter()
        .zip(&direction)
        .map(|(p, d)| std::array::from_fn(|a| p[a] - h * d[a]))
        .collect();
    let a = model.evaluate(&initial, &plus, &active, None)?;
    let b = model.evaluate(&initial, &minus, &active, None)?;
    let difference: f64 = kd
        .iter()
        .enumerate()
        .map(|(g, v)| (v - (a.internal_kn[g] - b.internal_kn[g]) / (2. * h)).powi(2))
        .sum();
    let norm: f64 = kd.iter().map(|v| v * v).sum();
    assert!((difference / norm).sqrt() < 1e-8);
    let work: f64 = e
        .internal_kn
        .iter()
        .enumerate()
        .map(|(g, f)| f * direction[g / 3][g % 3])
        .sum();
    let energy = (a.report.anchors[0].elastic_energy_kn_m
        - b.report.anchors[0].elastic_energy_kn_m)
        / (2. * h);
    assert!((energy - work).abs() < 1e-6);
    Ok(())
}

#[test]
fn installation_reference_and_lost_host_are_atomic_and_explicit() -> Result<()> {
    let (mesh, model) = fixture(false)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let u: Vec<_> = mesh.nodes.iter().map(|p| p.map(|x| 0.003 * x)).collect();
    let initial = model.initial_state();
    let names = vec!["tie".into()];
    let installed = model.change(&initial, &u, &names, true)?;
    let evaluated = model.evaluate(&installed, &u, &active, None)?;
    assert!(
        evaluated.state.anchors[0]
            .segment_forces_kn
            .iter()
            .all(|&f| f == 0.)
    );
    assert!(model.change(&installed, &u, &names, true).is_err());
    assert!(
        model
            .evaluate(&installed, &u, &vec![false; active.len()], None)
            .is_err()
    );
    assert!(!initial.anchors[0].active && installed.anchors[0].active);
    let removed = model.change(&installed, &u, &names, false)?;
    assert!(
        model
            .evaluate(&removed, &u, &vec![false; active.len()], None)?
            .report
            .anchors
            .is_empty()
    );
    assert_eq!(
        removed.anchors[0].reference_extensions_m,
        installed.anchors[0].reference_extensions_m
    );
    Ok(())
}

fn coupled(mcc: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let top = soil_fem_3d::boundary::surface_owners(&mesh)
        .into_iter()
        .filter(|(f, o)| o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][2] == 1.))
        .map(|(f, _)| f)
        .collect();
    mesh.face_sets.insert("top".into(), top);
    let mesh = mesh.elevate_straight()?;
    let mut ground: Ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-geostatic-mcc-column.json")
    } else {
        include_str!("../examples/ground-geostatic-mc-column.json")
    })?;
    if !mcc {
        for m in &mut ground.materials {
            if let soil::Model::MohrCoulomb(p) = &mut m.model {
                p.cohesion_kpa = 1000.;
                p.tension_cutoff_kpa = None;
            }
        }
    }
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.settings.initial_increment = 0.25;
    case.anchor_system = Some(serde_json::from_value(
        json!({"schema":"soil-fem-anchors/1","embedding_tolerance_m":1e-10,
        "anchors":[{"id":"tie","young_kpa":2e8,"area_m2":0.001,"head":{"position_m":[0.5,0.5,1.]},
        "bonded_path":[{"position_m":[0.5,0.5,0.3]},{"position_m":[0.5,0.5,0.15]},{"position_m":[0.5,0.5,0.05]}],"initially_active":false}]}),
    )?);
    case.stages = serde_json::from_value(json!([
        {"kind":"install_anchors","name":"install","anchors":["tie"]},
        {"kind":"jack_anchors","name":"jack","anchors":["tie"],"target_force_kn":2.},
        {"kind":"pressure","name":"preload","target_pressure_kpa":200.},
        {"kind":"pressure","name":"unload","target_pressure_kpa":0.},
        {"kind":"jack_anchors","name":"detension","anchors":["tie"],"target_force_kn":0.},
        {"kind":"remove_anchors","name":"remove","anchors":["tie"]}]))?;
    case.verification_reject_once_stages = vec!["jack".into(), "preload".into()];
    Ok((mesh, ground, case))
}

#[test]
fn anchor_jack_lock_slack_and_plastic_soil_restore_exactly_on_cutback_and_restart() -> Result<()> {
    for mcc in [false, true] {
        let (mesh, ground, case) = coupled(mcc)?;
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete, "{:?}", full.failure_reason);
        assert_eq!(full.schema, "soil-fem-stage-result/8");
        let locked = &full.frames[2].state.anchors.as_ref().unwrap().anchors[0];
        assert_eq!(locked.force_control_kn, None);
        assert!((locked.segment_forces_kn[0] - 2.).abs() < 1e-8);
        assert!(!full.frames[3].state.anchors.as_ref().unwrap().anchors[0].segment_taut[0]);
        for f in &full.frames {
            assert!(f.force_balance_error_kn.iter().all(|x| x.abs() < 1e-6));
            assert!(
                f.global_moment_balance_error_kn_m
                    .unwrap()
                    .iter()
                    .all(|x| x.abs() < 1e-6)
            );
        }
        assert!(full.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_structure_changed
            && a.rejected_trial_material_changed
            && a.committed_state_unchanged_on_rejection));
        if mcc {
            assert!(
                full.attempts
                    .iter()
                    .any(|a| !a.accepted && a.rejected_trial_plastic_changed)
            );
        }
        let mut clean = case.clone();
        clean.verification_reject_once_stages.clear();
        clean.stages = serde_json::from_value({
            let mut v = serde_json::to_value(&case.stages)?;
            for stage in v.as_array_mut().unwrap() {
                if ["jack", "preload"].contains(&stage["name"].as_str().unwrap()) {
                    stage["initial_increment"] = json!(0.125);
                }
            }
            v
        })?;
        let fresh = stages::solve(&mesh, ground.clone(), clean)?;
        assert!(fresh.complete);
        assert_eq!(
            serde_json::to_value(&full.last_committed)?,
            serde_json::to_value(&fresh.last_committed)?
        );
        let geo = full
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .count();
        let install = full
            .attempts
            .iter()
            .filter(|a| a.accepted && a.stage == "install")
            .count();
        let stop = geo + install + 1;
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), stop)?;
        assert!(
            partial.last_committed.anchors.as_ref().unwrap().anchors[0]
                .force_control_kn
                .is_some()
        );
        let cp = Checkpoint::seal(&mesh, partial)?;
        let cp: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&cp)?)?;
        let resumed = cp.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert_eq!(
            serde_json::to_value(&full.last_committed)?,
            serde_json::to_value(&resumed.last_committed)?
        );
        assert_eq!(
            serde_json::to_value(&full.frames)?,
            serde_json::to_value(&resumed.frames)?
        );
        let mut partial = stages::solve_until(&mesh, ground.clone(), case.clone(), stop)?;
        partial.last_committed.anchors.as_mut().unwrap().anchors[0].segment_forces_kn[0] += 1.;
        assert!(
            Checkpoint::seal(&mesh, partial)?
                .resume(&mesh, ground, case, None)
                .is_err()
        );
    }
    Ok(())
}

#[test]
fn excavation_with_retained_anchor_and_beam_releases_only_removed_soil_stress() -> Result<()> {
    let (mut mesh, ground, mut case) = coupled(false)?;
    let node = mesh.nodes.iter().position(|p| *p == [0., 0., 0.5]).unwrap();
    mesh.node_sets.insert("beam_base".into(), vec![node]);
    let anchor = &mut case.anchor_system.as_mut().unwrap().anchors[0];
    anchor.head.position_m = [0.5, 0.5, 0.4];
    case.beam_system = Some(serde_json::from_value(
        json!({"schema":"soil-fem-beams/1","rotation_length_scale_m":1.,"maximum_condensed_dofs":100,
        "nodes":[{"id":"base","position_m":[0.,0.,0.5],"soil_node_set":"beam_base"},{"id":"tip","position_m":[1.5,0.,0.5]}],
        "sections":[{"id":"s","young_kpa":2e8,"poisson":0.3,"area_m2":0.02,"iy_m4":0.0001,"iz_m4":0.0001,"torsion_j_m4":0.0001,"unit_weight_kn_m3":0.,"formulation":{"kind":"euler_bernoulli"}}],
        "members":[{"id":"beam","nodes":["base","tip"],"section_id":"s","local_y":[0,1,0],"initially_active":true}],
        "rotation_constraints":[{"node_id":"base","direction":[1,0,0],"value":0.},{"node_id":"base","direction":[0,1,0],"value":0.},{"node_id":"base","direction":[0,0,1],"value":0.}],
        "nodal_loads":[{"node_id":"tip","constant_force_kn":[0,0,-1]}]}),
    )?);
    case.stages.truncate(2);
    case.stages.push(serde_json::from_value(json!({"kind":"excavate","name":"cut_upper","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.}))?);
    case.verification_reject_once_stages = vec!["cut_upper".into()];
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.excavation_audits[0].zero_release_relative_force_residual < 1e-9);
    assert_eq!(
        run.frames
            .last()
            .unwrap()
            .beam_results
            .as_ref()
            .unwrap()
            .active_members,
        1
    );
    assert_eq!(
        run.frames
            .last()
            .unwrap()
            .anchor_results
            .as_ref()
            .unwrap()
            .anchors
            .len(),
        1
    );
    assert!(
        run.frames
            .last()
            .unwrap()
            .global_moment_balance_error_kn_m
            .unwrap()
            .iter()
            .all(|v| v.abs() < 1e-6)
    );
    Ok(())
}
