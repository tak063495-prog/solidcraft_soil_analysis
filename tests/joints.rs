use nalgebra::{DVector, Vector3};
use soil_fem_3d::{
    Result, boundary,
    checkpoint::Checkpoint,
    ground::Ground,
    joints,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool) -> Result<(Mesh, joints::Model)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/joint-column-mesh.json"))?;
    let mesh = if quadratic {
        mesh.elevate_straight()?
    } else {
        mesh
    };
    let case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
    let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
    Ok((mesh, model))
}

#[test]
fn opposing_tri3_tri6_transfer_three_direction_traction_work_and_moment() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, model) = fixture(quadratic)?;
        let upper: std::collections::BTreeSet<_> = mesh
            .region_ids
            .iter()
            .enumerate()
            .filter(|(_, r)| **r == 2)
            .flat_map(|(ei, _)| mesh.cell_nodes(ei))
            .collect();
        let u: Vec<_> = mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let rigid = [0.003 - 0.002 * p[1], -0.004 + 0.002 * p[0], 0.001];
                if upper.contains(&i) {
                    std::array::from_fn(|a| rigid[a] + [0.00005, 0.00004, -0.0001][a])
                } else {
                    rigid
                }
            })
            .collect();
        let active = vec![true; mesh.tetrahedra.len()];
        let evaluated = model.evaluate(&model.initial_state()?, &u, &active, None, true)?;
        assert!((model.layout().iter().map(|p| p.weight_m2).sum::<f64>() - 1.).abs() < 1e-12);
        assert!(
            evaluated
                .report
                .resultant_internal_force_kn
                .iter()
                .all(|x| x.abs() < 1e-10)
        );
        assert!(
            evaluated
                .report
                .resultant_internal_moment_kn_m
                .iter()
                .all(|x| x.abs() < 1e-10)
        );
        let work: f64 = evaluated
            .internal_kn
            .iter()
            .zip(u.iter().flatten())
            .map(|(f, u)| f * u)
            .sum();
        let expected = 1e6 * 1e-8 + 2e5 * (0.00005_f64.powi(2) + 0.00004_f64.powi(2));
        assert!((work - expected).abs() < 1e-12);
        assert!((evaluated.report.elastic_energy_kn_m - expected / 2.).abs() < 1e-12);
        for p in &evaluated.report.points {
            assert!((p.normal_pressure_kpa - 100.).abs() < 1e-9);
        }
        model.validate_state(&evaluated.state, &u, &active)?;
        let archived =
            model.evaluate(&evaluated.state, &u, &vec![false; active.len()], None, true)?;
        assert_eq!(archived.state, evaluated.state);
        assert!(archived.internal_kn.iter().all(|x| *x == 0.));
        assert!(archived.tangent.is_empty());
    }
    Ok(())
}

#[test]
fn interface_assembly_tangent_matches_independent_direction_and_rejects_bad_layouts() -> Result<()>
{
    let (mesh, model) = fixture(true)?;
    let mut u = vec![[0.; 3]; mesh.nodes.len()];
    for &n in &mesh.node_sets["Seam:b"] {
        u[n] = [0.001, 0.0004, -0.0001];
    }
    let initial = model.initial_state()?;
    let active = vec![true; mesh.tetrahedra.len()];
    let evaluation = model.evaluate(&initial, &u, &active, None, true)?;
    let direction: Vec<[f64; 3]> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, p)| [0.2 + (i as f64).sin(), -0.1 + p[0], 0.3 - p[1]])
        .collect();
    let h = 2e-8;
    let mut plus = u.clone();
    let mut minus = u.clone();
    for n in 0..u.len() {
        for a in 0..3 {
            plus[n][a] += h * direction[n][a];
            minus[n][a] -= h * direction[n][a];
        }
    }
    let p = model.evaluate(&initial, &plus, &active, None, false)?;
    let m = model.evaluate(&initial, &minus, &active, None, false)?;
    let fd = DVector::from_iterator(
        3 * u.len(),
        p.internal_kn
            .iter()
            .zip(&m.internal_kn)
            .map(|(a, b)| (a - b) / (2. * h)),
    );
    let mut predicted = DVector::zeros(fd.len());
    for (g, h, k) in evaluation.tangent {
        predicted[g] += k * direction[h / 3][h % 3];
    }
    assert!((&predicted - &fd).norm() / fd.norm() < 1e-7);
    assert!(model.evaluate(&initial, &u, &[], None, false).is_err());
    let mut bad = initial.clone();
    bad.points[0].normal = [1., 0., 0.];
    assert!(model.validate_state(&bad, &u, &active).is_err());
    let mut bad_mesh = mesh.clone();
    bad_mesh
        .face_sets
        .insert("Seam:b".into(), bad_mesh.face_sets["Seam:a"].clone());
    assert!(joints::Model::new(&bad_mesh, model.configuration.clone()).is_err());
    // Shared DOFs are not a contact interface.
    let uncut: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    assert!(joints::Model::new(&uncut, model.configuration.clone()).is_err());
    Ok(())
}

fn coupled(mcc: bool) -> Result<(Mesh, Ground, StageCase)> {
    let (mesh, _) = fixture(true)?;
    let ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-joint-mcc.json")
    } else {
        include_str!("../examples/ground-joint-mc.json")
    })?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
    case.verification_reject_once_stages = vec!["shear".into(), "open".into(), "reverse".into()];
    Ok((mesh, ground, case))
}

#[test]
fn global_joint_closure_sliding_opening_recontact_cutback_and_restart() -> Result<()> {
    for mcc in [false, true] {
        let (mesh, ground, case) = coupled(mcc)?;
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete, "mcc={mcc}: {:?}", full.failure_reason);
        eprintln!(
            "mcc={mcc}: {} attempts; rejected {:?}",
            full.attempts.len(),
            full.attempts
                .iter()
                .filter(|a| !a.accepted)
                .map(|a| (
                    &a.stage,
                    a.progress_from,
                    a.progress_to,
                    a.rejected_trial_joint_slip_changed,
                    &a.reason
                ))
                .collect::<Vec<_>>()
        );
        assert_eq!(full.schema, "soil-fem-stage-result/9");
        let slip = full
            .frames
            .iter()
            .flat_map(|f| &f.state.joints.as_ref().unwrap().points)
            .map(|p| p.accumulated_plastic_slip_m)
            .fold(0., f64::max);
        assert!(slip > 1e-5, "mcc={mcc}: slip={slip}");
        let opened = &full
            .frames
            .iter()
            .find(|f| f.name == "open")
            .unwrap()
            .joint_results
            .as_ref()
            .unwrap();
        assert!(
            opened
                .points
                .iter()
                .all(|p| p.opening_m > 1e-6 && p.normal_pressure_kpa == 0.),
            "mcc={mcc}: open gaps {:?}",
            opened.points.iter().map(|p| p.gap_m).collect::<Vec<_>>()
        );
        let closed = full
            .frames
            .iter()
            .find(|f| f.name == "recontact")
            .unwrap()
            .joint_results
            .as_ref()
            .unwrap();
        assert!(
            closed
                .points
                .iter()
                .all(|p| p.normal_pressure_kpa > 1. && p.opening_m == 0.)
        );
        assert!(
            full.frames
                .iter()
                .all(|f| f.force_balance_error_kn.iter().all(|x| x.abs() < 1e-6))
        );
        assert!(full.frames.iter().all(|f| {
            f.global_moment_balance_error_kn_m
                .unwrap()
                .iter()
                .all(|x| x.abs() < 1e-6)
        }));
        assert!(full.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_joint_slip_changed
            && a.committed_state_unchanged_on_rejection));
        let mut clean = case.clone();
        clean.verification_reject_once_stages.clear();
        let mut stages = serde_json::to_value(&clean.stages)?;
        for s in stages.as_array_mut().unwrap() {
            if ["shear", "open", "reverse"].contains(&s["name"].as_str().unwrap()) {
                let step = s["initial_increment"]
                    .as_f64()
                    .unwrap_or(case.settings.initial_increment);
                s["initial_increment"] = serde_json::json!(step / 2.);
            }
        }
        clean.stages = serde_json::from_value(stages)?;
        let fresh = stages::solve(&mesh, ground.clone(), clean)?;
        assert!(fresh.complete);
        assert_eq!(
            serde_json::to_value(&fresh.last_committed)?,
            serde_json::to_value(&full.last_committed)?
        );
        let geo = full
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .count();
        let shear = full
            .attempts
            .iter()
            .filter(|a| a.stage == "shear" && a.accepted)
            .count();
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), geo + shear + 1)?;
        let sealed = Checkpoint::seal(&mesh, partial)?;
        let restored: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&sealed)?)?;
        let resumed = restored.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert_eq!(
            serde_json::to_value(&resumed.frames)?,
            serde_json::to_value(&full.frames)?
        );
        let mut partial =
            stages::solve_until(&mesh, ground.clone(), case.clone(), geo + shear + 1)?;
        partial.last_committed.joints.as_mut().unwrap().points[0].traction_kpa[2] += 1.;
        assert!(
            Checkpoint::seal(&mesh, partial)?
                .resume(&mesh, ground, case, None)
                .is_err()
        );
    }
    Ok(())
}

#[test]
fn contact_normal_points_from_side_a_into_side_b() -> Result<()> {
    let (mesh, model) = fixture(true)?;
    let owners = boundary::surface_owners(&mesh);
    for p in model.layout() {
        assert!((Vector3::from(p.normal) - Vector3::new(0., 0., 1.)).norm() < 1e-12);
        assert_ne!(p.owner_a, p.owner_b);
    }
    assert!(
        mesh.face_sets["Seam:a"]
            .iter()
            .all(|f| owners[f].len() == 1)
    );
    Ok(())
}

#[test]
fn excavation_of_contact_side_archives_history_and_releases_only_its_old_force() -> Result<()> {
    let (mesh, ground, mut case) = coupled(false)?;
    case.verification_reject_once_stages = vec!["excavate_upper".into()];
    case.stages = serde_json::from_value(
        serde_json::json!([{"kind":"excavate","name":"excavate_upper","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.}]),
    )?;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
    assert!(
        run.frames[0]
            .joint_results
            .as_ref()
            .unwrap()
            .points
            .iter()
            .any(|p| p.normal_pressure_kpa > 10.)
    );
    assert!(
        run.frames[1]
            .joint_results
            .as_ref()
            .unwrap()
            .points
            .iter()
            .all(|p| !p.active)
    );
    assert_eq!(run.frames[0].state.joints, run.frames[1].state.joints);
    assert!(
        run.last_committed
            .excavation_release_load_kn
            .iter()
            .all(|x| *x == 0.)
    );
    assert!(run.frames[1].relative_force_residual < 1e-8);
    Ok(())
}
