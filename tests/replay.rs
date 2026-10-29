use soil_fem_3d::{
    Result,
    ground::Ground,
    mesh::Mesh,
    replay::{self, Outcome, PsiPolicy, Specification, TensionPolicy},
    soil,
    stages::{self, StageCase},
};

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let top = soil_fem_3d::boundary::surface_owners(&mesh)
        .into_iter()
        .filter(|(f, o)| o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][2] == 1.))
        .map(|(f, _)| f)
        .collect();
    mesh.face_sets.insert("top".into(), top);
    let ground = serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let case = serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    Ok((mesh.elevate_straight()?, ground, case))
}
fn spec() -> Specification {
    Specification {
        schema: "soil-fem-replay/1".into(),
        material_ids: vec![1, 2],
        factors: vec![1., 1.1, 1.2],
        psi_policy: PsiPolicy::Keep,
        tension_policy: TensionPolicy::Keep,
    }
}

#[test]
fn unit_factor_is_exact_and_plastic_histories_do_not_leak_between_factor_orders() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let baseline = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(baseline.complete);
    let forward = replay::run(&mesh, ground.clone(), case.clone(), spec())?;
    assert_eq!(forward.stability_classification, "not_evaluated");
    assert!(
        forward
            .trials
            .iter()
            .all(|t| t.outcome == Outcome::Completed)
    );
    let unit = forward.trials[0].run.as_ref().unwrap();
    assert_eq!(
        serde_json::to_vec(&baseline.last_committed)?,
        serde_json::to_vec(&unit.last_committed)?
    );
    assert_eq!(
        serde_json::to_vec(&baseline.frames)?,
        serde_json::to_vec(&unit.frames)?
    );
    assert!(
        unit.last_committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|x| x.abs() > 1e-7))
    );
    let mut reverse = spec();
    reverse.factors.reverse();
    let reverse = replay::run(&mesh, ground, case, reverse)?;
    for trial in &forward.trials {
        let other = reverse
            .trials
            .iter()
            .find(|t| t.factor == trial.factor)
            .unwrap();
        assert_eq!(
            serde_json::to_vec(&trial.run.as_ref().unwrap().last_committed)?,
            serde_json::to_vec(&other.run.as_ref().unwrap().last_committed)?
        );
        assert_eq!(
            serde_json::to_vec(&trial.run.as_ref().unwrap().frames)?,
            serde_json::to_vec(&other.run.as_ref().unwrap().frames)?
        );
    }
    assert_ne!(
        serde_json::to_vec(
            &forward.trials[0]
                .run
                .as_ref()
                .unwrap()
                .last_committed
                .material_states
        )?,
        serde_json::to_vec(
            &forward.trials[2]
                .run
                .as_ref()
                .unwrap()
                .last_committed
                .material_states
        )?
    );
    Ok(())
}

#[test]
fn selected_strength_and_policies_are_explicit_and_mcc_is_rejected() -> Result<()> {
    let (mesh, mut ground, case) = fixture()?;
    let mut request = spec();
    request.factors = vec![2.];
    request.material_ids = vec![2];
    request.psi_policy = PsiPolicy::Cap;
    request.tension_policy = TensionPolicy::DivideByFactor;
    let reduced = request.adjusted(&ground, 2.)?;
    assert_eq!(
        serde_json::to_vec(&ground.materials[0])?,
        serde_json::to_vec(&reduced.materials[0])?
    );
    if let soil::Model::MohrCoulomb(p) = &reduced.materials[1].model {
        assert_eq!(p.cohesion_kpa, 5.);
        assert!((p.friction_deg - 16.102113751986014).abs() < 1e-12);
    }
    let mcc: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    ground.materials[0].model = mcc.materials[0].model.clone();
    assert!(
        replay::run(&mesh, ground, case, request)
            .unwrap_err()
            .to_string()
            .contains("MCC")
    );
    Ok(())
}

#[test]
fn policy_errors_and_numerical_incompleteness_do_not_produce_stability_bounds() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let mut request = spec();
    request.factors = vec![1., 100.];
    let report = replay::run(&mesh, ground.clone(), case.clone(), request)?;
    assert_eq!(report.trials[0].outcome, Outcome::Completed);
    assert_eq!(report.trials[1].outcome, Outcome::RejectedInput);
    assert!(report.trials[1].adjusted_ground.is_none());
    assert!(report.trials[1].run.is_none());
    let mut hard = case;
    hard.settings.maximum_newton_iterations = 1;
    hard.settings.initial_increment = 1.;
    hard.settings.minimum_increment = 1.;
    let mut ground = ground;
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    let mut request = spec();
    request.factors = vec![1.];
    let report = replay::run(&mesh, ground, hard, request)?;
    assert_eq!(report.trials[0].outcome, Outcome::Incomplete);
    let failed = report.trials[0].run.as_ref().unwrap();
    assert!(!failed.complete && failed.frames.is_empty());
    assert!(
        failed
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .all(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    assert_eq!(report.stability_classification, "not_evaluated");
    Ok(())
}

#[test]
fn stage_gate_prefix_factor_one_and_factor_order_preserve_independent_histories() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete);
    for through in [0, 1, case.stages.len()] {
        let gate = replay::run_through_stage(&mesh, ground.clone(), case.clone(), spec(), through)?;
        assert_eq!(gate.requested_through_construction_stages, through);
        assert_eq!(gate.original_case.stages.len(), case.stages.len());
        assert_eq!(gate.replay.case.stages.len(), through);
        assert_eq!(
            gate.strength_applied_from,
            "original_inputs_before_initial_equilibration"
        );
        assert_eq!(gate.stability_classification, "not_evaluated");
        assert!(
            gate.replay
                .trials
                .iter()
                .all(|t| t.outcome == Outcome::Completed)
        );
        let unit = gate.replay.trials[0].run.as_ref().unwrap();
        assert_eq!(
            serde_json::to_vec(&unit.frames)?,
            serde_json::to_vec(&full.frames[..through + 1])?
        );
        assert_eq!(
            serde_json::to_vec(&unit.last_committed)?,
            serde_json::to_vec(&full.frames[through].state)?
        );
        let mut reverse = spec();
        reverse.factors.reverse();
        let other =
            replay::run_through_stage(&mesh, ground.clone(), case.clone(), reverse, through)?;
        for a in &gate.replay.trials {
            let b = other
                .replay
                .trials
                .iter()
                .find(|b| b.factor == a.factor)
                .unwrap();
            assert_eq!(
                serde_json::to_vec(&a.run.as_ref().unwrap().frames)?,
                serde_json::to_vec(&b.run.as_ref().unwrap().frames)?
            );
        }
    }
    Ok(())
}

#[test]
fn stage_gate_validates_bounds_and_keeps_original_suffix_diagnostics_and_policy_errors()
-> Result<()> {
    let (mesh, ground, mut case) = fixture()?;
    assert!(
        replay::run_through_stage(
            &mesh,
            ground.clone(),
            case.clone(),
            spec(),
            case.stages.len() + 1
        )
        .is_err()
    );
    assert!(
        replay::run_through_stage(&mesh, ground.clone(), case.clone(), spec(), usize::MAX).is_err()
    );
    case.verification_reject_once_stages = vec!["small_unload".into()];
    let report = replay::run_through_stage(&mesh, ground.clone(), case.clone(), spec(), 1)?;
    assert_eq!(
        report.original_case.verification_reject_once_stages,
        case.verification_reject_once_stages
    );
    assert!(
        report
            .replay
            .case
            .verification_reject_once_stages
            .is_empty()
    );
    assert_eq!(
        report.source_identity.mesh_sha256,
        soil_fem_3d::checkpoint::hash_json(&mesh)?
    );
    assert_eq!(
        report.source_identity.ground_sha256,
        soil_fem_3d::checkpoint::hash_json(&ground)?
    );
    assert_eq!(
        report.source_identity.case_sha256,
        soil_fem_3d::checkpoint::hash_json(&case)?
    );
    assert_ne!(
        report.source_identity.case_sha256,
        report.replay.source_identity.case_sha256
    );
    let mut request = spec();
    request.factors = vec![100.];
    let bad = replay::run_through_stage(&mesh, ground.clone(), case.clone(), request, 1)?;
    assert_eq!(bad.replay.trials[0].outcome, Outcome::RejectedInput);
    assert!(bad.replay.trials[0].run.is_none());
    let mcc: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    assert!(
        replay::run_through_stage(&mesh, mcc, case, spec(), 0)
            .unwrap_err()
            .to_string()
            .contains("MCC")
    );
    Ok(())
}

#[test]
fn stage_gate_numeric_incompleteness_preserves_rejected_history_without_stability_bounds()
-> Result<()> {
    let (mesh, mut ground, mut case) = fixture()?;
    ground
        .initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    case.settings.maximum_newton_iterations = 1;
    case.settings.initial_increment = 1.;
    case.settings.minimum_increment = 1.;
    let mut request = spec();
    request.factors = vec![1.];
    let gate = replay::run_through_stage(&mesh, ground, case, request, 1)?;
    assert_eq!(gate.replay.trials[0].outcome, Outcome::Incomplete);
    let run = gate.replay.trials[0].run.as_ref().unwrap();
    assert!(!run.complete && run.frames.is_empty());
    assert!(
        run.geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .all(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    assert_eq!(gate.stability_classification, "not_evaluated");
    assert_eq!(gate.replay.stability_classification, "not_evaluated");
    Ok(())
}
