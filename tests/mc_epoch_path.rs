use soil_fem_3d::{
    Result,
    checkpoint::hash_json,
    ground::Ground,
    mc::{Flow, Parameters},
    mc_epoch_path::{Advance, Checkpoint, Controller, FailureKind, Specification, StopReason},
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};
fn fixture(tension: bool) -> Result<(Mesh, stages::Run)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    for m in &mut ground.materials {
        m.unit_weight_kn_m3 = 1e-10;
        m.model = soil::Model::MohrCoulomb(Parameters {
            young_kpa: 30_000.,
            poisson: 0.25,
            cohesion_kpa: if tension { 20. } else { 40. },
            friction_deg: if tension { 30. } else { 0. },
            dilation_deg: 0.,
            tension_cutoff_kpa: None,
            flow: Flow::NonAssociated,
        });
    }
    if tension {
        // Prescribed affine extension produces exactly sigma_i=-10 kPa.
        // Pure MC with phi=30, psi=0, c=0 has no hydrostatic tensile return.
        let constraints = (0..3)
            .map(|axis| {
                let mut direction = [0.; 3];
                direction[axis] = 1.;
                let mut gradient = [0.; 3];
                gradient[axis] = 1. / 6000.;
                soil_fem_3d::boundary::Constraint {
                    node_set: "all".into(),
                    direction,
                    gradient,
                    value_m: 0.,
                    reference_point_m: [0.; 3],
                }
            })
            .collect();
        case.boundary = Some(soil_fem_3d::boundary::Conditions {
            constraints,
            nodal_loads: vec![],
            surface_loads: vec![],
        });
        case.stages =
            serde_json::from_str(r#"[{"kind":"boundary","name":"extend","target_factor":1.0}]"#)?;
    } else {
        case.stages = serde_json::from_str(
            r#"[{"kind":"pressure","name":"load","target_pressure_kpa":100.0}]"#,
        )?;
    }
    let origin = stages::solve(&mesh, ground, case)?;
    assert!(origin.complete, "{:?}", origin.failure_reason);
    Ok((mesh, origin))
}

fn specification(origin: &stages::Run, c: f64) -> Specification {
    let target_parameters = origin
        .ground
        .materials
        .iter()
        .map(|m| {
            let soil::Model::MohrCoulomb(mut p) = m.model.clone() else {
                panic!("MC")
            };
            p.cohesion_kpa = c;
            p
        })
        .collect();
    Specification {
        schema: "soil-fem-mc-epoch-path/1".into(),
        target_parameters,
        initial_increment: 0.25,
        minimum_increment: 1. / 1024.,
        maximum_attempts: 100,
    }
}
fn finish(controller: &mut Controller<'_>) -> Result<()> {
    for _ in 0..1000 {
        if controller.advance()? == Advance::Finished {
            return Ok(());
        }
    }
    panic!("path did not terminate within declared test budget")
}
fn clone_checkpoint(cp: &Checkpoint) -> Result<Checkpoint> {
    Ok(serde_json::from_slice(&serde_json::to_vec(cp)?)?)
}
fn reseal(cp: &mut Checkpoint) -> Result<()> {
    cp.payload_sha256 = hash_json(&cp.payload)?;
    Ok(())
}

#[test]
fn feasible_plastic_path_and_resume_match_exactly_and_closed_form() -> Result<()> {
    let (mesh, origin) = fixture(false)?;
    let spec = specification(&origin, 20.);
    let mut a = Controller::new(&mesh, &origin, spec)?;
    assert_eq!(a.advance()?, Advance::Accepted);
    assert_eq!(a.advance()?, Advance::Accepted);
    assert_eq!(a.continuation().progress, 0.5);
    let checkpoint = a.checkpoint()?;
    let cp = clone_checkpoint(&checkpoint)?;
    let mut b = Controller::resume(&mesh, &origin, &cp)?;
    finish(&mut a)?;
    finish(&mut b)?;
    assert!(a.complete() && b.complete());
    assert_eq!(a.continuation().accepted.len(), 4);
    assert_eq!(hash_json(&a.checkpoint()?)?, hash_json(&b.checkpoint()?)?);
    assert_eq!(a.stability_classification(), "not_evaluated");
    for (p, u) in mesh
        .nodes
        .iter()
        .zip(&a.mechanics().committed().displacement_m)
    {
        assert!((u[2] + 110. / 30000. * p[2]).abs() < 1e-9);
    }
    assert_eq!(hash_json(&checkpoint)?, hash_json(&cp)?);
    Ok(())
}

#[test]
fn actual_material_failures_cut_back_restore_parameters_and_resume_prefix_exactly() -> Result<()> {
    let (mesh, origin) = fixture(true)?;
    let mut spec = specification(&origin, 0.);
    spec.initial_increment = 1.;
    let mut a = Controller::new(&mesh, &origin, spec)?;
    let old_state = hash_json(a.mechanics().committed())?;
    let old_parameters = hash_json(&a.mechanics().parameters())?;
    assert_eq!(a.advance()?, Advance::Rejected);
    let attempt = &a.continuation().attempts[0];
    assert_eq!(attempt.failure_kind, Some(FailureKind::Material));
    assert_eq!(hash_json(a.mechanics().committed())?, old_state);
    assert_eq!(hash_json(&a.mechanics().parameters())?, old_parameters);
    assert_eq!(a.continuation().next_increment, 0.5);
    let cp = clone_checkpoint(&a.checkpoint()?)?;
    let mut b = Controller::resume(&mesh, &origin, &cp)?;
    finish(&mut a)?;
    finish(&mut b)?;
    assert!(!a.complete());
    assert_eq!(
        a.continuation().stop_reason,
        Some(StopReason::CutbackExhausted)
    );
    assert!(a.continuation().accepted.len() > 2);
    assert_eq!(hash_json(&a.checkpoint()?)?, hash_json(&b.checkpoint()?)?);
    for trial in a.continuation().attempts.iter().filter(|x| !x.accepted) {
        assert_eq!(trial.before_state_sha256, trial.after_state_sha256);
        assert_eq!(
            trial.before_parameters_sha256,
            trial.after_parameters_sha256
        );
    }
    // Independent hydrostatic tensile admissibility: c/tan(phi)>=10kPa.
    let limit = 1. - 10. * 30_f64.to_radians().tan() / 20.;
    let progress = a.continuation().progress;
    assert!(progress <= limit + 1e-10 && limit - progress < 2. / 1024.);
    assert_eq!(a.stability_classification(), "not_evaluated");
    Ok(())
}

#[test]
fn resealed_progress_increment_target_and_history_chain_tampering_are_rejected() -> Result<()> {
    let (mesh, origin) = fixture(false)?;
    let mut a = Controller::new(&mesh, &origin, specification(&origin, 20.))?;
    a.advance()?;
    let original = a.checkpoint()?;
    for option in 0..6 {
        let mut cp = clone_checkpoint(&original)?;
        match option {
            0 => cp.payload.continuation.progress = 0.5,
            1 => cp.payload.continuation.next_increment = 0.4,
            2 => cp.payload.specification.target_parameters[0].cohesion_kpa = 25.,
            3 => {
                cp.payload.continuation.accepted.clear();
            }
            4 => cp.payload.continuation.attempts[0].before_state_sha256 = "changed".into(),
            _ => cp.payload.continuation.stop_reason = Some(StopReason::CutbackExhausted),
        }
        reseal(&mut cp)?;
        assert!(
            Controller::resume(&mesh, &origin, &cp).is_err(),
            "tamper{option}"
        );
    }
    Ok(())
}

#[test]
fn attempt_limit_and_invalid_controls_do_not_produce_stability_bounds() -> Result<()> {
    let (mesh, origin) = fixture(false)?;
    let mut spec = specification(&origin, 20.);
    spec.maximum_attempts = 1;
    let mut a = Controller::new(&mesh, &origin, spec)?;
    a.advance()?;
    assert!(!a.complete());
    assert_eq!(a.continuation().stop_reason, Some(StopReason::AttemptLimit));
    let cp = a.checkpoint()?;
    let mut b = Controller::resume(&mesh, &origin, &cp)?;
    assert_eq!(b.advance()?, Advance::Finished);
    assert_eq!(b.stability_classification(), "not_evaluated");
    for option in 0..5 {
        let mut s = specification(&origin, 20.);
        match option {
            0 => s.minimum_increment = f64::NAN,
            1 => s.maximum_attempts = 0,
            2 => s.target_parameters[0].young_kpa *= 2.,
            3 => s.target_parameters[0].tension_cutoff_kpa = Some(0.),
            _ => s.target_parameters[0].cohesion_kpa = -1.,
        };
        assert!(Controller::new(&mesh, &origin, s).is_err());
    }
    Ok(())
}
