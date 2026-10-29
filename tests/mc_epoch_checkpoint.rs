use soil_fem_3d::{
    Error, Result,
    checkpoint::hash_json,
    ground::Ground,
    mc::{Flow, Parameters},
    mesh::Mesh,
    soil,
    stages::{self, McEpochCheckpoint, McParameterMechanics, StageCase},
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
fn cohesion(controller: &McParameterMechanics<'_>, c: f64) -> Vec<Parameters> {
    controller
        .parameters()
        .iter()
        .cloned()
        .map(|mut p| {
            p.cohesion_kpa = c;
            p
        })
        .collect()
}
fn clone_checkpoint(cp: &McEpochCheckpoint) -> Result<McEpochCheckpoint> {
    Ok(serde_json::from_slice(&serde_json::to_vec(cp)?)?)
}
fn reseal(cp: &mut McEpochCheckpoint) -> Result<()> {
    cp.payload_sha256 = hash_json(&cp.payload)?;
    Ok(())
}
fn stable(v: serde_json::Value) -> serde_json::Value {
    match v {
        serde_json::Value::Object(map) => serde_json::Value::Object(
            map.into_iter()
                .filter(|(k, _)| k != "elapsed_seconds")
                .map(|(k, v)| (k, stable(v)))
                .collect(),
        ),
        serde_json::Value::Array(a) => {
            serde_json::Value::Array(a.into_iter().map(stable).collect())
        }
        v => v,
    }
}

#[test]
fn accepted_plastic_epoch_roundtrip_and_later_continuation_are_exact() -> Result<()> {
    let (mesh, origin) = fixture(false)?;
    let mut continuous = McParameterMechanics::new(&mesh, &origin)?;
    continuous.try_parameters(cohesion(&continuous, 30.))?;
    let checkpoint = continuous.checkpoint()?;
    assert_eq!(checkpoint.payload.parameters[0].cohesion_kpa, 30.);
    assert!(
        checkpoint
            .payload
            .committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|v| v.abs() > 1e-6))
    );
    let cp = clone_checkpoint(&checkpoint)?;
    let mut resumed = McParameterMechanics::resume(&mesh, &origin, &cp)?;
    assert_eq!(
        serde_json::to_vec(resumed.parameters())?,
        serde_json::to_vec(continuous.parameters())?
    );
    assert_eq!(
        hash_json(resumed.committed())?,
        hash_json(continuous.committed())?
    );
    let final_parameters = cohesion(&continuous, 20.);
    let expected = continuous.try_parameters(final_parameters.clone())?;
    let actual = resumed.try_parameters(final_parameters)?;
    assert_eq!(
        stable(serde_json::to_value(expected)?),
        stable(serde_json::to_value(actual)?)
    );
    assert_eq!(
        hash_json(resumed.committed())?,
        hash_json(continuous.committed())?
    );
    assert_eq!(hash_json(&checkpoint)?, hash_json(&cp)?);
    Ok(())
}

#[test]
fn tampering_and_resealed_wrong_material_load_origin_or_constraints_are_rejected() -> Result<()> {
    let (mesh, origin) = fixture(false)?;
    let mut controller = McParameterMechanics::new(&mesh, &origin)?;
    controller.try_parameters(cohesion(&controller, 30.))?;
    let original = controller.checkpoint()?;
    let before = hash_json(controller.committed())?;
    for option in 0..8 {
        let mut bad = clone_checkpoint(&original)?;
        match option {
            0 => bad.payload.parameters[0].cohesion_kpa = 31.,
            1 => bad.payload.parameters[0].young_kpa *= 2.,
            2 => bad.payload.committed.applied_nodal_load_kn[0] += 1.,
            3 => bad.payload.origin_run_sha256 = "changed".into(),
            4 => bad.payload.identity.solver_binary_sha256 = "changed".into(),
            5 => bad.payload.committed.displacement_m[0][0] += 0.1,
            6 => {
                let soil::State::MohrCoulomb(s) = &mut bad.payload.committed.material_states[0]
                else {
                    panic!("MC")
                };
                s.void_ratio = -1.;
            }
            _ => bad
                .payload
                .committed
                .material_states
                .pop()
                .map(|_| ())
                .unwrap(),
        }
        if option != 0 {
            reseal(&mut bad)?;
        }
        assert!(
            McParameterMechanics::resume(&mesh, &origin, &bad).is_err(),
            "accepted tamper option{option}"
        );
    }
    assert_eq!(hash_json(controller.committed())?, before);
    Ok(())
}

#[test]
fn genuine_global_material_failure_preserves_both_history_and_parameter_epoch() -> Result<()> {
    let (mesh, origin) = fixture(true)?;
    for state in &origin.last_committed.material_states {
        for s in &state.stress()[0..3] {
            assert!((s + 10.).abs() < 1e-7);
        }
    }
    let mut controller = McParameterMechanics::new(&mesh, &origin)?;
    let before = controller.checkpoint()?;
    let bad = cohesion(&controller, 0.);
    let failure = controller.try_parameters(bad);
    assert!(matches!(failure, Err(Error::Material(_))), "{failure:?}");
    assert_eq!(hash_json(&controller.checkpoint()?)?, hash_json(&before)?);
    // A feasible smaller parameter jump is still usable after the failed trial.
    let accepted = controller.try_parameters(cohesion(&controller, 10.))?;
    assert_eq!(accepted.iterations, 0);
    let snapshot = controller.checkpoint()?;
    let resumed = McParameterMechanics::resume(&mesh, &origin, &snapshot)?;
    assert_eq!(
        hash_json(resumed.committed())?,
        hash_json(controller.committed())?
    );
    assert_eq!(resumed.parameters()[0].cohesion_kpa, 10.);
    Ok(())
}
