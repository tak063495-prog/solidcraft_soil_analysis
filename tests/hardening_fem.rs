use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::{self, Ground},
    hardening,
    linalg::LinearSolver,
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(if quadratic {
        include_str!("../examples/hardening-column-tet10-mesh.json")
    } else {
        include_str!("../examples/mc-column-mesh.json")
    })?;
    Ok((
        mesh,
        serde_json::from_str(include_str!("../examples/ground-hardening-column.json"))?,
        serde_json::from_str(include_str!("../examples/stages-hardening-column.json"))?,
    ))
}

#[test]
fn idh_preload_unload_reload_and_excavation_satisfy_current_force_and_sparse_dense_gates()
-> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, case) = fixture(quadratic)?;
        let sparse = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(
            sparse.complete,
            "quadratic={quadratic}: {:?}",
            sparse.failure_reason
        );
        assert_eq!(sparse.schema, "soil-fem-stage-result/10");
        assert!(
            sparse
                .frames
                .iter()
                .all(|f| f.relative_force_residual <= case.settings.relative_force_tolerance),
            "unloading cannot normalize its endpoint residual by a larger previous load"
        );
        let top = |i: usize| {
            sparse.frames[i].monitors["ground_surface"]
                .mean_displacement_z_m
                .unwrap()
        };
        assert!(top(1) < 0. && top(2) > top(1) && top(2) < 0. && top(3) < top(2));
        let first = 100. / (-top(1));
        let unloading = 100. / (top(2) - top(1));
        assert!(
            unloading > 2. * first,
            "first loading must be softer than unloading"
        );
        let floor = &sparse.frames[4].monitors["excavation_floor"];
        assert!(floor.stage_increment_z_m.unwrap() > 0.);
        for (before, after) in sparse.frames[1]
            .state
            .material_states
            .iter()
            .zip(&sparse.frames[2].state.material_states)
        {
            let (
                soil::State::InvariantDoubleHardening(a),
                soil::State::InvariantDoubleHardening(b),
            ) = (before, after)
            else {
                panic!("wrong history type")
            };
            assert_eq!(a.cap_pressure_kpa, b.cap_pressure_kpa);
            assert_eq!(a.shear_plastic_strain, b.shear_plastic_strain);
        }
        let mut dense_case = case.clone();
        dense_case.settings.linear_solver = LinearSolver::DenseLu;
        dense_case.settings.maximum_free_dofs = 2000;
        let dense = stages::solve(&mesh, ground, dense_case)?;
        assert!(dense.complete, "{:?}", dense.failure_reason);
        for (a, b) in sparse.frames.iter().zip(&dense.frames) {
            let delta: f64 = a
                .state
                .displacement_m
                .iter()
                .flatten()
                .zip(b.state.displacement_m.iter().flatten())
                .map(|(a, b)| (a - b).powi(2))
                .sum();
            let norm: f64 = b.state.displacement_m.iter().flatten().map(|b| b * b).sum();
            assert!(delta.sqrt() / norm.sqrt().max(1e-12) < 1e-7);
            for (a, b) in a.state.material_states.iter().zip(&b.state.material_states) {
                let (
                    soil::State::InvariantDoubleHardening(a),
                    soil::State::InvariantDoubleHardening(b),
                ) = (a, b)
                else {
                    panic!("wrong history type")
                };
                assert!((a.cap_pressure_kpa - b.cap_pressure_kpa).abs() < 1e-6);
                assert!((a.shear_plastic_strain - b.shear_plastic_strain).abs() < 1e-8);
            }
        }
        assert!(sparse.excavation_audits[0].zero_release_relative_force_residual < 1e-8);
        for (i, point) in sparse.integration_points.iter().enumerate() {
            if !sparse.frames[4].state.active_elements[point.element_index] {
                assert_eq!(
                    serde_json::to_value(&sparse.frames[3].state.material_states[i])?,
                    serde_json::to_value(&sparse.frames[4].state.material_states[i])?
                );
            }
        }
    }
    Ok(())
}

#[test]
fn idh_cap_shear_void_ratio_and_release_histories_restore_after_cutback_and_all_phase_restarts()
-> Result<()> {
    let (mesh, ground, case) = fixture(true)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    for name in ["preload", "excavate_upper"] {
        assert!(full.attempts.iter().any(|a| a.stage == name
            && !a.accepted
            && a.rejected_trial_plastic_changed
            && a.committed_state_unchanged_on_rejection));
    }
    let mut clean = serde_json::to_value(&case)?;
    clean["verification_reject_once_stages"] = serde_json::json!([]);
    for stage in clean["stages"].as_array_mut().unwrap() {
        let name = stage["name"].as_str().unwrap();
        stage["initial_increment"] = serde_json::json!(
            full.attempts
                .iter()
                .find(|a| a.accepted && a.stage == name)
                .unwrap()
                .progress_to
        );
    }
    let clean = stages::solve(&mesh, ground.clone(), serde_json::from_value(clean)?)?;
    assert!(clean.complete);
    assert!(
        serde_json::to_value(&clean.frames)? == serde_json::to_value(&full.frames)?,
        "rejected material history leaked into a subsequent trial"
    );
    assert!(
        serde_json::to_value(&clean.last_committed)? == serde_json::to_value(&full.last_committed)?
    );
    let mut count = 0;
    for name in ["preload", "unload", "reload", "excavate_upper"] {
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count + 1)?;
        assert!(
            partial.last_committed.stage_progress > 0.
                && partial.last_committed.stage_progress < 1.
        );
        let sealed = Checkpoint::seal(&mesh, partial)?;
        let saved: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&sealed)?)?;
        let resumed = saved.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert!(resumed.complete);
        assert!(
            serde_json::to_value(&resumed.frames)? == serde_json::to_value(&full.frames)?,
            "phase restart changed IDH states"
        );
        assert!(
            serde_json::to_value(&resumed.last_committed)?
                == serde_json::to_value(&full.last_committed)?
        );
        count += full
            .attempts
            .iter()
            .filter(|a| a.accepted && a.stage == name)
            .count();
    }
    Ok(())
}

#[test]
fn idh_initial_cap_is_distinct_from_mcc_ocr_and_tampered_history_is_rejected() -> Result<()> {
    let (mesh, ground, case) = fixture(false)?;
    let initial = ground::initialize(&mesh, ground.clone())?;
    assert!(
        initial
            .point_states
            .iter()
            .all(|p| p.mcc_pc_over_mean_effective.is_none())
    );
    let mut invalid = serde_json::to_value(&ground)?;
    invalid["layers"][0]["consolidation"] =
        serde_json::json!({"definition":"mcc_pc_over_mean_effective","ratio":2.});
    assert!(
        serde_json::from_value::<Ground>(invalid)?
            .validate()
            .is_err()
    );
    let mut invalid = ground.clone();
    invalid.layers[0].k0_x = 0.4;
    invalid.layers[0].k0_y = 0.4;
    assert!(ground::initialize(&mesh, invalid).is_err());
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let mut partial = partial;
    let soil::State::InvariantDoubleHardening(s) = &mut partial.last_committed.material_states[0]
    else {
        panic!("wrong history")
    };
    s.cap_pressure_kpa += 1.;
    let resealed = Checkpoint::seal(&mesh, partial)?;
    assert!(resealed.resume(&mesh, ground.clone(), case, None).is_err());
    let parameters: hardening::Parameters =
        serde_json::from_str(include_str!("../examples/hardening-parameters.json"))?;
    let state = soil::State::InvariantDoubleHardening(parameters.initial_state(
        [100., 100., 100., 0., 0., 0.],
        0.7,
        200.,
    )?);
    let decoded: soil::State = serde_json::from_value(serde_json::to_value(&state)?)?;
    assert!(matches!(decoded, soil::State::InvariantDoubleHardening(_)));
    let mcc = soil::Model::Mcc(soil_fem_3d::mcc::Parameters {
        critical_state_m: 1.2,
        lambda: 0.2,
        kappa: 0.04,
        poisson: 0.3,
        reference_pressure_kpa: 1.,
        specific_volume_ncl_at_reference: 3.5,
        integration: soil_fem_3d::mcc::Integration::BackwardEuler,
    });
    assert!(mcc.validate_state(&state).is_err());
    Ok(())
}
