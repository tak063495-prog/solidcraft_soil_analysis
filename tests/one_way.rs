use soil_fem_3d::{
    Result, checkpoint,
    ground::Ground,
    mesh::Mesh,
    one_way::{self, Case},
    stages,
};

fn fixture(quadratic: bool, mcc: bool) -> Result<(Mesh, Ground, Case)> {
    let mesh = serde_json::from_str(if quadratic {
        include_str!("../examples/seepage-column-tet10-mesh.json")
    } else {
        include_str!("../examples/seepage-column-tet4-mesh.json")
    })?;
    let ground = serde_json::from_str(if mcc {
        include_str!("../examples/ground-one-way-mcc.json")
    } else {
        include_str!("../examples/ground-one-way-mc.json")
    })?;
    let case = serde_json::from_str(if mcc {
        include_str!("../examples/one-way-mcc.json")
    } else {
        include_str!("../examples/one-way-mc.json")
    })?;
    Ok((mesh, ground, case))
}
fn sum_force(f: &[f64]) -> [f64; 3] {
    let mut result = [0.; 3];
    for (g, x) in f.iter().enumerate() {
        result[g % 3] += x;
    }
    result
}

#[test]
fn raw_hydraulic_quadrature_transfers_water_mass_and_symmetric_pore_force_without_double_counting()
-> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, case) = fixture(quadratic, false)?;
        let gamma = ground.initialization.water_unit_weight_kn_m3;
        let run = one_way::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(
            run.last_committed.fluid.hydraulic_points.len(),
            27 * mesh.tetrahedra.len()
        );
        if !quadratic {
            assert_eq!(
                run.last_committed.fluid.mechanical_points.len(),
                mesh.tetrahedra.len()
            );
        }
        for f in &run.frames {
            let fluid = &f.state.fluid;
            let total = sum_force(&fluid.water_body_nodal_load_kn);
            assert!((total[2] + gamma * fluid.water_equivalent_mass_volume_m3).abs() < 1e-12);
            assert!(
                sum_force(&fluid.fluid_internal_nodal_force_kn)
                    .iter()
                    .all(|v| v.abs() < 1e-12)
            );
            let mut moment = [0.; 3];
            for (p, q) in mesh
                .nodes
                .iter()
                .zip(fluid.fluid_internal_nodal_force_kn.as_chunks::<3>().0)
            {
                moment[0] += p[1] * q[2] - p[2] * q[1];
                moment[1] += p[2] * q[0] - p[0] * q[2];
                moment[2] += p[0] * q[1] - p[1] * q[0];
            }
            assert!(moment.iter().all(|v| v.abs() < 1e-12));
            assert!(
                f.state
                    .mechanical
                    .pore_pressure_kpa
                    .as_ref()
                    .unwrap()
                    .iter()
                    .all(|p| *p == 0.)
            );
            for ((physical, pore), rhs) in fluid
                .physical_nodal_load_kn
                .iter()
                .zip(&fluid.fluid_internal_nodal_force_kn)
                .zip(&fluid.equivalent_mechanical_rhs_kn)
            {
                assert_eq!(*rhs, physical - pore);
            }
            if f.initialization_progress == 1. {
                assert_eq!(
                    f.state.mechanical.applied_nodal_load_kn,
                    fluid.equivalent_mechanical_rhs_kn
                );
            }
        }
        assert!(run.last_committed.fluid.hydraulic_points.iter().any(|p| {
            p.physical_pore_pressure_kpa < 0.
                && p.bishop_chi > 0.
                && p.bishop_chi < 1.
                && p.effective_pore_stress_kpa != p.physical_pore_pressure_kpa
        }));
    }
    Ok(())
}

#[test]
fn initial_unsaturated_column_matches_independent_continuous_elastic_integrals() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, mut case) = fixture(quadratic, false)?;
        case.flow.analysis = serde_json::from_value(
            serde_json::json!({"kind":"transient","periods":[{"name":"held","end_time_s":1.,"boundaries":[{"face_set":"bottom","condition":{"kind":"head","value":{"head_m":0.}}}]}]}),
        )?;
        case.exterior_periods = vec![one_way::ExteriorPeriod {
            name: "held".into(),
            reservoirs: vec![],
        }];
        let run = one_way::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        let initial = &run
            .frames
            .iter()
            .find(|f| f.initialization_progress == 1.)
            .unwrap()
            .state;
        // sigma'_zz - seed = gamma_w int_z^1 theta(s) ds - chi(z)u(z).
        // Integrate this stress from the fixed bottom; direct VG expressions
        // and independent 1D composite Simpson, not native quadrature/loads.
        let n = 16000;
        let law = |z: f64| {
            let m = 1. - 1. / 1.6;
            let se = (1. + (1.2 * z).powf(1.6)).powf(-m);
            (0.08 + (0.7 / 1.7 - 0.08) * se, -10. * z * se)
        };
        for (i, p) in mesh.nodes.iter().enumerate() {
            let z = p[2];
            let mut body = 0.;
            let mut suction = 0.;
            for j in 0..=n {
                let t = j as f64 / n as f64;
                let weight = if j == 0 || j == n {
                    1.
                } else if j % 2 == 0 {
                    2.
                } else {
                    4.
                };
                body += weight * z.min(t) * law(t).0;
                let s = z * t;
                suction += weight * law(s).1 * z;
            }
            let expected = -(10. * body - suction) / (3. * n as f64 * 12000.);
            assert!(
                (initial.mechanical.displacement_m[i][2] - expected).abs() < 1e-6,
                "Tet10={quadratic} z={z} actual={} expected={expected}",
                initial.mechanical.displacement_m[i][2]
            );
        }
    }
    Ok(())
}

#[test]
fn saturated_limit_recovers_existing_pore_pressure_displacement_increments() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, mut case) = fixture(quadratic, false)?;
        case.flow.initial_head.head_m = 1.;
        case.flow.analysis = serde_json::from_value(
            serde_json::json!({"kind":"transient","periods":[{"name":"pressurize","end_time_s":1.,"boundaries":[{"face_set":"bottom","condition":{"kind":"head","value":{"head_m":1.5}}}]}]}),
        )?;
        case.initial_exterior_water = vec![one_way::Reservoir {
            face_set: "top".into(),
            level_m: 1.,
        }];
        case.exterior_periods = vec![one_way::ExteriorPeriod {
            name: "pressurize".into(),
            reservoirs: vec![one_way::Reservoir {
                face_set: "top".into(),
                level_m: 1.,
            }],
        }];
        let mut legacy = ground.clone();
        legacy.initialization.water_level_m = Some(1.);
        for m in &mut legacy.materials {
            m.unit_weight_kn_m3 += 10. * (0.7 / 1.7);
        }
        let mut stages = case.mechanics.clone();
        let mut flow: serde_json::Value =
            serde_json::from_str(include_str!("../examples/stages-seepage-column.json"))?;
        flow = flow["stages"][0]["flow"].clone();
        flow["heads"] = serde_json::json!([{"face_set":"bottom","head_m":1.5}]);
        stages.stages = serde_json::from_value(
            serde_json::json!([{"kind":"steady_seepage","name":"pressurize","flow":flow,"exterior_water_level_m":1.,"initial_increment":1.}]),
        )?;
        let old = stages::solve(&mesh, legacy, stages)?;
        assert!(old.complete);
        let run = one_way::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert!(
            run.last_committed
                .fluid
                .hydraulic_points
                .iter()
                .all(|p| p.bishop_chi == 1.)
        );
        for i in 0..mesh.nodes.len() {
            for a in 0..3 {
                let increment = run.last_committed.mechanical.displacement_m[i][a]
                    - run.initial_reference_displacement_m[i][a];
                let expected = old.frames.last().unwrap().state.displacement_m[i][a]
                    - old.frames[0].state.displacement_m[i][a];
                assert!(
                    (increment - expected).abs() < 5e-11,
                    "Tet10={quadratic} node={i} actual={increment} expected={expected}"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn bishop_zero_one_and_separate_exterior_pressure_match_exact_elastic_limits() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, mut zero) = fixture(quadratic, false)?;
        zero.flow.analysis = serde_json::from_value(
            serde_json::json!({"kind":"transient","periods":[{"name":"held","end_time_s":1.,"boundaries":[{"face_set":"bottom","condition":{"kind":"head","value":{"head_m":0.}}}]}]}),
        )?;
        zero.exterior_periods = vec![one_way::ExteriorPeriod {
            name: "held".into(),
            reservoirs: vec![],
        }];
        zero.bishop = one_way::Bishop::Constant { value: 0. };
        let mut one = zero.clone();
        one.bishop = one_way::Bishop::Constant { value: 1. };
        let a = one_way::solve(&mesh, ground.clone(), zero.clone())?;
        let b = one_way::solve(&mesh, ground.clone(), one)?;
        assert!(a.complete && b.complete);
        assert!(
            a.last_committed
                .fluid
                .fluid_internal_nodal_force_kn
                .iter()
                .all(|v| *v == 0.)
        );
        assert_eq!(
            a.last_committed.fluid.water_body_nodal_load_kn,
            b.last_committed.fluid.water_body_nodal_load_kn
        );
        for (i, p) in mesh.nodes.iter().enumerate() {
            let actual = b.last_committed.mechanical.displacement_m[i][2]
                - a.last_committed.mechanical.displacement_m[i][2];
            assert!(
                (actual + 5. * p[2] * p[2] / 12000.).abs() < 5e-11,
                "Tet10={quadratic} node={i} chi difference={actual}"
            );
        }
        // A separately prescribed exterior pressure adds its own traction;
        // it does not alter the stored internal physical water field.
        zero.initial_exterior_water = vec![one_way::Reservoir {
            face_set: "top".into(),
            level_m: 1.5,
        }];
        zero.exterior_periods[0].reservoirs = zero.initial_exterior_water.clone();
        let c = one_way::solve(&mesh, ground, zero)?;
        assert!(c.complete);
        for (i, p) in mesh.nodes.iter().enumerate() {
            let actual = c.last_committed.mechanical.displacement_m[i][2]
                - a.last_committed.mechanical.displacement_m[i][2];
            assert!(
                (actual + 5. * p[2] / 12000.).abs() < 5e-11,
                "Tet10={quadratic} exterior node={i} difference={actual}"
            );
        }
        assert!((sum_force(&c.last_committed.fluid.exterior_nodal_load_kn)[2] + 5.).abs() < 1e-12);
    }
    Ok(())
}

#[test]
fn declared_compressible_water_mass_includes_positive_storage_and_closed_head_is_preserved()
-> Result<()> {
    let (mesh, ground, mut case) = fixture(false, false)?;
    case.flow.initial_head.head_m = 2.;
    for m in &mut case.flow.materials {
        m.retention.saturated_specific_storage_per_m = 0.1;
    }
    case.flow.analysis = serde_json::from_value(
        serde_json::json!({"kind":"transient","periods":[{"name":"closed","end_time_s":1.,"boundaries":[]}]}),
    )?;
    case.exterior_periods = vec![one_way::ExteriorPeriod {
        name: "closed".into(),
        reservoirs: vec![],
    }];
    let run = one_way::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.last_committed
            .hydraulic
            .relative_cell_head_m
            .iter()
            .all(|h| *h == Some(2.))
    );
    let mass = 0.7 / 1.7 + 0.1 * 1.5;
    assert!((run.last_committed.fluid.water_equivalent_mass_volume_m3 - mass).abs() < 1e-12);
    assert!(
        (sum_force(&run.last_committed.fluid.water_body_nodal_load_kn)[2] + 10. * mass).abs()
            < 1e-12
    );
    Ok(())
}

#[test]
fn genuine_mechanical_exhaustion_does_not_commit_the_successful_hydraulic_trial() -> Result<()> {
    let (mesh, ground, mut case) = fixture(true, true)?;
    let soil_fem_3d::richards::Analysis::Transient { periods } = &mut case.flow.analysis else {
        unreachable!()
    };
    periods.truncate(1);
    periods[0].end_time_s = 100.;
    case.exterior_periods.truncate(1);
    case.mechanics.settings.maximum_newton_iterations = 4;
    case.flow.settings.minimum_time_step_s = 100.;
    let run = one_way::solve(&mesh, ground, case)?;
    assert!(!run.complete);
    let failed = run.attempts.last().unwrap();
    assert_eq!(failed.phase, "transient");
    assert!(failed.hydraulic_candidate_changed);
    assert!(!failed.accepted);
    assert!(
        failed
            .reason
            .as_ref()
            .unwrap()
            .contains("Newton iteration limit")
    );
    assert_eq!(
        failed.committed_pair_sha256_before,
        failed.committed_pair_sha256_after
    );
    assert_eq!(run.last_committed.hydraulic.time_s, 0.);
    assert_eq!(
        serde_json::to_value(&run.last_committed)?,
        serde_json::to_value(&run.frames.last().unwrap().state)?
    );
    Ok(())
}

#[test]
fn plastic_coupled_cutbacks_restore_both_histories_and_same_accepted_times_replay_exactly()
-> Result<()> {
    let (mesh, ground, mut case) = fixture(true, true)?;
    let soil_fem_3d::richards::Analysis::Transient { periods } = &mut case.flow.analysis else {
        unreachable!()
    };
    periods.truncate(1);
    periods[0].end_time_s = 100.;
    case.exterior_periods.truncate(1);
    case.verification_reject_once_at_times_s = vec![100.];
    let run = one_way::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(run.attempts.iter().any(|a| !a.accepted
        && a.hydraulic_candidate_changed
        && a.mechanical_candidate_changed
        && a.plastic_candidate_changed));
    for a in run.attempts.iter().filter(|a| !a.accepted) {
        assert_eq!(
            a.committed_pair_sha256_before,
            a.committed_pair_sha256_after
        );
    }
    let mut clean = case;
    clean.verification_reject_once_at_times_s.clear();
    let soil_fem_3d::richards::Analysis::Transient { periods } = &mut clean.flow.analysis else {
        unreachable!()
    };
    let original = periods[0].clone();
    *periods = run
        .frames
        .iter()
        .filter(|f| f.phase == "transient")
        .enumerate()
        .map(|(i, f)| soil_fem_3d::richards::Period {
            name: format!("clean-{i}"),
            end_time_s: f.state.hydraulic.time_s,
            boundaries: original.boundaries.clone(),
        })
        .collect();
    clean.exterior_periods = periods
        .iter()
        .map(|p| one_way::ExteriorPeriod {
            name: p.name.clone(),
            reservoirs: vec![],
        })
        .collect();
    let direct = one_way::solve(&mesh, ground, clean)?;
    assert!(direct.complete);
    assert_eq!(
        serde_json::to_value(run.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?,
        serde_json::to_value(direct.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?
    );
    Ok(())
}

#[test]
fn checkpoints_replay_initialization_and_physical_periods_and_reject_resealed_changes() -> Result<()>
{
    let (mesh, ground, case) = fixture(false, true)?;
    let direct = one_way::solve(&mesh, ground.clone(), case.clone())?;
    assert!(direct.complete);
    for count in [1, 2, 3, 12, 22, 31] {
        let partial = one_way::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        let saved = one_way::Checkpoint::seal(partial)?;
        if count == 3 {
            for mode in 0..3 {
                let mut value = serde_json::to_value(&saved)?;
                match mode {
                    0 => value["run"]["runtime"]["next_time_step_s"] = serde_json::json!(0.),
                    1 => {
                        value["run"]["solver_binary_sha256"] =
                            serde_json::json!("changed executable")
                    }
                    _ => {
                        value["run"]["last_committed"]["hydraulic"]["cell_storage_offset_m3"][0] =
                            serde_json::json!(999.)
                    }
                }
                let mut forged: one_way::Checkpoint = serde_json::from_value(value)?;
                forged.run_sha256 = checkpoint::hash_json(&forged.run)?;
                assert!(
                    forged
                        .resume(&mesh, ground.clone(), case.clone(), None)
                        .is_err()
                );
            }
        }
        let mut bad = serde_json::to_value(&saved)?;
        bad["run"]["last_committed"]["fluid"]["mechanical_points"][0]["physical_pore_pressure_kpa"] =
            serde_json::json!(999.);
        let mut bad: one_way::Checkpoint = serde_json::from_value(bad)?;
        bad.run_sha256 = checkpoint::hash_json(&bad.run)?;
        assert!(
            bad.resume(&mesh, ground.clone(), case.clone(), None)
                .is_err()
        );
        let resumed = saved.resume(&mesh, ground.clone(), case.clone(), None)?;
        assert_eq!(
            one_way::normalized(&resumed)?,
            one_way::normalized(&direct)?
        );
    }
    Ok(())
}

#[test]
fn invalid_closure_porosity_exterior_and_mixed_stage_policies_are_rejected() -> Result<()> {
    let (mesh, ground, case) = fixture(false, false)?;
    let mut bad = case.clone();
    bad.initial_porosity_tolerance = 0.;
    assert!(one_way::solve(&mesh, ground.clone(), bad).is_err());
    let mut bad = case.clone();
    bad.flow.materials[0].retention.saturated_water_content = 0.9;
    assert!(one_way::solve(&mesh, ground.clone(), bad).is_err());
    let mut bad = case.clone();
    bad.exterior_periods[1].reservoirs = vec![one_way::Reservoir {
        face_set: "missing".into(),
        level_m: 0.,
    }];
    assert!(one_way::solve(&mesh, ground.clone(), bad).is_err());
    let mut bad = case.clone();
    bad.bishop = one_way::Bishop::Constant { value: 1.1 };
    assert!(one_way::solve(&mesh, ground.clone(), bad).is_err());
    let mut bad = case;
    bad.mechanics.stages = vec![stages::Stage::HydrostaticWater {
        name: "unsupported mixing".into(),
        target_level_m: 0.,
        initial_increment: None,
    }];
    assert!(one_way::solve(&mesh, ground, bad).is_err());
    Ok(())
}
