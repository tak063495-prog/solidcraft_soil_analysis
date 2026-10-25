use soil_fem_3d::{
    Result, checkpoint,
    ground::Ground,
    mesh::Mesh,
    one_way,
    one_way_construction::{self as construction, Case},
    richards,
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
    let mut base: one_way::Case = serde_json::from_str(if mcc {
        include_str!("../examples/one-way-mcc.json")
    } else {
        include_str!("../examples/one-way-mc.json")
    })?;
    base.flow.analysis = serde_json::from_value(serde_json::json!({"kind":"transient","periods":[
        {"name":"before","end_time_s":200.,"boundaries":[{"face_set":"bottom","condition":{"kind":"head","value":{"head_m":0.}}},{"face_set":"top","condition":{"kind":"atmospheric_no_ponding","rainfall_vertical_m_s":0.}}]},
        {"name":"after","end_time_s":400.,"boundaries":[{"face_set":"bottom","condition":{"kind":"head","value":{"head_m":0.}}},{"face_set":"floor","condition":{"kind":"atmospheric_no_ponding","rainfall_vertical_m_s":0.}}]}
    ]}))?;
    base.exterior_periods = vec![
        one_way::ExteriorPeriod {
            name: "before".into(),
            reservoirs: vec![],
        },
        one_way::ExteriorPeriod {
            name: "after".into(),
            reservoirs: vec![],
        },
    ];
    let case = serde_json::from_value(
        serde_json::json!({"schema":"soil-fem-one-way-construction/1","one_way":base,"excavations":[
            {"name":"cut","before_period":"after","bodies":["UpperClay"],"water_boundary_policy":"explicit_period_boundaries","cut_traction_policy":"explicit_period_reservoirs"}
        ]}),
    )?;
    Ok((mesh, ground, case))
}

#[test]
fn saturated_flooded_cut_matches_exact_oedometer_unloading_and_mass_theory() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, mut case) = fixture(quadratic, false)?;
        case.one_way.flow.initial_head.head_m = 1.;
        case.one_way.mechanics.pressure_on_original_surface_only = true;
        let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
            unreachable!()
        };
        for p in periods.iter_mut() {
            p.boundaries[0].condition = richards::Condition::Head {
                value: richards::AffineHead {
                    head_m: 1.,
                    gradient: [0.; 3],
                    reference_point_m: [0.; 3],
                },
            };
            p.boundaries[1].condition = richards::Condition::OutwardFlux { flux_m_s: 0. };
        }
        case.one_way.exterior_periods[1].reservoirs = vec![one_way::Reservoir {
            face_set: "floor".into(),
            level_m: 1.,
        }];
        let run = construction::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        let record = &run.excavations[0];
        assert!((record.removal.water_removed_m3 - 0.5 * 0.7 / 1.7).abs() < 1e-13);
        assert!(record.zero_release_relative_force_residual < 1e-8);
        let origin = &run
            .frames
            .iter()
            .rfind(|f| f.phase == "transient" && f.state.hydraulic.time_s <= 200.)
            .unwrap()
            .state;
        let change = 100. + 0.5 * (16. + 10. * 0.7 / 1.7) - 5.;
        for (i, p) in mesh.nodes.iter().enumerate().filter(|(_, p)| p[2] <= 0.5) {
            let actual = run.last_committed.mechanical.displacement_m[i][2]
                - origin.mechanical.displacement_m[i][2];
            assert!(
                (actual - change * p[2] / 12000.).abs() < 2e-10,
                "Tet10={quadratic} actual={actual} z={}",
                p[2]
            );
        }
        assert!(
            run.frames
                .iter()
                .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-11)
        );
        for f in run
            .frames
            .iter()
            .filter(|f| f.phase == "excavation_release")
        {
            assert_eq!(f.state.hydraulic, record.removal.state);
            assert_eq!(f.state.hydraulic.time_s, 200.);
            assert!(f.hydraulic_balance.is_none());
        }
        assert!(
            run.frames
                .iter()
                .rfind(|f| f.phase == "excavation_release")
                .unwrap()
                .holding_load_kn
                .iter()
                .all(|h| *h == 0.)
        );
    }
    Ok(())
}

#[test]
fn domain_death_retains_remaining_heads_storage_fluxes_traces_and_archives_removed_water()
-> Result<()> {
    let (mesh, ground, case) = fixture(true, true)?;
    let run = construction::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let record = &run.excavations[0];
    let transfer = &record.removal;
    assert!(!transfer.newly_exposed_corners.is_empty());
    for e in 0..mesh.tetrahedra.len() {
        if transfer.state.relative_cell_head_m[e].is_some() {
            assert_eq!(
                transfer.state.relative_cell_head_m[e],
                record.origin_hydraulic.relative_cell_head_m[e]
            );
            assert_eq!(
                transfer.state.cell_storage_offset_m3[e],
                record.origin_hydraulic.cell_storage_offset_m3[e]
            );
            assert_eq!(
                transfer.state.cell_outflow_m3_s[e],
                record.origin_hydraulic.cell_outflow_m3_s[e]
            );
        } else {
            assert!(transfer.state.cell_storage_offset_m3[e].is_none());
            assert!(transfer.state.cell_outflow_m3_s[e].is_none());
        }
    }
    for (i, j) in transfer.new_to_old_face.iter().enumerate() {
        assert_eq!(
            transfer.new_faces[i].corners,
            transfer.old_faces[*j].corners
        );
        assert_eq!(
            transfer.state.relative_face_head_m[i],
            record.origin_hydraulic.relative_face_head_m[*j]
        );
    }
    for cell in &transfer.removed_cells {
        assert_eq!(
            Some(cell.relative_head_m),
            record.origin_hydraulic.relative_cell_head_m[cell.element_index]
        );
        assert_eq!(
            Some(cell.storage_offset_m3),
            record.origin_hydraulic.cell_storage_offset_m3[cell.element_index]
        );
        assert!(cell.water_equivalent_volume_m3 > 0.);
    }
    let origin = run
        .frames
        .iter()
        .rfind(|f| f.phase == "transient" && f.state.hydraulic.time_s <= 200.)
        .unwrap();
    for &i in &record.removed_mechanical_points {
        assert_eq!(
            serde_json::to_value(&origin.state.mechanical.material_states[i])?,
            serde_json::to_value(&run.last_committed.mechanical.material_states[i])?
        );
    }
    assert!(
        run.frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-9)
    );
    assert!(
        run.initialization
            .last_committed
            .mechanical
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|v| *v != 0.))
            || origin
                .state
                .mechanical
                .material_states
                .iter()
                .any(|s| s.plastic_strain().iter().any(|v| *v != 0.))
    );
    Ok(())
}

#[test]
fn rejected_domain_candidate_restores_history_and_clean_release_path_is_exact() -> Result<()> {
    let (mesh, ground, mut case) = fixture(false, true)?;
    case.verification_reject_once_excavations = vec!["cut".into()];
    let run = construction::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.attempts
            .iter()
            .any(|a| !a.accepted && a.candidate_changed)
    );
    for a in run.attempts.iter().filter(|a| !a.accepted) {
        assert_eq!(
            a.committed_pair_sha256_before,
            a.committed_pair_sha256_after
        );
    }
    case.verification_reject_once_excavations.clear();
    case.excavations[0].initial_increment = Some(0.25);
    let clean = construction::solve(&mesh, ground, case)?;
    assert!(clean.complete, "{:?}", clean.failure_reason);
    assert_eq!(
        serde_json::to_value(run.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?,
        serde_json::to_value(clean.frames.iter().map(|f| &f.state).collect::<Vec<_>>())?
    );
    Ok(())
}

#[test]
fn construction_checkpoints_restore_nonzero_holding_domain_archive_and_next_physical_time()
-> Result<()> {
    let (mesh, ground, mut case) = fixture(false, false)?;
    case.verification_reject_once_excavations = vec!["cut".into()];
    let direct = construction::solve(&mesh, ground.clone(), case.clone())?;
    assert!(direct.complete);
    for count in [1, 2, 3, 4, 5, 6] {
        let partial = construction::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        assert!(!partial.complete);
        let saved = construction::Checkpoint::seal(partial)?;
        if saved.run.runtime.as_ref().unwrap().pending.is_some() {
            for mode in 0..4 {
                let mut data = serde_json::to_value(&saved)?;
                match mode {
                    0 => {
                        data["run"]["last_committed"]["hydraulic"]["relative_face_head_m"][0] =
                            serde_json::json!(999.)
                    }
                    1 => {
                        data["run"]["runtime"]["pending"]["record"]["removal"]["water_removed_m3"] =
                            serde_json::json!(999.)
                    }
                    2 => {
                        data["run"]["runtime"]["pending"]["record"]["initial_holding_load_kn"][0] =
                            serde_json::json!(999.)
                    }
                    _ => {
                        data["run"]["last_committed"]["mechanical"]["active_elements"][0] =
                            serde_json::json!(false)
                    }
                }
                let mut forged: construction::Checkpoint = serde_json::from_value(data)?;
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
            construction::normalized(&direct)?,
            construction::normalized(&resumed)?
        );
    }
    Ok(())
}

#[test]
fn missing_cut_boundaries_overlap_reactivation_and_corrupt_origin_are_rejected() -> Result<()> {
    let (mesh, ground, case) = fixture(false, false)?;
    let mut bad = case.clone();
    let richards::Analysis::Transient { periods } = &mut bad.one_way.flow.analysis else {
        unreachable!()
    };
    periods[1].boundaries.pop();
    assert!(construction::solve(&mesh, ground.clone(), bad).is_err());
    let mut bad = case.clone();
    bad.excavations[0].bodies = vec!["UpperClay".into(), "UpperClay".into()];
    assert!(construction::solve(&mesh, ground.clone(), bad).is_err());
    let run = construction::solve(&mesh, ground, case.clone())?;
    assert!(run.complete);
    let r = &run.excavations[0];
    let before = vec![true; mesh.tetrahedra.len()];
    let after = run.last_committed.mechanical.active_elements.clone();
    assert!(
        richards::remove_domain(&mesh, &case.one_way.flow, &after, &before, &r.removal.state)
            .is_err()
    );
    let mut wrong = r.origin_hydraulic.clone();
    wrong.cell_storage_offset_m3[0] = Some(999.);
    assert!(richards::remove_domain(&mesh, &case.one_way.flow, &before, &after, &wrong).is_err());
    Ok(())
}

#[test]
fn compressible_water_mass_removed_at_fixed_time_includes_positive_storage() -> Result<()> {
    let (mesh, ground, mut case) = fixture(true, false)?;
    case.one_way.flow.initial_head.head_m = 1.;
    for material in &mut case.one_way.flow.materials {
        material.retention.saturated_specific_storage_per_m = 0.1;
    }
    let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
        unreachable!()
    };
    for p in periods {
        p.boundaries[0].condition = richards::Condition::Head {
            value: richards::AffineHead {
                head_m: 1.,
                gradient: [0.; 3],
                reference_point_m: [0.; 3],
            },
        };
        p.boundaries[1].condition = richards::Condition::OutwardFlux { flux_m_s: 0. };
    }
    case.one_way.exterior_periods[1].reservoirs = vec![one_way::Reservoir {
        face_set: "floor".into(),
        level_m: 1.,
    }];
    let run = construction::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let r = &run.excavations[0].removal;
    assert!((r.water_removed_m3 - (0.5 * 0.7 / 1.7 + 0.1 * 0.125)).abs() < 1e-13);
    assert!((r.water_remaining_m3 - (0.5 * 0.7 / 1.7 + 0.1 * 0.375)).abs() < 1e-13);
    assert!(
        run.frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-10)
    );
    Ok(())
}

#[test]
fn two_local_removals_preserve_separate_archives_and_cumulative_water_budget() -> Result<()> {
    let (mut mesh, mut ground, mut case) = fixture(true, false)?;
    ground.excavation_bodies = vec!["CutA".into(), "CutB".into()];
    case.one_way.mechanics.settings.kinematics = soil_fem_3d::stages::Kinematics::SideRollers;
    mesh.excavation_sets.insert("CutA".into(), vec![9, 10, 11]);
    mesh.excavation_sets.insert("CutB".into(), vec![6, 7, 8]);
    let owners = soil_fem_3d::boundary::surface_owners(&mesh);
    let active_a = (0..12).map(|e| e < 9).collect::<Vec<_>>();
    mesh.face_sets.insert(
        "cut_a".into(),
        owners
            .iter()
            .filter(|(_, o)| o.len() == 2 && o.iter().filter(|&&e| active_a[e]).count() == 1)
            .map(|(k, _)| *k)
            .collect(),
    );
    let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
        unreachable!()
    };
    periods[0].end_time_s = 100.;
    let mut middle = periods[0].clone();
    middle.name = "middle".into();
    middle.end_time_s = 200.;
    middle.boundaries.truncate(1);
    middle.boundaries.push(richards::Boundary {
        face_set: "cut_a".into(),
        condition: richards::Condition::AtmosphericNoPonding {
            rainfall_vertical_m_s: 0.,
        },
    });
    periods[1].end_time_s = 300.;
    periods.insert(1, middle);
    let mut event = case.excavations[0].clone();
    event.name = "first".into();
    event.before_period = "middle".into();
    event.bodies = vec!["CutA".into()];
    case.excavations[0].bodies = vec!["CutB".into()];
    case.excavations.insert(0, event);
    case.one_way.exterior_periods.insert(
        1,
        one_way::ExteriorPeriod {
            name: "middle".into(),
            reservoirs: vec![],
        },
    );
    let run = construction::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(run.excavations.len(), 2);
    let removed = run
        .excavations
        .iter()
        .map(|r| r.removal.water_removed_m3)
        .sum::<f64>();
    assert_eq!(
        run.frames.last().unwrap().cumulative_removed_water_m3,
        removed
    );
    assert!(
        run.frames
            .iter()
            .all(|f| f.signed_global_water_budget_error_m3.abs() < 1e-10)
    );
    let first = &run.excavations[0].removal;
    let second = &run.excavations[1].origin_hydraulic;
    for cell in &first.removed_cells {
        assert!(second.relative_cell_head_m[cell.element_index].is_none());
    }
    let count = run
        .attempts
        .iter()
        .take_while(|a| a.name != "cut")
        .filter(|a| a.accepted)
        .count()
        + 1;
    let partial = construction::solve_until(&mesh, ground.clone(), case.clone(), count)?;
    assert!(partial.runtime.as_ref().unwrap().pending.is_some());
    let resumed = construction::Checkpoint::seal(partial)?.resume(&mesh, ground, case, None)?;
    assert_eq!(
        construction::normalized(&run)?,
        construction::normalized(&resumed)?
    );
    Ok(())
}

#[test]
fn genuine_post_cut_hydraulic_exhaustion_preserves_the_last_accepted_domain_and_mass_archive()
-> Result<()> {
    let (mesh, ground, mut case) = fixture(false, false)?;
    case.one_way.flow.settings.maximum_iterations = 1;
    case.one_way.flow.settings.minimum_time_step_s = 100.;
    let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
        unreachable!()
    };
    periods[1].boundaries[1].condition = richards::Condition::AtmosphericNoPonding {
        rainfall_vertical_m_s: 1e-3,
    };
    let run = construction::solve(&mesh, ground, case)?;
    assert!(!run.complete);
    assert!(
        run.failure_reason.as_ref().unwrap().contains("iteration"),
        "{:?}",
        run.failure_reason
    );
    assert_eq!(run.excavations.len(), 1);
    assert_eq!(run.last_committed.hydraulic.time_s, 200.);
    assert_eq!(
        run.last_committed.hydraulic,
        run.excavations[0].removal.state
    );
    assert_eq!(
        serde_json::to_value(&run.last_committed)?,
        serde_json::to_value(&run.frames.last().unwrap().state)?
    );
    assert!(
        run.frames
            .last()
            .unwrap()
            .holding_load_kn
            .iter()
            .all(|h| *h == 0.)
    );
    let failed = run.attempts.last().unwrap();
    assert!(!failed.accepted && failed.phase == "transient");
    assert_eq!(
        failed.committed_pair_sha256_before,
        failed.committed_pair_sha256_after
    );
    Ok(())
}

#[test]
fn unchanged_domain_reuses_geometry_but_applies_each_new_period_boundary() -> Result<()> {
    let (mesh, ground, mut case) = fixture(false, false)?;
    let richards::Analysis::Transient { periods } = &mut case.one_way.flow.analysis else {
        unreachable!()
    };
    let mut first = periods[0].clone();
    first.name = "first".into();
    first.end_time_s = 100.;
    periods[0].boundaries[0].condition = richards::Condition::Head {
        value: richards::AffineHead {
            head_m: 0.2,
            gradient: [0.; 3],
            reference_point_m: [0.; 3],
        },
    };
    periods.insert(0, first);
    case.one_way.exterior_periods.insert(
        0,
        one_way::ExteriorPeriod {
            name: "first".into(),
            reservoirs: vec![],
        },
    );
    let mut standalone = case.one_way.flow.clone();
    let richards::Analysis::Transient { periods } = &mut standalone.analysis else {
        unreachable!()
    };
    periods.truncate(2);
    let reference = richards::solve(&mesh, standalone, &vec![true; mesh.tetrahedra.len()], None)?;
    assert!(reference.complete);
    let run = construction::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let precut = run
        .frames
        .iter()
        .find(|f| f.phase == "transient" && f.state.hydraulic.time_s == 200.)
        .unwrap();
    assert_eq!(precut.state.hydraulic, reference.last_committed);
    Ok(())
}
