use nalgebra::{SVector, Vector3};
use serde_json::json;
use soil_fem_3d::{
    Result,
    beams::{self, Configuration, Formulation, Kernel, Section},
    checkpoint::Checkpoint,
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};

fn section(timoshenko: bool) -> Section {
    Section {
        id: "steel".into(),
        young_kpa: 2e8,
        poisson: 0.3,
        area_m2: 0.02,
        iy_m4: 0.0001,
        iz_m4: 0.0002,
        torsion_j_m4: 0.00005,
        unit_weight_kn_m3: 0.,
        formulation: if timoshenko {
            Formulation::Timoshenko {
                shear_area_y_m2: 0.016,
                shear_area_z_m2: 0.014,
            }
        } else {
            Formulation::EulerBernoulli
        },
    }
}
fn configuration(base: [f64; 3], tip: [f64; 3], timo: bool) -> Result<Configuration> {
    Ok(serde_json::from_value(
        json!({"schema":"soil-fem-beams/1","nodes":[{"id":"base","position_m":base,"soil_node_set":"BeamBase"},{"id":"tip","position_m":tip,"soil_node_set":null}],
        "sections":[section(timo)],"members":[{"id":"beam","nodes":["base","tip"],"section_id":"steel","local_y":[0,1,0],"initially_active":true}],
        "rotation_length_scale_m":2.,"maximum_condensed_dofs":100,
        "rotation_constraints":[{"node_id":"base","direction":[1,0,0],"value":0.},{"node_id":"base","direction":[0,1,0],"value":0.},{"node_id":"base","direction":[0,0,1],"value":0.}],
        "nodal_loads":[{"node_id":"tip","constant_force_kn":[1.,2.,3.],"constant_moment_kn_m":[4.,0.,0.]}]}),
    )?)
}
fn mesh() -> Result<Mesh> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    mesh.node_sets.insert("BeamBase".into(), vec![0]);
    Ok(mesh)
}

#[test]
fn kernels_pass_six_rigid_modes_axial_patch_and_arbitrary_orientation() -> Result<()> {
    for timo in [false, true] {
        let a = [0.3, -0.8, 2.];
        let b = [2.1, 1.3, 4.7];
        let s = section(timo);
        let kernel = Kernel::new(a, b, [0.1, 1., 0.2], &s, [[false; 3]; 2])?;
        for mode in 0..6 {
            let mut d = SVector::<f64, 12>::zeros();
            for (end, x) in [a, b].iter().enumerate() {
                if mode < 3 {
                    d[6 * end + mode] = 1.;
                } else {
                    let mut w = Vector3::zeros();
                    w[mode - 3] = 1.;
                    let t = w.cross(&Vector3::from(*x));
                    for q in 0..3 {
                        d[6 * end + q] = t[q];
                        d[6 * end + 3 + q] = w[q];
                    }
                }
            }
            assert!((kernel.stiffness * d).norm() < 1e-12 * kernel.stiffness.norm() * d.norm());
        }
        let mut axial = SVector::<f64, 12>::zeros();
        for q in 0..3 {
            axial[6 + q] = (b[q] - a[q]) * 1e-5;
        }
        let (_, action) = kernel.recover(&axial, [0.; 3], 0.);
        assert!((action[6] - s.young_kpa * s.area_m2 * 1e-5).abs() < 1e-8);
        assert!((action[0] + action[6]).abs() < 1e-8);
    }
    Ok(())
}

#[test]
fn cantilever_axial_biaxial_bending_torsion_and_scaling_match_closed_forms() -> Result<()> {
    let mesh = mesh()?;
    let active = vec![true; mesh.nodes.len()];
    let u = vec![[0.; 3]; mesh.nodes.len()];
    for timo in [false, true] {
        for scaling in [0.5, 2., 10.] {
            let mut c = configuration([0.; 3], [2., 0., 0.], timo)?;
            c.rotation_length_scale_m = scaling;
            let model = beams::Model::new(&mesh, c)?;
            let result = model.evaluate(&model.initial_state(), &u, &active, 0.)?;
            let s = section(timo);
            let g = s.young_kpa / (2. * (1. + s.poisson));
            let d = result.state.node_displacement_m[1];
            let r = result.state.node_rotation_rad[1];
            let shear = if timo {
                [2. * 2. / (g * 0.016), 3. * 2. / (g * 0.014)]
            } else {
                [0., 0.]
            };
            let expected = [
                1. * 2. / (s.young_kpa * s.area_m2),
                2. * 8. / (3. * s.young_kpa * s.iz_m4) + shear[0],
                3. * 8. / (3. * s.young_kpa * s.iy_m4) + shear[1],
            ];
            for a in 0..3 {
                assert!(
                    (d[a] - expected[a]).abs() < 1e-10 * expected[a].abs().max(1e-10),
                    "{d:?} != {expected:?}"
                );
            }
            assert!((r[0] - 4. * 2. / (g * s.torsion_j_m4)).abs() < 1e-12);
            assert!((r[1] + 3. * 4. / (2. * s.young_kpa * s.iy_m4)).abs() < 1e-12);
            assert!((r[2] - 2. * 4. / (2. * s.young_kpa * s.iz_m4)).abs() < 1e-12);
            assert!(result.report.true_relative_interior_residual < 1e-11);
            for a in 0..3 {
                assert!((result.interface_internal_kn[a] + [1., 2., 3.][a]).abs() < 1e-8);
            }
            let applied = Vector3::from(result.report.external_moment_origin_kn_m)
                + Vector3::from(result.report.support_moment_origin_kn_m);
            assert!(applied.norm() < 1e-8);
        }
    }
    Ok(())
}

#[test]
fn distributed_loading_and_end_releases_recover_fixed_and_simply_supported_moments() -> Result<()> {
    let mesh = mesh()?;
    let mut c = configuration([0.; 3], [2., 0., 0.], false)?;
    c.nodal_loads.clear();
    c.members[0].constant_line_load_kn_m = [0., -10., 0.];
    let model = beams::Model::new(&mesh, c.clone())?;
    let result = model.evaluate(
        &model.initial_state(),
        &vec![[0.; 3]; mesh.nodes.len()],
        &vec![true; mesh.nodes.len()],
        0.,
    )?;
    let expected = -10. * 2_f64.powi(4) / (8. * section(false).young_kpa * section(false).iz_m4);
    assert!((result.state.node_displacement_m[1][1] - expected).abs() < 1e-12);
    c.members[0].rotational_releases = [[false, false, true]; 2];
    c.translation_constraints = vec![serde_json::from_value(
        json!({"node_id":"tip","direction":[0,1,0],"value":0.}),
    )?];
    c.rotation_constraints.push(serde_json::from_value(
        json!({"node_id":"tip","direction":[0,0,1],"value":0.}),
    )?);
    let model = beams::Model::new(&mesh, c)?;
    let result = model.evaluate(
        &model.initial_state(),
        &vec![[0.; 3]; mesh.nodes.len()],
        &vec![true; mesh.nodes.len()],
        0.,
    )?;
    let m = &result.report.members[0];
    assert!(m.local_end_actions[5].abs() < 1e-8 && m.local_end_actions[11].abs() < 1e-8);
    assert!((m.section_forces_at_start_mid_end[1][5] - 10. * 4. / 8.).abs() < 1e-8);
    assert!(
        (m.local_end_actions[1] - 10.).abs() < 1e-8 && (m.local_end_actions[7] - 10.).abs() < 1e-8
    );
    assert!(m.local_end_relative_rotation_rad[0][2].abs() > 1e-6);
    Ok(())
}

#[test]
fn unsupported_torsion_and_invalid_geometry_are_rejected_without_stabilization() -> Result<()> {
    let mesh = mesh()?;
    let mut c = configuration([0.; 3], [2., 0., 0.], false)?;
    c.rotation_constraints
        .retain(|r| r.direction != [1., 0., 0.]);
    let model = beams::Model::new(&mesh, c)?;
    let state = model.initial_state();
    assert!(
        model
            .evaluate(
                &state,
                &vec![[0.; 3]; mesh.nodes.len()],
                &vec![true; mesh.nodes.len()],
                0.
            )
            .is_err()
    );
    assert_eq!(state, model.initial_state());
    assert!(
        Kernel::new(
            [0.; 3],
            [0.; 3],
            [0., 1., 0.],
            &section(false),
            [[false; 3]; 2]
        )
        .is_err()
    );
    assert!(
        Kernel::new(
            [0.; 3],
            [2., 0., 0.],
            [1., 0., 0.],
            &section(false),
            [[false; 3]; 2]
        )
        .is_err()
    );
    Ok(())
}

fn coupled(mcc: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh = mesh()?;
    let node = mesh.nodes.iter().position(|p| *p == [0., 0., 1.]).unwrap();
    mesh.node_sets.insert("BeamBase".into(), vec![node]);
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
    let mut beam = configuration([0., 0., 1.], [2., 0., 1.], false)?;
    beam.members[0].initially_active = false;
    beam.nodal_loads[0].constant_force_kn = [0.; 3];
    beam.nodal_loads[0].constant_moment_kn_m = [0.; 3];
    beam.nodal_loads[0].boundary_force_kn = [0., 0., -5.];
    case.beam_system = Some(beam);
    case.settings.initial_increment = 0.25;
    case.stages = serde_json::from_value(json!([
        {"kind":"install_beams","name":"install","members":["beam"],"placement":"cad_new_nodes"},
        {"kind":"boundary","name":"load_beam","target_factor":1.},
        {"kind":"boundary","name":"unload_beam","target_factor":0.},
        {"kind":"remove_beams","name":"remove","members":["beam"]}]))?;
    case.verification_reject_once_stages = vec!["install".into(), "load_beam".into()];
    Ok((mesh, ground, case))
}

#[test]
fn condensed_beam_is_coupled_to_mc_and_mcc_force_moment_and_exact_restart() -> Result<()> {
    for mcc in [false, true] {
        let (mesh, ground, case) = coupled(mcc)?;
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete, "{:?}", full.failure_reason);
        assert_eq!(full.schema, "soil-fem-stage-result/7");
        if !mcc {
            assert!(full.frames.iter().flat_map(|f|f.state.material_states.iter()).all(|s|matches!(s,soil::State::MohrCoulomb(s) if s.equivalent_plastic_strain==0.)));
        }
        if mcc {
            assert!(full.attempts.iter().any(|a| !a.accepted
                && a.rejected_trial_structure_changed
                && a.rejected_trial_plastic_changed));
        }
        for frame in &full.frames {
            assert!(frame.force_balance_error_kn.iter().all(|x| x.abs() < 1e-6));
            assert!(
                frame
                    .global_moment_balance_error_kn_m
                    .unwrap()
                    .iter()
                    .all(|x| x.abs() < 1e-6)
            );
        }
        assert!(full.attempts.iter().any(|a| !a.accepted
            && a.rejected_trial_structure_changed
            && a.committed_state_unchanged_on_rejection));
        let initial = full
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
        for count in [initial + 1, initial + install + 1] {
            let cp = Checkpoint::seal(
                &mesh,
                stages::solve_until(&mesh, ground.clone(), case.clone(), count)?,
            )?;
            let resumed = cp.resume(&mesh, ground.clone(), case.clone(), None)?;
            assert_eq!(
                serde_json::to_vec(&full.last_committed)?,
                serde_json::to_vec(&resumed.last_committed)?
            );
            assert_eq!(
                serde_json::to_vec(&full.frames)?,
                serde_json::to_vec(&resumed.frames)?
            );
        }
        if !mcc {
            let mut reference = case.clone();
            reference.beam_system = None;
            reference.stages = serde_json::from_value(
                json!([{"kind":"boundary","name":"load","target_factor":1.},{"kind":"boundary","name":"unload","target_factor":0.}]),
            )?;
            reference.verification_reject_once_stages.clear();
            reference.boundary.as_mut().unwrap().nodal_loads =
                serde_json::from_value(json!([{"node_set":"BeamBase","total_force_kn":[0,0,-5]}]))?;
            let direct = stages::solve(&mesh, ground, reference)?;
            assert!(direct.complete);
            for (a, b) in full.frames[2]
                .state
                .displacement_m
                .iter()
                .flatten()
                .zip(direct.frames[1].state.displacement_m.iter().flatten())
            {
                assert!(
                    (a - b).abs() < 1e-8,
                    "beam soil u={a:e}, direct u={b:e}, beam report={:?}",
                    full.frames[2].beam_results
                );
            }
        }
    }
    Ok(())
}

#[test]
fn timoshenko_uniform_load_and_nonzero_constraint_lifting_match_theory() -> Result<()> {
    let mesh = mesh()?;
    let active = vec![true; mesh.nodes.len()];
    let zero = vec![[0.; 3]; mesh.nodes.len()];
    for shear_area in [0.00001, 0.016, 1e9] {
        let mut c = configuration([0.; 3], [2., 0., 0.], true)?;
        c.sections[0].formulation = Formulation::Timoshenko {
            shear_area_y_m2: shear_area,
            shear_area_z_m2: shear_area,
        };
        c.nodal_loads.clear();
        c.members[0].constant_line_load_kn_m = [0., -10., 0.];
        let model = beams::Model::new(&mesh, c)?;
        let result = model.evaluate(&model.initial_state(), &zero, &active, 0.)?;
        let s = section(true);
        let g = s.young_kpa / (2. * (1. + s.poisson));
        let expected = -10.
            * (2_f64.powi(4) / (8. * s.young_kpa * s.iz_m4)
                + 2_f64.powi(2) / (2. * g * shear_area));
        assert!((result.state.node_displacement_m[1][1] - expected).abs() < 1e-10);
    }
    let mut c = configuration([0.; 3], [2., 0., 0.], true)?;
    c.nodal_loads.clear();
    let theta = Vector3::new(0.01, 0.02, -0.03);
    for (i, condition) in c.rotation_constraints.iter_mut().enumerate() {
        condition.value = theta[i];
    }
    let mut displacement = zero;
    displacement[0] = [0.004, -0.003, 0.002];
    let model = beams::Model::new(&mesh, c)?;
    let result = model.evaluate(&model.initial_state(), &displacement, &active, 0.)?;
    let expected = Vector3::from(displacement[0]) + theta.cross(&Vector3::new(2., 0., 0.));
    assert!((Vector3::from(result.state.node_displacement_m[1]) - expected).norm() < 1e-12);
    assert!((Vector3::from(result.state.node_rotation_rad[1]) - theta).norm() < 1e-12);
    assert!(result.interface_internal_kn.norm() < 1e-8);
    Ok(())
}

#[test]
fn rehashed_checkpoint_rejects_corrupted_recovered_beam_kinematics() -> Result<()> {
    let (mesh, ground, case) = coupled(true)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete);
    let geostatic = full
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let installed = full
        .attempts
        .iter()
        .filter(|a| a.accepted && a.stage == "install")
        .count();
    let mut partial = stages::solve_until(
        &mesh,
        ground.clone(),
        case.clone(),
        geostatic + installed + 1,
    )?;
    partial
        .last_committed
        .beams
        .as_mut()
        .unwrap()
        .node_rotation_rad[1][0] += 0.001;
    // Re-seal the altered payload, so this exercises physical cache validation
    // beyond the ordinary JSON checksum guard.
    let checkpoint = Checkpoint::seal(&mesh, partial)?;
    assert!(checkpoint.resume(&mesh, ground, case, None).is_err());
    Ok(())
}
