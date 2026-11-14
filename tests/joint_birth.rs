use serde_json::json;
use soil_fem_3d::{
    Result,
    checkpoint::{self, Checkpoint},
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};

fn fixture(policy: &str, quadratic: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/joint-column-mesh.json"))?;
    let mesh = if quadratic {
        mesh.elevate_straight()?
    } else {
        mesh
    };
    let mut ground: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-joint-mc.json"))?;
    ground["initialization"]["geostatic"]["base_loads"] = json!({});
    for body in ["LowerClay", "UpperClay"] {
        ground["initialization"]["geostatic"]["seed"]["states"][body]["effective_stress_compression_positive_kpa"] =
            json!([0., 0., 0., 0., 0., 0.]);
    }
    for m in ground["materials"].as_array_mut().unwrap() {
        m["unit_weight_kn_m3"] = json!(1e-20);
    }
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-normal-column.json"))?;
    case["joint_system"]["joints"][0]["birth_policy"] = json!(policy);
    case["verification_reject_once_stages"] = json!([]);
    case["stages"] = json!([
        {"kind":"pressure","name":"preload","target_pressure_kpa":50.},
        {"kind":"joint_birth","name":"contact","joints":["Seam"]},
        {"kind":"pressure","name":"reload","target_pressure_kpa":80.},
        {"kind":"pressure","name":"unload","target_pressure_kpa":20.}
    ]);
    Ok((
        mesh,
        serde_json::from_value(ground)?,
        serde_json::from_value(case)?,
    ))
}
#[test]
fn mode_a_column_preserves_displacements_all_nodal_forces_and_fresh_energy() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, ground, case) = fixture("transfer_force", quadratic)?;
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        let before = &run.frames[1];
        let born = &run.frames[2];
        for (a, b) in before
            .state
            .displacement_m
            .iter()
            .flatten()
            .zip(born.state.displacement_m.iter().flatten())
        {
            assert!((a - b).abs() < 1e-10);
        }
        let prepared = &run.transition_audits[0].prepared_state;
        assert!(
            run.transition_audits[0]
                .initial_holding_load_kn
                .iter()
                .all(|f| *f == 0.)
        );
        assert_eq!(
            born.state.applied_nodal_load_kn,
            before.state.applied_nodal_load_kn
        );
        assert_eq!(prepared.displacement_m, before.state.displacement_m);
        for p in &prepared.joints.as_ref().unwrap().points {
            assert!(
                (p.traction_kpa[2] + 50.).abs() < 1e-8,
                "{:?}",
                p.traction_kpa
            );
            assert_eq!(p.plastic_slip_m, [0.; 3]);
            assert_eq!(p.accumulated_plastic_slip_m, 0.);
            assert_eq!(p.plastic_dissipation_kpa_m, 0.);
            assert_eq!(p.numerical_dissipation_kpa_m, 0.);
            assert_eq!(p.endpoint_work_kpa_m, 0.);
        }
        for frame in &run.frames {
            assert!(frame.relative_force_residual < 1e-8);
            assert!(frame.maximum_constraint_error_m < 1e-10);
        }
        assert!(run.transition_audits[0].zero_release_relative_force_residual < 1e-10);
    }
    Ok(())
}
#[test]
fn mode_b_birth_reference_closes_by_pressure_over_kn() -> Result<()> {
    let (mesh, ground, case) = fixture("stress_free", true)?;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for p in &run.transition_audits[0]
        .prepared_state
        .joints
        .as_ref()
        .unwrap()
        .points
    {
        assert_eq!(p.traction_kpa, [0.; 3]);
        assert_eq!(p.birth.as_ref().unwrap().jump_m, p.jump_m);
    }
    for p in &run.frames[2].state.joints.as_ref().unwrap().points {
        let reference = p.birth.as_ref().unwrap();
        assert!((p.jump_m[2] - reference.jump_m[2] + 50. / 1e6).abs() < 1e-10);
    }
    Ok(())
}
#[test]
fn dry_new_mc_body_and_contact_have_separate_births() -> Result<()> {
    let (mesh, mut ground, case) = fixture("stress_free", true)?;
    for m in &mut ground.materials {
        m.unit_weight_kn_m3 = 20.;
    }
    let mut value = serde_json::to_value(case)?;
    let seed = soil::State::MohrCoulomb(soil_fem_3d::mc::State {
        effective_stress_compression_positive_kpa: [0.; 6],
        plastic_strain_compression_positive: [0.; 6],
        void_ratio: 0.7,
        equivalent_plastic_strain: 0.,
        plastic_work_kpa: 0.,
        last_return: soil_fem_3d::mc::Return::Elastic,
    });
    value["initially_inactive_bodies"] = json!(["UpperClay"]);
    value["stages"] = json!([
        {"kind":"birth","name":"fill","bodies":["UpperClay"],"states":{"UpperClay":seed},"placement":"cad_new_nodes"},
        {"kind":"joint_birth","name":"contact","joints":["Seam"]},
        {"kind":"pressure","name":"load","target_pressure_kpa":10.}
    ]);
    let run = stages::solve(&mesh, ground, serde_json::from_value(value)?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    assert!(
        run.frames[0]
            .state
            .displacement_m
            .iter()
            .flatten()
            .any(|v| v.abs() > 1e-8)
    );
    for p in &run.transition_audits[1]
        .prepared_state
        .joints
        .as_ref()
        .unwrap()
        .points
    {
        assert_eq!(p.traction_kpa, [0.; 3]);
        assert_eq!(p.plastic_slip_m, [0.; 3]);
    }
    for p in &run
        .frames
        .last()
        .unwrap()
        .state
        .joints
        .as_ref()
        .unwrap()
        .points
    {
        assert!((p.elastic_jump()[2] + 20. / 1e6).abs() < 1e-10, "{:?}", p);
    }
    Ok(())
}
fn without_clock(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::Object(m) => {
            m.remove("elapsed_seconds");
            for v in m.values_mut() {
                without_clock(v);
            }
        }
        serde_json::Value::Array(a) => {
            for v in a {
                without_clock(v);
            }
        }
        _ => {}
    }
}
#[test]
fn every_accepted_checkpoint_replays_and_rehashed_birth_fields_are_rejected() -> Result<()> {
    let (mesh, ground, mut case) = fixture("transfer_force", false)?;
    case.verification_reject_once_stages = vec!["contact".into(), "reload".into()];
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete, "{:?}", full.failure_reason);
    let count = full.attempts.iter().filter(|a| a.accepted).count()
        + full
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .count();
    let mut expected = serde_json::to_value(&full)?;
    without_clock(&mut expected);
    for limit in 1..=count {
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), limit)?;
        if partial.last_committed.stage_index == 2 {
            assert_eq!(
                partial.last_committed.applied_nodal_load_kn,
                full.frames[1].state.applied_nodal_load_kn,
                "mode A must not change any external load, even at intermediate progress"
            );
            assert!(
                partial
                    .last_committed
                    .transition_release_load_kn
                    .iter()
                    .all(|f| *f == 0.)
            );
        }
        let checkpoint = Checkpoint::seal(&mesh, partial)?;
        let encoded = serde_json::to_vec(&checkpoint)?;
        let resumed = serde_json::from_slice::<Checkpoint>(&encoded)?.resume(
            &mesh,
            ground.clone(),
            case.clone(),
            None,
        )?;
        let mut actual = serde_json::to_value(&resumed)?;
        without_clock(&mut actual);
        assert_eq!(actual, expected, "limit {limit}");
        if checkpoint
            .run
            .last_committed
            .joints
            .as_ref()
            .unwrap()
            .lifecycle[0]
            .generation
            > 0
        {
            for kind in 0..4 {
                let mut changed: Checkpoint = serde_json::from_slice(&encoded)?;
                let s = changed.run.last_committed.joints.as_mut().unwrap();
                match kind {
                    0 => s.lifecycle[0].generation += 1,
                    1 => s.points[0].birth.as_mut().unwrap().jump_m[0] += 0.01,
                    2 => s.points[0].birth.as_mut().unwrap().initial_traction_kpa[2] -= 1.,
                    _ => s.lifecycle[0].policy = Some(soil_fem_3d::joints::BirthPolicy::StressFree),
                }
                changed.run_sha256 = checkpoint::hash_json(&changed.run)?;
                assert!(
                    changed
                        .resume(&mesh, ground.clone(), case.clone(), None)
                        .is_err()
                );
            }
        }
    }
    Ok(())
}
#[test]
fn every_nodal_projection_and_cone_rejection_is_checked_before_birth() -> Result<()> {
    for quadratic in [false, true] {
        let (mesh, _, case) = fixture("transfer_force", quadratic)?;
        let j = soil_fem_3d::joints::Model::new(&mesh, case.joint_system.unwrap())?;
        let active = vec![true; mesh.tetrahedra.len()];
        let u = vec![[0.; 3]; mesh.nodes.len()];
        let rows: Vec<_> = (0..3 * u.len()).map(|g| vec![(g, 1.)]).collect();
        let force = |q: [f64; 3]| {
            let mut f = vec![0.; 3 * u.len()];
            for p in j.layout() {
                for ((a, b), shape) in p.nodes_a.iter().zip(&p.nodes_b).zip(&p.shape) {
                    for d in 0..3 {
                        f[3 * a + d] -= p.weight_m2 * shape * q[d];
                        f[3 * b + d] += p.weight_m2 * shape * q[d];
                    }
                }
            }
            f
        };
        let initial = j.initial_state()?;
        let ids = vec!["Seam".to_owned()];
        let good = force([12., 5., -100.]);
        let born = j.birth(&initial, &u, &active, &ids, &good, &rows)?;
        let assembled = j.evaluate(&born, &u, &active, None, false)?;
        for (a, b) in assembled.internal_kn.iter().zip(&good) {
            assert!((a - b).abs() < 1e-10);
        }
        assert!(
            j.birth(&initial, &u, &active, &ids, &force([0., 0., 100.]), &rows)
                .is_err()
        );
        assert!(
            j.birth(&initial, &u, &active, &ids, &force([31., 0., -100.]), &rows)
                .is_err()
        );
        // Keep both resultant and moment unchanged: a resultant-only gate
        // cannot detect this incompatible same-side nodal force pattern.
        let [a, b] = [j.layout()[0].nodes_a[0], j.layout()[0].nodes_a[1]];
        let direction =
            nalgebra::Vector3::from(mesh.nodes[b]) - nalgebra::Vector3::from(mesh.nodes[a]);
        let mut bad = good.clone();
        for axis in 0..3 {
            bad[3 * a + axis] += direction[axis];
            bad[3 * b + axis] -= direction[axis];
        }
        let resultant_and_moment = |f: &[f64]| {
            mesh.nodes.iter().enumerate().fold(
                (nalgebra::Vector3::zeros(), nalgebra::Vector3::zeros()),
                |(sum, moment), (n, &x)| {
                    let force = nalgebra::Vector3::from_row_slice(&f[3 * n..3 * n + 3]);
                    (
                        sum + force,
                        moment + nalgebra::Vector3::from(x).cross(&force),
                    )
                },
            )
        };
        let (before_force, before_moment) = resultant_and_moment(&good);
        let (after_force, after_moment) = resultant_and_moment(&bad);
        assert!((before_force - after_force).norm() < 1e-12);
        assert!((before_moment - after_moment).norm() < 1e-12);
        assert!(j.birth(&initial, &u, &active, &ids, &bad, &rows).is_err());
    }
    Ok(())
}
#[test]
fn death_freezes_slip_and_dissipation_and_rebirth_has_no_old_reference() -> Result<()> {
    let (mesh, _, case) = fixture("stress_free", true)?;
    let j = soil_fem_3d::joints::Model::new(&mesh, case.joint_system.unwrap())?;
    let active = vec![true; mesh.tetrahedra.len()];
    let mut u = vec![[0.; 3]; mesh.nodes.len()];
    let rows: Vec<_> = (0..3 * u.len()).map(|g| vec![(g, 1.)]).collect();
    let zero = vec![0.; 3 * u.len()];
    let ids = vec!["Seam".to_owned()];
    let first = j.birth(&j.initial_state()?, &u, &active, &ids, &zero, &rows)?;
    for (ei, &r) in mesh.region_ids.iter().enumerate() {
        if r == 2 {
            for n in mesh.cell_nodes(ei) {
                u[n] = [0.003, 0.001, -0.0001];
            }
        }
    }
    let slipped = j.evaluate(&first, &u, &active, None, false)?.state;
    assert!(
        slipped
            .points
            .iter()
            .all(|p| p.accumulated_plastic_slip_m > 0. && p.plastic_dissipation_kpa_m > 0.)
    );
    let removed = j.death(&slipped, &ids)?;
    let archive = removed.lifecycle[0].archived.clone();
    for x in &mut u {
        for v in x {
            *v += 0.02;
        }
    }
    let dormant = j.evaluate(&removed, &u, &active, None, false)?;
    assert!(dormant.internal_kn.iter().all(|v| *v == 0.));
    assert_eq!(dormant.state, removed);
    // The exact pre-birth tie is re-established independently of archived slip.
    u.fill([0.02; 3]);
    let born = j.birth(&removed, &u, &active, &ids, &zero, &rows)?;
    assert_eq!(born.lifecycle[0].generation, 2);
    assert_eq!(born.lifecycle[0].archived, archive);
    for (i, p) in &archive[0].points {
        assert_eq!(*p, slipped.points[*i]);
    }
    for p in &born.points {
        assert_eq!(p.traction_kpa, [0.; 3]);
        assert_eq!(p.plastic_slip_m, [0.; 3]);
        assert_eq!(p.accumulated_plastic_slip_m, 0.);
        assert_eq!(p.plastic_dissipation_kpa_m, 0.);
        assert_eq!(p.numerical_dissipation_kpa_m, 0.);
    }
    Ok(())
}

#[test]
fn water_birth_rebuilds_face_loads_and_rehashed_tampering_is_refused() -> Result<()> {
    let mesh: Mesh = serde_json::from_str(include_str!(
        "../examples/joint-water-column-tet10-mesh.json"
    ))?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-joint-water-mc.json"))?;
    ground.initialization.water_level_m = Some(1.5);
    for policy in ["column", "flow", "soil", "prescribed"] {
        let source = std::fs::read_to_string(format!(
            "{}/examples/stages-joint-water-{policy}.json",
            env!("CARGO_MANIFEST_DIR")
        ))?;
        let mut case: serde_json::Value = serde_json::from_str(&source)?;
        case["joint_system"]["joints"][0]["birth_policy"] = json!("transfer_force");
        let first = if policy == "column" {
            json!({"kind":"pressure","name":"preload","target_pressure_kpa":50.})
        } else {
            let mut s = case["stages"][0].clone();
            s["exterior_water_level_m"] = json!(1.5);
            s
        };
        case["stages"] = json!([first,{"kind":"joint_birth","name":"contact","joints":["Seam"]}]);
        let case: StageCase = serde_json::from_value(case)?;
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete, "{policy}: {:?}", full.failure_reason);
        let count = full
            .geostatic
            .as_ref()
            .unwrap()
            .attempts
            .iter()
            .filter(|a| a.accepted)
            .count()
            + full
                .attempts
                .iter()
                .filter(|a| a.accepted && a.stage != "contact")
                .count()
            + 1;
        let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), count)?;
        assert_eq!(
            partial.last_committed.applied_nodal_load_kn,
            full.frames[1].state.applied_nodal_load_kn,
            "{policy}: mode A physical RHS changed"
        );
        assert!(
            partial
                .last_committed
                .transition_release_load_kn
                .iter()
                .all(|f| *f == 0.)
        );
        let saved = Checkpoint::seal(&mesh, partial)?;
        let encoded = serde_json::to_vec(&saved)?;
        saved.verify(&mesh, &ground, &case)?;
        for kind in 0..3 {
            let mut changed: Checkpoint = serde_json::from_slice(&encoded)?;
            let field = changed.run.last_committed.joint_water.as_mut().unwrap();
            match kind {
                0 => field.face_nodal_force_a_kn[0][0][2] += 0.1,
                1 => field.nodal_load_kn[2] += 0.1,
                _ => field.physical_pressure_kpa[0] += 1.,
            }
            changed.run_sha256 = checkpoint::hash_json(&changed.run)?;
            assert!(
                changed
                    .resume(&mesh, ground.clone(), case.clone(), None)
                    .is_err(),
                "{policy} field {kind}"
            );
        }
    }
    Ok(())
}

#[test]
fn body_death_then_dry_rebirth_keeps_the_removed_slip_generation_frozen() -> Result<()> {
    let mesh: Mesh =
        serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
    let ground: Ground = serde_json::from_str(include_str!("../examples/ground-joint-mc.json"))?;
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
    case["joint_system"]["joints"][0]["birth_policy"] = json!("stress_free");
    case["verification_reject_once_stages"] = json!(["fill", "rebirth"]);
    let seed = json!({"effective_stress_compression_positive_kpa":[0.,0.,0.,0.,0.,0.],"plastic_strain_compression_positive":[0.,0.,0.,0.,0.,0.],"void_ratio":0.7,"equivalent_plastic_strain":0.,"plastic_work_kpa":0.,"last_return":"elastic"});
    case["stages"] = json!([
        {"kind":"joint_birth","name":"contact","joints":["Seam"]},
        {"kind":"boundary","name":"slide","target_factor":-0.1},
        {"kind":"excavate","name":"remove","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.},
        {"kind":"birth","name":"fill","bodies":["UpperClay"],"states":{"UpperClay":seed},"placement":"cad_new_nodes"},
        {"kind":"joint_birth","name":"rebirth","joints":["Seam"]}
    ]);
    let case: StageCase = serde_json::from_value(case)?;
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let slipped = run.frames[2].state.joints.as_ref().unwrap();
    assert!(
        slipped
            .points
            .iter()
            .any(|p| p.accumulated_plastic_slip_m > 0. && p.plastic_dissipation_kpa_m > 0.)
    );
    assert_eq!(
        run.frames[3].state.joints.as_ref().unwrap().points,
        slipped.points
    );
    let born = run
        .transition_audits
        .last()
        .unwrap()
        .prepared_state
        .joints
        .as_ref()
        .unwrap();
    assert_eq!(born.lifecycle[0].generation, 2);
    for (i, p) in &born.lifecycle[0].archived[0].points {
        assert_eq!(*p, slipped.points[*i]);
    }
    for p in &born.points {
        assert_eq!(p.traction_kpa, [0.; 3]);
        assert_eq!(p.plastic_slip_m, [0.; 3]);
        assert_eq!(p.plastic_dissipation_kpa_m, 0.);
        assert_eq!(p.numerical_dissipation_kpa_m, 0.);
    }
    let accepted = run
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let count = accepted
        + run
            .attempts
            .iter()
            .take_while(|a| a.stage != "rebirth")
            .filter(|a| a.accepted)
            .count()
        + 1;
    let saved = Checkpoint::seal(
        &mesh,
        stages::solve_until(&mesh, ground.clone(), case.clone(), count)?,
    )?;
    let resumed = saved.resume(&mesh, ground, case, None)?;
    let mut a = serde_json::to_value(&run)?;
    let mut b = serde_json::to_value(&resumed)?;
    without_clock(&mut a);
    without_clock(&mut b);
    assert_eq!(a, b);
    Ok(())
}
#[test]
fn omitted_policy_active_soil_birth_and_wet_crack_tip_combinations_stay_rejected() -> Result<()> {
    let (mesh, ground, case) = fixture("stress_free", true)?;
    let mut missing = case.clone();
    missing.joint_system.as_mut().unwrap().joints[0].birth_policy = None;
    assert!(missing.validate().is_err());
    let mut tip_case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-joint-crack-tip.json"))?;
    let tip_mesh: Mesh =
        serde_json::from_str(include_str!("../examples/joint-crack-tip-mesh.json"))?;
    tip_case.joint_system.as_mut().unwrap().joints[0].birth_policy =
        Some(soil_fem_3d::joints::BirthPolicy::StressFree);
    tip_case.joint_system.as_mut().unwrap().joints[0].water =
        Some(soil_fem_3d::joint_water::Policy::Continuous);
    assert!(soil_fem_3d::joints::Model::new(&tip_mesh, tip_case.joint_system.unwrap()).is_err());
    let mut value = serde_json::to_value(case)?;
    value["stages"] = json!([{"kind":"birth","name":"bad","bodies":["UpperClay"],"states":{},"placement":"cad_new_nodes"}]);
    assert!(stages::solve(&mesh, ground, serde_json::from_value(value)?).is_err());
    Ok(())
}
