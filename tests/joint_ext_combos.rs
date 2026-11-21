use nalgebra::Vector3;
use serde_json::json;
use soil_fem_3d::{
    Result,
    ground::Ground,
    joints,
    mesh::Mesh,
    stages::{self, StageCase},
};

fn tip(policy: &str) -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(include_str!("../examples/joint-crack-tip-mesh.json"))?;
    let ground = serde_json::from_str(include_str!("../examples/ground-joint-crack-tip.json"))?;
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-crack-tip.json"))?;
    case["joint_system"]["joints"][0]["birth_policy"] = json!(policy);
    case["stages"] = json!([
        {"kind":"joint_birth","name":"contact","joints":["Seam"]},
        {"kind":"pressure","name":"reload","target_pressure_kpa":20.},
        {"kind":"joint_death","name":"remove","joints":["Seam"]},
        {"kind":"joint_birth","name":"rebirth","joints":["Seam"]}
    ]);
    Ok((mesh, ground, serde_json::from_value(case)?))
}

#[test]
fn tip_force_projection_keeps_shared_nodes_zero_and_rejects_invalid_transfer() -> Result<()> {
    let (mesh, _, case) = tip("transfer_force")?;
    let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
    let active = vec![true; mesh.tetrahedra.len()];
    let u = vec![[0.; 3]; mesh.nodes.len()];
    let rows: Vec<_> = (0..3 * u.len()).map(|g| vec![(g, 1.)]).collect();
    let layout = model.layout();
    let force = |traction: [f64; 3]| {
        let mut f = vec![0.; rows.len()];
        for p in &layout {
            for ((&a, &b), &w) in p.nodes_a.iter().zip(&p.nodes_b).zip(&p.shape) {
                for d in 0..3 {
                    f[3 * a + d] -= p.weight_m2 * w * traction[d];
                    f[3 * b + d] += p.weight_m2 * w * traction[d];
                }
            }
        }
        f
    };
    let initial = model.initial_state()?;
    let ids = vec!["Seam".into()];
    let good = force([12., 5., -100.]);
    let born = model.birth(&initial, &u, &active, &ids, &good, &rows)?;
    let eval = model.evaluate(&born, &u, &active, None, false)?;
    for (actual, expected) in eval.internal_kn.iter().zip(&good) {
        assert!((actual - expected).abs() < 1e-10);
    }
    for n in model.shared_tip_nodes() {
        assert_eq!(&eval.internal_kn[3 * n..3 * n + 3], &[0.; 3]);
        assert!(
            model
                .tied_pairs(&initial, &active)?
                .iter()
                .all(|[a, b]| a != b && *a != n && *b != n)
        );
    }
    // A balanced force pair on the shared front cannot be carried by the
    // contact jump operator. Zero tip coefficients must remain zero ROWS in
    // force recovery, not disappear from its residual check.
    let tip: Vec<_> = model.shared_tip_nodes().into_iter().collect();
    let direction = Vector3::from(mesh.nodes[tip[1]]) - Vector3::from(mesh.nodes[tip[0]]);
    let mut bad_tip = good.clone();
    for axis in 0..3 {
        bad_tip[3 * tip[0] + axis] += direction[axis];
        bad_tip[3 * tip[1] + axis] -= direction[axis];
    }
    let mut resultant = Vector3::zeros();
    let mut moment = Vector3::zeros();
    for (n, &x) in mesh.nodes.iter().enumerate() {
        let delta = Vector3::from_row_slice(&bad_tip[3 * n..3 * n + 3])
            - Vector3::from_row_slice(&good[3 * n..3 * n + 3]);
        resultant += delta;
        moment += Vector3::from(x).cross(&delta);
    }
    assert!(resultant.norm() < 1e-10 && moment.norm() < 1e-10);
    assert!(
        model
            .birth(&initial, &u, &active, &ids, &bad_tip, &rows)
            .is_err(),
        "unreproducible shared-tip residual must not be discarded"
    );
    for q in [[0., 0., 100.], [40., 0., -100.]] {
        assert!(
            model
                .birth(&initial, &u, &active, &ids, &force(q), &rows)
                .is_err()
        );
    }
    let mut incompatible = good;
    let p = &layout[0];
    let n = p
        .nodes_a
        .iter()
        .zip(&p.nodes_b)
        .find(|(a, b)| a != b)
        .unwrap()
        .0;
    incompatible[3 * n] += 1.;
    assert!(
        model
            .birth(&initial, &u, &active, &ids, &incompatible, &rows)
            .is_err()
    );
    let mut partial = active.clone();
    partial[p.owner_b] = false;
    assert!(
        model
            .birth(&initial, &u, &partial, &ids, &force([0., 0., -100.]), &rows)
            .is_err()
    );
    Ok(())
}

#[test]
fn tip_mode_b_stores_the_actual_near_tip_reference_and_freezes_sliding_generation() -> Result<()> {
    let (mesh, _, case) = tip("stress_free")?;
    let model = joints::Model::new(&mesh, case.joint_system.unwrap())?;
    let active = vec![true; mesh.tetrahedra.len()];
    let mut u: Vec<_> = mesh
        .nodes
        .iter()
        .map(|x| [1e-4 * x[0], -2e-4 * x[1], 3e-4 * x[2]])
        .collect();
    // A nonzero, admissible roundoff-sized tie discrepancy must be stored,
    // including near the tip; do not replace the B reference by zero.
    let p = &model.layout()[0];
    let n = *p
        .nodes_b
        .iter()
        .zip(&p.nodes_a)
        .find(|(a, b)| a != b)
        .unwrap()
        .0;
    u[n][2] += 1e-12;
    let rows: Vec<_> = (0..3 * u.len()).map(|g| vec![(g, 1.)]).collect();
    let zero = vec![0.; rows.len()];
    let ids = vec!["Seam".into()];
    let born = model.birth(&model.initial_state()?, &u, &active, &ids, &zero, &rows)?;
    assert!(
        born.points
            .iter()
            .any(|p| Vector3::from(p.jump_m).norm() > 0.)
    );
    for p in &born.points {
        assert_eq!(p.birth.as_ref().unwrap().jump_m, p.jump_m);
        assert_eq!(p.traction_kpa, [0.; 3]);
    }
    for p in model.layout() {
        for (&a, &b) in p.nodes_a.iter().zip(&p.nodes_b) {
            if a != b {
                u[b] = [0.003, 0.001, -0.0001];
            }
        }
    }
    let slipped = model.evaluate(&born, &u, &active, None, false)?.state;
    assert!(
        slipped
            .points
            .iter()
            .any(|p| p.accumulated_plastic_slip_m > 0.)
    );
    let removed = model.death(&slipped, &ids)?;
    let archive = removed.lifecycle[0].archived.clone();
    u.fill([0.02; 3]);
    let reborn = model.birth(&removed, &u, &active, &ids, &zero, &rows)?;
    assert_eq!(reborn.lifecycle[0].generation, 2);
    assert_eq!(reborn.lifecycle[0].archived, archive);
    for p in reborn.points {
        assert_eq!(p.traction_kpa, [0.; 3]);
        assert_eq!(p.plastic_slip_m, [0.; 3]);
        assert_eq!(p.accumulated_plastic_slip_m, 0.);
        assert_eq!(p.plastic_dissipation_kpa_m, 0.);
        assert_eq!(p.numerical_dissipation_kpa_m, 0.);
        assert_eq!(p.endpoint_work_kpa_m, 0.);
    }
    Ok(())
}

#[test]
fn tip_global_modes_have_equilibrium_and_fresh_generations() -> Result<()> {
    for policy in ["transfer_force", "stress_free"] {
        let (mesh, ground, case) = tip(policy)?;
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{policy}: {:?}", run.failure_reason);
        for f in &run.frames {
            assert!(f.relative_force_residual < 1e-8);
            assert!(f.maximum_constraint_error_m < 1e-10);
        }
        for audit in &run.transition_audits {
            if audit.stage == "contact" || audit.stage == "rebirth" {
                for p in &audit.prepared_state.joints.as_ref().unwrap().points {
                    assert_eq!(p.plastic_slip_m, [0.; 3]);
                    assert_eq!(p.plastic_dissipation_kpa_m, 0.);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn hydrostatic_birth_requires_fresh_explicit_state_and_regenerates_covered_interface() -> Result<()>
{
    let mesh: Mesh = serde_json::from_str(include_str!(
        "../examples/joint-water-column-tet10-mesh.json"
    ))?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-joint-water-mc.json"))?;
    ground.initialization.water_level_m = Some(1.5);
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-water-column.json"))?;
    case["joint_system"]["joints"][0]["birth_policy"] = json!("stress_free");
    case["initially_inactive_bodies"] = json!(["UpperClay"]);
    let seed = json!({"effective_stress_compression_positive_kpa":[5.,5.,5.,0.,0.,0.],"plastic_strain_compression_positive":[0.,0.,0.,0.,0.,0.],"void_ratio":0.7,"equivalent_plastic_strain":0.,"plastic_work_kpa":0.,"last_return":"elastic"});
    case["stages"] = json!([
        {"kind":"birth","name":"fill","bodies":["UpperClay"],"states":{"UpperClay":seed},"placement":"cad_new_nodes","hydraulic_initialization":{"kind":"hydrostatic","level_m":1.5}},
        {"kind":"joint_birth","name":"contact","joints":["Seam"]},
        {"kind":"pressure","name":"load","target_pressure_kpa":20.}
    ]);
    case["verification_reject_once_stages"] = json!([]);
    let run = stages::solve(&mesh, ground.clone(), serde_json::from_value(case.clone())?)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let audit = &run.transition_audits[0];
    assert!(
        run.frames[0]
            .state
            .joint_water
            .as_ref()
            .unwrap()
            .nodal_load_kn
            .iter()
            .all(|v| *v == 0.)
    );
    assert!(
        run.frames[1]
            .state
            .joint_water
            .as_ref()
            .unwrap()
            .nodal_load_kn
            .iter()
            .any(|v| *v != 0.)
    );
    let mut point = 0;
    for (e, &region) in mesh.region_ids.iter().enumerate() {
        let count = run
            .integration_points
            .iter()
            .filter(|p| p.element_index == e)
            .count();
        if region == 2 {
            for p in &audit.prepared_state.material_states[point..point + count] {
                assert_eq!(p.stress(), &[5., 5., 5., 0., 0., 0.]);
            }
        }
        point += count;
    }
    for mode in 0..5 {
        let mut bad = case.clone();
        match mode {
            0 => {
                bad["stages"][0]
                    .as_object_mut()
                    .unwrap()
                    .remove("hydraulic_initialization");
            }
            1 => bad["stages"][0]["hydraulic_initialization"]["level_m"] = json!(1.4),
            2 => bad["stages"][0]["states"]["UpperClay"]["equivalent_plastic_strain"] = json!(0.1),
            3 => bad["initially_inactive_bodies"] = json!([]),
            _ => bad["stages"][0]["hydraulic_initialization"]["kind"] = json!("richards"),
        }
        let rejected = match serde_json::from_value::<StageCase>(bad) {
            Ok(c) => stages::solve(&mesh, ground.clone(), c).is_err(),
            Err(_) => true,
        };
        assert!(rejected);
    }
    let parsed: StageCase = serde_json::from_value(case)?;
    let init = run
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let saved = soil_fem_3d::checkpoint::Checkpoint::seal(
        &mesh,
        stages::solve_until(&mesh, ground.clone(), parsed.clone(), init + 1)?,
    )?;
    let mut changed = saved;
    changed
        .run
        .last_committed
        .joint_water
        .as_mut()
        .unwrap()
        .physical_pressure_kpa[0] += 1.;
    changed.run_sha256 = soil_fem_3d::checkpoint::hash_json(&changed.run)?;
    assert!(changed.resume(&mesh, ground, parsed, None).is_err());
    Ok(())
}

#[test]
fn tip_birth_checkpoints_replay_every_position_and_refuse_rehashed_force_reference_and_generation()
-> Result<()> {
    use soil_fem_3d::checkpoint::{Checkpoint, hash_json};
    fn strip(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Object(m) => {
                m.remove("elapsed_seconds");
                for v in m.values_mut() {
                    strip(v);
                }
            }
            serde_json::Value::Array(a) => {
                for v in a {
                    strip(v);
                }
            }
            _ => {}
        }
    }
    for policy in ["transfer_force", "stress_free"] {
        let (mesh, ground, mut case) = tip(policy)?;
        case.verification_reject_once_stages = vec!["contact".into(), "rebirth".into()];
        let full = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(full.complete);
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
        strip(&mut expected);
        let expected = serde_json::to_string(&expected)?;
        for limit in 1..=count {
            let part = stages::solve_until(&mesh, ground.clone(), case.clone(), limit)?;
            let saved = serde_json::to_vec(&Checkpoint::seal(&mesh, part)?)?;
            let resumed = serde_json::from_slice::<Checkpoint>(&saved)?.resume(
                &mesh,
                ground.clone(),
                case.clone(),
                None,
            )?;
            let mut actual = serde_json::to_value(&resumed)?;
            strip(&mut actual);
            assert_eq!(
                serde_json::to_string(&actual)?,
                expected,
                "{policy} at {limit}"
            );
            let cp: Checkpoint = serde_json::from_slice(&saved)?;
            if cp.run.last_committed.joints.as_ref().unwrap().lifecycle[0].generation > 0 {
                for mode in 0..4 {
                    let mut changed: Checkpoint = serde_json::from_slice(&saved)?;
                    let state = changed.run.last_committed.joints.as_mut().unwrap();
                    match mode {
                        0 => state.lifecycle[0].generation += 1,
                        1 => state.points[0].birth.as_mut().unwrap().jump_m[0] += 0.001,
                        2 => state.points[0].birth.as_mut().unwrap().initial_traction_kpa[2] -= 1.,
                        _ => state.points[0].birth.as_mut().unwrap().initial_energy_kpa_m += 1.,
                    }
                    changed.run_sha256 = hash_json(&changed.run)?;
                    assert!(
                        changed
                            .resume(&mesh, ground.clone(), case.clone(), None)
                            .is_err()
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn richards_wet_joint_infiltration_drainage_and_two_way_remain_rejected() -> Result<()> {
    use soil_fem_3d::one_way;
    let mesh: Mesh = serde_json::from_str(include_str!(
        "../examples/joint-water-column-tet10-mesh.json"
    ))?;
    let ground: Ground = serde_json::from_str(include_str!("../examples/ground-one-way-mc.json"))?;
    let mut input: serde_json::Value =
        serde_json::from_str(include_str!("../examples/one-way-mc.json"))?;
    let mut mechanics: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-water-column.json"))?;
    mechanics["stages"] = json!([]);
    mechanics["verification_reject_once_stages"] = json!([]);
    input["mechanics"] = mechanics;
    for head in [2., -1.] {
        input["flow"]["initial_head"]["head_m"] = json!(head);
        let err = one_way::solve(
            &mesh,
            ground.clone(),
            serde_json::from_value(input.clone())?,
        )
        .unwrap_err();
        assert!(
            err.to_string()
                .contains("construction/structures need their coupled policies"),
            "{err}"
        );
    }
    input["coupling"] = json!("two_way");
    assert!(serde_json::from_value::<one_way::Case>(input).is_err());
    Ok(())
}
