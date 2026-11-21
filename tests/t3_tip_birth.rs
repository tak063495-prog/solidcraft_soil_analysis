use nalgebra::Vector3;
use soil_fem_3d::{
    Result,
    joint_law::Material,
    joints::{BirthPolicy, Model},
    mesh::Mesh,
    stages::StageCase,
};

struct Patch {
    model: Model,
    u: Vec<[f64; 3]>,
    active: Vec<bool>,
    rows: Vec<Vec<(usize, f64)>>,
    ids: Vec<String>,
}

impl Patch {
    fn new(policy: BirthPolicy) -> Result<Self> {
        let mesh: Mesh =
            serde_json::from_str(include_str!("../examples/joint-crack-tip-mesh.json"))?;
        let case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-joint-crack-tip.json"))?;
        let mut config = case.joint_system.unwrap();
        config.joints[0].birth_policy = Some(policy);
        Ok(Self {
            model: Model::new(&mesh, config)?,
            u: vec![[0.; 3]; mesh.nodes.len()],
            active: vec![true; mesh.tetrahedra.len()],
            rows: (0..3 * mesh.nodes.len()).map(|g| vec![(g, 1.)]).collect(),
            ids: vec!["Seam".into()],
        })
    }
    fn force(&self, q: [f64; 3]) -> Vec<f64> {
        let mut f = vec![0.; self.rows.len()];
        for p in self.model.layout() {
            for ((a, b), shape) in p.nodes_a.iter().zip(&p.nodes_b).zip(p.shape) {
                for (axis, &traction) in q.iter().enumerate() {
                    let value = p.weight_m2 * shape * traction;
                    f[3 * a + axis] -= value;
                    f[3 * b + axis] += value;
                }
            }
        }
        f
    }
}

#[test]
fn partial_face_history_area_work_and_checkpoint_resume() -> Result<()> {
    let mut patch = Patch::new(BirthPolicy::StressFree)?;
    let j = &patch.model;
    let geometry = j.layout();
    let zero = vec![0.; patch.rows.len()];
    let born = j.birth(
        &j.initial_state()?,
        &patch.u,
        &patch.active,
        &patch.ids,
        &zero,
        &patch.rows,
    )?;
    for p in &geometry {
        for (&a, &b) in p.nodes_a.iter().zip(&p.nodes_b) {
            if a != b {
                patch.u[b] = [0.003, 0., -0.0001];
            }
        }
    }
    let slipped = j
        .evaluate(&born, &patch.u, &patch.active, None, false)?
        .state;
    let mut cut = patch.active.clone();
    cut[geometry[0].owner_b] = false;
    let frozen = j.freeze_removed_faces(&slipped, &patch.active, &cut)?;
    let archive = &frozen.lifecycle[0].archived[0];
    let expected_work: f64 = geometry
        .iter()
        .zip(&slipped.points)
        .filter(|(p, _)| p.face_index == 0)
        .map(|(p, s)| p.weight_m2 * s.endpoint_work_kpa_m)
        .sum();
    assert!(expected_work > 0.);
    assert_eq!(archive.endpoint_work_kn_m, Some(expected_work));
    let region = frozen.lifecycle[0].region.as_ref().unwrap();
    assert!((region.area_m2 - 0.5).abs() < 1e-12);
    assert!(!frozen.lifecycle[0].active && frozen.lifecycle[1].active);
    let eval = j.evaluate(&frozen, &patch.u, &cut, None, false)?;
    let area: f64 = geometry
        .iter()
        .zip(&eval.report.points)
        .filter(|(_, r)| r.active)
        .map(|(p, _)| p.weight_m2)
        .sum();
    assert!((area - 0.5).abs() < 1e-12);
    let mut both = cut.clone();
    both[region.owner_a] = false;
    let still_frozen = j.freeze_removed_faces(&eval.state, &cut, &both)?;
    assert_eq!(still_frozen.lifecycle[0], frozen.lifecycle[0]);
    patch.u.fill([0.; 3]);
    let at_placement = j.evaluate(&frozen, &patch.u, &cut, None, false)?.state;
    let tied = j.prepare_body_birth(&at_placement, &cut, &patch.active)?;
    let restart = serde_json::from_slice(&serde_json::to_vec(&tied)?)?;
    let next = j.birth(
        &restart,
        &patch.u,
        &patch.active,
        &patch.ids,
        &zero,
        &patch.rows,
    )?;
    assert_eq!(
        (next.lifecycle[0].generation, next.lifecycle[1].generation),
        (2, 1)
    );
    assert_eq!(next.lifecycle[0].archived, frozen.lifecycle[0].archived);
    for (i, p) in geometry.iter().enumerate() {
        if p.face_index == 1 {
            assert_eq!(next.points[i], at_placement.points[i]);
        } else {
            assert_eq!(next.points[i].birth.as_ref().unwrap().jump_m, [0.; 3]);
            assert_eq!(next.points[i].endpoint_work_kpa_m, 0.);
            assert_eq!(next.points[i].plastic_dissipation_kpa_m, 0.);
        }
    }
    assert!(
        j.birth(
            &next,
            &patch.u,
            &patch.active,
            &patch.ids,
            &zero,
            &patch.rows
        )
        .is_err()
    );
    let mut tampered = frozen.clone();
    tampered.lifecycle[0].region.as_mut().unwrap().area_m2 += 1.;
    assert!(j.validate_lifecycle(&tampered).is_err());
    Ok(())
}

#[test]
fn strict_point_cone_and_one_ulp_excess() {
    let mat = Material {
        normal_stiffness_kpa_m: 1e5,
        tangential_stiffness_kpa_m: 1e5,
        friction_coefficient: 0.3,
    };
    for (t, ok) in [
        (0., true),
        (30., true),
        (30_f64.next_up(), false),
        (40., false),
    ] {
        for (normal, q) in [
            ([0., 0., 1.], [t, 0., -100.]),
            ([1., 0., 0.], [-100., t, 0.]),
        ] {
            assert_eq!(mat.birth(normal, [0.; 3], q).is_ok(), ok);
        }
    }
    for q in [[0., 0., 1.], [1., 0., 0.]] {
        assert!(mat.birth([0., 0., 1.], [0.; 3], q).is_err());
    }
}

#[test]
fn rejected_transfer_records_trial_without_committing_work() -> Result<()> {
    let patch = Patch::new(BirthPolicy::TransferForce)?;
    let state = patch.model.initial_state()?;
    let checkpoint = serde_json::to_vec(&state)?;
    for q in [[40., 0., -100.], [0., 0., 100.]] {
        let f = patch.force(q);
        let err = patch
            .model
            .birth(&state, &patch.u, &patch.active, &patch.ids, &f, &patch.rows)
            .unwrap_err()
            .to_string();
        for token in [
            "GP=",
            "p_kpa=",
            "friction_excess_kpa=",
            "requested_kn=",
            "reproduced_kn=",
        ] {
            assert!(err.contains(token));
        }
        assert_eq!(serde_json::to_vec(&state)?, checkpoint);
    }
    let next = patch.model.birth(
        &state,
        &patch.u,
        &patch.active,
        &patch.ids,
        &patch.force([0., 0., -100.]),
        &patch.rows,
    )?;
    assert_eq!(next.lifecycle[0].generation, 1);
    assert!(
        next.points
            .iter()
            .all(|p| p.endpoint_work_kpa_m == 0. && p.plastic_dissipation_kpa_m == 0.)
    );
    Ok(())
}

#[test]
fn free_and_prescribed_tip_loads_use_the_same_zero_contact_row() -> Result<()> {
    let mut patch = Patch::new(BirthPolicy::TransferForce)?;
    let initial = patch.model.initial_state()?;
    let tip = *patch.model.shared_tip_nodes().first().unwrap();
    let mut f = patch.force([0., 0., -100.]);
    f[3 * tip] += 1.;
    let error = patch
        .model
        .birth(
            &initial,
            &patch.u,
            &patch.active,
            &patch.ids,
            &f,
            &patch.rows,
        )
        .unwrap_err()
        .to_string();
    assert!(error.contains("every free nodal residual") && error.contains("rank="));
    patch.rows[3 * tip].clear();
    let born = patch.model.birth(
        &initial,
        &patch.u,
        &patch.active,
        &patch.ids,
        &f,
        &patch.rows,
    )?;
    let value = patch
        .model
        .evaluate(&born, &patch.u, &patch.active, None, false)?;
    assert_eq!(&value.internal_kn[3 * tip..3 * tip + 3], &[0.; 3]);
    assert_eq!(value.internal_kn[3 * tip] - f[3 * tip], -1.);
    Ok(())
}

#[test]
fn projected_field_minimizes_area_norm_at_fixed_all_node_load() -> Result<()> {
    let patch = Patch::new(BirthPolicy::TransferForce)?;
    let target = [12., 5., -100.];
    let f = patch.force(target);
    let born = patch.model.birth(
        &patch.model.initial_state()?,
        &patch.u,
        &patch.active,
        &patch.ids,
        &f,
        &patch.rows,
    )?;
    let eval = patch
        .model
        .evaluate(&born, &patch.u, &patch.active, None, false)?;
    assert!(
        eval.internal_kn
            .iter()
            .zip(f)
            .all(|(a, b)| (a - b).abs() <= 1e-10)
    );
    // The forward constant field and the fitted field have identical nodal forces.
    // Their difference is thus a known nonzero null direction.
    let (mut cross, mut distance) = (0., 0.);
    for (p, state) in patch.model.layout().iter().zip(born.points) {
        let q = Vector3::from(state.birth.unwrap().initial_traction_kpa);
        let delta = Vector3::from(target) - q;
        cross += p.weight_m2 * q.dot(&delta);
        distance += p.weight_m2 * delta.norm_squared();
    }
    assert!(distance > 1.);
    assert!(cross.abs() < 1e-10);
    Ok(())
}

#[test]
fn active_tip_excavation_freezes_only_removed_face_and_balances_holding() -> Result<()> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/joint-crack-tip-mesh.json"))?;
    let mut mesh = serde_json::to_value(mesh)?;
    mesh["excavation_sets"] = serde_json::json!({"PartialFace":[13,14,15]});
    let mut ground: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-joint-crack-tip.json"))?;
    ground["excavation_bodies"] = serde_json::json!(["PartialFace"]);
    ground["initialization"]["geostatic"]["initial_boundary_factor"] = serde_json::json!(-0.001);
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-joint-crack-tip.json"))?;
    case["joint_system"]["joints"][0]["birth_policy"] = serde_json::json!("stress_free");
    case["verification_reject_once_stages"] = serde_json::json!([]);
    case["stages"] = serde_json::json!([
        {"kind":"joint_birth","name":"birth","joints":["Seam"]},
        {"kind":"excavate","name":"cut","bodies":["PartialFace"],"remaining_surface_pressure_kpa":0.}
    ]);
    let run = soil_fem_3d::stages::solve(
        &serde_json::from_value(mesh)?,
        serde_json::from_value(ground)?,
        serde_json::from_value(case)?,
    )?;
    assert!(run.complete, "{:?}", run.failure_reason);
    let old = run.frames[1].state.joints.as_ref().unwrap();
    let next = run.last_committed.joints.as_ref().unwrap();
    assert!(run.excavation_audits[0].zero_release_relative_force_residual <= 1e-10);
    for face in &next.lifecycle {
        let info = face.region.as_ref().unwrap();
        assert_eq!(face.active, info.owner_b != 13);
        assert_eq!(face.generation, 1);
        if face.active {
            assert!(face.archived.is_empty());
        } else {
            for (i, s) in &face.archived[0].points {
                assert_eq!(s, &old.points[*i]);
                assert_eq!(s, &next.points[*i]);
            }
        }
    }
    Ok(())
}
