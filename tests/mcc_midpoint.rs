use nalgebra::{Matrix3, Rotation3, SVector};
use soil_fem_3d::{
    Result,
    mcc::{self, Integration, Parameters, State},
};
fn parameters() -> Parameters {
    Parameters {
        critical_state_m: 1.2,
        lambda: 0.2,
        kappa: 0.04,
        poisson: 0.3,
        reference_pressure_kpa: 1.,
        specific_volume_ncl_at_reference: 3.,
        integration: Integration::Midpoint,
    }
}
fn state(stress: [f64; 6], ratio: f64) -> State {
    let (p, q) = mcc::invariants(&stress);
    State {
        effective_stress_compression_positive_kpa: stress,
        void_ratio: 1.,
        preconsolidation_pressure_kpa: ratio * (p + q * q / (1.44 * p)),
        plastic_strain_compression_positive: [0.; 6],
    }
}
#[test]
fn optional_scheme_preserves_legacy_json_and_rejects_unknown_integration() -> Result<()> {
    let mut p = parameters();
    p.integration = Integration::BackwardEuler;
    let value = serde_json::to_value(&p)?;
    assert!(value.get("integration").is_none());
    let restored: Parameters = serde_json::from_value(value.clone())?;
    assert_eq!(restored.integration, Integration::BackwardEuler);
    let mut explicit = value.clone();
    explicit["integration"] = serde_json::json!("backward_euler");
    let explicit: Parameters = serde_json::from_value(explicit)?;
    let old = state([100., 100., 100., 0., 0., 0.], 1.1);
    let de = [0.001, 0.002, 0.003, 0.001, -0.0003, 0.0002];
    assert_eq!(
        serde_json::to_value(restored.integrate_partition(&old, de, 8)?.state)?,
        serde_json::to_value(explicit.integrate_partition(&old, de, 8)?.state)?
    );
    let mut value = value;
    value["integration"] = serde_json::json!("midpoint");
    assert_eq!(
        serde_json::from_value::<Parameters>(value.clone())?.integration,
        Integration::Midpoint
    );
    value["integration"] = serde_json::json!("other");
    assert!(serde_json::from_value::<Parameters>(value).is_err());
    Ok(())
}
#[test]
fn ocr_first_yield_and_unloading_converge_quadratically_to_closed_form() -> Result<()> {
    let p = parameters();
    let old = state([100., 100., 100., 0., 0., 0.], 1.5);
    let ev: f64 = 0.02;
    let crossing = -(-p.kappa / 2. * 1.5_f64.ln()).ln_1p();
    let vy = 2. * (-crossing).exp();
    let v = 2. * (-ev).exp();
    let pressure = 150. * (vy / p.lambda * (-(-(ev - crossing)).exp_m1())).exp();
    let mut errors = Vec::new();
    for count in [8, 16, 32, 64] {
        let loaded = p
            .integrate_partition(&old, [ev / 3., ev / 3., ev / 3., 0., 0., 0.], count)?
            .state;
        let error = (loaded.effective_stress_compression_positive_kpa[0] - pressure).abs();
        errors.push(error);
        assert!((1. + loaded.void_ratio - v).abs() < 1e-13);
        assert!(
            (loaded.preconsolidation_pressure_kpa
                - loaded.effective_stress_compression_positive_kpa[0])
                .abs()
                < 1e-9
        );
        let unloaded = p
            .integrate_partition(
                &loaded,
                [-0.005 / 3., -0.005 / 3., -0.005 / 3., 0., 0., 0.],
                count,
            )?
            .state;
        let target = pressure * (v / p.kappa * (-0.005_f64.exp_m1())).exp();
        assert!((unloaded.effective_stress_compression_positive_kpa[0] - target).abs() < error);
        assert_eq!(
            unloaded.preconsolidation_pressure_kpa,
            loaded.preconsolidation_pressure_kpa
        );
        assert_eq!(
            unloaded.plastic_strain_compression_positive,
            loaded.plastic_strain_compression_positive
        );
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.5, "{errors:?}");
    }
    assert!(errors[3] < 4e-6);
    Ok(())
}
#[test]
fn full_partition_tangent_includes_first_yield_hardening_and_opposite_reyield() -> Result<()> {
    let p = parameters();
    for (stress, de, ratio) in [
        (
            [100., 100., 100., 0., 0., 0.],
            [0.004, 0.004, 0.004, 0., 0., 0.],
            1.,
        ),
        (
            [130., 85., 85., 8., -5., 3.],
            [0.002, 0.001, 0.003, 0.004, -0.002, 0.001],
            1.,
        ),
        (
            [80., 20., 20., 0., 0., 0.],
            [0.03, -0.015, -0.015, 0.004, -0.002, 0.003],
            1.,
        ),
        (
            [100., 100., 100., 0., 0., 0.],
            [0.02 / 3., 0.02 / 3., 0.02 / 3., 0., 0., 0.],
            1.5,
        ),
        (
            [
                150.59644256269408,
                74.70177871865296,
                74.70177871865296,
                0.,
                0.,
                0.,
            ],
            [-0.04, 0.02, 0.02, 0., 0., 0.],
            1.,
        ),
        (
            [110., 95., 95., 2., -1., 0.5],
            [0.001, -0.0002, 0.0004, 0.002, -0.001, 0.0007],
            5.,
        ),
        (
            [110., 95., 95., 2., -1., 0.5],
            [0.001, -0.001, 0., 0.002, -0.001, 0.0007],
            5.,
        ),
    ] {
        let old = state(stress, ratio);
        for count in [1, 8, 16, 32] {
            let (a, d) = p.implicit_stress_tangent_partition(&old, de, count)?;
            let direct = p.integrate_partition(&old, de, count)?;
            assert_eq!(
                serde_json::to_value(&a.state)?,
                serde_json::to_value(&direct.state)?
            );
            let mut numerator = 0.;
            let mut denominator = 0.;
            for j in 0..6 {
                let h = 2e-7;
                let mut plus = de;
                let mut minus = de;
                plus[j] += h;
                minus[j] -= h;
                let x = p.integrate_partition(&old, plus, count)?.state;
                let y = p.integrate_partition(&old, minus, count)?.state;
                for i in 0..6 {
                    let fd = (x.effective_stress_compression_positive_kpa[i]
                        - y.effective_stress_compression_positive_kpa[i])
                        / (2. * h);
                    assert!((fd - d[(i, j)]).abs() / (1. + fd.abs()) < 1e-6);
                    numerator += (fd - d[(i, j)]).powi(2);
                    denominator += fd * fd;
                }
            }
            assert!((numerator / denominator).sqrt() < 2e-8);
        }
    }
    Ok(())
}
fn matrix(s: [f64; 6], strain: bool) -> Matrix3<f64> {
    let k = if strain { 0.5 } else { 1. };
    Matrix3::from_row_slice(&[
        s[0],
        k * s[3],
        k * s[5],
        k * s[3],
        s[1],
        k * s[4],
        k * s[5],
        k * s[4],
        s[2],
    ])
}
fn vector(s: Matrix3<f64>, strain: bool) -> [f64; 6] {
    let k = if strain { 2. } else { 1. };
    [
        s[(0, 0)],
        s[(1, 1)],
        s[(2, 2)],
        k * s[(0, 1)],
        k * s[(1, 2)],
        k * s[(0, 2)],
    ]
}
#[test]
fn arbitrary_rotation_preserves_stress_plastic_history_and_directional_tangent() -> Result<()> {
    let p = parameters();
    let q = Rotation3::from_euler_angles(0.29, -0.31, 0.43).into_inner();
    for (stress, de, ratio) in [
        (
            [130., 85., 85., 8., -5., 3.],
            [0.002, 0.001, 0.003, 0.004, -0.002, 0.001],
            1.,
        ),
        (
            [80., 20., 20., 0., 0., 0.],
            [0.03, -0.015, -0.015, 0.004, -0.002, 0.003],
            1.,
        ),
        (
            [100., 100., 100., 0., 0., 0.],
            [0.02 / 3., 0.02 / 3., 0.02 / 3., 0., 0., 0.],
            1.5,
        ),
    ] {
        let mut old = state(stress, ratio);
        old.plastic_strain_compression_positive = [0.001, 0.002, -0.001, 0.0007, -0.0003, 0.0004];
        let mut rotated = old.clone();
        rotated.effective_stress_compression_positive_kpa =
            vector(q * matrix(stress, false) * q.transpose(), false);
        rotated.plastic_strain_compression_positive = vector(
            q * matrix(old.plastic_strain_compression_positive, true) * q.transpose(),
            true,
        );
        let (a, d) = p.implicit_stress_tangent_partition(&old, de, 8)?;
        let (b, c) = p.implicit_stress_tangent_partition(
            &rotated,
            vector(q * matrix(de, true) * q.transpose(), true),
            8,
        )?;
        let back = vector(
            q.transpose() * matrix(b.state.effective_stress_compression_positive_kpa, false) * q,
            false,
        );
        let history = vector(
            q.transpose() * matrix(b.state.plastic_strain_compression_positive, true) * q,
            true,
        );
        for i in 0..6 {
            assert!((back[i] - a.state.effective_stress_compression_positive_kpa[i]).abs() < 1e-8);
            assert!((history[i] - a.state.plastic_strain_compression_positive[i]).abs() < 1e-11);
        }
        let direction = [0.3, -0.2, 0.1, 0.4, -0.1, 0.2];
        let local = d * SVector::<f64, 6>::from(direction);
        let result =
            c * SVector::<f64, 6>::from(vector(q * matrix(direction, true) * q.transpose(), true));
        let back = SVector::<f64, 6>::from(vector(
            q.transpose() * matrix(result.into(), false) * q,
            false,
        ));
        assert!((local - back).norm() / local.norm() < 1e-9);
    }
    Ok(())
}
#[test]
fn elastic_reversal_and_permitted_positive_drift_unloading_retain_pc_without_projection()
-> Result<()> {
    let p = parameters();
    let old = state([110., 95., 95., 2., -1., 0.5], 5.);
    let de = [0.001, -0.0002, 0.0004, 0.002, -0.001, 0.0007];
    let loaded = p.integrate_partition(&old, de, 8)?.state;
    let back = p.integrate_partition(&loaded, de.map(|x| -x), 8)?.state;
    for i in 0..6 {
        assert!(
            (back.effective_stress_compression_positive_kpa[i]
                - old.effective_stress_compression_positive_kpa[i])
                .abs()
                < 1e-9
        );
    }
    assert_eq!(
        back.preconsolidation_pressure_kpa,
        old.preconsolidation_pressure_kpa
    );
    let mut old = state([100., 100., 100., 0., 0., 0.], 1.);
    old.effective_stress_compression_positive_kpa[..3].fill(100. * (1. + 1e-9));
    let initial = p.normalized_yield(&old);
    let unloaded =
        p.integrate_partition(&old, [-1e-12 / 3., -1e-12 / 3., -1e-12 / 3., 0., 0., 0.], 1)?;
    assert_eq!(unloaded.plastic_substeps, 0);
    assert_eq!(
        unloaded.state.preconsolidation_pressure_kpa,
        old.preconsolidation_pressure_kpa
    );
    assert_eq!(
        unloaded.state.plastic_strain_compression_positive,
        old.plastic_strain_compression_positive
    );
    assert!(
        1e-11 < p.normalized_yield(&unloaded.state)
            && p.normalized_yield(&unloaded.state) < initial
    );
    Ok(())
}
#[test]
fn invalid_midpoint_trials_leave_committed_material_history_unchanged() -> Result<()> {
    let p = parameters();
    let old = state([80., 20., 20., 0., 0., 0.], 1.);
    let before = serde_json::to_value(&old)?;
    assert!(p.integrate_partition(&old, [1.; 6], 1).is_err());
    assert!(p.integrate_partition(&old, [0.; 6], 0).is_err());
    assert!(p.integrate_partition(&old, [f64::NAN; 6], 1).is_err());
    assert_eq!(before, serde_json::to_value(&old)?);
    Ok(())
}

#[test]
fn adaptive_partition_and_its_consistent_tangent_use_the_identical_material_map() -> Result<()> {
    let p = parameters();
    let old = state([130., 85., 85., 8., -5., 3.], 1.);
    let strain = [0.002, 0.001, 0.003, 0.004, -0.002, 0.001];
    let adaptive = p.integrate(&old, strain, 1e-3)?;
    let fixed = p.integrate_partition(&old, strain, adaptive.substeps)?;
    let (implicit, _) = p.implicit_stress_tangent_partition(&old, strain, adaptive.substeps)?;
    assert_eq!(
        serde_json::to_value(&adaptive.state)?,
        serde_json::to_value(&fixed.state)?
    );
    assert_eq!(
        serde_json::to_value(&adaptive.state)?,
        serde_json::to_value(&implicit.state)?
    );
    Ok(())
}

fn hydrostatic_tet10(
    schema: &str,
) -> Result<(
    soil_fem_3d::mesh::Mesh,
    soil_fem_3d::ground::Ground,
    soil_fem_3d::stages::StageCase,
)> {
    use serde_json::json;
    let mesh: soil_fem_3d::mesh::Mesh = serde_json::from_value(json!({
      "schema":"soil-fem-mesh/2","units":"m-kN-kPa","axis":"z_up","crs":"LOCAL","origin_m":[0.,0.,0.],
      "nodes":[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],"tetrahedra":[[0,1,2,3]],"material_ids":[1],
      "node_sets":{"n0":[0],"n1":[1],"n2":[2],"all":[0,1,2,3]},
      "face_sets":{"all":[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]},
      "regions":[{"id":1,"body":"soil","material_id":1,"bounds_m":[[0.,0.,0.],[1.,1.,1.]],"volume_m3":1./6.,"source_sha256":"0".repeat(64)}],"region_ids":[1]
    }))?;
    let mut mesh = mesh.elevate_straight()?;
    // Manufactured hydrostatic problem: exact P2 counterloads cancel gravity.
    // Shape integrals are -V/20 at vertices and V/5 at the six edge nodes.
    let countergravity: Vec<_> = (0..10)
        .map(|i| {
            mesh.node_sets.insert(format!("n{i}"), vec![i]);
            let weight = if i < 4 { -1. / 120. } else { 1. / 30. };
            json!({"node_set":format!("n{i}"),"total_force_kn":[0.,0.,20.*weight]})
        })
        .collect();
    let mut p = parameters();
    p.integration = Integration::BackwardEuler;
    p.specific_volume_ncl_at_reference = 2. + p.lambda * 100_f64.ln();
    let mut ground: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    ground["materials"] = json!([{"id":1,"name":"clay","unit_weight_kn_m3":20.,"mcc":p}]);
    ground["bodies"] = json!([{"body":"soil","material_id":1}]);
    ground["initialization"]["geostatic"]["seed"] =
        json!({"kind":"body_states","states":{"soil":state([100.,100.,100.,0.,0.,0.],1.)}});
    ground["initialization"]["geostatic"]["base_loads"] = json!({"nodal_loads":countergravity,"surface_loads":[{"kind":"pressure","face_set":"all","pressure_kpa":100.}]});
    let ground = serde_json::from_value(ground)?;
    let mut case: serde_json::Value =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case["schema"] = json!(schema);
    case["settings"]["mcc_integration"] = json!("midpoint");
    case["settings"]["initial_increment"] = json!(0.25);
    case["settings"]["minimum_increment"] = json!(1. / 256.);
    case["settings"]["material_maximum_substep"] = json!(5e-5);
    case["settings"]["relative_force_tolerance"] = json!(1e-10);
    case["settings"]["relative_linear_tolerance"] = json!(1e-12);
    case["settings"]["maximum_free_dofs"] = json!(2000);
    case["settings"]["linear_solver"] = json!({"method":"dense_lu"});
    case["boundary"] = json!({"constraints":[
      {"node_set":"n0","direction":[1.,0.,0.],"value_m":0.},{"node_set":"n0","direction":[0.,1.,0.],"value_m":0.},{"node_set":"n0","direction":[0.,0.,1.],"value_m":0.},
      {"node_set":"n1","direction":[0.,1.,0.],"value_m":0.},{"node_set":"n1","direction":[0.,0.,1.],"value_m":0.},{"node_set":"n2","direction":[0.,0.,1.],"value_m":0.}],
      "surface_loads":[{"kind":"pressure","face_set":"all","pressure_kpa":20.}]});
    case["monitors"] = json!({});
    case["monitor_sets"] = json!({});
    case["verification_reject_once_stages"] = json!(["load"]);
    case["stages"] = json!([{"kind":"boundary","name":"load","target_factor":1.},{"kind":"boundary","name":"unload","target_factor":0.}]);
    Ok((mesh, ground, serde_json::from_value(case)?))
}
#[test]
fn both_global_mixed_forms_follow_hydrostatic_nc_and_unloading_closed_forms() -> Result<()> {
    use soil_fem_3d::stages;
    for schema in [
        "soil-fem-mixed-mcc-stages/1",
        "soil-fem-projected-mcc-stages/1",
    ] {
        let (mesh, ground, case) = hydrostatic_tet10(schema)?;
        assert_eq!(
            ground.materials[0].model.as_mcc().unwrap().integration,
            Integration::BackwardEuler
        );
        let run = stages::solve(&mesh, ground, case)?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(run.frames.len(), 3);
        assert!(
            run.attempts
                .iter()
                .any(|a| a.stage == "load" && !a.accepted)
        );
        for (frame, pressure, volume) in [
            (&run.frames[1], 120., 2. - 0.2 * 1.2_f64.ln()),
            (&run.frames[2], 100., 2. - (0.2 - 0.04) * 1.2_f64.ln()),
        ] {
            let strain = (2. / volume).ln() / 3.;
            for (point, u) in mesh.nodes.iter().zip(&frame.state.displacement_m) {
                for a in 0..3 {
                    assert!(
                        (u[a] + point[a] * strain).abs() < 2e-9,
                        "{schema}: {u:?} vs {point:?}"
                    );
                }
            }
            for material in &frame.state.material_states {
                let s = material.as_mcc().unwrap();
                for stress in &s.effective_stress_compression_positive_kpa[..3] {
                    assert!((*stress - pressure).abs() < 1e-7);
                }
                assert!((s.preconsolidation_pressure_kpa - 120.).abs() < 1e-6);
                assert!((1. + s.void_ratio - volume).abs() < 1e-8);
            }
        }
    }
    Ok(())
}
#[test]
fn global_midpoint_cutback_checkpoint_restores_plastic_history_and_binds_scheme() -> Result<()> {
    use soil_fem_3d::{checkpoint::Checkpoint, stages};
    let (mesh, ground, case) = hydrostatic_tet10("soil-fem-projected-mcc-stages/1")?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete);
    let initialization = full
        .geostatic
        .as_ref()
        .unwrap()
        .attempts
        .iter()
        .filter(|a| a.accepted)
        .count();
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), initialization + 1)?;
    assert!(!partial.complete);
    assert!(
        partial
            .last_committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|x| x.abs() > 1e-10))
    );
    let saved = Checkpoint::seal(&mesh, partial)?;
    let encoded = serde_json::to_vec(&saved)?;
    let saved: Checkpoint = serde_json::from_slice(&encoded)?;
    let resumed = saved.resume(&mesh, ground.clone(), case.clone(), None)?;
    assert_eq!(
        serde_json::to_value(&resumed.last_committed)?,
        serde_json::to_value(&full.last_committed)?
    );
    assert_eq!(
        serde_json::to_value(&resumed.frames)?,
        serde_json::to_value(&full.frames)?
    );
    let mut changed = case;
    changed.settings.mcc_integration = Some(Integration::BackwardEuler);
    let saved: Checkpoint = serde_json::from_slice(&encoded)?;
    assert!(saved.resume(&mesh, ground, changed, None).is_err());
    Ok(())
}

#[test]
fn both_global_schemes_start_nearly_hydrostatic_elastic_newton_steps() -> Result<()> {
    use soil_fem_3d::{geostatic::Seed, stages};
    for (label, scheme) in [
        ("be", Integration::BackwardEuler),
        ("mid", Integration::Midpoint),
    ] {
        for projected in [false, true] {
            let schema = if projected {
                "soil-fem-projected-mcc-stages/1"
            } else {
                "soil-fem-mixed-mcc-stages/1"
            };
            let (mesh, mut ground, mut case) = hydrostatic_tet10(schema)?;
            case.settings.mcc_integration = Some(scheme);
            let Seed::BodyStates { states } =
                &mut ground.initialization.geostatic.as_mut().unwrap().seed
            else {
                panic!("hydrostatic fixture requires body states")
            };
            let soil_fem_3d::soil::State::Mcc(state) = states.get_mut("soil").unwrap() else {
                panic!("hydrostatic fixture requires MCC")
            };
            state.preconsolidation_pressure_kpa *= 2.;
            let run = stages::solve(&mesh, ground, case)?;
            assert!(
                run.complete,
                "{label} projected={projected}: {:?}",
                run.failure_reason
            );
            assert_eq!(run.frames.len(), 3);
            assert!(
                run.last_committed
                    .material_states
                    .iter()
                    .all(|s| s.plastic_strain().iter().all(|x| x.abs() < 1e-12))
            );
        }
    }
    Ok(())
}
