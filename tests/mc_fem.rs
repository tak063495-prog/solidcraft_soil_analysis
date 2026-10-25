use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    ground::{self, Ground},
    mesh::Mesh,
    soil,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    Ok((
        if quadratic {
            mesh.elevate_straight()?
        } else {
            mesh
        },
        serde_json::from_str(include_str!("../examples/ground-mc-column.json"))?,
        serde_json::from_str(include_str!("../examples/stages-mc-column.json"))?,
    ))
}
// Independent scalar uniaxial-strain perfect-plastic branch, assembled in series.
fn column_stress(total_strain: f64) -> (f64, [f64; 2]) {
    let phi = 30f64.to_radians().sin();
    let psi = 5f64.to_radians().sin();
    let c = 10.;
    let nu = 0.1;
    let a = (1. - nu) / ((1. + nu) * (1. - 2. * nu));
    let b = nu / ((1. + nu) * (1. - 2. * nu));
    let m = [-(1. + psi) / 2., -(1. + psi) / 2., 1. - psi];
    let n = [-(1. + phi) / 2., -(1. + phi) / 2., 1. - phi];
    let dm = std::array::from_fn::<_, 3, _>(|i| a * m[i] + b * (m.iter().sum::<f64>() - m[i]));
    let denominator = n.iter().zip(dm).map(|(x, y)| x * y).sum::<f64>();
    let drive = (1. - phi) * a - (1. + phi) * b;
    let initial = (1. - phi) * 100. - (1. + phi) * 50. - 2. * c * 30f64.to_radians().cos();
    let strain = |stress: f64, e: f64| {
        let ey = -initial / (drive * e);
        let sy = 100. + a * e * ey;
        if stress <= sy {
            (stress - 100.) / (a * e)
        } else {
            ey + (stress - sy) / (e * (a - dm[2] * drive / denominator))
        }
    };
    let mut lo = 100.;
    let mut hi = 1000.;
    for _ in 0..100 {
        let mid = (lo + hi) / 2.;
        if (strain(mid, 10000.) + strain(mid, 20000.)) / 2. < total_strain {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let stress = (lo + hi) / 2.;
    (stress, [strain(stress, 10000.), strain(stress, 20000.)])
}

#[test]
fn mc_tet4_and_tet10_global_newton_match_series_column_theory() -> Result<()> {
    let (expected, eps) = column_stress(0.02);
    for quadratic in [false, true] {
        let (mesh, ground, case) = fixture(quadratic)?;
        let run = stages::solve(&mesh, ground.clone(), case.clone())?;
        assert!(run.complete, "{:?}", run.failure_reason);
        assert_eq!(run.schema, "soil-fem-stage-result/3");
        for (i, p) in run.integration_points.iter().enumerate() {
            let state = &run.frames[1].state.material_states[i];
            assert!(
                (state.stress()[2] - expected).abs() < 2e-5,
                "computed {:?},expected {expected}",
                state.stress()
            );
            let soil::State::MohrCoulomb(s) = state else {
                panic!("MC model lost its history type")
            };
            assert!(s.equivalent_plastic_strain > 0.);
            assert!(s.plastic_work_kpa >= 0.);
            assert!(p.weight_m3 > 0.);
        }
        for (n, x) in mesh.nodes.iter().enumerate() {
            if (x[2] - 0.5).abs() < 1e-10 {
                assert!((run.frames[1].state.displacement_m[n][2] + eps[0] * 0.5).abs() < 2e-8);
            }
        }
        assert!(
            run.attempts
                .iter()
                .filter(|a| a.accepted)
                .any(|a| a.newton_iterations.unwrap_or(0) > 1)
        );
        assert!(run.attempts.iter().any(|a| !a.accepted
            && a.committed_state_unchanged_on_rejection
            && a.rejected_trial_plastic_changed
            && a.rejected_trial_material_changed));
        assert!(
            run.frames
                .iter()
                .all(|f| f.relative_force_residual < 1e-8 && f.maximum_constraint_error_m < 1e-10)
        );
        if quadratic {
            let mut dense = case;
            dense.settings.linear_solver = soil_fem_3d::linalg::LinearSolver::DenseLu;
            dense.settings.maximum_free_dofs = 2000;
            let dense = stages::solve(&mesh, ground, dense)?;
            assert!(dense.complete, "{:?}", dense.failure_reason);
            for (a, b) in run
                .last_committed
                .displacement_m
                .iter()
                .flatten()
                .zip(dense.last_committed.displacement_m.iter().flatten())
            {
                assert!((a - b).abs() < 1e-9);
            }
            for (a, b) in run
                .last_committed
                .material_states
                .iter()
                .zip(&dense.last_committed.material_states)
            {
                for (a, b) in a.stress().iter().zip(b.stress()) {
                    assert!((a - b).abs() < 2e-5);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn mcc_and_mc_share_one_mesh_without_sharing_history_fields() -> Result<()> {
    let (mesh, mut ground, mut case) = fixture(true)?;
    let other: Ground = serde_json::from_str(include_str!("../examples/ground-local.json"))?;
    ground.materials[0].model = other.materials[0].model.clone();
    ground.layers[0].consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 2. };
    ground.layers[0].void_ratio = ground::VoidRatio::FromMccReference;
    case.boundary
        .as_mut()
        .unwrap()
        .constraints
        .iter_mut()
        .find(|c| c.node_set == "top")
        .unwrap()
        .value_m = -0.002;
    case.verification_reject_once_stages.clear();
    case.settings.initial_increment = 1.;
    let run = stages::solve(&mesh, ground, case)?;
    assert!(run.complete, "{:?}", run.failure_reason);
    for (i, p) in run.integration_points.iter().enumerate() {
        assert_eq!(
            matches!(run.last_committed.material_states[i], soil::State::Mcc(_)),
            mesh.material_ids[p.element_index] == 1
        );
    }
    let json = serde_json::to_value(&run)?;
    let states = json["last_committed"]["material_states"]
        .as_array()
        .unwrap();
    assert!(
        states
            .iter()
            .any(|s| s.get("preconsolidation_pressure_kpa").is_some()
                && s.get("equivalent_plastic_strain").is_none())
    );
    assert!(
        states
            .iter()
            .any(|s| s.get("preconsolidation_pressure_kpa").is_none()
                && s.get("equivalent_plastic_strain").is_some())
    );
    Ok(())
}
#[test]
fn mc_checkpoint_restores_nonzero_plastic_history_and_rejects_state_model_mismatch() -> Result<()> {
    let (mesh, ground, case) = fixture(true)?;
    let full = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(full.complete);
    let partial = stages::solve_until(&mesh, ground.clone(), case.clone(), 2)?;
    assert!(!partial.complete);
    assert!(
        partial
            .last_committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|v| v.abs() > 1e-7))
    );
    let sealed = Checkpoint::seal(&mesh, partial)?;
    let saved: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&sealed)?)?;
    let resumed = saved.resume(&mesh, ground.clone(), case.clone(), None)?;
    assert_eq!(
        serde_json::to_vec(&full.last_committed)?,
        serde_json::to_vec(&resumed.last_committed)?
    );
    assert_eq!(
        serde_json::to_vec(&full.frames)?,
        serde_json::to_vec(&resumed.frames)?
    );
    let mcc = soil::Model::Mcc(serde_json::from_value(
        serde_json::json!({"critical_state_m":1.2,"lambda":0.1,"kappa":0.02,"poisson":0.3,"reference_pressure_kpa":1.,"specific_volume_ncl_at_reference":2.7}),
    )?);
    assert!(
        mcc.validate_state(&full.last_committed.material_states[0])
            .is_err()
    );
    Ok(())
}
#[test]
fn mc_initial_k0_has_no_fictitious_ocr_and_mcc_json_remains_compatible() -> Result<()> {
    let (mesh, mut ground, _) = fixture(false)?;
    let result = ground::initialize(&mesh, ground.clone())?;
    assert!(
        result
            .point_states
            .iter()
            .all(|p| p.mcc_pc_over_mean_effective.is_none())
    );
    let raw = serde_json::to_value(result)?;
    assert!(raw["point_states"][0].get("mohr_coulomb_state").is_some());
    assert!(raw["point_states"][0].get("mcc_state").is_none());
    ground.layers[0].consolidation = ground::Consolidation::MccPcOverMeanEffective { ratio: 2. };
    assert!(ground.validate().is_err());
    let raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-local.json"))?;
    let parsed: Ground = serde_json::from_value(raw.clone())?;
    assert_eq!(raw, serde_json::to_value(parsed)?);
    let mut bad = raw;
    bad["materials"][0]["mohr_coulomb"] = serde_json::json!({"young_kpa":1.,"poisson":0.2,"cohesion_kpa":1.,"friction_deg":0.,"dilation_deg":0.,"flow":"non_associated","tension_cutoff_kpa":null});
    assert!(serde_json::from_value::<Ground>(bad).is_err());
    Ok(())
}
