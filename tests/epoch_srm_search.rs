use serde_json::json;
use soil_fem_3d::{
    Result, checkpoint,
    epoch_srm_search::{self as epoch, Specification},
    ground::Ground,
    mesh::Mesh,
    soil,
    srm_bound::Problem,
    srm_search::{Outcome, Requirement, StopReason},
    stages::{self, Run, StageCase},
};

fn inputs() -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let p: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
    raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    raw["materials"] = json!([{"id":1,"name":"unconfined epoch SRM","unit_weight_kn_m3":1e-10,"mohr_coulomb":p.material}]);
    raw["bodies"] = json!([{"body":"TestCube","material_id":1}]);
    let ground = serde_json::from_value(raw)?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(p.boundary);
    case.stages = serde_json::from_value(
        json!([{"kind":"boundary","name":"load","target_factor":1.,"initial_increment":1.}]),
    )?;
    case.verification_reject_once_stages.clear();
    case.monitors.clear();
    case.monitor_sets.clear();
    let search = serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    let spec = Specification {
        schema: "soil-fem-epoch-srm-search/1".into(),
        search,
        initial_increment: 0.25,
        minimum_increment: 1. / 4096.,
        maximum_attempts: 128,
    };
    Ok((mesh, ground, case, spec))
}
fn copy(run: &Run) -> Result<Run> {
    Ok(serde_json::from_slice(&serde_json::to_vec(run)?)?)
}

#[test]
fn accepted_origin_search_brackets_four_over_three_and_keeps_every_factor_independent() -> Result<()>
{
    let (mesh, g, c, spec) = inputs()?;
    let origin = stages::solve(&mesh, g, c)?;
    let before = checkpoint::hash_json(&origin)?;
    let report = epoch::run(&mesh, copy(&origin)?, spec)?;
    assert_eq!(before, report.origin_run_sha256);
    assert_eq!(before, checkpoint::hash_json(&origin)?);
    assert_eq!(
        before,
        checkpoint::hash_json(&report.original_unreduced_run)?
    );
    assert_eq!(report.requirement, Requirement::EquilibriumRequirementMet);
    assert_eq!(report.stop_reason, StopReason::FactorTolerance);
    let b = report.interval.as_ref().unwrap();
    assert!(b.attained_equilibrium_factor < 4. / 3.);
    assert!(b.independently_excluded_factor > 4. / 3. && b.width <= 0.002);
    assert!(!report.continuum_stability_or_design_safety_certified);
    for trial in &report.trials {
        assert!(trial.evidence.run.is_none());
        if trial.evidence.outcome == Outcome::EquilibriumAttained {
            let saved = trial.epoch.as_ref().unwrap();
            assert_eq!(saved.payload.continuation.progress, 1.);
            let inner = &saved.payload.material_checkpoint.payload;
            assert_eq!(inner.origin_run_sha256, before);
            assert_eq!(inner.origin_state_sha256, report.origin_state_sha256);
            assert_eq!(inner.identity, report.source_identity);
            assert_eq!(
                inner.parameters[0].cohesion_kpa,
                10. / trial.evidence.factor
            );
            let frame = &trial.evidence.endpoint_audit.as_ref().unwrap().audit;
            assert_eq!(
                checkpoint::hash_json(&frame.state)?,
                checkpoint::hash_json(&inner.committed)?
            );
            assert!(frame.relative_force_residual <= origin.case.settings.relative_force_tolerance);
            assert!(
                frame.maximum_normalized_yield <= 1e-8 && frame.maximum_constraint_error_m <= 1e-10
            );
        } else {
            assert_eq!(trial.evidence.outcome, Outcome::PhysicallyExcluded);
            assert!(trial.epoch.is_none());
            assert!(trial.evidence.factor > 4. / 3.);
        }
    }
    Ok(())
}

#[test]
fn required_factor_above_limit_has_independent_exclusion_without_an_epoch() -> Result<()> {
    let (mesh, g, c, mut spec) = inputs()?;
    spec.search.required_factor = 1.6;
    let report = epoch::run(&mesh, stages::solve(&mesh, g, c)?, spec)?;
    assert_eq!(report.requirement, Requirement::IndependentlyExcluded);
    let t = report
        .trials
        .iter()
        .find(|t| t.evidence.factor == 1.6)
        .unwrap();
    assert_eq!(t.evidence.outcome, Outcome::PhysicallyExcluded);
    assert!(t.epoch.is_none());
    Ok(())
}

#[test]
fn incomplete_epoch_paths_stay_unknown_and_cannot_create_an_attained_side() -> Result<()> {
    let (mesh, g, c, mut spec) = inputs()?;
    spec.maximum_attempts = 1;
    let report = epoch::run(&mesh, stages::solve(&mesh, g, c)?, spec)?;
    assert_eq!(report.requirement, Requirement::Indeterminate);
    assert_eq!(report.stop_reason, StopReason::IndeterminateTrial);
    assert!(report.interval.is_none());
    for t in &report.trials {
        if t.evidence.outcome == Outcome::Indeterminate {
            assert_eq!(
                t.evidence.failure_kind.as_deref(),
                Some("epoch_path_incomplete")
            );
            assert_eq!(
                t.epoch.as_ref().unwrap().payload.continuation.progress,
                0.25
            );
        } else {
            assert_eq!(t.evidence.outcome, Outcome::PhysicallyExcluded);
        }
    }
    Ok(())
}

#[test]
fn factor_trial_limit_is_reported_without_claiming_the_requested_interval_width() -> Result<()> {
    let (mesh, g, c, mut spec) = inputs()?;
    spec.search.maximum_trials = 3;
    let report = epoch::run(&mesh, stages::solve(&mesh, g, c)?, spec)?;
    assert_eq!(report.trials.len(), 3);
    assert_eq!(report.stop_reason, StopReason::TrialLimit);
    assert!(
        report.interval.as_ref().unwrap().width
            > report.specification.search.absolute_factor_tolerance
    );
    Ok(())
}

#[test]
fn malformed_schema_prefix_path_controls_and_virtual_dimensions_are_refused() -> Result<()> {
    let (mesh, g, c, spec) = inputs()?;
    let origin = stages::solve(&mesh, g, c)?;
    let mut cases = vec![];
    let mut s = spec.clone();
    s.schema = "wrong".into();
    cases.push(s);
    let mut s = spec.clone();
    s.search.through_stages = 0;
    cases.push(s);
    let mut s = spec.clone();
    s.initial_increment = 0.;
    cases.push(s);
    let mut s = spec.clone();
    s.minimum_increment = 0.5;
    cases.push(s);
    let mut s = spec.clone();
    s.maximum_attempts = 0;
    cases.push(s);
    let mut s = spec.clone();
    s.search.schema = "soil-fem-srm-search/2".into();
    s.search.virtual_displacements_m = vec![vec![[0.; 3]; mesh.nodes.len() - 1]];
    cases.push(s);
    for bad in cases {
        assert!(epoch::run(&mesh, copy(&origin)?, bad).is_err());
    }
    let mut raw = serde_json::to_value(spec)?;
    raw["unexpected"] = json!(true);
    assert!(serde_json::from_value::<Specification>(raw).is_err());
    Ok(())
}

#[test]
fn changed_producer_inputs_and_unbound_loads_cannot_define_an_origin() -> Result<()> {
    let (mesh, g, c, spec) = inputs()?;
    let origin = stages::solve(&mesh, g, c)?;
    let mut bad = copy(&origin)?;
    bad.case.boundary.as_mut().unwrap().constraints[0].value_m = 1e-4;
    assert!(epoch::run(&mesh, bad, spec.clone()).is_err());
    let mut bad = copy(&origin)?;
    bad.last_committed.applied_nodal_load_kn[0] += 0.1;
    assert!(epoch::run(&mesh, bad, spec.clone()).is_err());
    let mut bad = copy(&origin)?;
    bad.complete = false;
    assert!(epoch::run(&mesh, bad, spec).is_err());
    Ok(())
}

fn plastic_origin(cut: bool) -> Result<(Mesh, Run, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut g: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    for m in &mut g.materials {
        m.unit_weight_kn_m3 = 1e-10;
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            panic!("MC")
        };
        p.young_kpa = 100000.;
        p.poisson = 0.3;
        p.cohesion_kpa = 40.;
        p.friction_deg = 0.;
        p.dilation_deg = 0.;
        p.tension_cutoff_kpa = None;
    }
    let mut c: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    let mut steps = vec![
        json!({"kind":"pressure","name":"plastic_preload","target_pressure_kpa":200.}),
        json!({"kind":"pressure","name":"unload","target_pressure_kpa":150.}),
    ];
    if cut {
        steps.push(json!({"kind":"excavate","name":"cut_top","bodies":["UpperClay"],"remaining_surface_pressure_kpa":0.0}));
    }
    steps.push(json!({"kind":"reset_display","name":"reference"}));
    c.stages = serde_json::from_value(json!(steps))?;
    c.verification_reject_once_stages.clear();
    let (_, _, _, mut spec) = inputs()?;
    spec.search.material_ids = g.materials.iter().map(|m| m.id).collect();
    spec.search.required_factor = 2.;
    spec.search.maximum_factor = 3.;
    spec.search.through_stages = c.stages.len();
    let origin = stages::solve(&mesh, g, c)?;
    assert!(origin.complete, "{:?}", origin.failure_reason);
    Ok((mesh, origin, spec))
}

#[test]
fn plastic_preload_unload_epoch_fields_match_confined_column_compliance_and_separate_path()
-> Result<()> {
    let (mesh, origin, spec) = plastic_origin(false)?;
    let hash = checkpoint::hash_json(&origin)?;
    let report = epoch::run(&mesh, copy(&origin)?, spec)?;
    assert_eq!(report.requirement, Requirement::EquilibriumRequirementMet);
    assert_eq!(report.stop_reason, StopReason::NoPhysicalUpper);
    assert!(report.interval.is_none());
    let zmin = mesh
        .nodes
        .iter()
        .map(|p| p[2])
        .fold(f64::INFINITY, f64::min);
    for trial in &report.trials {
        assert_eq!(trial.evidence.outcome, Outcome::EquilibriumAttained);
        let checkpoint = trial.epoch.as_ref().unwrap();
        let state = &checkpoint.payload.material_checkpoint.payload.committed;
        let lateral =
            (120_f64 + (0.3 / 0.7) * (150. - 200.)).max(150. - 80. / trial.evidence.factor);
        let axial = 0.4 * (150. + 2. * lateral) / 100000.;
        for point in &state.material_states {
            let soil::State::MohrCoulomb(s) = point else {
                panic!("MC")
            };
            for (actual, expected) in s
                .effective_stress_compression_positive_kpa
                .iter()
                .zip([lateral, lateral, 150., 0., 0., 0.])
            {
                assert!((actual - expected).abs() < 1e-6);
            }
            assert!(s.equivalent_plastic_strain > 1e-5);
            assert!((s.void_ratio - (1.7 * (-axial).exp() - 1.)).abs() < 1e-9);
        }
        for (point, u) in mesh.nodes.iter().zip(&state.displacement_m) {
            assert!((u[2] + axial * (point[2] - zmin)).abs() < 1e-9);
        }
        assert_eq!(
            state.display_reference_displacement_m,
            origin.last_committed.display_reference_displacement_m
        );
    }
    let required = report
        .trials
        .iter()
        .find(|t| t.evidence.factor == 2.)
        .unwrap();
    let saved = required.epoch.as_ref().unwrap();
    let mut independent = soil_fem_3d::mc_epoch_path::Controller::new(
        &mesh,
        &origin,
        saved.payload.specification.clone(),
    )?;
    while !independent.complete() && independent.continuation().stop_reason.is_none() {
        independent.advance()?;
    }
    assert_eq!(
        serde_json::to_vec(saved)?,
        serde_json::to_vec(&independent.checkpoint()?)?
    );
    assert_eq!(hash, checkpoint::hash_json(&origin)?);
    Ok(())
}

#[test]
fn post_excavation_trials_keep_all_dormant_history_and_reject_resealed_frozen_metadata()
-> Result<()> {
    let (mesh, origin, spec) = plastic_origin(true)?;
    let report = epoch::run(&mesh, copy(&origin)?, spec)?;
    assert_eq!(
        report.trials[0].evidence.outcome,
        Outcome::EquilibriumAttained
    );
    assert_eq!(report.requirement, Requirement::EquilibriumRequirementMet);
    let dormant: Vec<_> = origin
        .integration_points
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (!origin.last_committed.active_elements[p.element_index]).then_some(i))
        .collect();
    assert!(!dormant.is_empty());
    for trial in &report.trials {
        if let Some(saved) = &trial.epoch {
            let current = &saved.payload.material_checkpoint.payload.committed;
            for &i in &dormant {
                assert_eq!(
                    checkpoint::hash_json(&current.material_states[i])?,
                    checkpoint::hash_json(&origin.last_committed.material_states[i])?
                );
            }
        }
    }
    let saved = report.trials.iter().find_map(|t| t.epoch.as_ref()).unwrap();
    let mut bad: stages::McEpochCheckpoint =
        serde_json::from_slice(&serde_json::to_vec(&saved.payload.material_checkpoint)?)?;
    bad.payload.committed.applied_nodal_load_kn[0] += 1e-12;
    bad.payload_sha256 = checkpoint::hash_json(&bad.payload)?;
    assert!(epoch::audit_epoch_endpoint(&mesh, &origin, &bad).is_err());
    Ok(())
}

#[test]
fn cli_creates_a_bound_epoch_report_refuses_foreign_producer_and_preserves_existing_output()
-> Result<()> {
    use std::{
        fs,
        path::PathBuf,
        process::Command,
        time::{SystemTime, UNIX_EPOCH},
    };
    let temp = fs::canonicalize(std::env::temp_dir())?;
    let directory = temp.join(format!(
        "soil-epoch-srm-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&directory)?;
    let (mesh, g, c, s) = inputs()?;
    for (name, value) in [
        ("mesh", serde_json::to_value(&mesh)?),
        ("ground", serde_json::to_value(&g)?),
        ("case", serde_json::to_value(&c)?),
        ("spec", serde_json::to_value(&s)?),
    ] {
        fs::write(
            directory.join(format!("{name}.json")),
            serde_json::to_vec(&value)?,
        )?;
    }
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_soil-fem-3d"));
    let path = |name: &str| directory.join(format!("{name}.json"));
    let solve = Command::new(&exe)
        .args(["solve-stages"])
        .args([path("mesh"), path("ground"), path("case"), path("origin")])
        .output()?;
    assert!(
        solve.status.success(),
        "{}",
        String::from_utf8_lossy(&solve.stderr)
    );
    let original_bytes = fs::read(path("origin"))?;
    let search = Command::new(&exe)
        .arg("search-mc-epoch-srm")
        .args([path("mesh"), path("origin"), path("spec"), path("report")])
        .output()?;
    assert!(
        search.status.success(),
        "{}",
        String::from_utf8_lossy(&search.stderr)
    );
    let report_bytes = fs::read(path("report"))?;
    let report: epoch::Report = serde_json::from_slice(&report_bytes)?;
    assert_eq!(report.requirement, Requirement::EquilibriumRequirementMet);
    assert_eq!(
        report.source_identity.solver_binary_sha256,
        checkpoint::hash_file(&exe)?
    );
    assert_eq!(fs::read(path("origin"))?, original_bytes);
    let repeat = Command::new(&exe)
        .arg("search-mc-epoch-srm")
        .args([
            path("missing"),
            path("missing"),
            path("missing"),
            path("report"),
        ])
        .output()?;
    assert!(!repeat.status.success());
    assert!(String::from_utf8_lossy(&repeat.stderr).contains("output already exists"));
    assert_eq!(fs::read(path("report"))?, report_bytes);
    fs::write(
        path("foreign"),
        serde_json::to_vec(&stages::solve(&mesh, g, c)?)?,
    )?;
    let foreign = Command::new(&exe)
        .arg("search-mc-epoch-srm")
        .args([path("mesh"), path("foreign"), path("spec"), path("bad")])
        .output()?;
    assert!(!foreign.status.success());
    assert!(!path("bad").exists());
    assert!(directory.is_absolute() && directory.parent() == Some(temp.as_path()));
    assert!(
        directory
            .file_name()
            .unwrap()
            .to_string_lossy()
            .starts_with("soil-epoch-srm-")
    );
    fs::remove_dir_all(directory)?;
    Ok(())
}
