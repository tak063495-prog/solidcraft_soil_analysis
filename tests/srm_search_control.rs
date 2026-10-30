use serde_json::json;
use soil_fem_3d::{
    Result,
    ground::Ground,
    mesh::Mesh,
    soil,
    srm_bound::Problem,
    srm_search::{self, Outcome, Requirement, Specification, StopReason},
    stages::StageCase,
};

fn fixture() -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let p: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
    raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    raw["materials"] = json!([{"id":1,"name":"Unconfined proof-of-search",
        "unit_weight_kn_m3":1e-10,"mohr_coulomb":p.material}]);
    raw["bodies"] = json!([{"body":"TestCube","material_id":1}]);
    let ground: Ground = serde_json::from_value(raw)?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(p.boundary);
    case.stages = serde_json::from_value(
        json!([{"kind":"boundary","name":"load","target_factor":1.0,"initial_increment":1.0}]),
    )?;
    case.verification_reject_once_stages = vec![];
    case.monitor_sets.clear();
    case.monitors.clear();
    let spec: Specification = serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    Ok((mesh, ground, case, spec))
}

fn plastic_column() -> Result<(Mesh, Ground, StageCase, Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut ground: Ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    for m in &mut ground.materials {
        m.unit_weight_kn_m3 = 1e-10;
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            panic!("MC fixture")
        };
        p.cohesion_kpa = 40.;
        p.friction_deg = 0.;
        p.dilation_deg = 0.;
        p.tension_cutoff_kpa = None;
    }
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.stages = serde_json::from_value(
        json!([{"kind":"pressure","name":"load","target_pressure_kpa":100.0}]),
    )?;
    case.verification_reject_once_stages = vec![];
    let (_, _, _, mut spec) = fixture()?;
    spec.material_ids = ground.materials.iter().map(|m| m.id).collect();
    spec.maximum_factor = 4.;
    spec.required_factor = 2.;
    Ok((mesh, ground, case, spec))
}

fn scientific(mut value: serde_json::Value) -> serde_json::Value {
    fn clean(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Object(map) => {
                map.remove("elapsed_seconds");
                for child in map.values_mut() {
                    clean(child);
                }
            }
            serde_json::Value::Array(values) => {
                for child in values {
                    clean(child);
                }
            }
            _ => {}
        }
    }
    clean(&mut value);
    value
}

#[test]
fn distributed_frictional_virtual_upper_survives_all_restarts_without_history_transfer()
-> Result<()> {
    use soil_fem_3d::{
        boundary::{Conditions, SurfaceLoad},
        srm_virtual_work as work,
    };
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-virtual-work-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let mut p: work::Problem =
        serde_json::from_str(include_str!("../examples/mc-virtual-work-problem.json"))?;
    p.materials[0].parameters.cohesion_kpa = 10.;
    p.materials[0].parameters.friction_deg = 30.;
    let b = Conditions {
        constraints: p.constraints.clone(),
        nodal_loads: vec![],
        surface_loads: vec![SurfaceLoad::Traction {
            face_set: "face-2-1".into(),
            traction_kpa: [0., 0., -30.],
            gradient_kpa_per_m: [[0.; 3]; 3],
            reference_point_m: [0.; 3],
        }],
    };
    let (_, _, mut case, mut spec) = fixture()?;
    case.boundary = Some(b);
    let mut g: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = g["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    g["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube":seed});
    g["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    g["materials"] = json!([{"id":1,"name":"Frictional reference","unit_weight_kn_m3":1e-10,"mohr_coulomb":p.materials[0].parameters}]);
    g["bodies"] = json!([{"body":"TestCube","material_id":1}]);
    let ground: Ground = serde_json::from_value(g)?;
    spec.schema = "soil-fem-srm-search/2".into();
    spec.maximum_factor = 1.5;
    spec.required_factor = 1.05;
    let phi = (30_f64.to_radians().tan() / 1.105).atan();
    let ratio = (1. + phi.sin()) / (1. - phi.sin());
    spec.virtual_displacements_m = vec![
        mesh.nodes
            .iter()
            .map(|x| [(ratio + 1e-4) * x[0], 0., -x[2]])
            .collect(),
    ];
    let expected = srm_search::run(&mesh, ground.clone(), case.clone(), spec.clone())?;
    let (saved, count) = advance_all(&mesh, &ground, &case, &spec)?;
    assert!(count > 0);
    let result = saved.report().unwrap();
    let bounds = result.interval.as_ref().unwrap();
    let critical = 2. * (10. * (30. * 30_f64.to_radians().tan() + 10.)).sqrt() / 30.;
    assert!(
        bounds.attained_equilibrium_factor < critical
            && bounds.independently_excluded_factor > critical
    );
    for trial in &result.trials {
        if trial.outcome == Outcome::PhysicallyExcluded {
            assert!(matches!(
                trial.physical_exclusion,
                soil_fem_3d::srm_bound::Assessment::NotApplicable { .. }
            ));
            assert!(trial.virtual_work.iter().any(|e|matches!(&e.assessment,work::Assessment::Evaluated{certificate} if certificate.excludes_admissible_discrete_equilibrium)));
            assert!(trial.run.is_none());
        }
    }
    assert_eq!(
        scientific(serde_json::to_value(result)?),
        scientific(serde_json::to_value(expected)?)
    );
    Ok(())
}

fn advance_all(
    mesh: &Mesh,
    g: &Ground,
    c: &StageCase,
    spec: &Specification,
) -> Result<(soil_fem_3d::srm_search_control::Snapshot, usize)> {
    use soil_fem_3d::srm_search_control::Snapshot;
    let mut saved = Snapshot::start(mesh, g, c, spec)?;
    let mut accepted_saves = 0;
    for _ in 0..1000 {
        if saved.is_complete() {
            return Ok((saved, accepted_saves));
        }
        // A real disk/JSON boundary: restore every single accepted chunk.
        let bytes = serde_json::to_vec(&saved)?;
        saved = serde_json::from_slice::<Snapshot>(&bytes)?.advance(mesh, g, c, spec, 1, None)?;
        if saved.has_accepted_mechanical_checkpoint() {
            accepted_saves += 1;
        }
    }
    panic!("bounded tiny mechanical reference did not finish")
}

#[test]
fn restart_every_increment_matches_uninterrupted_closed_form_search_and_all_scientific_history()
-> Result<()> {
    let (mesh, g, mut c, mut spec) = fixture()?;
    c.verification_reject_once_stages = vec!["load".into()];
    spec.absolute_factor_tolerance = 0.02;
    let expected = srm_search::run(&mesh, g.clone(), c.clone(), spec.clone())?;
    let (saved, count) = advance_all(&mesh, &g, &c, &spec)?;
    assert!(
        count > 4,
        "actual accepted mechanical checkpoints must occur across factors"
    );
    saved.verify(&mesh, &g, &c, &spec)?;
    let actual = saved.report().unwrap();
    let bounds = actual.interval.as_ref().unwrap();
    assert!(bounds.attained_equilibrium_factor < 20. / 15.);
    assert!(bounds.independently_excluded_factor > 20. / 15.);
    assert_eq!(actual.requirement, Requirement::EquilibriumRequirementMet);
    assert_eq!(
        scientific(serde_json::to_value(actual)?),
        scientific(serde_json::to_value(expected)?)
    );
    assert!(
        actual
            .original_unreduced_run
            .attempts
            .iter()
            .any(|a| !a.accepted)
    );
    assert!(!actual.continuum_stability_or_design_safety_certified);
    Ok(())
}

#[test]
fn original_geostatic_and_each_factor_keep_their_own_accepted_checkpoint() -> Result<()> {
    use soil_fem_3d::srm_search_control::Snapshot;
    let (mesh, mut g, mut c, mut spec) = fixture()?;
    for m in &mut g.materials {
        m.unit_weight_kn_m3 = 20.;
        let soil::Model::MohrCoulomb(p) = &mut m.model else {
            unreachable!()
        };
        p.cohesion_kpa = 1000.;
    }
    c.settings.initial_increment = 0.25;
    spec.required_factor = 1.;
    spec.maximum_factor = 1.1;
    spec.absolute_factor_tolerance = 0.1;
    let mut saved = Snapshot::start(&mesh, &g, &c, &spec)?;
    saved = saved.advance(&mesh, &g, &c, &spec, 1, None)?;
    assert!(saved.has_accepted_mechanical_checkpoint());
    assert_eq!(saved.current_factor(), Some(1.));
    let raw = serde_json::to_value(&saved)?;
    assert_eq!(
        raw["payload"]["inflight"]["run"]["last_committed"]["stage_index"],
        0
    );
    assert!(
        raw["payload"]["inflight"]["run"]["last_committed"]["stage_progress"]
            .as_f64()
            .unwrap()
            < 1.
    );
    let expected = srm_search::run(&mesh, g.clone(), c.clone(), spec.clone())?;
    let (completed, n) = advance_all(&mesh, &g, &c, &spec)?;
    assert!(n > 2);
    assert_eq!(
        scientific(serde_json::to_value(completed.report().unwrap())?),
        scientific(serde_json::to_value(expected)?)
    );
    Ok(())
}

#[test]
fn changed_material_mesh_full_case_specification_and_saved_history_are_rejected() -> Result<()> {
    use soil_fem_3d::srm_search_control::Snapshot;
    let (mesh, g, c, spec) = fixture()?;
    let saved = Snapshot::start(&mesh, &g, &c, &spec)?.advance(&mesh, &g, &c, &spec, 1, None)?;
    let mut changed_g = g.clone();
    changed_g.materials[0].unit_weight_kn_m3 *= 2.;
    assert!(saved.verify(&mesh, &changed_g, &c, &spec).is_err());
    let mut changed_mesh = mesh.clone();
    changed_mesh.nodes[0][0] += 0.001;
    assert!(saved.verify(&changed_mesh, &g, &c, &spec).is_err());
    let mut changed_case = c.clone();
    changed_case.monitor_sets.insert("new".into(), "all".into());
    assert!(saved.verify(&mesh, &g, &changed_case, &spec).is_err());
    let mut changed_spec = spec.clone();
    changed_spec.required_factor = 1.2;
    assert!(saved.verify(&mesh, &g, &c, &changed_spec).is_err());
    for key in ["payload_sha256", "specification_sha256"] {
        let mut corrupt = serde_json::to_value(&saved)?;
        corrupt[key] = json!("bad");
        assert!(
            serde_json::from_value::<Snapshot>(corrupt)?
                .verify(&mesh, &g, &c, &spec)
                .is_err()
        );
    }
    let mut corrupt = serde_json::to_value(&saved)?;
    corrupt["identity"]["solver_binary_sha256"] = json!("another-binary");
    assert!(
        serde_json::from_value::<Snapshot>(corrupt)?
            .verify(&mesh, &g, &c, &spec)
            .is_err()
    );
    let mut corrupt = serde_json::to_value(&saved)?;
    corrupt["payload"]["trials"] = json!([{}]);
    assert!(serde_json::from_value::<Snapshot>(corrupt).is_err());
    Ok(())
}

#[test]
fn genuine_nonlinear_failure_after_restarts_is_unknown_and_never_an_upper() -> Result<()> {
    let (mesh, mut g, mut c, mut spec) = plastic_column()?;
    g.initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    c.settings.maximum_newton_iterations = 3;
    c.settings.initial_increment = 1.;
    c.settings.minimum_increment = 1.;
    spec.required_factor = 4.;
    spec.maximum_factor = 4.;
    let expected = srm_search::run(&mesh, g.clone(), c.clone(), spec.clone())?;
    let (saved, count) = advance_all(&mesh, &g, &c, &spec)?;
    assert!(count > 0);
    let result = saved.report().unwrap();
    assert_eq!(result.stop_reason, StopReason::IndeterminateTrial);
    assert!(result.interval.is_none());
    let trial = result.trials.iter().find(|t| t.factor == 4.).unwrap();
    assert_eq!(trial.outcome, Outcome::Indeterminate);
    assert_eq!(trial.failure_kind.as_deref(), Some("incomplete"));
    let failed_run = trial
        .run
        .as_ref()
        .expect("incomplete result keeps its accepted history");
    assert!(!failed_run.complete);
    let reason = failed_run.failure_reason.as_deref().unwrap();
    assert!(reason.contains("Newton iteration limit") && reason.contains("cutback exhausted"));
    assert!(
        failed_run
            .attempts
            .iter()
            .any(|a| !a.accepted && a.committed_state_unchanged_on_rejection)
    );
    assert_eq!(
        scientific(serde_json::to_value(result)?),
        scientific(serde_json::to_value(expected)?)
    );
    Ok(())
}

#[test]
fn observer_failure_leaves_last_saved_snapshot_valid_and_cannot_create_a_trial_failure()
-> Result<()> {
    use soil_fem_3d::{
        Error,
        srm_search_control::Snapshot,
        stage_control::{Observer, Progress},
    };
    struct Broken;
    impl Observer for Broken {
        fn on_progress(&self, _: &Progress) -> Result<()> {
            Err(Error::Control("intentional observer IO refusal".into()))
        }
    }
    let (mesh, g, c, spec) = fixture()?;
    let saved = Snapshot::start(&mesh, &g, &c, &spec)?;
    let bytes = serde_json::to_vec(&saved)?;
    let attempt =
        serde_json::from_slice::<Snapshot>(&bytes)?.advance(&mesh, &g, &c, &spec, 1, Some(&Broken));
    assert!(matches!(attempt, Err(Error::Control(_))));
    saved.verify(&mesh, &g, &c, &spec)?;
    assert_eq!(serde_json::to_vec(&saved)?, bytes);
    assert!(saved.trials().is_empty() && saved.report().is_none());
    assert!(
        serde_json::from_slice::<Snapshot>(&bytes)?
            .advance(&mesh, &g, &c, &spec, 0, None)
            .is_err()
    );
    Ok(())
}
