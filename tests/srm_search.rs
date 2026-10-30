use serde_json::json;
use soil_fem_3d::{
    Result,
    checkpoint::hash_json,
    ground::Ground,
    mesh::Mesh,
    soil,
    srm_bound::{self, Assessment, NodalProblem, Problem, WaterPolicy},
    srm_search::{self, Outcome, Requirement, Specification, StopReason},
    stages::{self, StageCase},
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

#[test]
fn real_fem_search_brackets_independent_unconfined_strength_and_evaluates_requested_factor()
-> Result<()> {
    let (mesh, g, c, spec) = fixture()?;
    let result = srm_search::run(&mesh, g, c, spec)?;
    let interval = result
        .interval
        .as_ref()
        .expect("independent physical bracket");
    let critical = 20. / 15.;
    assert!(interval.attained_equilibrium_factor < critical);
    assert!(interval.independently_excluded_factor > critical);
    assert!(interval.width <= result.specification.absolute_factor_tolerance);
    assert_eq!(result.stop_reason, StopReason::FactorTolerance);
    assert_eq!(result.requirement, Requirement::EquilibriumRequirementMet);
    assert!(
        result
            .trials
            .iter()
            .any(|t| t.factor == 1.1 && t.outcome == Outcome::EquilibriumAttained)
    );
    assert!(result.unreduced_source_history_unchanged);
    assert!(!result.continuum_stability_or_design_safety_certified);
    for t in &result.trials {
        if t.outcome == Outcome::PhysicallyExcluded {
            assert!(t.run.is_none());
            assert!(
                matches!(&t.physical_exclusion,Assessment::Evaluated{certificate}
                if certificate.mean_stress_excludes_admissible_discrete_equilibrium)
            );
        } else {
            let audit = t
                .endpoint_audit
                .as_ref()
                .expect("actual complete FEM audited");
            assert!(audit.audit.relative_force_residual < 1e-8);
            assert!(audit.audit.maximum_normalized_yield <= 1e-8);
            assert!(audit.audit.maximum_constraint_error_m <= 1e-10);
        }
    }
    Ok(())
}

#[test]
fn requested_factor_above_the_closed_form_limit_is_independently_excluded() -> Result<()> {
    let (mesh, g, c, mut spec) = fixture()?;
    spec.required_factor = 1.6;
    let result = srm_search::run(&mesh, g, c, spec)?;
    assert_eq!(result.requirement, Requirement::IndependentlyExcluded);
    let t = result.trials.iter().find(|t| t.factor == 1.6).unwrap();
    assert_eq!(t.outcome, Outcome::PhysicallyExcluded);
    assert!(t.failure_kind.is_none() && t.run.is_none());
    Ok(())
}

#[test]
fn distributed_support_and_actual_plastic_history_never_invent_a_physical_upper() -> Result<()> {
    let (mesh, g, c, spec) = plastic_column()?;
    let result = srm_search::run(&mesh, g, c, spec)?;
    assert!(result.interval.is_none());
    assert_eq!(result.stop_reason, StopReason::NoPhysicalUpper);
    assert_eq!(result.requirement, Requirement::EquilibriumRequirementMet);
    let t = result.trials.iter().find(|t| t.factor == 4.).unwrap();
    assert!(matches!(
        t.physical_exclusion,
        Assessment::NotApplicable { .. }
    ));
    assert_eq!(t.outcome, Outcome::EquilibriumAttained);
    assert!(
        t.run
            .as_ref()
            .unwrap()
            .last_committed
            .material_states
            .iter()
            .any(|s| matches!(s,soil::State::MohrCoulomb(m) if m.equivalent_plastic_strain>0.))
    );
    Ok(())
}

#[test]
fn different_trial_order_reproduces_the_full_same_factor_plastic_state() -> Result<()> {
    let (mesh, g, c, mut spec) = plastic_column()?;
    let first = srm_search::run(&mesh, g.clone(), c.clone(), spec.clone())?;
    spec.required_factor = 3.;
    let second = srm_search::run(&mesh, g, c, spec)?;
    let last = |r: &srm_search::Report| {
        hash_json(
            &r.trials
                .iter()
                .find(|t| t.factor == 4.)
                .unwrap()
                .run
                .as_ref()
                .unwrap()
                .last_committed,
        )
    };
    assert_eq!(last(&first)?, last(&second)?);
    assert_ne!(first.trials[1].factor, second.trials[1].factor);
    assert_eq!(
        hash_json(&first.original_unreduced_run.last_committed)?,
        hash_json(&second.original_unreduced_run.last_committed)?
    );
    Ok(())
}

#[test]
fn actual_rejected_loading_trial_restores_history_for_every_fresh_replay() -> Result<()> {
    let (mesh, g, mut c, spec) = fixture()?;
    c.verification_reject_once_stages = vec!["load".into()];
    let result = srm_search::run(&mesh, g, c, spec)?;
    for run in std::iter::once(&result.original_unreduced_run)
        .chain(result.trials.iter().filter_map(|t| t.run.as_ref()))
    {
        assert!(run.attempts.iter().any(|a| !a.accepted));
        assert!(
            run.attempts
                .iter()
                .filter(|a| !a.accepted)
                .all(|a| a.committed_state_unchanged_on_rejection)
        );
    }
    Ok(())
}

#[test]
fn actual_nonlinear_failure_is_indeterminate_without_a_fabricated_upper() -> Result<()> {
    let (mesh, mut g, mut c, mut spec) = plastic_column()?;
    // A separate actual A4 diagnostic retained the original failing origin:
    // its geostatic verification rejection exhausted the one-step budget.
    // With that unrelated injection disabled, cap3 accepts F1 and rejects F4.
    g.initialization
        .geostatic
        .as_mut()
        .unwrap()
        .verification_reject_once = false;
    c.settings.maximum_newton_iterations = 3;
    c.settings.initial_increment = 1.;
    c.settings.minimum_increment = 1.;
    spec.required_factor = 4.;
    let result = srm_search::run(&mesh, g, c, spec)?;
    let failed = result.trials.iter().find(|t| t.factor == 4.).unwrap();
    assert_eq!(failed.outcome, Outcome::Indeterminate);
    assert!(failed.reason.is_some());
    assert!(result.interval.is_none());
    assert_eq!(result.requirement, Requirement::Indeterminate);
    assert_eq!(result.stop_reason, StopReason::IndeterminateTrial);
    assert!(
        result
            .trials
            .iter()
            .all(|t| t.outcome != Outcome::PhysicallyExcluded)
    );
    Ok(())
}

#[test]
fn saved_endpoint_reconstruction_rejects_changed_load_stress_and_mesh() -> Result<()> {
    let (mesh, g, c, spec) = fixture()?;
    let result = srm_search::run(&mesh, g, c, spec)?;
    let mut changed: stages::Run =
        serde_json::from_value(serde_json::to_value(&result.original_unreduced_run)?)?;
    changed.last_committed.applied_nodal_load_kn[3] += 1.;
    assert!(srm_search::audit_unreduced_endpoint(&mesh, &changed).is_err());
    changed = serde_json::from_value(serde_json::to_value(&result.original_unreduced_run)?)?;
    let soil::State::MohrCoulomb(s) = &mut changed.last_committed.material_states[0] else {
        panic!("MC")
    };
    s.effective_stress_compression_positive_kpa[0] += 100.;
    assert!(srm_search::audit_unreduced_endpoint(&mesh, &changed).is_err());
    let mut moved = mesh.clone();
    moved.nodes[1][0] += 0.01;
    assert!(srm_search::audit_unreduced_endpoint(&moved, &result.original_unreduced_run).is_err());
    Ok(())
}

#[test]
fn assembled_nodal_witness_matches_surface_loads_and_rejects_inactive_node_forces() -> Result<()> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let p: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut nodal = NodalProblem {
        schema: "soil-fem-mc-nodal-mean-stress-bound/1".into(),
        material_id: 1,
        water_policy: WaterPolicy::Dry,
        material: p.material.clone(),
        constraints: p.boundary.constraints.clone(),
        physical_nodal_force_kn: p.boundary.loads(&mesh, &[true; 6])?,
        active_elements: vec![true; 6],
    };
    let Assessment::Evaluated {
        certificate: original,
    } = srm_bound::analyze(&mesh, &p)?
    else {
        panic!("witness")
    };
    let Assessment::Evaluated {
        certificate: assembled,
    } = srm_bound::analyze_nodal(&mesh, &nodal)?
    else {
        panic!("witness")
    };
    assert_eq!(
        original.mean_effective_stress_compression_positive_kpa,
        assembled.mean_effective_stress_compression_positive_kpa
    );
    nodal.active_elements = vec![true, false, false, false, false, false];
    nodal.physical_nodal_force_kn[6] = 1.;
    assert!(srm_bound::analyze_nodal(&mesh, &nodal).is_err());
    Ok(())
}

#[test]
fn invalid_search_controls_and_unsupported_water_davis_or_unknown_policies_reject() -> Result<()> {
    let (_, g, c, spec) = fixture()?;
    for f in [f64::NAN, 0.5, f64::INFINITY, 101.] {
        let mut bad = spec.clone();
        bad.maximum_factor = f;
        assert!(bad.validate(&g, &c).is_err());
    }
    let mut raw = serde_json::to_value(&spec)?;
    raw["psi_policy"] = json!("implicit");
    assert!(serde_json::from_value::<Specification>(raw).is_err());
    let mut wet = g.clone();
    wet.initialization.water_level_m = Some(0.5);
    assert!(spec.validate(&wet, &c).is_err());
    let mut davis = g;
    let soil::Model::MohrCoulomb(p) = &mut davis.materials[0].model else {
        panic!("MC")
    };
    p.flow = soil_fem_3d::mc::Flow::Davis;
    assert!(spec.validate(&davis, &c).is_err());
    Ok(())
}

#[test]
fn intermediate_prefix_keeps_the_original_unreduced_continuation_state() -> Result<()> {
    let (mesh, g, mut c, spec) = fixture()?;
    c.stages.push(serde_json::from_value(
        json!({"kind":"boundary","name":"unload","target_factor":0.5}),
    )?);
    let source = hash_json(&(&g, &c))?;
    let search = srm_search::run(&mesh, g.clone(), c.clone(), spec)?;
    let reference = stages::solve(&mesh, g.clone(), c.clone())?;
    let after_search = stages::solve(&mesh, g.clone(), c.clone())?;
    assert_eq!(
        hash_json(&reference.last_committed)?,
        hash_json(&after_search.last_committed)?
    );
    assert_eq!(hash_json(&(&g, &c))?, source);
    assert!(reference.complete && after_search.complete);
    assert_ne!(
        hash_json(&search.original_unreduced_run.last_committed)?,
        hash_json(&reference.last_committed)?
    );
    Ok(())
}

#[test]
fn reconstructed_endpoint_restores_the_actual_nonzero_geostatic_reference() -> Result<()> {
    let (mesh, mut g, c, spec) = plastic_column()?;
    for m in &mut g.materials {
        m.unit_weight_kn_m3 = 20.;
    }
    let result = srm_search::run(&mesh, g, c, spec)?;
    assert!(
        result
            .original_unreduced_run
            .reference_displacement_m
            .iter()
            .flatten()
            .any(|v| v.abs() > 1e-10)
    );
    let original = result.original_unreduced_run.frames.last().unwrap();
    let reconstructed = &result.original_endpoint_audit.audit;
    assert_eq!(
        reconstructed.maximum_absolute_engineering_strain_from_initial_equilibrium,
        original.maximum_absolute_engineering_strain_from_initial_equilibrium,
    );
    assert!(
        original
            .maximum_absolute_engineering_strain_from_initial_equilibrium
            .unwrap()
            > 1e-6
    );
    Ok(())
}
