use serde_json::{Value, json};
use soil_fem_3d::{
    Result, epoch_srm_search as epoch, exact_fs,
    ground::Ground,
    mesh::Mesh,
    srm_search as fresh,
    stages::{self, StageCase},
};
fn inputs() -> Result<(Mesh, Ground, StageCase, epoch::Specification)> {
    Ok((
        serde_json::from_str(include_str!("../examples/exact-fs-uniaxial/mesh.json"))?,
        serde_json::from_str(include_str!("../examples/exact-fs-uniaxial/ground.json"))?,
        serde_json::from_str(include_str!("../examples/exact-fs-uniaxial/case.json"))?,
        serde_json::from_str(include_str!(
            "../examples/exact-fs-uniaxial/epoch-spec.json"
        ))?,
    ))
}
fn report() -> Result<(Mesh, epoch::Specification, epoch::Report)> {
    let (m, g, c, s) = inputs()?;
    let r = epoch::run(&m, stages::solve(&m, g, c)?, s.clone())?;
    Ok((m, s, r))
}
#[test]
fn completed_material_path_and_reconstructed_mean_stress_certify_only_the_uniaxial_interval()
-> Result<()> {
    let (m, s, r) = report()?;
    let a = exact_fs::epoch(&m, &r.original_unreduced_run, &s, &r)?;
    assert!(a.exact_fs_for_this_case);
    assert!(a.attained_equilibrium_factor.unwrap() < 4. / 3.);
    assert!(a.independently_excluded_factor.unwrap() > 4. / 3.);
    assert!(a.width.unwrap() <= 0.002);
    assert!(!a.continuum_stability_or_design_safety_certified);
    assert!(serde_json::to_value(a)?.get("midpoint").is_none());
    Ok(())
}
#[test]
fn fresh_history_report_is_reaudited_without_replaying_a_solve() -> Result<()> {
    let (m, g, c, s) = inputs()?;
    let r = fresh::run(&m, g.clone(), c.clone(), s.search.clone())?;
    assert!(exact_fs::fresh(&m, &g, &c, &s.search, &r)?.exact_fs_for_this_case);
    let raw = serde_json::to_value(&r)?;
    let lo = r.interval.as_ref().unwrap().attained_trial_index;
    let hi = r.interval.as_ref().unwrap().excluded_trial_index;
    let mut changes = vec![];
    let mut b = raw.clone();
    b["trials"][lo]["endpoint_audit"] = Value::Null;
    changes.push(b);
    let mut b = raw.clone();
    b["trials"][hi]["physical_exclusion"]["certificate"]["mean_stress_excludes_admissible_discrete_equilibrium"] =
        json!(false);
    changes.push(b);
    let mut b = raw.clone();
    b["trials"][lo]["run"]["complete"] = json!(false);
    changes.push(b);
    let mut b = raw.clone();
    b["trials"].as_array_mut().unwrap().remove(1);
    changes.push(b);
    let mut b = raw;
    b["interval"]["width"] = json!(0.);
    changes.push(b);
    for b in changes {
        let r = serde_json::from_value(b)?;
        assert!(exact_fs::fresh(&m, &g, &c, &s.search, &r).is_err());
    }
    Ok(())
}
#[test]
fn missing_changed_or_fabricated_epoch_evidence_is_rejected() -> Result<()> {
    let (m, s, r) = report()?;
    let raw = serde_json::to_value(&r)?;
    let lo = r.interval.as_ref().unwrap().attained_trial_index;
    let hi = r.interval.as_ref().unwrap().excluded_trial_index;
    let mut mutations = vec![];
    for path in ["epoch", "evidence"] {
        let mut b = raw.clone();
        b["trials"][lo][path] = Value::Null;
        mutations.push(b);
    }
    let mut b = raw.clone();
    b["trials"].as_array_mut().unwrap().remove(1);
    mutations.push(b);
    let mut b = raw.clone();
    b["trials"][hi]["evidence"]["physical_exclusion"] =
        json!({"status":"not_applicable","reason":"nonconvergence"});
    mutations.push(b);
    let mut b = raw.clone();
    b["interval"]["width"] = json!(0.00001);
    mutations.push(b);
    let mut b = raw.clone();
    b["trials"][lo]["epoch"]["payload"]["continuation"]["progress"] = json!(0.5);
    mutations.push(b);
    let mut b = raw.clone();
    b["trials"][lo]["evidence"]["endpoint_audit"]["physical_force_kn"][0] = json!(123.);
    mutations.push(b);
    let mut b = raw.clone();
    b["trials"][lo]["evidence"]["factor"] = json!(2.);
    mutations.push(b);
    let mut b = raw.clone();
    b["continuum_stability_or_design_safety_certified"] = json!(true);
    mutations.push(b);
    let mut b = raw.clone();
    b["specification"]["search"]["absolute_factor_tolerance"] = json!(1.);
    mutations.push(b);
    for b in mutations {
        if let Ok(r) = serde_json::from_value::<epoch::Report>(b) {
            assert!(exact_fs::epoch(&m, &r.original_unreduced_run, &s, &r).is_err());
        }
    }
    Ok(())
}
#[test]
fn trial_limit_and_incomplete_paths_never_become_exact_or_an_upper_bound() -> Result<()> {
    for incomplete in [false, true] {
        let (m, g, c, mut s) = inputs()?;
        if incomplete {
            s.maximum_attempts = 1;
        } else {
            s.search.maximum_trials = 3;
        }
        let r = epoch::run(&m, stages::solve(&m, g, c)?, s.clone())?;
        assert!(!exact_fs::epoch(&m, &r.original_unreduced_run, &s, &r)?.exact_fs_for_this_case);
        let mut b = serde_json::to_value(&r)?;
        b["stop_reason"] = json!("factor_tolerance");
        assert!(
            exact_fs::epoch(
                &m,
                &r.original_unreduced_run,
                &s,
                &serde_json::from_value(b)?
            )
            .is_err()
        );
        let mut b = serde_json::to_value(&r)?;
        b["trials"][0]["evidence"]["outcome"] = json!("physically_excluded");
        assert!(
            exact_fs::epoch(
                &m,
                &r.original_unreduced_run,
                &s,
                &serde_json::from_value(b)?
            )
            .is_err()
        );
    }
    Ok(())
}
#[test]
fn no_physical_upper_may_continue_requirement_gate_but_does_not_certify_exact_fs() -> Result<()> {
    let (m, g, c, mut s) = inputs()?;
    s.search.maximum_factor = 1.2;
    let r = epoch::run(&m, stages::solve(&m, g, c)?, s.clone())?;
    assert_eq!(
        soil_fem_3d::epoch_srm_gate::effect(&r),
        soil_fem_3d::epoch_srm_gate::GateEffect::Continue
    );
    assert!(!exact_fs::epoch(&m, &r.original_unreduced_run, &s, &r)?.exact_fs_for_this_case);
    assert!(r.interval.is_none());
    Ok(())
}
#[test]
fn fresh_failed_trial_cannot_be_relabelled_as_physical_exclusion_or_equilibrium() -> Result<()> {
    let (m, g, c, s) = inputs()?;
    let r = fresh::run(&m, g.clone(), c.clone(), s.search.clone())?;
    for outcome in ["physically_excluded", "equilibrium_attained"] {
        let mut b = serde_json::to_value(&r)?;
        b["trials"][0]["run"] = Value::Null;
        b["trials"][0]["endpoint_audit"] = Value::Null;
        b["trials"][0]["failure_kind"] = json!("nonlinear");
        b["trials"][0]["outcome"] = json!(outcome);
        assert!(exact_fs::fresh(&m, &g, &c, &s.search, &serde_json::from_value(b)?).is_err());
    }
    Ok(())
}
