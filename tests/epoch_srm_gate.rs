use serde_json::json;
use soil_fem_3d::{
    Result,
    epoch_srm_gate::{self, Decision, GateEffect, Specification},
    epoch_srm_search::Specification as EpochSearch,
    ground::Ground,
    mesh::Mesh,
    srm_bound::Problem,
    srm_search::{self, Requirement, StopReason},
    stage_control::{Observer, Progress},
    stages::{self, Run, StageCase},
};

fn column(
    stages: serde_json::Value,
) -> Result<(Mesh, Ground, StageCase, srm_search::Specification)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/srm-mean-stress-mesh.json"))?;
    let mesh = mesh.elevate_straight()?;
    let problem: Problem =
        serde_json::from_str(include_str!("../examples/srm-mean-stress-problem.json"))?;
    let mut raw: serde_json::Value =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let seed = raw["initialization"]["geostatic"]["seed"]["states"]["LowerClay"].clone();
    raw["initialization"]["geostatic"]["seed"]["states"] = json!({"TestCube": seed});
    raw["initialization"]["geostatic"]["verification_reject_once"] = json!(false);
    raw["materials"] = json!([{
        "id": 1,
        "name": "unconfined epoch gate",
        "unit_weight_kn_m3": 1e-10,
        "mohr_coulomb": problem.material
    }]);
    raw["bodies"] = json!([{"body": "TestCube", "material_id": 1}]);
    let ground = serde_json::from_value(raw)?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.boundary = Some(problem.boundary);
    case.stages = serde_json::from_value(stages)?;
    case.verification_reject_once_stages.clear();
    case.monitors.clear();
    case.monitor_sets.clear();
    let mut search: srm_search::Specification =
        serde_json::from_str(include_str!("../examples/srm-search.json"))?;
    search.through_stages = 0;
    Ok((mesh, ground, case, search))
}

fn two_stages() -> serde_json::Value {
    json!([
        {"kind": "boundary", "name": "load", "target_factor": 1.0, "initial_increment": 1.0},
        {"kind": "reset_display", "name": "reference"}
    ])
}

fn gate(
    mut search: srm_search::Specification,
    gates: Vec<usize>,
    attempts: usize,
    trials: usize,
) -> Specification {
    search.maximum_trials = trials;
    Specification {
        schema: "soil-fem-epoch-srm-stage-gate/1".into(),
        search: EpochSearch {
            schema: "soil-fem-epoch-srm-search/1".into(),
            search,
            initial_increment: 0.25,
            minimum_increment: 1. / 4096.,
            maximum_attempts: attempts,
        },
        gate_after_stages: gates,
    }
}

fn bytes<T: serde::Serialize>(value: &T) -> Result<Vec<u8>> {
    Ok(serde_json::to_vec(value)?)
}
fn scientific(value: &(impl serde::Serialize + ?Sized)) -> Result<serde_json::Value> {
    let mut value = serde_json::to_value(value)?;
    fn clean(value: &mut serde_json::Value) {
        match value {
            serde_json::Value::Object(map) => {
                // Linear-solve timings and hashes of the whole Run change between
                // executions. The timed payload itself is still compared.
                map.remove("elapsed_seconds");
                map.remove("origin_run_sha256");
                map.remove("payload_sha256");
                for child in map.values_mut() {
                    clean(child);
                }
            }
            serde_json::Value::Array(items) => {
                for child in items {
                    clean(child);
                }
            }
            _ => {}
        }
    }
    clean(&mut value);
    Ok(value)
}
fn assert_same_history(
    left: &(impl serde::Serialize + ?Sized),
    right: &(impl serde::Serialize + ?Sized),
) -> Result<()> {
    let left = scientific(left)?;
    let right = scientific(right)?;
    if left != right {
        panic!(
            "{}",
            first_diff(&left, &right, "").unwrap_or_else(|| "equal".into())
        );
    }
    Ok(())
}
fn first_diff(left: &serde_json::Value, right: &serde_json::Value, path: &str) -> Option<String> {
    use serde_json::Value;
    match (left, right) {
        (Value::Object(a), Value::Object(b)) => {
            let mut keys: Vec<_> = a.keys().chain(b.keys()).collect();
            keys.sort();
            keys.dedup();
            for key in keys {
                match (a.get(key), b.get(key)) {
                    (None, Some(_)) => return Some(format!("{path}/{key} missing on left")),
                    (Some(_), None) => return Some(format!("{path}/{key} missing on right")),
                    (Some(x), Some(y)) if x != y => {
                        return first_diff(x, y, &format!("{path}/{key}"));
                    }
                    _ => {}
                }
            }
            None
        }
        (Value::Array(a), Value::Array(b)) if a.len() != b.len() => {
            Some(format!("{path} array length {} != {}", a.len(), b.len()))
        }
        (Value::Array(a), Value::Array(b)) => {
            a.iter().zip(b).enumerate().find_map(|(i, (x, y))| {
                (x != y)
                    .then(|| first_diff(x, y, &format!("{path}/{i}")))
                    .flatten()
            })
        }
        _ => Some(format!(
            "{path}: {} != {}",
            &left.to_string()[..left.to_string().len().min(180)],
            &right.to_string()[..right.to_string().len().min(180)]
        )),
    }
}
fn reload(snapshot: &epoch_srm_gate::Snapshot) -> Result<epoch_srm_gate::Snapshot> {
    Ok(serde_json::from_slice(&serde_json::to_vec(snapshot)?)?)
}
fn drive(
    mesh: &Mesh,
    ground: &Ground,
    case: &StageCase,
    spec: &Specification,
) -> Result<epoch_srm_gate::Report> {
    let mut snapshot = epoch_srm_gate::Snapshot::new(mesh, ground, case, spec)?;
    for _ in 0..100_000 {
        snapshot.verify(mesh, ground, case, spec)?;
        if snapshot.is_complete() {
            return Ok(serde_json::from_slice(&bytes(snapshot.report().unwrap())?)?);
        }
        snapshot = reload(&snapshot)?.advance(mesh, ground, case, spec, 1, None)?;
    }
    Err(soil_fem_3d::Error::Control(
        "gate unit drive did not finish".into(),
    ))
}

fn same_prefix(prefix: &Run, full: &Run) -> Result<()> {
    let count = prefix.frames.len();
    assert!(count > 0 && count < full.frames.len());
    let mut earlier = scientific(&full.frames[..count])?;
    let mut prefix_frames = scientific(&prefix.frames)?;
    // A later ResetDisplay makes the full case retain material IDs from the
    // first increment. The truncated prefix does not. Compare the mechanics.
    fn drop_ids(value: &mut serde_json::Value) {
        match value {
            serde_json::Value::Object(map) => {
                map.remove("element_material_ids");
                for child in map.values_mut() {
                    drop_ids(child);
                }
            }
            serde_json::Value::Array(items) => {
                for child in items {
                    drop_ids(child);
                }
            }
            _ => {}
        }
    }
    drop_ids(&mut earlier);
    drop_ids(&mut prefix_frames);
    if prefix_frames != earlier {
        panic!(
            "{}",
            first_diff(&prefix_frames, &earlier, "frames").unwrap_or_else(|| "equal".into())
        );
    }
    assert_same_history(&prefix.attempts, &full.attempts[..prefix.attempts.len()])?;
    let mut committed = scientific(&prefix.last_committed)?;
    let mut frame_state = scientific(&full.frames[count - 1].state)?;
    drop_ids(&mut committed);
    drop_ids(&mut frame_state);
    if committed != frame_state {
        panic!(
            "{}",
            first_diff(&committed, &frame_state, "committed").unwrap_or_else(|| "equal".into())
        );
    }
    let mut prefix_geo = scientific(&prefix.geostatic)?;
    let mut full_geo = scientific(&full.geostatic)?;
    drop_ids(&mut prefix_geo);
    drop_ids(&mut full_geo);
    if prefix_geo != full_geo {
        panic!(
            "{}",
            first_diff(&prefix_geo, &full_geo, "geostatic").unwrap_or_else(|| "equal".into())
        );
    }
    assert_same_history(&prefix.integration_points, &full.integration_points)?;
    Ok(())
}

#[test]
fn stepped_full_case_keeps_material_ids() -> Result<()> {
    let (mesh, ground, case, _) = column(two_stages())?;
    let direct = stages::solve(&mesh, ground.clone(), case.clone())?;
    let mut run = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let mut steps = 1;
    while !run.complete && steps < 100 {
        let saved = soil_fem_3d::checkpoint::Checkpoint::seal(&mesh, run)?;
        run = saved.resume(&mesh, ground.clone(), case.clone(), Some(1))?;
        steps += 1;
    }
    assert!(run.complete, "{:?}", run.failure_reason);
    assert_eq!(
        direct.last_committed.element_material_ids.len(),
        run.last_committed.element_material_ids.len(),
        "steps {steps}"
    );
    assert_same_history(&direct, &run)?;
    Ok(())
}

#[test]
fn fresh_history_srm_bracket_remains_four_thirds() -> Result<()> {
    let (mesh, ground, case, mut search) = column(json!([
        {"kind": "boundary", "name": "load", "target_factor": 1.0, "initial_increment": 1.0}
    ]))?;
    search.through_stages = 1;
    let result = srm_search::run(&mesh, ground, case, search)?;
    let interval = result.interval.as_ref().unwrap();
    let critical = 20. / 15.;
    assert!(interval.attained_equilibrium_factor < critical);
    assert!(interval.independently_excluded_factor > critical);
    assert!(interval.width <= result.specification.absolute_factor_tolerance);
    assert_eq!(result.requirement, Requirement::EquilibriumRequirementMet);
    assert_eq!(
        result.replay_origin,
        "original_inputs_before_initial_equilibration_for_every_factor"
    );
    Ok(())
}

#[test]
fn passed_gates_match_the_continuous_run_after_every_accepted_unit() -> Result<()> {
    let (mesh, ground, case, search) = column(two_stages())?;
    let spec = gate(search, vec![1, 2], 128, 20);
    let direct = epoch_srm_gate::run(&mesh, &ground, &case, &spec)?;
    let stepped = drive(&mesh, &ground, &case, &spec)?;
    let left = scientific(&direct)?;
    let right = scientific(&stepped)?;
    if left != right {
        panic!(
            "diff {}",
            first_diff(&left, &right, "ROOT").unwrap_or_else(|| "equal".into())
        );
    }
    assert_eq!(direct.decision, Decision::Passed);
    assert!(
        direct
            .gates
            .iter()
            .all(|gate| gate.effect == GateEffect::Continue)
    );
    assert_eq!(
        direct.gates[0].requirement,
        Requirement::EquilibriumRequirementMet
    );
    assert_eq!(direct.gates[0].stop_reason, StopReason::FactorTolerance);
    assert!(!direct.stopped_before_later_stages);
    assert!(!direct.nonconvergence_used_as_exclusion_or_collapse);
    assert!(!direct.continuum_stability_or_design_safety_certified);
    let continuous = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert_eq!(
        direct
            .construction
            .last_committed
            .element_material_ids
            .len(),
        continuous.last_committed.element_material_ids.len(),
        "stages {} vs {} decision {:?}",
        direct.construction.case.stages.len(),
        continuous.case.stages.len(),
        direct.decision
    );
    assert_same_history(&direct.construction, &continuous)?;
    same_prefix(&direct.gates[0].search.original_unreduced_run, &continuous)?;
    assert_eq!(direct.construction.case.stages.len(), 2);
    assert!(direct.construction.complete && direct.construction.failure_reason.is_none());
    Ok(())
}

#[test]
fn independent_exclusion_stops_before_the_later_stage_and_keeps_the_unreduced_prefix() -> Result<()>
{
    let (mesh, ground, case, mut search) = column(two_stages())?;
    search.required_factor = 1.6;
    let spec = gate(search, vec![1, 2], 128, 20);
    let report = epoch_srm_gate::run(&mesh, &ground, &case, &spec)?;
    assert_same_history(&report, &drive(&mesh, &ground, &case, &spec)?)?;
    assert_eq!(report.decision, Decision::StoppedExcluded);
    assert_eq!(report.gates.len(), 1);
    assert_eq!(report.gates[0].effect, GateEffect::StopExcluded);
    assert_eq!(
        report.gates[0].requirement,
        Requirement::IndependentlyExcluded
    );
    assert_eq!(report.construction.case.stages.len(), 1);
    assert_eq!(
        serde_json::to_value(&report.construction.case.stages[0])?["name"],
        "load"
    );
    assert!(report.construction.complete && report.construction.failure_reason.is_none());
    assert!(report.stopped_before_later_stages);
    assert!(!report.nonconvergence_used_as_exclusion_or_collapse);
    let continuous = stages::solve(&mesh, ground, case)?;
    same_prefix(&report.construction, &continuous)?;
    assert_eq!(continuous.case.stages.len(), 2);
    Ok(())
}

#[test]
fn nonconvergence_is_indeterminate_and_does_not_exclude_or_continue() -> Result<()> {
    let (mesh, ground, case, search) = column(two_stages())?;
    let spec = gate(search, vec![1, 2], 1, 20);
    let report = epoch_srm_gate::run(&mesh, &ground, &case, &spec)?;
    assert_eq!(report.decision, Decision::StoppedIndeterminate);
    assert_eq!(report.gates.len(), 1);
    assert_eq!(report.gates[0].effect, GateEffect::StopIndeterminate);
    assert_eq!(report.gates[0].requirement, Requirement::Indeterminate);
    assert_eq!(report.gates[0].stop_reason, StopReason::IndeterminateTrial);
    assert_ne!(report.decision, Decision::StoppedExcluded);
    assert_eq!(report.construction.case.stages.len(), 1);
    assert!(report.construction.complete && report.construction.failure_reason.is_none());
    assert!(report.stopped_before_later_stages);
    assert!(!report.nonconvergence_used_as_exclusion_or_collapse);
    assert!(!report.continuum_stability_or_design_safety_certified);
    Ok(())
}

#[test]
fn trial_limit_does_not_authorize_the_next_stage() -> Result<()> {
    let (mesh, ground, case, search) = column(two_stages())?;
    let spec = gate(search, vec![1, 2], 128, 3);
    let report = epoch_srm_gate::run(&mesh, &ground, &case, &spec)?;
    assert_eq!(
        report.gates[0].requirement,
        Requirement::EquilibriumRequirementMet
    );
    assert_eq!(report.gates[0].stop_reason, StopReason::TrialLimit);
    assert_eq!(report.gates[0].effect, GateEffect::StopIndeterminate);
    assert_eq!(report.decision, Decision::StoppedIndeterminate);
    assert_eq!(report.construction.case.stages.len(), 1);
    assert!(report.stopped_before_later_stages);
    Ok(())
}

struct AlwaysStop;
impl Observer for AlwaysStop {
    fn on_progress(&self, _: &Progress) -> Result<()> {
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        Ok(true)
    }
}

#[test]
fn stop_at_an_accepted_unit_keeps_that_snapshot_and_resumes_to_the_direct_report() -> Result<()> {
    let (mesh, ground, case, search) = column(two_stages())?;
    let spec = gate(search, vec![1, 2], 128, 20);
    let expected = epoch_srm_gate::run(&mesh, &ground, &case, &spec)?;
    let mut snapshot = epoch_srm_gate::Snapshot::new(&mesh, &ground, &case, &spec)?;
    snapshot = reload(&snapshot)?.advance(&mesh, &ground, &case, &spec, 1, None)?;
    let saved = bytes(&snapshot)?;
    let stopped = reload(&snapshot)?.advance(&mesh, &ground, &case, &spec, 8, Some(&AlwaysStop))?;
    assert_eq!(saved, bytes(&stopped)?);
    assert!(!stopped.is_complete());
    let mut snapshot = reload(&stopped)?;
    for _ in 0..100_000 {
        if snapshot.is_complete() {
            break;
        }
        snapshot = reload(&snapshot)?.advance(&mesh, &ground, &case, &spec, 1, None)?;
    }
    assert_same_history(snapshot.report().unwrap(), &expected)?;
    let mut changed = spec.clone();
    changed.gate_after_stages.push(3);
    assert!(stopped.verify(&mesh, &ground, &case, &changed).is_err());
    Ok(())
}

#[test]
fn malformed_stage_list_and_fixed_through_stages_are_refused_before_solving() -> Result<()> {
    let (mesh, ground, case, search) = column(two_stages())?;
    let mut spec = gate(search, vec![2, 1], 128, 20);
    assert!(epoch_srm_gate::Snapshot::new(&mesh, &ground, &case, &spec).is_err());
    spec.gate_after_stages = vec![1, 2];
    spec.search.search.through_stages = 1;
    assert!(epoch_srm_gate::Snapshot::new(&mesh, &ground, &case, &spec).is_err());
    spec.search.search.through_stages = 0;
    spec.gate_after_stages.clear();
    assert!(epoch_srm_gate::run(&mesh, &ground, &case, &spec).is_err());
    Ok(())
}
