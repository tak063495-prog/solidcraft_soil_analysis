//! Replay the bounded native route and retain the step proposals around the reload stall.
use serde_json::json;
use soil_fem_3d::{Result, ground::Ground, mesh::Mesh, stages};
use std::{env, fs, path::PathBuf};

fn main() -> Result<()> {
    let paths: Vec<PathBuf> = env::args_os().skip(1).map(PathBuf::from).collect();
    let mesh: Mesh = serde_json::from_slice(&fs::read(&paths[0])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&paths[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&paths[2])?)?;
    let limit = case.settings.attempt_budget_policy.stage_limit(&mesh)?;
    assert_eq!((mesh.tetrahedra10.len(), limit), (12, 2065));
    let run = stages::solve(&mesh, ground, case)?;
    let attempts: Vec<_> = run
        .attempts
        .iter()
        .filter(|a| a.stage == "reload")
        .collect();
    assert!(attempts.len() <= limit);
    let mut records = Vec::new();
    for a in &attempts {
        if a.progress_from < 926. / 1024. || a.progress_from >= 927. / 1024. {
            continue;
        }
        let e = a.joint_path_error.as_ref().expect("step-doubling audit");
        let target = e.friction_event_target.map(|(t, i)| {
            json!({"load": t, "point": i, "old": e.old_status[i], "mid": e.midpoint_status[i],
                   "coarse": e.coarse_status[i], "fine": e.fine_status[i]})
        });
        records.push(json!({"from": a.progress_from, "to": a.progress_to, "accepted": a.accepted,
                            "reason": a.reason, "event": target, "slip_errors_m": e.slip_errors_m, "coarse": e.coarse_status, "mid": e.midpoint_status, "fine": e.fine_status, "old": e.old_status, "mid_force": e.midpoint_relative_force_residual, "fine_force": e.fine_relative_force_residual, "proposed": e.proposed_increment,
                            "slip_ratio": e.maximum_slip_error_m/e.slip_budget_m,
                            "dissipation_ratio": e.maximum_dissipation_error_kpa_m/e.dissipation_budget_kpa_m,
                            "midpoint_linear_solves": e.midpoint_linear_solves.len(),
                            "fine_linear_solves": e.fine_linear_solves.len(),
                            "coarse_linear_solves": e.coarse_linear_solves.as_ref().map(Vec::len)}));
    }
    println!(
        "{}",
        json!({"attempts": attempts.len(), "complete": run.complete,
                          "last": run.last_committed.stage_progress, "records": records})
    );
    Ok(())
}
