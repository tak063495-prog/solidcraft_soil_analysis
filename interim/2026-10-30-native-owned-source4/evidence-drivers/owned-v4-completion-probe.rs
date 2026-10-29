//! Public owned-process protocol regression; the test worker performs no FEM.
use serde_json::{Value, json};
use solidcraft_engine::Session;
use std::{fs, path::Path, thread, time::{Duration, Instant}};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let document = Path::new(args.first().ok_or("document required")?);
    let host = Path::new(args.get(1).ok_or("test host required")?);
    let output = Path::new(args.get(2).ok_or("new output folder required")?);
    let fixed = args.get(3).ok_or("expected fixed flag required")? == "fixed";
    if output.exists() { return Err("new output required".into()); }
    fs::create_dir(output)?;
    let mut session = Session::new(solidcraft_engine::doc::Document::from_json(&fs::read_to_string(document)?)?);
    session.set_soil_job_root(host)?;
    let snapshot = session.execute("analysis.soil.snapshot", &json!({}))?;
    let input = snapshot.get("input_sha256").ok_or("snapshot identity missing")?;
    if snapshot.pointer("/snapshot/stages/schema").and_then(Value::as_str) != Some("soil-fem-projected-mcc-stages/1") {
        return Err("original projected CAD snapshot required".into());
    }
    let mut rows = Vec::new();
    for (mode, expected) in [("legacy", true), ("source4", fixed), ("wrong_family", false),
        ("wrong_descriptor", false), ("wrong_storage", false), ("missing_semantics", false), ("future_view", false)] {
        fs::write(host.join("mode.json"), serde_json::to_vec(&json!({"mode":mode}))?)?;
        let directory = output.join(mode);
        session.execute("analysis.soil.job.start", &json!({"directory":directory,
            "expected_input_sha256":input,"mesh_size_m":1.,"element_order":2,"quadrature_order":3}))?;
        let began = Instant::now();
        let status = loop {
            let status = session.execute("analysis.soil.job.status", &json!({}))?;
            if status.get("process_alive") == Some(&Value::Bool(false)) { break status; }
            if began.elapsed() > Duration::from_secs(60) { return Err("test worker observation deadline".into()); }
            thread::sleep(Duration::from_millis(25));
        };
        if status.get("result_usable").and_then(Value::as_bool) != Some(expected) {
            return Err(format!("unexpected owned completion for {mode}: {status}").into());
        }
        if status.get("state").and_then(Value::as_str) != Some(if expected {"succeeded"} else {"failed"}) {
            return Err(format!("unexpected terminal state for {mode}: {status}").into());
        }
        rows.push(json!({"mode":mode,"expected_result_usable":expected,"observed_status":status}));
    }
    let receipt = json!({"passed":true,"original_projected_CAD_snapshot":input,
        "fixed_behavior_expected":fixed,"real_owned_child_handle_and_public_command_path":true,
        "hash_resealed_receipts_for_every_variant":true,"protocol_test_worker_does_not_compute_FEM":true,
        "actual_large_result_reader_GUI_and_physics_not_tested_here":true,"cases":rows,"p6_complete":false});
    fs::write(output.join("receipt.json"), serde_json::to_vec_pretty(&receipt)?)?;
    println!("PASS public owned-process completion protocol, fixed={fixed}");
    Ok(())
}
