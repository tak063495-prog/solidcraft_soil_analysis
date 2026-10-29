//! Fresh public owned MCC workflow. This worker performs the real FEM and audits.
use serde_json::{Value, json};
use solidcraft_engine::Session;
use std::{fs, path::Path, thread, time::Duration};

fn save(path: &Path, value: &Value) -> Result<(), Box<dyn std::error::Error>> {
    let pending = path.with_extension("pending");
    fs::write(&pending, serde_json::to_vec_pretty(value)?)?;
    fs::rename(pending, path)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let document = Path::new(args.first().ok_or("document required")?);
    let host = Path::new(args.get(1).ok_or("real FEM host required")?);
    let output = Path::new(args.get(2).ok_or("new run directory required")?);
    let pause = Path::new(args.get(3).ok_or("pause marker required")?);
    if output.exists() || pause.exists() {
        return Err("new directory and active execution policy required".into());
    }
    fs::create_dir(output)?;
    let mut session = Session::new(solidcraft_engine::doc::Document::from_json(
        &fs::read_to_string(document)?,
    )?);
    session.set_soil_job_root(host)?;
    let snapshot = session.execute("analysis.soil.snapshot", &json!({}))?;
    if snapshot
        .pointer("/snapshot/stages/schema")
        .and_then(Value::as_str)
        != Some("soil-fem-projected-mcc-stages/1")
    {
        return Err("original projected MCC CAD input required".into());
    }
    let input = snapshot
        .get("input_sha256")
        .ok_or("snapshot identity missing")?;
    save(&output.join("snapshot.json"), &snapshot)?;
    let started = session.execute(
        "analysis.soil.job.start",
        &json!({"directory":output.join("owned-case"), "expected_input_sha256":input,
            "mesh_size_m":2.0, "element_order":2, "quadrature_order":3, "cooperative":true}),
    )?;
    save(&output.join("started.json"), &started)?;
    let mut cancellation_requested = false;
    let mut phase = Value::Null;
    let status = loop {
        let status = session.execute("analysis.soil.job.status", &json!({}))?;
        save(&output.join("latest-status.json"), &status)?;
        let next_phase = status
            .pointer("/worker_report/phase")
            .cloned()
            .unwrap_or(Value::Null);
        if next_phase != phase {
            println!("real owned workflow phase: {next_phase}");
            phase = next_phase;
        }
        if status.get("process_alive") == Some(&Value::Bool(false)) {
            break status;
        }
        if pause.exists() && !cancellation_requested {
            let token = status
                .get("job")
                .and_then(Value::as_str)
                .ok_or("owned job token missing")?;
            let cancelled = session.execute("analysis.soil.job.cancel", &json!({"job":token}))?;
            save(&output.join("cooperative-stop-request.json"), &cancelled)?;
            cancellation_requested = true;
        }
        thread::sleep(Duration::from_secs(2));
    };
    save(&output.join("terminal-status.json"), &status)?;
    if matches!(
        status.get("state").and_then(Value::as_str),
        Some("checkpointed" | "cancelled")
    ) {
        save(
            &output.join("stopped.json"),
            &json!({"passed":false,
            "cooperative_stop":true,"status":status,"full_solve_not_accepted":true}),
        )?;
        std::process::exit(3);
    }
    if status.get("state").and_then(Value::as_str) != Some("succeeded")
        || status.get("result_usable") != Some(&Value::Bool(true))
    {
        return Err(format!("real owned solve not accepted: {status}").into());
    }
    let manifest = status
        .get("result")
        .and_then(Value::as_str)
        .ok_or("result manifest missing")?;
    let opened = session.execute(
        "analysis.soil.results.open",
        &json!({
        "manifest":manifest,"expected_input_sha256":input}),
    )?;
    save(&output.join("opened.json"), &opened)?;
    let metadata: Value = serde_json::from_slice(&fs::read(manifest)?)?;
    if metadata.get("schema").and_then(Value::as_str) != Some("soil-fem-solidcraft-view/4")
        || metadata
            .get("frames")
            .and_then(Value::as_array)
            .map(Vec::len)
            != Some(5)
    {
        return Err("five-stage source-point view required".into());
    }
    let mut settings = serde_json::to_value(
        session
            .soil_result_settings()
            .ok_or("result settings missing")?,
    )?;
    let mut frames = Vec::new();
    for index in 0..5 {
        let object = settings.as_object_mut().ok_or("settings object required")?;
        object.insert("frame".into(), json!(index));
        object.insert("quantity".into(), json!("vertical_displacement"));
        object.insert("displacement_reference".into(), json!("stage"));
        let configured = session.execute("analysis.soil.results.configure", &settings)?;
        let node = session.execute("analysis.soil.results.query", &json!({"node":0}))?;
        frames.push(json!({"frame":index,"configured":configured,"node0":node}));
    }
    save(
        &output.join("receipt.json"),
        &json!({"passed":true,
        "fresh_real_owned_MCC_solve_and_audits":true,"public_start_status_and_open_configure_query":true,
        "original_projected_CAD_input_sha256":input,"view4_all5_frames":true,"frames":frames,
        "small_workflow_discretization_only":true,"actual_GUI_and_large14909_and_spatial_convergence_separate":true,
        "p6_complete":false}),
    )?;
    println!("PASS fresh real owned MCC source-point workflow and five public result frames");
    Ok(())
}
