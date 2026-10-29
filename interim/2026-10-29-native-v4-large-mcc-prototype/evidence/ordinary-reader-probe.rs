//! Actual existing MCC result regression: no source rewriting or GUI claim.
use serde_json::{Value,json};
use solidcraft_engine::Session;
use std::{fs,path::Path};
fn main()->Result<(),Box<dyn std::error::Error>> {
    let a=std::env::args().skip(1).collect::<Vec<_>>();
    let document=Path::new(a.first().ok_or("document required")?);
    let manifest=Path::new(a.get(1).ok_or("manifest required")?);
    let output=Path::new(a.get(2).ok_or("new output required")?);
    if output.exists(){return Err("output already exists".into());}
    let mut session=Session::new(solidcraft_engine::doc::Document::from_json(&fs::read_to_string(document)?)?);
    let snapshot=session.execute("analysis.soil.snapshot",&json!({}))?;
    let opened=session.execute("analysis.soil.results.open",&json!({"manifest":manifest,"expected_input_sha256":snapshot.get("input_sha256")}))?;
    let mut settings=serde_json::to_value(session.soil_result_settings().ok_or("missing result settings")?)?;
    let mut frames=Vec::new();
    for frame in 0..5 {
        let object=settings.as_object_mut().ok_or("missing settings object")?;
        object.insert("frame".into(),json!(frame));
        object.insert("quantity".into(),json!("mean_effective_stress"));
        object.insert("displacement_reference".into(),json!("stage"));
        let configured=session.execute("analysis.soil.results.configure",&settings)?;
        let mut points=Vec::new();
        for point in [0usize,123,75788] {
            points.push(json!({"id":point,"value":session.execute("analysis.soil.results.query",&json!({"integration_point":point}))?}));
        }
        let mut nodes=Vec::new();
        for node in [0usize,1,100] {
            nodes.push(json!({"id":node,"value":session.execute("analysis.soil.results.query",&json!({"node":node}))?}));
        }
        frames.push(json!({"frame":frame,"configured":configured,"points":points,"nodes":nodes}));
        println!("ordinary MCC frame {frame} reviewed");
    }
    let result=json!({"opened":opened,"frames":frames,"actual_existing_MCC_source":true,"GUI_acceptance_separate":true});
    let file=fs::File::create_new(output)?;
    serde_json::to_writer_pretty(file,&result)?;
    Ok(())
}
