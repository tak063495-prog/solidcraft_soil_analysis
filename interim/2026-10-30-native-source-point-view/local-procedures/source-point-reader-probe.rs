//! Public API exercise of the real version4 packet; no solver or imported binary runs.
use serde_json::{Value,json};
use solidcraft_engine::Session;
use std::{fs,path::Path};

fn configure(session:&mut Session, settings:&mut Value, frame:usize, quantity:&str, section:Value, scale:f64)->Result<Value,Box<dyn std::error::Error>> {
    let object=settings.as_object_mut().ok_or("settings object required")?;
    object.insert("frame".into(),json!(frame));
    object.insert("quantity".into(),json!(quantity));
    object.insert("section".into(),section);
    object.insert("deformation_scale".into(),json!(scale));
    object.insert("displacement_reference".into(),json!("stage"));
    Ok(session.execute("analysis.soil.results.configure",settings)?)
}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let a=std::env::args().skip(1).collect::<Vec<_>>();
    let document=Path::new(a.first().ok_or("document required")?);
    let manifest=Path::new(a.get(1).ok_or("manifest required")?);
    let plan:Value=serde_json::from_slice(&fs::read(a.get(2).ok_or("query plan required")?)?)?;
    let output=Path::new(a.get(3).ok_or("new output required")?);
    if output.exists(){return Err("output already exists".into());}
    let mut session=Session::new(solidcraft_engine::doc::Document::from_json(&fs::read_to_string(document)?)?);
    let snapshot=session.execute("analysis.soil.snapshot",&json!({}))?;
    let opened=session.execute("analysis.soil.results.open",&json!({"manifest":manifest,"expected_input_sha256":snapshot.get("input_sha256")}))?;
    let mut settings=serde_json::to_value(session.soil_result_settings().ok_or("missing result settings")?)?;
    let mut frames=Vec::new();
    for request in plan.get("frames").and_then(Value::as_array).ok_or("frame plan required")? {
        let frame=request.get("index").and_then(Value::as_u64).ok_or("frame index required")? as usize;
        let configured=configure(&mut session,&mut settings,frame,"vertical_displacement",Value::Null,1.)?;
        let mut points=Vec::new();
        for point in request.get("point_ids").and_then(Value::as_array).ok_or("point plan required")? {
            points.push(json!({"id":point,"value":session.execute("analysis.soil.results.query",&json!({"integration_point":point}))?}));
        }
        let mut nodes=Vec::new();
        for node in request.get("node_ids").and_then(Value::as_array).ok_or("node plan required")? {
            nodes.push(json!({"id":node,"value":session.execute("analysis.soil.results.query",&json!({"node":node}))?}));
        }
        let section=configure(&mut session,&mut settings,frame,"primary_mean_effective_stress",json!({"normal":[1.,0.,0.],"offset_m":5.,"slice_only":true}),12.5)?;
        let point=request.get("point_ids").and_then(|v|v.get(0)).ok_or("point required")?;
        let after_settings=session.execute("analysis.soil.results.query",&json!({"integration_point":point}))?;
        let inactive_rejected=if let Some(id)=request.get("inactive_point").filter(|v|!v.is_null()) {
            let rejected=session.execute("analysis.soil.results.query",&json!({"integration_point":id})).is_err();
            if !rejected{return Err("inactive source point accepted".into());}
            Some(rejected)
        }else{None};
        let mut plates=Vec::new();
        for plate in request.get("plate_ids").and_then(Value::as_array).ok_or("plate plan required")? {
            plates.push(json!({"id":plate,"value":session.execute("analysis.soil.results.query",&json!({"plate":plate}))?}));
        }
        let wall=configure(&mut session,&mut settings,frame,"plate_bending_xx",Value::Null,1.)?;
        frames.push(json!({"frame":frame,"configured":configured,"points":points,"nodes":nodes,"section":section,
            "raw_after_settings":after_settings,"inactive_point_rejected":inactive_rejected,"plates":plates,"wall":wall}));
        println!("version4 authoritative-source frame {frame} reviewed");
    }
    serde_json::to_writer_pretty(fs::File::create_new(output)?,&json!({"opened":opened,"frames":frames,"actual_public_reader":true,"GUI_acceptance_separate":true}))?;
    Ok(())
}
