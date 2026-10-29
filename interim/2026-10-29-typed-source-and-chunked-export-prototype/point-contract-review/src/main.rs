#![forbid(unsafe_code)]
#![allow(dead_code)]
type Result<T> = std::result::Result<T, String>;
mod params {
    pub fn bad(_: &str, message: impl ToString) -> String { message.to_string() }
}
#[path = "../../pending/source_points.rs"]
mod source_points;
fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value,json};
    fn fixture() -> Value { serde_json::from_str(include_str!("../python-descriptor.json")).unwrap() }
    fn valid(value: Value) -> bool {
        let Ok(value) = serde_json::from_value::<source_points::Descriptor>(value) else { return false };
        value.validate(&source_points::Context { source_schema:"soil-fem-projected-mcc-stage-result/2",source_sha256:&"a".repeat(64),
            frame_index:0,frame_name:"excavate",point_count:3,active_point_count:2 }).is_ok()
    }
    #[test]
    fn actual_python_descriptor_decodes_and_matches_rust_context() { assert!(valid(fixture())); }
    #[test]
    fn bindings_counts_units_sign_component_order_and_unknown_fields_rejected() {
        for (key,value) in [("source_result_sha256",json!("b".repeat(64))),("frame_name",json!("preload")),
            ("frame_index",json!(1)),("frame_index",json!(true)),("integration_point_count",json!(2)),
            ("active_integration_point_count",json!(3)),("length_unit",json!("mm")),("stress_unit",json!("Pa")),
            ("stress_sign",json!("tension_positive")),("stress_components",json!(["zz","yy","xx","xy","yz","xz"])),
            ("schema",json!("soil-fem-view-integration-points/3")),("unknown",json!(1))] {
            let mut changed=fixture();changed[key]=value;assert!(!valid(changed),"{key}");
        }
    }
    #[test]
    fn duplicate_descriptor_keys_are_rejected() {
        let raw=serde_json::to_string(&fixture()).unwrap().replacen("\"frame_index\":0","\"frame_index\":0,\"frame_index\":0",1);
        assert!(serde_json::from_str::<source_points::Descriptor>(&raw).is_err());
    }
}
