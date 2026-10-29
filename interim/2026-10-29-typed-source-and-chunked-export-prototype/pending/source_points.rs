//! Versioned references to raw projected-MCC points in the authoritative Run.
use crate::{Result, params::bad};
use serde::Deserialize;

const CMD: &str = "analysis.soil.results";
pub(super) const SCHEMA: &str = "soil-fem-view-source-points/1";
pub(super) const VIEW_SCHEMA: &str = "soil-fem-solidcraft-view/4";

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Descriptor {
    pub schema: String,
    pub source_result_schema: String,
    pub source_result_sha256: String,
    pub frame_index: usize,
    pub frame_name: String,
    pub integration_point_count: usize,
    pub active_integration_point_count: usize,
    pub length_unit: String,
    pub stress_unit: String,
    pub stress_sign: String,
    pub stress_components: [String; 6],
}

pub(super) struct Context<'a> {
    pub source_schema: &'a str,
    pub source_sha256: &'a str,
    pub frame_index: usize,
    pub frame_name: &'a str,
    pub point_count: usize,
    pub active_point_count: usize,
}

impl Descriptor {
    pub fn validate(&self, expected: &Context<'_>) -> Result<()> {
        let supported = matches!(expected.source_schema,
            "soil-fem-projected-mcc-stage-result/1" | "soil-fem-projected-mcc-stage-result/2");
        let hash = expected.source_sha256.len() == 64
            && expected.source_sha256.bytes().all(|v| v.is_ascii_digit() || (b'a'..=b'f').contains(&v));
        if self.schema != SCHEMA || !supported || !hash
            || self.source_result_schema != expected.source_schema || self.source_result_sha256 != expected.source_sha256
            || self.frame_index != expected.frame_index || self.frame_name != expected.frame_name || self.frame_name.is_empty()
            || expected.point_count == 0 || expected.point_count > 5_000_000
            || expected.active_point_count == 0 || expected.active_point_count > expected.point_count
            || self.integration_point_count != expected.point_count || self.active_integration_point_count != expected.active_point_count
            || self.length_unit != "m" || self.stress_unit != "kPa" || self.stress_sign != "compression_positive"
            || self.stress_components != ["xx", "yy", "zz", "xy", "yz", "xz"]
        {
            return Err(bad(CMD, "source-point descriptor differs from authoritative frame"));
        }
        Ok(())
    }
}
