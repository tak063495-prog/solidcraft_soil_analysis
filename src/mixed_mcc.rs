//! Primary mixed-MCC state. Mean effective stress is not pore pressure.
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PressureDof {
    pub node_index: usize,
    pub material_id: u32,
    pub compression_positive_kpa: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    #[serde(default = "legacy_tangent_method")]
    pub tangent_method: String,
    pub pressure_dofs: Vec<PressureDof>,
    /// Immutable GP indexing, including frozen removed points.
    pub cumulative_volume_strain_compression_positive: Vec<f64>,
    pub reference_void_ratio: Vec<f64>,
    pub last_increment_partitions: Vec<usize>,
    pub relative_volume_equation_residual: f64,
    pub relative_pressure_consistency_residual: f64,
    /// Present only for the explicitly identified P1-volume formulation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub projected_volume: Option<ProjectedVolume>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VolumeDof {
    pub node_index: usize,
    pub material_id: u32,
    pub compression_positive: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompatibilityOffset {
    pub node_index: usize,
    pub material_id: u32,
    pub integral_m3: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectedVolume {
    pub formulation: String,
    pub volume_dofs: Vec<VolumeDof>,
    /// Integral of N (div(u)+theta) frozen when the active domain changes.
    pub compatibility_offsets: Vec<CompatibilityOffset>,
    /// Diagnostic only: consistency is weak for this identified formulation.
    pub maximum_gp_pressure_projection_difference_kpa: f64,
}

pub(crate) const PROJECTED_FORMULATION: &str =
    "tet10_u_material_p1_p_p1_volume_frozen_compatibility";

fn legacy_tangent_method() -> String {
    "finite_difference_fixed_partition".into()
}
pub(crate) const TANGENT_METHOD: &str = "implicit_branch_fixed_partition";
