//! Water bookkeeping policies for BIRTH and MATSET at fixed physical time.
use crate::{Result, richards};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum WaterPolicy {
    Inactive,
    InitializeHead {
        value: richards::AffineHead,
    },
    PreserveHead,
    PreserveWater {
        reference_head: richards::AffineHead,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChangedCell {
    pub element_index: usize,
    pub old_material_id: Option<u32>,
    pub new_material_id: u32,
    pub before: Option<richards::RemovedCell>,
    pub after: richards::RemovedCell,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WaterChange {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub internal_transfers: Vec<CellTransfer>,
    pub state: richards::State,
    pub old_faces: Vec<richards::Face>,
    pub new_faces: Vec<richards::Face>,
    pub new_to_old_face: Vec<Option<usize>>,
    pub changed_cells: Vec<ChangedCell>,
    pub water_before_m3: f64,
    pub water_after_m3: f64,
    pub water_added_m3: f64,
    pub water_removed_m3: f64,
    pub absolute_mass_transfer_error_m3: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellTransfer {
    pub donor_element: usize,
    pub receiver_element: usize,
    pub volume_m3: f64,
}
pub struct Domain<'a> {
    pub active: &'a [bool],
    pub material_ids: &'a [u32],
}
pub(crate) struct TransitionInput<'a> {
    pub before: Domain<'a>,
    pub after: Domain<'a>,
    pub selected: &'a [usize],
    pub policy: &'a WaterPolicy,
    pub internal_donors: Option<&'a std::collections::BTreeMap<usize, usize>>,
}
pub fn change(
    mesh: &crate::mesh::Mesh,
    case: &richards::Case,
    before: Domain<'_>,
    after: Domain<'_>,
    old: &richards::State,
    selected: &[usize],
    policy: &WaterPolicy,
) -> Result<WaterChange> {
    richards::transition_domain(
        mesh,
        case,
        old,
        TransitionInput {
            before,
            after,
            selected,
            policy,
            internal_donors: None,
        },
        None,
    )
}

/// Physical party on the other side of an instantaneous soil inventory change.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExternalWaterParty {
    pub reservoir_id: String,
}
