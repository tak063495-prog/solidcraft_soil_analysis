//! General-geometry constitutive seeds and geostatic continuation records.
//! Equilibration uses the same global material-force Newton as construction.
use crate::{
    Error, Result,
    boundary::Loads,
    checkpoint,
    ground::Ground,
    mesh::Mesh,
    soil,
    stages::{Attempt, Committed},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Body {
    pub body: String,
    pub material_id: u32,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointSeed {
    pub element_index: usize,
    pub local_index: usize,
    pub material_id: u32,
    pub position_m: [f64; 3],
    pub state: soil::State,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Seed {
    BodyStates {
        states: BTreeMap<String, soil::State>,
    },
    IntegrationPoints {
        mesh_sha256: String,
        points: Vec<PointSeed>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    pub seed: Seed,
    pub initial_boundary_factor: f64,
    pub maximum_initial_engineering_strain: f64,
    pub relative_volume_tolerance: f64,
    #[serde(default, skip_serializing_if = "Loads::is_empty")]
    pub base_loads: Loads,
    /// Initial exterior faces exposed to the prescribed hydrostatic water body.
    /// Newly excavated faces retain the existing flooded-cut convention.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub hydrostatic_face_sets: Vec<String>,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub verification_reject_once: bool,
}
impl Configuration {
    pub fn validate(&self) -> Result<()> {
        if !self.initial_boundary_factor.is_finite()
            || !self.maximum_initial_engineering_strain.is_finite()
            || !(1e-8..=0.1).contains(&self.maximum_initial_engineering_strain)
            || !self.relative_volume_tolerance.is_finite()
            || !(1e-10..=0.05).contains(&self.relative_volume_tolerance)
            || self.hydrostatic_face_sets.iter().any(|n| n.is_empty())
            || self
                .hydrostatic_face_sets
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != self.hydrostatic_face_sets.len()
        {
            return Err(Error::Input(
                "invalid geostatic factor/strain/volume/water-face settings".into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub complete: bool,
    pub progress: f64,
    pub next_increment: f64,
    pub verification_rejection_done: bool,
    pub seed: Committed,
    pub last_committed: Committed,
    pub initial_holding_load_kn: Vec<f64>,
    pub attempts: Vec<Attempt>,
    pub maximum_absolute_engineering_strain_from_seed: f64,
    pub final_relative_physical_force_residual: Option<f64>,
    pub failure_reason: Option<String>,
}

pub fn states(mesh: &Mesh, ground: &Ground) -> Result<Vec<soil::State>> {
    mesh.validate()?;
    ground.validate()?;
    let cfg = ground
        .initialization
        .geostatic
        .as_ref()
        .ok_or_else(|| Error::Input("missing geostatic settings".into()))?;
    if mesh.regions.len() != ground.bodies.len()
        || !["soil-fem-mesh/2", "soil-fem-mesh/3"].contains(&mesh.schema.as_str())
    {
        return Err(Error::Input(
            "general initialization requires a body-resolved mesh and all body mappings".into(),
        ));
    }
    let bodies: BTreeMap<_, _> = ground
        .bodies
        .iter()
        .map(|b| (b.body.as_str(), b.material_id))
        .collect();
    let material: BTreeMap<_, _> = ground.materials.iter().map(|m| (m.id, &m.model)).collect();
    let regions: BTreeMap<_, _> = mesh.regions.iter().map(|r| (r.id, r)).collect();
    let tools: BTreeSet<_> = ground.excavation_bodies.iter().collect();
    if tools.len() != ground.excavation_bodies.len()
        || tools != mesh.excavation_sets.keys().collect()
    {
        return Err(Error::Input(
            "geostatic excavation-tool mapping mismatch".into(),
        ));
    }
    for r in &mesh.regions {
        if bodies.get(r.body.as_str()) != Some(&r.material_id) {
            return Err(Error::Input(
                "geostatic CAD body/material mapping mismatch".into(),
            ));
        }
    }
    cfg.base_loads.conditions().validate_definition(mesh)?;
    cfg.base_loads
        .conditions()
        .loads(mesh, &vec![true; mesh.tetrahedra.len()])?;
    let owners = crate::boundary::surface_owners(mesh);
    for name in &cfg.hydrostatic_face_sets {
        if mesh.face_sets.get(name).is_none_or(Vec::is_empty) {
            return Err(Error::Input(format!(
                "unknown/empty hydrostatic face set {name}"
            )));
        }
        if mesh.face_sets[name].iter().any(|f| {
            let mut f = *f;
            f.sort_unstable();
            owners.get(&f).is_none_or(|v| v.len() != 1)
        }) {
            return Err(Error::Input(
                "initial hydrostatic faces must be on the initial exterior".into(),
            ));
        }
    }
    let layout = mesh.integration_layout()?;
    let mut volumes = BTreeMap::<u32, f64>::new();
    for q in &layout {
        *volumes.entry(mesh.region_ids[q.element_index]).or_default() += q.weight_m3;
    }
    for r in &mesh.regions {
        if (volumes[&r.id] - r.volume_m3).abs()
            > 1e-10 + cfg.relative_volume_tolerance * r.volume_m3
        {
            return Err(Error::Input(format!(
                "geostatic mesh/CAD volume difference exceeds configured tolerance in {}",
                r.body
            )));
        }
    }
    let out = match &cfg.seed {
        Seed::BodyStates { states } => {
            if states.keys().map(String::as_str).collect::<BTreeSet<_>>()
                != bodies.keys().copied().collect()
            {
                return Err(Error::Input(
                    "seed must specify every body exactly once".into(),
                ));
            }
            layout
                .iter()
                .map(|q| states[&regions[&mesh.region_ids[q.element_index]].body].clone())
                .collect::<Vec<_>>()
        }
        Seed::IntegrationPoints {
            mesh_sha256,
            points,
        } => {
            if *mesh_sha256 != checkpoint::hash_json(mesh)? || points.len() != layout.len() {
                return Err(Error::Input(
                    "initial-state import mesh identity/point count mismatch".into(),
                ));
            }
            for (p, q) in points.iter().zip(&layout) {
                if p.element_index != q.element_index
                    || p.local_index != q.local_index
                    || p.material_id != mesh.material_ids[q.element_index]
                    || p.position_m
                        .iter()
                        .zip(q.position_m)
                        .any(|(a, b)| !a.is_finite() || (a - b).abs() > 1e-8)
                {
                    return Err(Error::Input(
                        "initial-state import point position/index/material mismatch".into(),
                    ));
                }
            }
            points.iter().map(|p| p.state.clone()).collect()
        }
    };
    for (state, q) in out.iter().zip(&layout) {
        material[&mesh.material_ids[q.element_index]].validate_state(state)?;
    }
    Ok(out)
}
