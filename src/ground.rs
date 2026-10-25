//! Admissible K0 initialization for conforming horizontal cuboid layers only.
//! Computes the FE residual of the prescribed state; it does not silently relax it.
use crate::{Error, Result, hardening, mc, mcc, mesh::Mesh, soil};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(try_from = "RawMaterial", into = "RawMaterial")]
pub struct Material {
    pub id: u32,
    pub name: String,
    pub unit_weight_kn_m3: f64,
    pub model: soil::Model,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawMaterial {
    id: u32,
    name: String,
    unit_weight_kn_m3: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    mcc: Option<mcc::Parameters>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    mohr_coulomb: Option<mc::Parameters>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    invariant_double_hardening: Option<hardening::Parameters>,
}
impl TryFrom<RawMaterial> for Material {
    type Error = String;
    fn try_from(r: RawMaterial) -> std::result::Result<Self, Self::Error> {
        let model =
            match (r.mcc, r.mohr_coulomb, r.invariant_double_hardening) {
                (Some(p), None, None) => soil::Model::Mcc(p),
                (None, Some(p), None) => soil::Model::MohrCoulomb(p),
                (None, None, Some(p)) => soil::Model::InvariantDoubleHardening(p),
                _ => return Err(
                    "exactly one of mcc, mohr_coulomb and invariant_double_hardening is required"
                        .into(),
                ),
            };
        Ok(Self {
            id: r.id,
            name: r.name,
            unit_weight_kn_m3: r.unit_weight_kn_m3,
            model,
        })
    }
}
impl From<Material> for RawMaterial {
    fn from(m: Material) -> Self {
        let (mcc, mohr_coulomb, invariant_double_hardening) = match m.model {
            soil::Model::Mcc(p) => (Some(p), None, None),
            soil::Model::MohrCoulomb(p) => (None, Some(p), None),
            soil::Model::InvariantDoubleHardening(p) => (None, None, Some(p)),
        };
        Self {
            id: m.id,
            name: m.name,
            unit_weight_kn_m3: m.unit_weight_kn_m3,
            mcc,
            mohr_coulomb,
            invariant_double_hardening,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "definition", rename_all = "snake_case", deny_unknown_fields)]
pub enum Consolidation {
    /// Required explicit acknowledgement for an MC material with no cap/OCR.
    NotApplicable,
    /// Isotropic MCC pressure ratio, not the oedometer vertical OCR.
    MccPcOverMeanEffective {
        ratio: f64,
    },
    /// Vertical OCR alone does not determine MCC pc: require its explicit calibration.
    VerticalOcr {
        ratio: f64,
        preconsolidation_pressure_kpa: f64,
    },
    ExplicitPc {
        preconsolidation_pressure_kpa: f64,
    },
    HardeningCapOverMeanEffective {
        ratio: f64,
    },
    HardeningExplicitCap {
        cap_pressure_kpa: f64,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub enum VoidRatio {
    FromMccReference,
    Explicit { value: f64 },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Layer {
    pub body: String,
    pub material_id: u32,
    pub k0_x: f64,
    pub k0_y: f64,
    pub consolidation: Consolidation,
    pub void_ratio: VoidRatio,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Initialization {
    pub method: String,
    pub surface_surcharge_kpa: f64,
    pub water_level_m: Option<f64>,
    pub water_unit_weight_kn_m3: f64,
    pub relative_equilibrium_tolerance: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub geostatic: Option<crate::geostatic::Configuration>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Ground {
    pub schema: String,
    pub units: String,
    pub initialization: Initialization,
    pub materials: Vec<Material>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub construction_material_ids: Vec<u32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub layers: Vec<Layer>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bodies: Vec<crate::geostatic::Body>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub excavation_bodies: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Point {
    pub element_index: usize,
    pub local_index: usize,
    pub weight_m3: f64,
    pub natural_coordinates: [f64; 3],
    pub region_id: u32,
    pub material_id: u32,
    pub centroid_m: [f64; 3],
    pub total_stress_compression_positive_kpa: [f64; 6],
    pub pore_pressure_kpa: f64,
    #[serde(flatten)]
    pub constitutive_state: soil::InitialPointState,
    pub porosity: f64,
    pub mcc_pc_over_mean_effective: Option<f64>,
    pub vertical_ocr: Option<f64>,
    pub past_maximum_vertical_effective_stress_kpa: Option<f64>,
    pub normalized_yield: f64,
    pub bulk_modulus_kpa: f64,
    pub shear_modulus_kpa: f64,
}

#[derive(Debug, Serialize)]
pub struct InitialState {
    pub schema: &'static str,
    pub mesh_schema: String,
    pub ground: Ground,
    pub point_states: Vec<Point>,
    pub body_volumes_m3: BTreeMap<String, f64>,
    pub shared_interface_triangles: usize,
    pub true_relative_free_force_residual: f64,
    pub applied_force_kn: [f64; 3],
    pub support_reaction_kn: [f64; 3],
    pub equilibrium_error_kn: [f64; 3],
    pub initialized_displacement_m: Vec<[f64; 3]>,
}

impl Ground {
    pub fn is_general(&self) -> bool {
        self.initialization.method == "equilibrate"
    }
    pub fn validate(&self) -> Result<()> {
        let init = &self.initialization;
        if self.schema != "soil-fem-ground/1"
            || self.units != "m-kN-kPa"
            || !["horizontal_k0", "equilibrate"].contains(&init.method.as_str())
            || !init.surface_surcharge_kpa.is_finite()
            || init.surface_surcharge_kpa < 0.
            || !init.water_unit_weight_kn_m3.is_finite()
            || init.water_unit_weight_kn_m3 <= 0.
            || init.water_level_m.is_some_and(|v| !v.is_finite())
            || !init.relative_equilibrium_tolerance.is_finite()
            || !(1e-12..=1e-6).contains(&init.relative_equilibrium_tolerance)
            || self.materials.is_empty()
            || (!self.is_general() && self.layers.is_empty())
        {
            return Err(Error::Input(
                "invalid ground initialization schema/settings".into(),
            ));
        }
        let mut ids = BTreeSet::new();
        for m in &self.materials {
            m.model.validate()?;
            if m.id == 0
                || !ids.insert(m.id)
                || m.name.is_empty()
                || !m.unit_weight_kn_m3.is_finite()
                || m.unit_weight_kn_m3 <= 0.
            {
                return Err(Error::Input("invalid/duplicate ground material".into()));
            }
        }
        let mut names = BTreeSet::new();
        let mut referenced: BTreeSet<_> = self.construction_material_ids.iter().copied().collect();
        if referenced.len() != self.construction_material_ids.len() || !referenced.is_subset(&ids) {
            return Err(Error::Input(
                "invalid/duplicate construction material IDs".into(),
            ));
        }
        if self.is_general() {
            if !self.layers.is_empty() || self.bodies.is_empty() {
                return Err(Error::Input(
                    "general initialization requires bodies, not K0 layers".into(),
                ));
            }
            self.initialization
                .geostatic
                .as_ref()
                .ok_or_else(|| Error::Input("equilibrate requires geostatic settings".into()))?
                .validate()?;
            for b in &self.bodies {
                if b.body.is_empty() || !names.insert(&b.body) || !ids.contains(&b.material_id) {
                    return Err(Error::Input(
                        "invalid/duplicate general body/material mapping".into(),
                    ));
                }
                referenced.insert(b.material_id);
            }
            if referenced != ids {
                return Err(Error::Input("unreferenced general material".into()));
            }
            return Ok(());
        }
        if self.initialization.geostatic.is_some() || !self.bodies.is_empty() {
            return Err(Error::Input(
                "horizontal K0 cannot also specify a general-geometry initialization".into(),
            ));
        }
        for layer in &self.layers {
            if layer.body.is_empty()
                || !names.insert(&layer.body)
                || !ids.contains(&layer.material_id)
                || !layer.k0_x.is_finite()
                || !layer.k0_y.is_finite()
                || layer.k0_x <= 0.
                || layer.k0_y <= 0.
            {
                return Err(Error::Input("invalid body/material/K0 mapping".into()));
            }
            referenced.insert(layer.material_id);
            let material = self
                .materials
                .iter()
                .find(|m| m.id == layer.material_id)
                .ok_or_else(|| Error::Input("unknown layer material".into()))?;
            if !material.model.is_mcc() {
                if matches!(material.model, soil::Model::InvariantDoubleHardening(_)) {
                    let cap = match layer.consolidation {
                        Consolidation::HardeningCapOverMeanEffective { ratio } => {
                            ratio.is_finite() && ratio >= 1.
                        }
                        Consolidation::HardeningExplicitCap { cap_pressure_kpa } => {
                            cap_pressure_kpa.is_finite() && cap_pressure_kpa > 0.
                        }
                        _ => false,
                    };
                    if !cap
                        || !matches!(layer.void_ratio,VoidRatio::Explicit{value} if value.is_finite()&&value>0.)
                    {
                        return Err(Error::Input("IDH layer requires its own cap state and explicit positive void ratio; MCC OCR definitions are not interchangeable".into()));
                    }
                    continue;
                }
                if !matches!(layer.consolidation, Consolidation::NotApplicable)
                    || !matches!(layer.void_ratio,VoidRatio::Explicit{value} if value.is_finite()&&value>0.)
                {
                    return Err(Error::Input("MC layer requires consolidation not_applicable and explicit positive initial void ratio".into()));
                }
                continue;
            }
            match layer.consolidation {
                Consolidation::MccPcOverMeanEffective { ratio }
                    if ratio.is_finite() && ratio >= 1. => {}
                Consolidation::VerticalOcr {
                    ratio,
                    preconsolidation_pressure_kpa,
                } if ratio.is_finite()
                    && ratio >= 1.
                    && preconsolidation_pressure_kpa.is_finite()
                    && preconsolidation_pressure_kpa > 0. => {}
                Consolidation::ExplicitPc {
                    preconsolidation_pressure_kpa,
                } if preconsolidation_pressure_kpa.is_finite()
                    && preconsolidation_pressure_kpa > 0. => {}
                _ => {
                    return Err(Error::Input(
                        "invalid consolidation pressure/OCR definition".into(),
                    ));
                }
            }
            if let VoidRatio::Explicit { value } = layer.void_ratio
                && (!value.is_finite() || value <= 0.)
            {
                return Err(Error::Input(
                    "explicit initial void ratio must be finite and positive".into(),
                ));
            }
        }
        if referenced != ids {
            return Err(Error::Input(
                "unreferenced material in ground configuration".into(),
            ));
        }
        Ok(())
    }
}

pub fn initialize(mesh: &Mesh, ground: Ground) -> Result<InitialState> {
    if ground.is_general() {
        return Err(Error::Input("general initialization needs displacement constraints and Newton settings; use solve-stages".into()));
    }
    let names: BTreeSet<_> = ground.excavation_bodies.iter().cloned().collect();
    if names.len() != ground.excavation_bodies.len()
        || names != mesh.excavation_sets.keys().cloned().collect()
    {
        return Err(Error::Input(
            "CAD excavation tools/ground mapping mismatch".into(),
        ));
    }
    mesh.validate()?;
    ground.validate()?;
    if !["soil-fem-mesh/2", "soil-fem-mesh/3"].contains(&mesh.schema.as_str())
        || mesh.regions.len() != ground.layers.len()
    {
        return Err(Error::Input(
            "body-resolved mesh/2 with one configured layer per CAD body required".into(),
        ));
    }
    let region_map: BTreeMap<_, _> = mesh.regions.iter().map(|r| (r.id, r)).collect();
    let materials: BTreeMap<_, _> = ground.materials.iter().map(|m| (m.id, m)).collect();
    let layers: BTreeMap<_, _> = ground.layers.iter().map(|l| (l.body.as_str(), l)).collect();
    let mut sorted: Vec<_> = mesh.regions.iter().collect();
    sorted.sort_by(|a, b| a.bounds_m[0][2].total_cmp(&b.bounds_m[0][2]));
    let base = sorted[0];
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-8 * (1. + a.abs().max(b.abs()));
    for (i, r) in sorted.iter().enumerate() {
        let layer = layers
            .get(r.body.as_str())
            .ok_or_else(|| Error::Input("CAD body missing in ground configuration".into()))?;
        if layer.material_id != r.material_id {
            return Err(Error::Input(
                "CAD region/material assignment changed".into(),
            ));
        }
        for side in 0..2 {
            for axis in 0..2 {
                if !near(r.bounds_m[side][axis], base.bounds_m[side][axis]) {
                    return Err(Error::Input("horizontal K0 requires layers with identical rectangular footprints; sloping/general geometry needs geostatic solve".into()));
                }
            }
        }
        if i > 0 && !near(sorted[i - 1].bounds_m[1][2], r.bounds_m[0][2]) {
            return Err(Error::Input(
                "horizontal layers have a gap or overlap".into(),
            ));
        }
        let cuboid = (0..3)
            .map(|a| r.bounds_m[1][a] - r.bounds_m[0][a])
            .product::<f64>();
        if !near(cuboid, r.volume_m3) {
            return Err(Error::Input(
                "K0 initializer accepts horizontal cuboid layers only".into(),
            ));
        }
    }
    let top = sorted
        .last()
        .ok_or_else(|| Error::Input("no layers".into()))?
        .bounds_m[1][2];
    if ground
        .initialization
        .water_level_m
        .is_some_and(|w| w > top + 1e-9)
    {
        return Err(Error::Input(
            "water above ground needs external water tractions; not supported by this initializer"
                .into(),
        ));
    }
    if let Some(h) = ground.initialization.water_level_m
        && h > base.bounds_m[0][2]
        && h < top
        && !sorted
            .iter()
            .any(|r| near(h, r.bounds_m[0][2]) || near(h, r.bounds_m[1][2]))
    {
        return Err(Error::Input("water table must coincide with a horizontal layer boundary for centroid stress quadrature".into()));
    }
    let n = mesh.nodes.len() * 3;
    let mut internal = vec![0.; n];
    let mut external = vec![0.; n];
    let mut points = Vec::new();
    let mut volumes = BTreeMap::<String, f64>::new();
    // Face identities also prove that adjacent materials share the same nodes.
    let mut faces = BTreeMap::<[usize; 3], Vec<(usize, u32)>>::new();
    for (ei, tet) in mesh.tetrahedra.iter().enumerate() {
        let r = region_map
            .get(&mesh.region_ids[ei])
            .ok_or_else(|| Error::Input("unknown region".into()))?;
        let layer = layers[r.body.as_str()];
        let material = materials[&r.material_id];
        let nodes = mesh.cell_nodes(ei);
        let p: Vec<_> = nodes.iter().map(|&i| mesh.nodes[i]).collect();
        if p.iter().any(|v| {
            (0..3).any(|a| v[a] < r.bounds_m[0][a] - 1e-8 || v[a] > r.bounds_m[1][a] + 1e-8)
        }) {
            return Err(Error::Input("element crosses its CAD layer bounds".into()));
        }
        for (qi, q) in mesh.integration(ei)?.into_iter().enumerate() {
            let centroid = q.position;
            let mut sv = ground.initialization.surface_surcharge_kpa;
            for rr in &sorted {
                let dz = (rr.bounds_m[1][2] - centroid[2].max(rr.bounds_m[0][2])).max(0.);
                sv += dz * materials[&rr.material_id].unit_weight_kn_m3;
            }
            let pore = ground.initialization.water_level_m.map_or(0., |h| {
                ground.initialization.water_unit_weight_kn_m3 * (h - centroid[2]).max(0.)
            });
            let sv_eff = sv - pore;
            if sv_eff <= 0. {
                return Err(Error::Input(
                    "non-positive initial vertical effective stress".into(),
                ));
            }
            let eff = [layer.k0_x * sv_eff, layer.k0_y * sv_eff, sv_eff, 0., 0., 0.];
            let mean = mcc::invariants(&eff).0;
            let (state, ocr_v, pc_ratio) = if let soil::Model::Mcc(parameters) = &material.model {
                let (pc, ocr_v) = match layer.consolidation {
                    Consolidation::MccPcOverMeanEffective { ratio } => (ratio * mean, None),
                    Consolidation::VerticalOcr {
                        ratio,
                        preconsolidation_pressure_kpa,
                    } => (preconsolidation_pressure_kpa, Some(ratio)),
                    Consolidation::ExplicitPc {
                        preconsolidation_pressure_kpa,
                    } => (preconsolidation_pressure_kpa, None),
                    _ => {
                        return Err(Error::Input("MCC requires a consolidation state".into()));
                    }
                };
                let expected_e = parameters.reference_void_ratio(mean, pc)?;
                let e = match layer.void_ratio {
                    VoidRatio::FromMccReference => expected_e,
                    VoidRatio::Explicit { value } => {
                        if !near(value, expected_e) {
                            return Err(Error::Input(format!(
                                "explicit e/pc/reference inconsistency in {} at element {ei}",
                                r.body
                            )));
                        }
                        value
                    }
                };
                let state = soil::State::Mcc(mcc::State {
                    effective_stress_compression_positive_kpa: eff,
                    void_ratio: e,
                    preconsolidation_pressure_kpa: pc,
                    plastic_strain_compression_positive: [0.; 6],
                });
                (state, ocr_v, Some(pc / mean))
            } else {
                let VoidRatio::Explicit { value: e } = layer.void_ratio else {
                    return Err(Error::Input("MC requires explicit void ratio".into()));
                };
                let state = match &material.model {
                    soil::Model::MohrCoulomb(_) => soil::State::MohrCoulomb(mc::State {
                        effective_stress_compression_positive_kpa: eff,
                        plastic_strain_compression_positive: [0.; 6],
                        void_ratio: e,
                        equivalent_plastic_strain: 0.,
                        plastic_work_kpa: 0.,
                        last_return: mc::Return::Elastic,
                    }),
                    soil::Model::InvariantDoubleHardening(parameters) => {
                        let cap = match layer.consolidation {
                            Consolidation::HardeningCapOverMeanEffective { ratio } => ratio * mean,
                            Consolidation::HardeningExplicitCap { cap_pressure_kpa } => {
                                cap_pressure_kpa
                            }
                            _ => return Err(Error::Input("IDH requires its own cap state".into())),
                        };
                        soil::State::InvariantDoubleHardening(
                            parameters.initial_state(eff, e, cap)?,
                        )
                    }
                    _ => return Err(Error::Input("unexpected initial material type".into())),
                };
                (state, None, None)
            };
            material.model.validate_state(&state)?;
            let e = state.void_ratio();
            let mut total = eff;
            for s in total.iter_mut().take(3) {
                *s += pore;
            }
            let volume = q.weight_m3;
            let b = q.b;
            *volumes.entry(r.body.clone()).or_default() += volume;
            for i in 0..3 * nodes.len() {
                let gi = 3 * nodes[i / 3] + i % 3;
                internal[gi] -= volume * (0..6).map(|j| b[(j, i)] * total[j]).sum::<f64>();
                if i % 3 == 2 {
                    external[gi] -= material.unit_weight_kn_m3 * volume * q.shape[i / 3];
                }
            }
            if qi == 0 {
                for mut f in [
                    [tet[0], tet[1], tet[2]],
                    [tet[0], tet[1], tet[3]],
                    [tet[0], tet[2], tet[3]],
                    [tet[1], tet[2], tet[3]],
                ] {
                    f.sort_unstable();
                    faces.entry(f).or_default().push((ei, r.id));
                }
            }
            let (bulk, shear) = material.model.moduli(&state)?;
            points.push(Point {
                element_index: ei,
                local_index: qi,
                weight_m3: volume,
                natural_coordinates: q.natural,
                region_id: r.id,
                material_id: r.material_id,
                centroid_m: centroid,
                total_stress_compression_positive_kpa: total,
                pore_pressure_kpa: pore,
                porosity: e / (1. + e),
                mcc_pc_over_mean_effective: pc_ratio,
                vertical_ocr: ocr_v,
                past_maximum_vertical_effective_stress_kpa: ocr_v.map(|ratio| ratio * sv_eff),
                normalized_yield: material.model.normalized_yield(&state)?,
                bulk_modulus_kpa: bulk,
                shear_modulus_kpa: shear,
                constitutive_state: state.into(),
            });
        }
    }
    for r in &mesh.regions {
        if !near(volumes[&r.body], r.volume_m3) {
            return Err(Error::Input("mesh/CAD region volume mismatch".into()));
        }
    }
    let mut neighbors = BTreeSet::new();
    let mut interfaces = 0;
    for (face, owners) in &faces {
        if owners.len() == 2 && owners[0].1 != owners[1].1 {
            let (a, b) = (owners[0].1, owners[1].1);
            neighbors.insert((a.min(b), a.max(b)));
            interfaces += 1;
        }
        if owners.len() == 1 && face.iter().all(|&i| near(mesh.nodes[i][2], top)) {
            let nodes = mesh.face_nodes(owners[0].0, face)?;
            let points: Vec<_> = nodes.iter().map(|&i| mesh.nodes[i]).collect();
            for q in crate::elements::face_points(&points)? {
                for (i, &node) in nodes.iter().enumerate() {
                    for a in 0..3 {
                        external[3 * node + a] -= ground.initialization.surface_surcharge_kpa
                            * q.shape[i]
                            * q.area_vector[a];
                    }
                }
            }
        }
    }
    for pair in sorted.windows(2) {
        if !neighbors.contains(&(pair[0].id.min(pair[1].id), pair[0].id.max(pair[1].id))) {
            return Err(Error::Input(
                "adjacent CAD layers do not share conforming interface faces".into(),
            ));
        }
    }
    let mut constrained = vec![false; n];
    for (i, p) in mesh.nodes.iter().enumerate() {
        constrained[3 * i] = near(p[0], base.bounds_m[0][0]) || near(p[0], base.bounds_m[1][0]);
        constrained[3 * i + 1] = near(p[1], base.bounds_m[0][1]) || near(p[1], base.bounds_m[1][1]);
        constrained[3 * i + 2] = near(p[2], base.bounds_m[0][2]);
    }
    let norm = external.iter().map(|v| v * v).sum::<f64>().sqrt();
    let free = (0..n)
        .filter(|&i| !constrained[i])
        .map(|i| (external[i] - internal[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    let residual = free / norm;
    if !residual.is_finite() || residual > ground.initialization.relative_equilibrium_tolerance {
        return Err(Error::Solve(format!(
            "K0 is not in FE equilibrium: true free residual {residual:e}; no state accepted"
        )));
    }
    let mut applied = [0.; 3];
    let mut reaction = [0.; 3];
    for i in 0..n {
        applied[i % 3] += external[i];
        if constrained[i] {
            reaction[i % 3] += internal[i] - external[i];
        }
    }
    Ok(InitialState {
        schema: if ground.materials.iter().all(|m| m.model.is_mcc()) {
            "soil-fem-initial-state/1"
        } else {
            "soil-fem-initial-state/2"
        },
        mesh_schema: mesh.schema.clone(),
        ground,
        point_states: points,
        body_volumes_m3: volumes,
        shared_interface_triangles: interfaces,
        true_relative_free_force_residual: residual,
        applied_force_kn: applied,
        support_reaction_kn: reaction,
        equilibrium_error_kn: std::array::from_fn(|i| applied[i] + reaction[i]),
        initialized_displacement_m: vec![[0.; 3]; mesh.nodes.len()],
    })
}
