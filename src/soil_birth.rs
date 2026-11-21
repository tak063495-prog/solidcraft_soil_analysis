//! Explicit soil placement states at a construction event.
use crate::{
    Error, Result, darcy, hydraulic_transition::ExternalWaterParty, mcc, mesh::Mesh,
    richards::AffineHead,
};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mixture {
    pub porosity: f64,
    pub upper_saturation: f64,
    pub particle_density_kg_m3: f64,
    pub water_density_kg_m3: f64,
    pub gravity_m_s2: f64,
}
impl Mixture {
    pub fn unit_weights(&self, gamma_w: f64) -> Result<[f64; 2]> {
        let water = self.gravity_m_s2 * self.water_density_kg_m3 / 1000.;
        if !(0. < self.porosity && self.porosity < 1.)
            || !(0. ..=1.).contains(&self.upper_saturation)
            || [
                self.gravity_m_s2,
                self.particle_density_kg_m3,
                self.water_density_kg_m3,
                water,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.)
            || !consistent(water, gamma_w)
        {
            return Err(Error::Input(
                "invalid two-region density, saturation or water unit weight".into(),
            ));
        }
        let dry = (1. - self.porosity) * self.gravity_m_s2 * self.particle_density_kg_m3 / 1000.;
        Ok([
            dry + self.porosity * self.upper_saturation * water,
            dry + self.porosity * water,
        ])
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Formation {
    pub history: String,
    pub hydraulic_porosity: f64,
    pub pc_over_mean_effective: Option<f64>,
}
pub(crate) fn consistent(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-8 * (1. + a.abs().max(b.abs()))
}
impl Formation {
    pub fn validate(&self, parameters: &mcc::Parameters, state: &mcc::State) -> Result<()> {
        parameters.validate_state(state)?;
        let p = mcc::invariants(&state.effective_stress_compression_positive_kpa).0;
        let pc = state.preconsolidation_pressure_kpa;
        if self.history.trim().is_empty()
            || !(0. < self.hydraulic_porosity && self.hydraulic_porosity < 1.)
            || !consistent(
                self.hydraulic_porosity,
                state.void_ratio / (1. + state.void_ratio),
            )
            || !consistent(state.void_ratio, parameters.reference_void_ratio(p, pc)?)
            || self
                .pc_over_mean_effective
                .is_some_and(|ocr| !consistent(ocr, pc / p))
        {
            return Err(Error::Input(
                "inconsistent saturated MCC formation history, porosity, e or mean-stress OCR"
                    .into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlacedWater {
    pub generation: usize,
    pub material_id: u32,
    pub element_index: usize,
    pub lower_volume_m3: f64,
    pub upper_volume_m3: f64,
    pub lower_water_m3: f64,
    pub upper_water_m3: f64,
    pub lower_weight_kn: f64,
    pub upper_weight_kn: f64,
    pub water_source: ExternalWaterParty,
}

/// A quadratic tet is inside the saturated half-space when all Bernstein
/// control values of its affine pressure-head function are nonnegative.
pub(crate) fn saturated_element(
    mesh: &Mesh,
    e: usize,
    pressure_head: impl Fn([f64; 3]) -> f64,
) -> bool {
    let depths: Vec<_> = mesh
        .cell_nodes(e)
        .iter()
        .map(|&n| pressure_head(mesh.nodes[n]))
        .collect();
    depths.iter().all(|&d| d.is_finite() && d >= 0.)
        && (depths.len() == 4
            || crate::elements::EDGES
                .iter()
                .enumerate()
                .all(|(k, &[i, j])| 2. * depths[k + 4] - (depths[i] + depths[j]) / 2. >= 0.))
}
/// Keep the old supported field immutable. Only freshly supported nodes can
/// acquire a head; shared nodes must satisfy the supplied placement condition.
pub(crate) fn extend_heads(
    mesh: &Mesh,
    old: &darcy::HeadField,
    selected: &[usize],
    head: &AffineHead,
) -> Result<darcy::HeadField> {
    head.validate()?;
    let mut field = old.clone();
    let datum = selected
        .iter()
        .flat_map(|&e| mesh.cell_nodes(e))
        .find(|&n| old.defined_nodes[n])
        .map_or(head.head_m, |n| old.datum_head_m[n]);
    for &e in selected {
        if !saturated_element(mesh, e, |p| head.relative(p, p[2])) {
            return Err(Error::Input(
                "post-Darcy placement requires full saturation".into(),
            ));
        }
        let nodes = mesh.cell_nodes(e);
        for n in nodes {
            let value = head.relative(mesh.nodes[n], datum);
            if !value.is_finite() || head.relative(mesh.nodes[n], mesh.nodes[n][2]) < 0. {
                return Err(Error::Input(
                    "post-Darcy placement requires finite saturated head".into(),
                ));
            }
            if old.defined_nodes[n] {
                if (old.datum_head_m[n] - datum + old.relative_head_m[n] - value).abs() > 1e-9 {
                    return Err(Error::Input(
                        "placement head conflicts with committed soil trace".into(),
                    ));
                }
            } else {
                field.defined_nodes[n] = true;
                field.datum_head_m[n] = datum;
                field.relative_head_m[n] = value;
                field.head_m[n] = datum + value;
            }
        }
    }
    Ok(field)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlacementWater {
    pub material_id: u32,
    pub water_source: ExternalWaterParty,
}
