//! A pointwise MC support-function witness for discrete FE equilibrium exclusion.
//! A virtual displacement is a separating test field, not an actual flow rule.
use crate::{Error, Result, boundary, checkpoint, mc, mesh::Mesh, srm_bound::WaterPolicy};
use nalgebra::{DVector, Matrix3};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub id: u32,
    pub parameters: mc::Parameters,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Problem {
    pub schema: String,
    pub water_policy: WaterPolicy,
    pub materials: Vec<Material>,
    pub constraints: Vec<boundary::Constraint>,
    pub physical_nodal_force_kn: Vec<f64>,
    pub active_elements: Vec<bool>,
    pub virtual_displacement_m: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointCheck {
    pub element_index: usize,
    pub local_index: usize,
    pub material_id: u32,
    pub weight_m3: f64,
    pub position_m: [f64; 3],
    pub natural_coordinates: [f64; 3],
    pub strain_tensor_scale: f64,
    /// xx, yy, zz, xy, yz, xz; off-diagonal values are tensor shear.
    pub normalized_strain_tensor_tension_positive: [f64; 6],
    pub normalized_principal_strains_descending: [f64; 3],
    pub effective_cohesion_kpa: f64,
    pub effective_friction_radians: f64,
    pub relative_dual_margin: Option<f64>,
    pub accepted_finite_support: bool,
    pub support_density_kpa: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Certificate {
    pub schema: String,
    pub mesh_sha256: String,
    pub problem_sha256: String,
    pub virtual_displacement_sha256: String,
    pub point_checks_sha256: String,
    pub point_hash_framing: String,
    pub active_elements: usize,
    pub integration_points: usize,
    pub nonzero_strain_points: usize,
    pub volume_m3: f64,
    pub constrained_components_checked: usize,
    pub all_constrained_components_exactly_zero: bool,
    pub minimum_relative_dual_margin: f64,
    pub required_relative_dual_margin: f64,
    pub external_virtual_work_kn_m: f64,
    pub maximum_admissible_internal_virtual_work_kn_m: f64,
    pub virtual_work_excess_kn_m: f64,
    pub numerical_work_clearance_kn_m: f64,
    pub excludes_admissible_discrete_equilibrium: bool,
    pub pure_mc_superset_used_when_rankine_is_present: bool,
    pub continuum_collapse_or_stability_certified: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Assessment {
    Evaluated { certificate: Box<Certificate> },
    NotApplicable { reason: String },
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Audit {
    pub assessment: Assessment,
    pub point_checks: Vec<PointCheck>,
    pub complete_point_audit: bool,
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let t = self.value + x;
        self.correction += if self.value.abs() >= x.abs() {
            (self.value - t) + x
        } else {
            (x - t) + self.value
        };
        self.value = t;
    }
    fn total(&self) -> f64 {
        self.value + self.correction
    }
}

fn unavailable(reason: impl Into<String>) -> Result<Assessment> {
    Ok(Assessment::NotApplicable {
        reason: reason.into(),
    })
}

pub fn analyze(mesh: &Mesh, problem: &Problem) -> Result<Assessment> {
    analyze_with_observer(mesh, problem, &mut |_| {})
}

/// Optional full diagnostic. Normal SRM uses analyze and retains only the
/// certificate and framed hash, avoiding a GP array for every search candidate.
pub fn audit(mesh: &Mesh, problem: &Problem) -> Result<Audit> {
    let mut points = Vec::new();
    let assessment = analyze_with_observer(mesh, problem, &mut |p| points.push(p.clone()))?;
    let complete = matches!(assessment, Assessment::Evaluated { .. });
    Ok(Audit {
        assessment,
        point_checks: points,
        complete_point_audit: complete,
    })
}

pub fn analyze_with_observer(
    mesh: &Mesh,
    problem: &Problem,
    observer: &mut impl FnMut(&PointCheck),
) -> Result<Assessment> {
    mesh.validate()?;
    if problem.schema != "soil-fem-mc-virtual-work-bound/1"
        || problem.active_elements.len() != mesh.tetrahedra.len()
        || !problem.active_elements.iter().any(|a| *a)
        || problem.physical_nodal_force_kn.len() != 3 * mesh.nodes.len()
        || problem.virtual_displacement_m.len() != mesh.nodes.len()
        || !problem
            .physical_nodal_force_kn
            .iter()
            .all(|x| x.is_finite())
        || !problem
            .virtual_displacement_m
            .iter()
            .flatten()
            .all(|x| x.is_finite())
    {
        return Err(Error::Input(
            "invalid virtual-work schema, dimensions or nonfinite data".into(),
        ));
    }
    let mut materials = BTreeMap::new();
    for material in &problem.materials {
        material.parameters.validate()?;
        if material.id == 0
            || materials
                .insert(material.id, &material.parameters)
                .is_some()
        {
            return Err(Error::Input(
                "invalid/duplicate virtual-work material ID".into(),
            ));
        }
    }
    let conditions = boundary::Conditions {
        constraints: problem.constraints.clone(),
        nodal_loads: vec![],
        surface_loads: vec![],
    };
    conditions.validate_definition(mesh)?;
    let mut used = vec![false; mesh.nodes.len()];
    for (element, enabled) in problem.active_elements.iter().enumerate() {
        if *enabled {
            if !materials.contains_key(&mesh.material_ids[element]) {
                return Err(Error::Input(
                    "active virtual-work element has no material".into(),
                ));
            }
            for n in mesh.cell_nodes(element) {
                used[n] = true;
            }
        }
    }
    for (n, active) in used.iter().enumerate() {
        if !active
            && (problem.virtual_displacement_m[n] != [0.; 3]
                || problem.physical_nodal_force_kn[3 * n..3 * n + 3]
                    .iter()
                    .any(|f| *f != 0.))
        {
            return Err(Error::Input(
                "virtual displacement/load on inactive node".into(),
            ));
        }
    }
    let mut constrained = 0;
    for c in &problem.constraints {
        if c.value_m != 0. || c.gradient != [0.; 3] {
            return unavailable("nonzero prescribed displacements need their own work policy");
        }
        let axes: Vec<_> = c
            .direction
            .iter()
            .enumerate()
            .filter_map(|(i, x)| (*x != 0.).then_some(i))
            .collect();
        if axes.len() != 1 {
            return unavailable(
                "initial virtual-work verifier requires axis-aligned support directions",
            );
        }
        let axis = axes[0];
        for &n in &mesh.node_sets[&c.node_set] {
            if used[n] {
                constrained += 1;
                if problem.virtual_displacement_m[n][axis] != 0. {
                    return unavailable(format!(
                        "virtual displacement violates supported component at node {n}, axis {axis}"
                    ));
                }
            }
        }
    }
    // Unsupported rigid-body components are invalid input, not a soil upper.
    let _supported =
        boundary::Transform::new(mesh, &problem.constraints, problem.active_elements.clone())?;
    let mut volume = Sum::default();
    let mut capacity = Sum::default();
    let mut external = Sum::default();
    let mut absolute_external = Sum::default();
    for (n, u) in problem.virtual_displacement_m.iter().enumerate() {
        for (axis, v) in u.iter().enumerate() {
            let work = problem.physical_nodal_force_kn[3 * n + axis] * v;
            if !work.is_finite() {
                return Err(Error::Input("virtual external work overflow".into()));
            }
            external.add(work);
            absolute_external.add(work.abs());
        }
    }
    let mut point_hash = Sha256::new();
    let mut count = 0;
    let mut nonzero = 0;
    let required = 1e-8;
    let mut minimum = f64::INFINITY;
    for (element, enabled) in problem.active_elements.iter().enumerate() {
        if !enabled {
            continue;
        }
        let id = mesh.material_ids[element];
        let (cohesion, phi, _) = materials[&id].effective_strength();
        if phi <= 1e-8 {
            return unavailable(
                "zero/near-zero friction requires independently established divergence-free virtual fields; no small trace is rounded to zero",
            );
        }
        let ratio = (1. + phi.sin()) / (1. - phi.sin());
        let apex = cohesion / phi.tan();
        if !ratio.is_finite() || !apex.is_finite() {
            return Err(Error::Input(
                "nonfinite MC virtual support parameters".into(),
            ));
        }
        let nodes = mesh.cell_nodes(element);
        let origin = problem.virtual_displacement_m[nodes[0]];
        let local = DVector::from_iterator(
            3 * nodes.len(),
            nodes.iter().flat_map(|&n| {
                (0..3).map(move |axis| problem.virtual_displacement_m[n][axis] - origin[axis])
            }),
        );
        for (index, q) in mesh.integration(element)?.iter().enumerate() {
            if !q.weight_m3.is_finite() || q.weight_m3 <= 0. {
                return Err(Error::Input(
                    "virtual-work convex bound requires positive quadrature weights".into(),
                ));
            }
            let e = &q.b * &local;
            if !e.iter().all(|x| x.is_finite()) {
                return Err(Error::Input("virtual strain overflow".into()));
            }
            let tensor = Matrix3::new(
                e[0],
                e[3] / 2.,
                e[5] / 2.,
                e[3] / 2.,
                e[1],
                e[4] / 2.,
                e[5] / 2.,
                e[4] / 2.,
                e[2],
            );
            let scale = tensor.iter().map(|x| x.abs()).fold(0., f64::max);
            let mut point = PointCheck {
                element_index: element,
                local_index: index,
                material_id: id,
                weight_m3: q.weight_m3,
                position_m: q.position,
                natural_coordinates: q.natural,
                strain_tensor_scale: scale,
                normalized_strain_tensor_tension_positive: [0.; 6],
                normalized_principal_strains_descending: [0.; 3],
                effective_cohesion_kpa: cohesion,
                effective_friction_radians: phi,
                relative_dual_margin: None,
                accepted_finite_support: true,
                support_density_kpa: Some(0.),
            };
            if scale > 0. {
                let normal = tensor / scale;
                let mut eigen: [f64; 3] = normal.symmetric_eigen().eigenvalues.into();
                eigen.sort_by(|a, b| b.total_cmp(a));
                let norm = normal.norm();
                let margin =
                    normal.trace() + (ratio - 1.) * eigen.iter().map(|x| x.min(0.)).sum::<f64>();
                let relative = margin / ((1. + ratio) * norm);
                if !relative.is_finite() {
                    return Err(Error::Input("nonfinite virtual dual evaluation".into()));
                }
                point.normalized_strain_tensor_tension_positive = [
                    normal[(0, 0)],
                    normal[(1, 1)],
                    normal[(2, 2)],
                    normal[(0, 1)],
                    normal[(1, 2)],
                    normal[(0, 2)],
                ];
                point.normalized_principal_strains_descending = eigen;
                point.relative_dual_margin = Some(relative);
                point.accepted_finite_support = relative > required;
                point.support_density_kpa = point
                    .accepted_finite_support
                    .then_some(apex * normal.trace() * scale);
                if !point.accepted_finite_support {
                    observer(&point);
                    return unavailable(format!(
                        "virtual strain outside strict MC dual domain at element {element}, point {index}"
                    ));
                }
                nonzero += 1;
                minimum = minimum.min(relative);
            }
            let density = point.support_density_kpa.unwrap();
            let work = density * q.weight_m3;
            if !work.is_finite() || density < 0. {
                return Err(Error::Input("nonfinite/negative virtual MC support".into()));
            }
            capacity.add(work);
            volume.add(q.weight_m3);
            count += 1;
            let raw = serde_json::to_vec(&point)?;
            point_hash.update((raw.len() as u64).to_le_bytes());
            point_hash.update(&raw);
            observer(&point);
        }
    }
    if nonzero == 0 {
        return unavailable(
            "zero strain field cannot establish a material collapse upper; check support mechanisms separately",
        );
    }
    let work = external.total();
    let cap = capacity.total();
    let v = volume.total();
    let scale = absolute_external
        .total()
        .max(cap.abs())
        .max(work.abs())
        .max(f64::MIN_POSITIVE);
    let clearance = (1e-8 + 4096. * f64::EPSILON * (count + mesh.nodes.len() + 1) as f64) * scale;
    let excess = work - cap;
    if ![work, cap, v, scale, clearance, excess]
        .iter()
        .all(|x| x.is_finite())
        || v <= 0.
    {
        return Err(Error::Input(
            "nonfinite accumulated virtual-work bound".into(),
        ));
    }
    Ok(Assessment::Evaluated {
        certificate: Box::new(Certificate {
            schema: "soil-fem-mc-virtual-work-bound-result/1".into(),
            mesh_sha256: checkpoint::hash_json(mesh)?,
            problem_sha256: checkpoint::hash_json(problem)?,
            virtual_displacement_sha256: checkpoint::hash_json(&problem.virtual_displacement_m)?,
            point_checks_sha256: format!("{:x}", point_hash.finalize()),
            point_hash_framing: "u64_le_length_then_serde_json_for_each_active_GP_in_order".into(),
            active_elements: problem.active_elements.iter().filter(|a| **a).count(),
            integration_points: count,
            nonzero_strain_points: nonzero,
            volume_m3: v,
            constrained_components_checked: constrained,
            all_constrained_components_exactly_zero: true,
            minimum_relative_dual_margin: minimum,
            required_relative_dual_margin: required,
            external_virtual_work_kn_m: work,
            maximum_admissible_internal_virtual_work_kn_m: cap,
            virtual_work_excess_kn_m: excess,
            numerical_work_clearance_kn_m: clearance,
            excludes_admissible_discrete_equilibrium: excess > clearance,
            pure_mc_superset_used_when_rankine_is_present: problem
                .materials
                .iter()
                .any(|m| m.parameters.tension_cutoff_kpa.is_some()),
            continuum_collapse_or_stability_certified: false,
        }),
    })
}
