//! A necessary yield-admissibility condition independent of Newton status.
//!
//! With exactly six independent rigid-body gauge reactions, global force and
//! moment balance determines every support reaction. Affine virtual work then
//! determines the volume-mean stress. MC and Rankine yield domains are convex:
//! an inadmissible mean excludes any admissible stress field at the positive
//! integration points. This is a discrete physical upper *candidate*, not a
//! continuum collapse certificate or a sufficiency/stability test.
use crate::{Error, Result, boundary, checkpoint, mc, mesh::Mesh};
use nalgebra::{Matrix3, SMatrix, SVector, Vector3};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WaterPolicy {
    Dry,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Problem {
    pub schema: String,
    pub material_id: u32,
    pub water_policy: WaterPolicy,
    /// The actual material at this trial, after the caller's explicit reduction.
    pub material: mc::Parameters,
    pub boundary: boundary::Conditions,
    pub body_force_kn_m3: [f64; 3],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub active_elements: Option<Vec<bool>>,
}

/// Already assembled physical nodal forces. Construction callers must bind and
/// independently reconstruct these loads; this API does not authenticate them.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NodalProblem {
    pub schema: String,
    pub material_id: u32,
    pub water_policy: WaterPolicy,
    pub material: mc::Parameters,
    pub constraints: Vec<boundary::Constraint>,
    pub physical_nodal_force_kn: Vec<f64>,
    pub active_elements: Vec<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GaugeReaction {
    pub node_index: usize,
    pub unit_direction: [f64; 3],
    pub magnitude_kn: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Certificate {
    pub schema: String,
    pub mesh_sha256: String,
    pub problem_sha256: String,
    pub volume_m3: f64,
    pub integration_points: usize,
    pub moment_origin_m: [f64; 3],
    pub moment_length_scale_m: f64,
    pub gauge_condition_number: f64,
    pub gauge_reactions: Vec<GaugeReaction>,
    pub relative_resultant_force_and_moment_residual: f64,
    pub mean_effective_stress_compression_positive_kpa: [f64; 6],
    pub principal_stresses_descending_kpa: [f64; 3],
    pub mean_mc_yield_margin_kpa: f64,
    pub mean_rankine_yield_margin_kpa: Option<f64>,
    pub numerical_clearance_kpa: f64,
    pub mean_stress_excludes_admissible_discrete_equilibrium: bool,
    pub positive_integration_weight_convexity_argument: bool,
    pub lower_bound_or_stability_proven: bool,
    pub continuum_collapse_factor_certified: bool,
    pub dry_effective_equals_total: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum Assessment {
    Evaluated { certificate: Box<Certificate> },
    NotApplicable { reason: String },
}

fn unavailable(reason: &str) -> Result<Assessment> {
    Ok(Assessment::NotApplicable {
        reason: reason.into(),
    })
}

pub fn analyze(mesh: &Mesh, problem: &Problem) -> Result<Assessment> {
    analyze_impl(mesh, problem, None, checkpoint::hash_json(problem)?)
}

pub fn analyze_nodal(mesh: &Mesh, nodal: &NodalProblem) -> Result<Assessment> {
    if nodal.schema != "soil-fem-mc-nodal-mean-stress-bound/1"
        || nodal.physical_nodal_force_kn.len() != 3 * mesh.nodes.len()
        || !nodal.physical_nodal_force_kn.iter().all(|x| x.is_finite())
    {
        return Err(Error::Input("invalid nodal exclusion problem".into()));
    }
    let problem = Problem {
        schema: "soil-fem-mc-mean-stress-bound/1".into(),
        material_id: nodal.material_id,
        water_policy: nodal.water_policy,
        material: nodal.material.clone(),
        boundary: boundary::Conditions {
            constraints: nodal.constraints.clone(),
            nodal_loads: vec![],
            surface_loads: vec![],
        },
        body_force_kn_m3: [0.; 3],
        active_elements: Some(nodal.active_elements.clone()),
    };
    analyze_impl(
        mesh,
        &problem,
        Some(&nodal.physical_nodal_force_kn),
        checkpoint::hash_json(nodal)?,
    )
}

fn analyze_impl(
    mesh: &Mesh,
    problem: &Problem,
    assembled_forces: Option<&[f64]>,
    problem_sha256: String,
) -> Result<Assessment> {
    mesh.validate()?;
    problem.material.validate()?;
    problem.boundary.validate_definition(mesh)?;
    if problem.schema != "soil-fem-mc-mean-stress-bound/1"
        || problem.material_id == 0
        || !problem.body_force_kn_m3.iter().all(|x| x.is_finite())
    {
        return Err(Error::Input("invalid mean-stress exclusion problem".into()));
    }
    let active = problem
        .active_elements
        .clone()
        .unwrap_or_else(|| vec![true; mesh.tetrahedra.len()]);
    if active.len() != mesh.tetrahedra.len() || !active.iter().any(|x| *x) {
        return Err(Error::Input("invalid exclusion active mask".into()));
    }
    if active
        .iter()
        .enumerate()
        .any(|(i, a)| *a && mesh.material_ids[i] != problem.material_id)
    {
        return unavailable(
            "a single homogeneous MC material ID is required; no mixture bound is inferred",
        );
    }
    let mut used = vec![false; mesh.nodes.len()];
    for (i, a) in active.iter().enumerate() {
        if *a {
            for n in mesh.cell_nodes(i) {
                used[n] = true;
            }
        }
    }
    if assembled_forces.is_some_and(|f| {
        used.iter()
            .enumerate()
            .any(|(n, active)| !active && f[3 * n..3 * n + 3].iter().any(|v| *v != 0.))
    }) {
        return Err(Error::Input("nodal exclusion load on inactive node".into()));
    }
    let mut gauges = Vec::new();
    for c in &problem.boundary.constraints {
        if c.value_m != 0. || c.gradient != [0.; 3] {
            return unavailable("only zero rigid-body gauges are supported by this witness");
        }
        let direction = Vector3::from(c.direction).normalize();
        for &node in &mesh.node_sets[&c.node_set] {
            if !used[node] {
                return unavailable("a gauge lies outside the active domain");
            }
            gauges.push((node, direction));
        }
    }
    if gauges.len() != 6 {
        return unavailable(
            "exactly six independent gauge reactions are required; distributed supports leave the mean stress undetermined",
        );
    }
    let mut lo = Vector3::repeat(f64::INFINITY);
    let mut hi = Vector3::repeat(f64::NEG_INFINITY);
    for (i, a) in used.iter().enumerate() {
        if *a {
            let x = Vector3::from(mesh.nodes[i]);
            lo = lo.inf(&x);
            hi = hi.sup(&x);
        }
    }
    let origin = lo + (hi - lo) * 0.5;
    let length = (hi - lo).amax();
    if !length.is_finite() || length <= 0. {
        return Err(Error::Input("invalid exclusion length scale".into()));
    }
    let mut a = SMatrix::<f64, 6, 6>::zeros();
    for (j, (node, direction)) in gauges.iter().enumerate() {
        let r = (Vector3::from(mesh.nodes[*node]) - origin) / length;
        a.fixed_view_mut::<3, 1>(0, j).copy_from(direction);
        a.fixed_view_mut::<3, 1>(3, j)
            .copy_from(&r.cross(direction));
    }
    let svd = a.svd(false, false);
    let largest = svd.singular_values.max();
    let smallest = svd.singular_values.min();
    if !smallest.is_finite() || smallest <= largest * 1e-8 {
        return unavailable("gauge force/moment map is singular or too ill-conditioned");
    }
    let condition = largest / smallest;
    // Also reject disconnected islands or constraints that restrict strain.
    let transform = boundary::Transform::new(mesh, &problem.boundary.constraints, active.clone())?;
    if transform.count + 6 != 3 * used.iter().filter(|x| **x).count() {
        return unavailable("constraints are not six independent rigid-body gauges");
    }
    let mut forces = match assembled_forces {
        Some(f) => f.to_vec(),
        None => problem.boundary.loads(mesh, &active)?,
    };
    let body = Vector3::from(problem.body_force_kn_m3);
    let mut volume = 0.;
    let mut points = 0;
    for (element, enabled) in active.iter().enumerate() {
        if !enabled {
            continue;
        }
        let nodes = mesh.cell_nodes(element);
        for q in mesh.integration(element)? {
            if !q.weight_m3.is_finite() || q.weight_m3 <= 0. {
                return Err(Error::Input(
                    "positive integration weights required for convexity witness".into(),
                ));
            }
            volume += q.weight_m3;
            points += 1;
            for (&n, shape) in nodes.iter().zip(&q.shape) {
                for axis in 0..3 {
                    forces[3 * n + axis] += shape * q.weight_m3 * body[axis];
                }
            }
        }
    }
    if !volume.is_finite() || volume <= 0. || !forces.iter().all(|x| x.is_finite()) {
        return Err(Error::Input(
            "nonfinite exclusion volume/load assembly".into(),
        ));
    }
    let mut resultant = SVector::<f64, 6>::zeros();
    for (n, force) in forces.as_chunks::<3>().0.iter().enumerate() {
        let f = Vector3::from_column_slice(force);
        let r = (Vector3::from(mesh.nodes[n]) - origin) / length;
        let torque = r.cross(&f);
        for axis in 0..3 {
            resultant[axis] += f[axis];
            resultant[axis + 3] += torque[axis];
        }
    }
    let lambda = a
        .lu()
        .solve(&(-resultant))
        .ok_or_else(|| Error::Input("gauge reaction solve failed".into()))?;
    let residual = (a * lambda + resultant).norm() / resultant.norm().max(1.);
    if !residual.is_finite() || residual > 1e-11 {
        return unavailable(
            "gauge reactions do not resolve force/moment balance accurately enough",
        );
    }
    let mut reactions = Vec::new();
    for (j, (n, direction)) in gauges.iter().enumerate() {
        for axis in 0..3 {
            forces[3 * n + axis] += direction[axis] * lambda[j];
        }
        reactions.push(GaugeReaction {
            node_index: *n,
            unit_direction: (*direction).into(),
            magnitude_kn: lambda[j],
        });
    }
    let mut moment = Matrix3::<f64>::zeros();
    let mut absolute_moment = 0.;
    for (n, force) in forces.as_chunks::<3>().0.iter().enumerate() {
        let r = Vector3::from(mesh.nodes[n]) - origin;
        let f = Vector3::from_column_slice(force);
        moment += r * f.transpose();
        absolute_moment += r.norm() * f.norm();
    }
    let stress = -(moment + moment.transpose()) / (2. * volume);
    if !stress.iter().all(|v| v.is_finite())
        || (moment - moment.transpose()).norm() > 1e-10 * absolute_moment.max(1.)
    {
        return unavailable(
            "load moment is nonfinite or lacks stress symmetry within the numerical guard",
        );
    }
    let mut principal: [f64; 3] = stress.symmetric_eigen().eigenvalues.into();
    principal.sort_by(|a, b| b.total_cmp(a));
    let (cohesion, phi, _) = problem.material.effective_strength();
    let margin = principal[0]
        - principal[2]
        - (principal[0] + principal[2]) * phi.sin()
        - 2. * cohesion * phi.cos();
    let rankine = problem
        .material
        .tension_cutoff_kpa
        .map(|t| -principal[2] - t);
    let scale = principal
        .iter()
        .map(|x| x.abs())
        .fold(cohesion.max(1.), f64::max);
    let clearance =
        1e-8 * scale + 4096. * f64::EPSILON * (1. + condition) * absolute_moment / volume;
    if !margin.is_finite() || rankine.is_some_and(|x| !x.is_finite()) || !clearance.is_finite() {
        return Err(Error::Input(
            "nonfinite mean-stress yield/clearance evaluation".into(),
        ));
    }
    let excluded = margin > clearance || rankine.is_some_and(|x| x > clearance);
    Ok(Assessment::Evaluated {
        certificate: Box::new(Certificate {
            schema: "soil-fem-mc-mean-stress-bound-result/1".into(),
            mesh_sha256: checkpoint::hash_json(mesh)?,
            problem_sha256,
            volume_m3: volume,
            integration_points: points,
            moment_origin_m: origin.into(),
            moment_length_scale_m: length,
            gauge_condition_number: condition,
            gauge_reactions: reactions,
            relative_resultant_force_and_moment_residual: residual,
            mean_effective_stress_compression_positive_kpa: [
                stress[(0, 0)],
                stress[(1, 1)],
                stress[(2, 2)],
                stress[(0, 1)],
                stress[(1, 2)],
                stress[(0, 2)],
            ],
            principal_stresses_descending_kpa: principal,
            mean_mc_yield_margin_kpa: margin,
            mean_rankine_yield_margin_kpa: rankine,
            numerical_clearance_kpa: clearance,
            mean_stress_excludes_admissible_discrete_equilibrium: excluded,
            positive_integration_weight_convexity_argument: true,
            lower_bound_or_stability_proven: false,
            continuum_collapse_factor_certified: false,
            dry_effective_equals_total: true,
        }),
    })
}
