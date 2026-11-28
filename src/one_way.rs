//! Frozen-porosity Richards -> small-strain soil mechanics at physical time.
//! Compression-positive sigma_total = sigma_constitutive + chi*u I.
//! Fluid stress is integrated on the hydraulic quadrature, then moved to the
//! mechanical RHS. The constitutive kernel receives zero additional GP pore
//! pressure. There is no extra Darcy body force or strain-to-storage feedback.
use crate::{
    Error, Result, boundary, checkpoint, elements, ground::Ground, hydrostatic,
    linalg::LinearReport, mesh::Mesh, retention, richards, stages,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Bishop {
    EffectiveSaturation,
    WaterSaturation,
    Constant { value: f64 },
}
impl Bishop {
    fn value(&self, s: &retention::State) -> f64 {
        match self {
            Self::EffectiveSaturation => s.effective_saturation,
            Self::WaterSaturation => s.water_saturation,
            Self::Constant { value } => *value,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WaterWeight {
    WaterEquivalentStorage,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolidWeightReference {
    DryBulkUnitWeight,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PorosityPolicy {
    FixedHydraulicMaterial,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reservoir {
    pub face_set: String,
    pub level_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExteriorPeriod {
    pub name: String,
    pub reservoirs: Vec<Reservoir>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// Fully water-filled aperture at zero mechanical jump, in metres.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_aperture_storage_m: Option<f64>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub joint_reservoir_heads: BTreeMap<String, crate::joint_water::rt0::HeadHistory>,
    pub schema: String,
    pub flow: richards::Case,
    pub mechanics: stages::StageCase,
    pub bishop: Bishop,
    pub water_weight: WaterWeight,
    pub solid_weight_reference: SolidWeightReference,
    pub porosity_policy: PorosityPolicy,
    pub initial_porosity_tolerance: f64,
    pub maximum_initial_engineering_strain: f64,
    pub initial_exterior_water: Vec<Reservoir>,
    pub exterior_periods: Vec<ExteriorPeriod>,
    #[serde(default = "default_store_points")]
    pub store_transfer_integration_points: bool,
    /// Include the face-local integration point index in JOINT trace evidence.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub verification_record_joint_trace_points: bool,
    #[serde(default)]
    pub verification_reject_once_at_times_s: Vec<f64>,
}
fn default_store_points() -> bool {
    true
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub element_index: usize,
    pub local_index: usize,
    pub position_m: [f64; 3],
    pub weight_m3: f64,
    pub pressure_head_m: f64,
    pub physical_pore_pressure_kpa: f64,
    pub bishop_chi: f64,
    pub effective_pore_stress_kpa: f64,
    pub effective_saturation: f64,
    pub water_saturation: f64,
    pub water_content: f64,
    pub water_equivalent_mass_per_volume: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fluid {
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub joint_water_volume_m3: BTreeMap<usize, f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub joint_traces: Vec<crate::joint_water::PointTraces>,
    pub hydraulic_points: Vec<Sample>,
    pub mechanical_points: Vec<Sample>,
    pub fluid_internal_nodal_force_kn: Vec<f64>,
    pub water_body_nodal_load_kn: Vec<f64>,
    pub exterior_nodal_load_kn: Vec<f64>,
    pub physical_nodal_load_kn: Vec<f64>,
    pub equivalent_mechanical_rhs_kn: Vec<f64>,
    pub water_equivalent_mass_volume_m3: f64,
    pub maximum_seed_porosity_difference: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub hydraulic: richards::State,
    /// Kernel state: applied_nodal_load_kn is F_physical - I_fluid, and its
    /// pore_pressure_kpa is zero. Physical u and chi*u are stored in fluid.
    pub mechanical: stages::Committed,
    pub fluid: Fluid,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MechanicalAudit {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_results: Option<crate::joints::Report>,
    pub newton_iterations: usize,
    pub relative_effective_rhs_residual: f64,
    pub maximum_constraint_error_m: f64,
    pub maximum_normalized_yield: f64,
    pub physical_applied_force_kn: [f64; 3],
    pub support_reaction_kn: [f64; 3],
    pub physical_force_balance_error_kn: [f64; 3],
    pub linear_solves: Vec<LinearReport>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub staggered_iterations: Option<usize>,
    pub name: String,
    pub phase: String,
    pub initialization_progress: f64,
    pub state: State,
    pub hydraulic_balance: Option<richards::Balance>,
    pub mechanical_audit: Option<MechanicalAudit>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Attempt {
    pub phase: String,
    pub period: String,
    pub start_time_s: f64,
    pub time_step_s: Option<f64>,
    pub initialization_from: f64,
    pub initialization_to: f64,
    pub accepted: bool,
    pub reason: Option<String>,
    pub hydraulic_candidate_changed: bool,
    pub mechanical_candidate_changed: bool,
    pub plastic_candidate_changed: bool,
    pub committed_pair_sha256_before: String,
    pub committed_pair_sha256_after: String,
    pub hydraulic_iterations: Vec<richards::Iteration>,
    pub mechanical_audit: Option<MechanicalAudit>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Runtime {
    pub initialization_progress: f64,
    pub initialization_increment: f64,
    pub seed_internal_kn: Vec<f64>,
    pub period_index: usize,
    pub next_time_step_s: f64,
    pub injected_times: BTreeSet<usize>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub complete: bool,
    pub failure_reason: Option<String>,
    pub mesh_sha256: String,
    pub ground_sha256: String,
    pub case_sha256: String,
    pub solver_binary_sha256: String,
    pub ground: Ground,
    pub case: Case,
    pub mechanical_integration_points: Vec<elements::PointInfo>,
    pub hydraulic_faces: Vec<richards::Face>,
    pub seed: State,
    pub initial_reference_displacement_m: Vec<[f64; 3]>,
    pub frames: Vec<Frame>,
    pub attempts: Vec<Attempt>,
    pub last_committed: State,
    pub runtime: Option<Runtime>,
}

struct Point {
    sample: elements::PointInfo,
    shape: Vec<f64>,
    gradient: Vec<[f64; 3]>,
}
pub(crate) struct Transfer<'a> {
    mesh: &'a Mesh,
    case: &'a Case,
    ground: &'a Ground,
    points: Vec<Point>,
    mechanical: Vec<elements::PointInfo>,
    laws: BTreeMap<u32, retention::VanGenuchten>,
    base: Vec<f64>,
    seed_porosity: Vec<f64>,
    material_ids: Vec<u32>,
}
impl<'a> Transfer<'a> {
    pub(crate) fn new(
        mesh: &'a Mesh,
        ground: &'a Ground,
        case: &'a Case,
        seed: &stages::Committed,
        mechanical: Vec<elements::PointInfo>,
    ) -> Result<Self> {
        let laws = case
            .flow
            .materials
            .iter()
            .map(|m| (m.id, m.retention))
            .collect();
        let mut points = Vec::new();
        for e in 0..mesh.tetrahedra.len() {
            let xyz = mesh
                .cell_nodes(e)
                .iter()
                .map(|&n| mesh.nodes[n])
                .collect::<Vec<_>>();
            for (j, &(r, w)) in elements::volume_rule(mesh.quadrature_order.unwrap_or(3))?
                .iter()
                .enumerate()
            {
                let (position, jac) = elements::tet_mapping(&xyz, r)?;
                let shape = if xyz.len() == 4 {
                    vec![1. - r.iter().sum::<f64>(), r[0], r[1], r[2]]
                } else {
                    elements::tet10_shape(r).0.to_vec()
                };
                points.push(Point {
                    sample: elements::PointInfo {
                        element_index: e,
                        local_index: j,
                        natural_coordinates: r,
                        position_m: position,
                        weight_m3: w * jac.determinant(),
                        shape_functions: shape.clone(),
                    },
                    shape,
                    gradient: elements::scalar_gradient(&xyz, r)?,
                });
            }
        }
        Ok(Self {
            mesh,
            case,
            ground,
            points,
            mechanical,
            laws,
            base: seed.applied_nodal_load_kn.clone(),
            seed_porosity: seed
                .material_states
                .iter()
                .map(|s| s.void_ratio() / (1. + s.void_ratio()))
                .collect(),
            material_ids: if seed.element_material_ids.is_empty() {
                mesh.material_ids.clone()
            } else {
                seed.element_material_ids.clone()
            },
        })
    }
    pub(crate) fn set_material_ids(&mut self, ids: &[u32]) -> Result<()> {
        if ids.len() != self.mesh.tetrahedra.len()
            || ids.iter().any(|id| !self.laws.contains_key(id))
        {
            return Err(Error::Input("missing lifecycle hydraulic material".into()));
        }
        self.material_ids = ids.to_vec();
        Ok(())
    }
    pub(crate) fn set_porosity_reference(&mut self, values: &[f64]) -> Result<()> {
        if values.len() != self.mechanical.len()
            || values
                .iter()
                .any(|p| !p.is_finite() || *p <= 0. || *p >= 1.)
        {
            return Err(Error::Input("invalid hydraulic porosity reference".into()));
        }
        self.seed_porosity = values.to_vec();
        Ok(())
    }
    fn sample(&self, p: &elements::PointInfo, state: &richards::State) -> Result<Sample> {
        let law = self
            .laws
            .get(&self.material_ids[p.element_index])
            .ok_or_else(|| Error::Input("ONE_WAY hydraulic material missing".into()))?;
        let h = state.relative_cell_head_m[p.element_index]
            .ok_or_else(|| Error::Input("missing active hydraulic head".into()))?;
        let psi = (self.case.flow.head_datum_m - p.position_m[2]) + h;
        let s = law.evaluate(psi)?;
        let u = psi * self.ground.initialization.water_unit_weight_kn_m3;
        let chi = self.case.bishop.value(&s);
        Ok(Sample {
            element_index: p.element_index,
            local_index: p.local_index,
            position_m: p.position_m,
            weight_m3: p.weight_m3,
            pressure_head_m: psi,
            physical_pore_pressure_kpa: u,
            bishop_chi: chi,
            effective_pore_stress_kpa: chi * u,
            effective_saturation: s.effective_saturation,
            water_saturation: s.water_saturation,
            water_content: s.water_content,
            water_equivalent_mass_per_volume: s.water_content
                + law.saturated_specific_storage_per_m * psi.max(0.),
        })
    }
    fn exterior(&self, reservoirs: &[Reservoir]) -> Result<Vec<f64>> {
        self.exterior_active(reservoirs, None)
    }
    fn exterior_active(
        &self,
        reservoirs: &[Reservoir],
        active_elements: Option<&[bool]>,
    ) -> Result<Vec<f64>> {
        let owners = boundary::surface_owners(self.mesh);
        let mut used = BTreeSet::new();
        let mut load = vec![0.; self.base.len()];
        for reservoir in reservoirs {
            let faces = self
                .mesh
                .face_sets
                .get(&reservoir.face_set)
                .ok_or_else(|| Error::Input("unknown exterior water face set".into()))?;
            if faces.is_empty() || !reservoir.level_m.is_finite() {
                return Err(Error::Input("empty/invalid exterior water boundary".into()));
            }
            for face in faces {
                let mut key = *face;
                key.sort_unstable();
                let active = owners
                    .get(&key)
                    .into_iter()
                    .flatten()
                    .copied()
                    .filter(|&e| active_elements.is_none_or(|s| s[e]))
                    .collect::<Vec<_>>();
                if active.len() != 1 {
                    return Err(Error::Input("exterior water requires exposed faces".into()));
                }
                let owner = active[0];
                if !used.insert(key) {
                    return Err(Error::Input("overlapping exterior water faces".into()));
                }
                let nodes = self.mesh.face_nodes(owner, &key)?;
                let xyz = nodes
                    .iter()
                    .map(|&n| self.mesh.nodes[n])
                    .collect::<Vec<_>>();
                let force = hydrostatic::face_load(
                    &xyz,
                    reservoir.level_m,
                    self.ground.initialization.water_unit_weight_kn_m3,
                    Default::default(),
                )?;
                for (j, n) in nodes.iter().enumerate() {
                    for a in 0..3 {
                        load[3 * n + a] += force.nodal_force_kn[j][a];
                    }
                }
            }
        }
        Ok(load)
    }
    pub(crate) fn exterior_on_domain(
        &self,
        reservoirs: &[Reservoir],
        active: &[bool],
    ) -> Result<Vec<f64>> {
        self.exterior_active(reservoirs, Some(active))
    }
    fn field(&self, state: &richards::State, reservoirs: &[Reservoir]) -> Result<Fluid> {
        self.field_with_base(state, reservoirs, &self.base)
    }
    pub(crate) fn field_with_base(
        &self,
        state: &richards::State,
        reservoirs: &[Reservoir],
        base: &[f64],
    ) -> Result<Fluid> {
        if base.len() != self.base.len() {
            return Err(Error::Input("invalid domain load length".into()));
        }
        let mut internal = vec![0.; self.base.len()];
        let mut body = internal.clone();
        let mut raw = Vec::new();
        let mut volume = 0.;
        for p in &self.points {
            if state.relative_cell_head_m[p.sample.element_index].is_none() {
                continue;
            }
            let s = self.sample(&p.sample, state)?;
            let w = p.sample.weight_m3;
            let nodes = self.mesh.cell_nodes(s.element_index);
            for (j, n) in nodes.iter().enumerate() {
                body[3 * n + 2] -= w
                    * p.shape[j]
                    * s.water_equivalent_mass_per_volume
                    * self.ground.initialization.water_unit_weight_kn_m3;
                for a in 0..3 {
                    internal[3 * n + a] -= w * p.gradient[j][a] * s.effective_pore_stress_kpa;
                }
            }
            volume += w * s.water_equivalent_mass_per_volume;
            if self.case.store_transfer_integration_points {
                raw.push(s);
            }
        }
        let mechanical = if self.case.store_transfer_integration_points {
            self.mechanical
                .iter()
                .filter(|p| state.relative_cell_head_m[p.element_index].is_some())
                .map(|p| self.sample(p, state))
                .collect::<Result<Vec<_>>>()?
        } else {
            vec![]
        };
        let mut porosity_difference: f64 = 0.;
        for (q, p) in self.mechanical.iter().zip(&self.seed_porosity) {
            if state.relative_cell_head_m[q.element_index].is_none() {
                continue;
            }
            porosity_difference = porosity_difference.max(
                (self.laws[&self.material_ids[q.element_index]].saturated_water_content - p).abs(),
            );
        }
        if porosity_difference > self.case.initial_porosity_tolerance {
            return Err(Error::Input(format!(
                "initial hydraulic porosity differs from constitutive seed by {porosity_difference}; exceeds explicit tolerance"
            )));
        }
        let active = state
            .relative_cell_head_m
            .iter()
            .map(Option::is_some)
            .collect::<Vec<_>>();
        let exterior = self.exterior_active(reservoirs, Some(&active))?;
        let physical = base
            .iter()
            .zip(&body)
            .zip(&exterior)
            .map(|((a, b), c)| a + b + c)
            .collect::<Vec<_>>();
        let rhs = physical
            .iter()
            .zip(&internal)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        if !rhs
            .iter()
            .chain(&physical)
            .chain(&body)
            .chain(&internal)
            .all(|v| v.is_finite())
            || !volume.is_finite()
        {
            return Err(Error::Solve("non-finite ONE_WAY transfer".into()));
        }
        Ok(Fluid {
            joint_water_volume_m3: BTreeMap::new(),
            joint_traces: vec![],
            hydraulic_points: raw,
            mechanical_points: mechanical,
            fluid_internal_nodal_force_kn: internal,
            water_body_nodal_load_kn: body,
            exterior_nodal_load_kn: exterior,
            physical_nodal_load_kn: physical,
            equivalent_mechanical_rhs_kn: rhs,
            water_equivalent_mass_volume_m3: volume,
            maximum_seed_porosity_difference: porosity_difference,
        })
    }
}

impl Fluid {
    pub(crate) fn apply_joint(
        &mut self,
        wet: Option<&crate::joint_water::Model>,
        active: &[bool],
        state: &richards::State,
        faces: &[richards::Face],
        datum: f64,
        time: f64,
    ) -> Result<Option<crate::joint_water::Field>> {
        let Some(wet) = wet else { return Ok(None) };
        let (field, traces) = wet.richards_field(active, state, faces, datum, time)?;
        for ((physical, rhs), wall) in self
            .physical_nodal_load_kn
            .iter_mut()
            .zip(&mut self.equivalent_mechanical_rhs_kn)
            .zip(&field.nodal_load_kn)
        {
            *physical += wall;
            *rhs += wall;
        }
        self.joint_traces = traces;
        Ok(Some(field))
    }
}

pub(crate) fn audit(t: &stages::ExternalTrial, fluid: &Fluid) -> MechanicalAudit {
    let mut physical = [0.; 3];
    for (g, f) in fluid.physical_nodal_load_kn.iter().enumerate() {
        physical[g % 3] += f;
    }
    MechanicalAudit {
        joint_results: t.frame.joint_results.clone(),
        newton_iterations: t.iterations,
        relative_effective_rhs_residual: t.frame.relative_force_residual,
        maximum_constraint_error_m: t.frame.maximum_constraint_error_m,
        maximum_normalized_yield: t.frame.maximum_normalized_yield,
        physical_applied_force_kn: physical,
        support_reaction_kn: t.frame.support_reaction_kn,
        physical_force_balance_error_kn: std::array::from_fn(|a| {
            physical[a] + t.frame.support_reaction_kn[a]
        }),
        linear_solves: t.linear_solves.clone(),
    }
}
fn changed(a: &impl Serialize, b: &impl Serialize) -> Result<bool> {
    Ok(checkpoint::hash_json(a)? != checkpoint::hash_json(b)?)
}
pub(crate) fn validate(case: &Case) -> Result<&[richards::Period]> {
    let richards::Analysis::Transient { periods } = &case.flow.analysis else {
        return Err(Error::Input(
            "ONE_WAY needs transient Richards periods".into(),
        ));
    };
    if let Some(aperture) = case.joint_aperture_storage_m
        && (!aperture.is_finite()
            || aperture <= 0.
            || case.schema != "soil-fem-staggered-transfer/1"
            || case.mechanics.joint_system.as_ref().is_none_or(|c| {
                c.joints.is_empty()
                    || c.joints.iter().any(|d| {
                        !matches!(d.water, Some(crate::joint_water::Policy::Continuous))
                            || d.birth_policy.is_some()
                    })
            }))
    {
        return Err(Error::Input("joint aperture storage requires positive reference aperture and active CONTINUOUS staggered joints".into()));
    }
    if !(case.schema == "soil-fem-one-way/1"
        || (case.schema == "soil-fem-staggered-transfer/1"
            && case.mechanics.joint_system.is_some()))
        || !case.initial_porosity_tolerance.is_finite()
        || !(1e-12..=0.2).contains(&case.initial_porosity_tolerance)
        || !case.maximum_initial_engineering_strain.is_finite()
        || !(1e-8..=0.1).contains(&case.maximum_initial_engineering_strain)
        || matches!(case.bishop,Bishop::Constant{value} if !value.is_finite() || !(0. ..=1.).contains(&value))
        || case.exterior_periods.len() != periods.len()
        || case
            .exterior_periods
            .iter()
            .zip(periods)
            .any(|(a, b)| a.name != b.name)
        || case
            .verification_reject_once_at_times_s
            .iter()
            .any(|t| !t.is_finite() || *t <= 0. || periods.last().is_none_or(|p| *t > p.end_time_s))
        || case
            .verification_reject_once_at_times_s
            .iter()
            .map(|t| t.to_bits())
            .collect::<BTreeSet<_>>()
            .len()
            != case.verification_reject_once_at_times_s.len()
    {
        return Err(Error::Input(
            "invalid ONE_WAY closure, initialization or exterior history".into(),
        ));
    }
    for (id, h) in &case.joint_reservoir_heads {
        if !case.mechanics.joint_system.as_ref().is_some_and(|j| {
            j.joints.iter().any(|d| {
                d.id == *id
                    && matches!(
                        d.water,
                        Some(crate::joint_water::Policy::Prescribed { level_m: None })
                    )
            })
        }) {
            return Err(Error::Input(
                "head history needs a declared PRESCRIBED wall without a fixed level".into(),
            ));
        }
        h.validate(
            periods
                .last()
                .ok_or_else(|| Error::Input("missing hydraulic period".into()))?
                .end_time_s,
        )?;
    }
    Ok(periods)
}

impl Case {
    pub(crate) fn clip_dt(&self, t: f64, dt: f64) -> f64 {
        self.joint_reservoir_heads
            .values()
            .flat_map(|h| &h.0)
            .map(|k| k[0])
            .filter(|s| *s > t)
            .map(|s| s - t)
            .fold(dt, f64::min)
    }
}

pub fn solve(mesh: &Mesh, ground: Ground, case: Case) -> Result<Report> {
    advance(mesh, ground, case, None, None, false, None)
}
pub(crate) fn initialize(mesh: &Mesh, ground: Ground, case: Case) -> Result<Report> {
    advance(mesh, ground, case, None, None, true, None)
}
pub(crate) fn initialize_on_domain(
    mesh: &Mesh,
    ground: Ground,
    case: Case,
    active: &[bool],
) -> Result<Report> {
    advance(mesh, ground, case, None, None, true, Some(active))
}
pub fn solve_until(mesh: &Mesh, ground: Ground, case: Case, count: usize) -> Result<Report> {
    if count == 0 {
        return Err(Error::Input("checkpoint count must be positive".into()));
    }
    advance(mesh, ground, case, None, Some(count), false, None)
}
pub fn resume(
    mesh: &Mesh,
    ground: Ground,
    case: Case,
    saved: Report,
    additional: Option<usize>,
) -> Result<Report> {
    if saved.complete || saved.failure_reason.is_some() || saved.runtime.is_none() {
        return Err(Error::Input(
            "ONE_WAY resume needs an accepted intermediate state".into(),
        ));
    }
    let count = saved.attempts.iter().filter(|a| a.accepted).count();
    let canonical = solve_until(mesh, ground.clone(), case.clone(), count)?;
    if normalized(&canonical)? != normalized(&saved)? {
        return Err(Error::Input(
            "ONE_WAY checkpoint inputs, executable, fields, runtime or history changed".into(),
        ));
    }
    drop(canonical);
    advance(mesh, ground, case, Some(saved), additional, false, None)
}
pub fn normalized(run: &Report) -> Result<serde_json::Value> {
    fn strip(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Object(m) => {
                for (k, v) in m {
                    if k == "elapsed_seconds" {
                        *v = serde_json::json!(0.);
                    } else {
                        strip(v);
                    }
                }
            }
            serde_json::Value::Array(a) => {
                for v in a {
                    strip(v);
                }
            }
            _ => {}
        }
    }
    let mut v = serde_json::to_value(run)?;
    strip(&mut v);
    Ok(v)
}

fn advance(
    mesh: &Mesh,
    ground: Ground,
    case: Case,
    saved: Option<Report>,
    limit: Option<usize>,
    initialization_only: bool,
    initial_active: Option<&[bool]>,
) -> Result<Report> {
    if limit == Some(0) {
        return Err(Error::Input("checkpoint count must be positive".into()));
    }
    let periods = validate(&case)?;
    let mut kernel = stages::ExternalMechanics::new(mesh, &ground, &case.mechanics)?;
    kernel.set_head_histories(&case.joint_reservoir_heads);
    kernel.set_joint_trace_point_recording(case.verification_record_joint_trace_points);
    let mut seed_mechanics = if let Some(active) = initial_active {
        kernel.seed_on_domain(active)?
    } else {
        kernel.seed()?
    };
    let active = seed_mechanics.active_elements.clone();
    let flow = if let Some(wet) = kernel.wet_joints() {
        let pairs = wet.richards_pairs(&active)?;
        let walls: BTreeSet<_> = pairs
            .iter()
            .flat_map(|p| [p.a.corners, p.b.corners])
            .collect();
        for reservoir in case
            .initial_exterior_water
            .iter()
            .chain(case.exterior_periods.iter().flat_map(|p| &p.reservoirs))
        {
            if mesh
                .face_sets
                .get(&reservoir.face_set)
                .into_iter()
                .flatten()
                .any(|f| {
                    let mut key = *f;
                    key.sort_unstable();
                    walls.contains(&key)
                })
            {
                return Err(Error::Input(
                    "same wall registered for exterior and JOINT pressure".into(),
                ));
            }
        }
        richards::Stepper::initial_with_joints(
            mesh,
            &case.flow,
            &wet.soil_domain(&active)?,
            &pairs,
        )?
    } else if let Some(active) = initial_active {
        richards::Stepper::initial_on_domain(mesh, &case.flow, active)?
    } else {
        richards::Stepper::new(mesh, &case.flow)?
    };
    let layout = kernel.layout();
    let transfer = Transfer::new(mesh, &ground, &case, &seed_mechanics, layout.clone())?;
    // Validate every exterior period before producing or advancing any state.
    transfer.exterior(&case.initial_exterior_water)?;
    for p in &case.exterior_periods {
        transfer.exterior(&p.reservoirs)?;
    }
    let initial = flow.initial()?;
    let faces = flow.faces();
    let mut fluid = transfer.field(&initial.state, &case.initial_exterior_water)?;
    seed_mechanics.joint_water = fluid.apply_joint(
        kernel.wet_joints(),
        &active,
        &initial.state,
        &faces,
        case.flow.head_datum_m,
        0.,
    )?;
    let volumes = |wet: Option<&crate::joint_water::Model>,
                   mechanics: &stages::Committed,
                   fluid: &Fluid|
     -> Result<BTreeMap<usize, f64>> {
        match case.joint_aperture_storage_m {
            Some(b0) => wet
                .ok_or_else(|| Error::Input("missing wet joints".into()))?
                .saturated_volumes(
                    b0,
                    mechanics
                        .joints
                        .as_ref()
                        .ok_or_else(|| Error::Input("missing joint state".into()))?,
                    &fluid.joint_traces,
                ),
            None => Ok(BTreeMap::new()),
        }
    };
    fluid.joint_water_volume_m3 = volumes(kernel.wet_joints(), &seed_mechanics, &fluid)?;
    let internal = kernel.internal(&seed_mechanics)?;
    let mut initial_mechanics = seed_mechanics;
    initial_mechanics.applied_nodal_load_kn = internal.clone();
    let initial_state = State {
        hydraulic: initial.state,
        mechanical: initial_mechanics,
        fluid,
    };
    let mut report = if let Some(saved) = saved {
        saved
    } else {
        Report {
            schema: "soil-fem-one-way-result/1".into(),
            complete: false,
            failure_reason: None,
            mesh_sha256: checkpoint::hash_json(mesh)?,
            ground_sha256: checkpoint::hash_json(&ground)?,
            case_sha256: checkpoint::hash_json(&case)?,
            solver_binary_sha256: checkpoint::hash_file(&std::env::current_exe()?)?,
            ground: ground.clone(),
            case: case.clone(),
            mechanical_integration_points: layout,
            hydraulic_faces: flow.faces(),
            seed: initial_state.clone(),
            initial_reference_displacement_m: vec![],
            frames: vec![Frame {
                staggered_iterations: None,
                name: "seed".into(),
                phase: "initialization_seed".into(),
                initialization_progress: 0.,
                state: initial_state.clone(),
                hydraulic_balance: None,
                mechanical_audit: None,
            }],
            attempts: vec![],
            last_committed: initial_state,
            runtime: Some(Runtime {
                initialization_progress: 0.,
                initialization_increment: case.mechanics.settings.initial_increment,
                seed_internal_kn: internal,
                period_index: 0,
                next_time_step_s: case.flow.settings.initial_time_step_s,
                injected_times: BTreeSet::new(),
            }),
        }
    };
    let mut rt = report
        .runtime
        .take()
        .ok_or_else(|| Error::Input("ONE_WAY runtime missing".into()))?;
    kernel.set_reference(report.initial_reference_displacement_m.clone());
    let mut accepted = 0;
    let initial_strain_limit = ground.initialization.geostatic.as_ref().map_or(
        case.maximum_initial_engineering_strain,
        |cfg| {
            case.maximum_initial_engineering_strain
                .min(cfg.maximum_initial_engineering_strain)
        },
    );
    while rt.initialization_progress < 1. || rt.period_index < periods.len() {
        if report.attempts.len() >= case.flow.settings.maximum_steps {
            report.failure_reason = Some("ONE_WAY attempt budget exhausted".into());
            break;
        }
        let initializing = rt.initialization_progress < 1.;
        let old = &report.last_committed;
        let hash = checkpoint::hash_json(old)?;
        let from = rt.initialization_progress;
        let to = if initializing {
            (from + rt.initialization_increment).min(1.)
        } else {
            1.
        };
        let dt = if initializing {
            None
        } else {
            Some(
                case.clip_dt(old.hydraulic.time_s, rt.next_time_step_s)
                    .min(periods[rt.period_index].end_time_s - old.hydraulic.time_s),
            )
        };
        let name = if initializing {
            "initial_equilibrium".to_string()
        } else {
            periods[rt.period_index].name.clone()
        };
        let mut candidate = None;
        let mut hydraulic_log = vec![];
        let mut mechanical_audit = None;
        let mut hchanged = false;
        let mut mchanged = false;
        let mut plastic = false;
        let mut invalid = false;
        let result = (|| -> Result<()> {
            let staggered = !initializing && case.schema == "soil-fem-staggered-transfer/1";
            let mut previous = None;
            let mut rates = BTreeMap::new();
            for iteration in 1..=if staggered { 2 } else { 1 } {
                let mut hydraulic = if let Some(dt) = dt {
                    let (trial, log) =
                        flow.trial_with_storage(&old.hydraulic, rt.period_index, dt, &rates);
                    hydraulic_log.extend(log);
                    trial?
                } else {
                    flow.initial()?
                };
                hchanged = changed(&old.hydraulic, &hydraulic.state)?;
                let reservoirs = if initializing {
                    &case.initial_exterior_water
                } else {
                    &case.exterior_periods[rt.period_index].reservoirs
                };
                let mut fluid = transfer.field(&hydraulic.state, reservoirs)?;
                let water = fluid.apply_joint(
                    kernel.wet_joints(),
                    &active,
                    &hydraulic.state,
                    &faces,
                    case.flow.head_datum_m,
                    old.hydraulic.time_s + dt.unwrap_or(0.),
                )?;
                let rhs = if initializing {
                    rt.seed_internal_kn
                        .iter()
                        .zip(&fluid.equivalent_mechanical_rhs_kn)
                        .map(|(a, b)| if to == 1. { *b } else { a + (b - a) * to })
                        .collect()
                } else {
                    fluid.equivalent_mechanical_rhs_kn.clone()
                };
                let trial = kernel.trial(&old.mechanical, rhs, name.clone(), water)?;
                mchanged = changed(&old.mechanical, &trial.frame.state)?;
                plastic = old
                    .mechanical
                    .material_states
                    .iter()
                    .zip(&trial.frame.state.material_states)
                    .any(|(a, b)| a.plastic_strain() != b.plastic_strain());
                let audit = audit(&trial, &fluid);
                mechanical_audit = Some(audit.clone());
                if initializing && kernel.maximum_strain(&trial.frame.state) > initial_strain_limit
                {
                    return Err(Error::Input("initial ONE_WAY engineering strain exceeds declared limit; no equilibrium accepted".into()));
                }
                fluid.joint_water_volume_m3 =
                    volumes(kernel.wet_joints(), &trial.frame.state, &fluid)?;
                if let Some(seconds) = dt
                    && case.joint_aperture_storage_m.is_some()
                {
                    rates = fluid
                        .joint_water_volume_m3
                        .iter()
                        .map(|(&fi, &v)| (fi, (v - old.fluid.joint_water_volume_m3[&fi]) / seconds))
                        .collect();
                    if iteration == 1 {
                        continue;
                    }
                    // Check the original water equations with the actual endpoint
                    // opening, rather than accepting the previous exchange iterate.
                    hydraulic.balance = flow.balance_with_storage(
                        &hydraulic.state,
                        &old.hydraulic,
                        rt.period_index,
                        seconds,
                        &rates,
                    )?;
                    if hydraulic.balance.maximum_scaled_residual > 1. {
                        return Err(Error::Solve(
                            "joint aperture storage exchange residual exceeds Richards criterion"
                                .into(),
                        ));
                    }
                }
                let state = State {
                    hydraulic: hydraulic.state,
                    mechanical: trial.frame.state,
                    fluid,
                };
                if staggered && case.joint_aperture_storage_m.is_none() {
                    // With no reverse constitutive law the exchange operator is constant.
                    // Re-solve from n and require exactly the same final state (including
                    // the two subsolver fields and contact history) before committing.
                    let exchange = checkpoint::hash_json(&state)?;
                    if iteration == 1 {
                        previous = Some(exchange);
                        continue;
                    }
                    if previous.as_ref() != Some(&exchange) {
                        return Err(Error::Solve(
                            "staggered transfer fields did not settle".into(),
                        ));
                    }
                }
                candidate = Some(Frame {
                    staggered_iterations: staggered.then_some(iteration),
                    name: name.clone(),
                    phase: if initializing {
                        "initialization"
                    } else {
                        "transient"
                    }
                    .into(),
                    initialization_progress: to,
                    state,
                    hydraulic_balance: (!initializing).then_some(hydraulic.balance),
                    mechanical_audit: Some(audit),
                });
            }
            if !initializing {
                for (i, t) in case.verification_reject_once_at_times_s.iter().enumerate() {
                    if !rt.injected_times.contains(&i)
                        && candidate
                            .as_ref()
                            .is_some_and(|f| f.state.hydraulic.time_s == *t)
                    {
                        rt.injected_times.insert(i);
                        return Err(Error::Solve("verification injection after both hydraulic and mechanical trials changed; restoring pair".into()));
                    }
                }
            }
            Ok(())
        })();
        let reason = result.as_ref().err().map(ToString::to_string);
        if matches!(result, Err(Error::Input(_))) {
            invalid = true;
        }
        let success = result.is_ok();
        let start_time = old.hydraulic.time_s;
        if success {
            let frame = candidate.ok_or_else(|| Error::Solve("missing coupled trial".into()))?;
            report.last_committed = frame.state.clone();
            report.frames.push(frame);
            accepted += 1;
            if initializing {
                rt.initialization_progress = to;
                rt.initialization_increment = (rt.initialization_increment * 1.5)
                    .min(case.mechanics.settings.initial_increment);
                if to == 1. {
                    report.initial_reference_displacement_m =
                        report.last_committed.mechanical.displacement_m.clone();
                    kernel.set_reference(report.initial_reference_displacement_m.clone());
                }
            } else {
                if report.last_committed.hydraulic.time_s >= periods[rt.period_index].end_time_s {
                    rt.period_index += 1;
                }
                rt.next_time_step_s =
                    (rt.next_time_step_s * 1.5).min(case.flow.settings.maximum_time_step_s);
            }
        }
        let after = checkpoint::hash_json(&report.last_committed)?;
        if !success && after != hash {
            return Err(Error::Solve(
                "rejected ONE_WAY trial changed committed pair".into(),
            ));
        }
        report.attempts.push(Attempt {
            phase: if initializing {
                "initialization"
            } else {
                "transient"
            }
            .into(),
            period: name,
            start_time_s: start_time,
            time_step_s: dt,
            initialization_from: from,
            initialization_to: to,
            accepted: success,
            reason: reason.clone(),
            hydraulic_candidate_changed: hchanged,
            mechanical_candidate_changed: mchanged,
            plastic_candidate_changed: plastic,
            committed_pair_sha256_before: hash,
            committed_pair_sha256_after: after,
            hydraulic_iterations: hydraulic_log,
            mechanical_audit,
        });
        if success {
            if limit == Some(accepted) || (initialization_only && rt.initialization_progress == 1.)
            {
                break;
            }
        } else if initializing
            && !invalid
            && rt.initialization_increment / 2. >= case.mechanics.settings.minimum_increment
        {
            rt.initialization_increment /= 2.;
        } else if !initializing
            && !invalid
            && dt.is_some_and(|d| d / 2. >= case.flow.settings.minimum_time_step_s)
        {
            rt.next_time_step_s = dt.unwrap_or(0.) / 2.;
        } else {
            report.failure_reason = reason;
            break;
        }
    }
    report.complete = report.failure_reason.is_none()
        && rt.initialization_progress == 1.
        && rt.period_index == periods.len();
    report.runtime = (!report.complete && report.failure_reason.is_none()).then_some(rt);
    Ok(report)
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Checkpoint {
    pub schema: String,
    pub run_sha256: String,
    pub run: Report,
}
impl Checkpoint {
    pub fn seal(run: Report) -> Result<Self> {
        if run.complete
            || run.failure_reason.is_some()
            || run.runtime.is_none()
            || !run.attempts.iter().any(|a| a.accepted)
        {
            return Err(Error::Input(
                "ONE_WAY checkpoint requires an accepted intermediate state".into(),
            ));
        }
        Ok(Self {
            schema: "soil-fem-one-way-checkpoint/1".into(),
            run_sha256: checkpoint::hash_json(&run)?,
            run,
        })
    }
    pub fn resume(
        self,
        mesh: &Mesh,
        ground: Ground,
        case: Case,
        additional: Option<usize>,
    ) -> Result<Report> {
        if self.schema != "soil-fem-one-way-checkpoint/1"
            || self.run_sha256 != checkpoint::hash_json(&self.run)?
        {
            return Err(Error::Input(
                "ONE_WAY checkpoint payload hash changed".into(),
            ));
        }
        resume(mesh, ground, case, self.run, additional)
    }
}
