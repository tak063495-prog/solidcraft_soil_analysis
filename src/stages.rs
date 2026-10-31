//! Small-strain FEM with hydrostatic or steady ONE_WAY pore pressure.
//! CSR GMRES/ILU(0) and an
//! independent bounded dense LU backend share the same material-force Newton.
use crate::{
    Error, Result,
    boundary::{Conditions, Transform},
    elements::{self, PointInfo, VolumePoint},
    geostatic,
    ground::{self, Ground},
    linalg::{self, CsrMatrix, CsrPattern, LinearReport, LinearSolver},
    mesh::Mesh,
    soil,
};
use nalgebra::{DMatrix, DVector, SVector, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::{sync::Arc, time::Instant};

#[path = "stages_mixed.rs"]
mod mixed;
pub(crate) mod srm_adapter;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Kinematics {
    SideRollers,
    Oedometer,
    Explicit,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Settings {
    pub kinematics: Kinematics,
    pub relative_force_tolerance: f64,
    pub relative_linear_tolerance: f64,
    pub maximum_newton_iterations: usize,
    pub initial_increment: f64,
    pub minimum_increment: f64,
    pub material_maximum_substep: f64,
    /// Numerical override for every MCC material; other constitutive laws retain
    /// their own integration. Ground input and its physical parameters stay intact.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mcc_integration: Option<crate::mcc::Integration>,
    pub tangent_perturbation: f64,
    pub maximum_free_dofs: usize,
    #[serde(default, skip_serializing_if = "LinearSolver::is_dense")]
    pub linear_solver: LinearSolver,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Stage {
    /// ONE_WAY: solve steady flow on the current active geometry, then apply
    /// its pore pressure to mechanics. Exterior water is specified separately.
    SteadySeepage {
        name: String,
        flow: crate::darcy::Case,
        exterior_water_level_m: f64,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    HydrostaticWater {
        name: String,
        target_level_m: f64,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    InstallAnchors {
        name: String,
        anchors: Vec<String>,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    RemoveAnchors {
        name: String,
        anchors: Vec<String>,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    JackAnchors {
        name: String,
        anchors: Vec<String>,
        target_force_kn: f64,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    #[serde(alias = "install_structures")]
    InstallBeams {
        name: String,
        members: Vec<String>,
        placement: BirthPlacement,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    #[serde(alias = "remove_structures")]
    RemoveBeams {
        name: String,
        members: Vec<String>,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    ResetDisplay {
        name: String,
    },
    ResetState {
        name: String,
        bodies: Vec<String>,
        states: BTreeMap<String, soil::State>,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    Birth {
        name: String,
        bodies: Vec<String>,
        states: BTreeMap<String, soil::State>,
        placement: BirthPlacement,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    MaterialSet {
        name: String,
        bodies: Vec<String>,
        material_id: u32,
        history: MaterialHistory,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    Boundary {
        name: String,
        target_factor: f64,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    Pressure {
        name: String,
        target_pressure_kpa: f64,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
    Excavate {
        name: String,
        bodies: Vec<String>,
        remaining_surface_pressure_kpa: f64,
        #[serde(default)]
        initial_increment: Option<f64>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub enum MaterialHistory {
    Preserve,
    Replace {
        states: BTreeMap<String, soil::State>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BirthPlacement {
    CurrentDisplacement,
    CadNewNodes,
}
impl Stage {
    fn is_transition(&self) -> bool {
        matches!(
            self,
            Self::Birth { .. }
                | Self::MaterialSet { .. }
                | Self::ResetState { .. }
                | Self::InstallBeams { .. }
                | Self::RemoveBeams { .. }
                | Self::InstallAnchors { .. }
                | Self::RemoveAnchors { .. }
        )
    }
    pub(crate) fn name(&self) -> &str {
        match self {
            Self::Pressure { name, .. }
            | Self::SteadySeepage { name, .. }
            | Self::HydrostaticWater { name, .. }
            | Self::Excavate { name, .. }
            | Self::Birth { name, .. }
            | Self::ResetState { name, .. }
            | Self::ResetDisplay { name }
            | Self::MaterialSet { name, .. }
            | Self::Boundary { name, .. } => name,
            Self::InstallBeams { name, .. } | Self::RemoveBeams { name, .. } => name,
            Self::InstallAnchors { name, .. }
            | Self::RemoveAnchors { name, .. }
            | Self::JackAnchors { name, .. } => name,
        }
    }
    fn step(&self, settings: &Settings) -> f64 {
        match self {
            Self::ResetDisplay { .. } => 1.,
            Self::Pressure {
                initial_increment, ..
            }
            | Self::HydrostaticWater {
                initial_increment, ..
            }
            | Self::SteadySeepage {
                initial_increment, ..
            }
            | Self::Excavate {
                initial_increment, ..
            }
            | Self::Boundary {
                initial_increment, ..
            }
            | Self::Birth {
                initial_increment, ..
            }
            | Self::MaterialSet {
                initial_increment, ..
            }
            | Self::ResetState {
                initial_increment, ..
            }
            | Self::InstallBeams {
                initial_increment, ..
            }
            | Self::RemoveBeams {
                initial_increment, ..
            }
            | Self::InstallAnchors {
                initial_increment, ..
            }
            | Self::RemoveAnchors {
                initial_increment, ..
            }
            | Self::JackAnchors {
                initial_increment, ..
            } => initial_increment.unwrap_or(settings.initial_increment),
        }
    }
    fn pressure(&self, previous: f64) -> f64 {
        match self {
            Self::Pressure {
                target_pressure_kpa,
                ..
            } => *target_pressure_kpa,
            Self::Excavate {
                remaining_surface_pressure_kpa,
                ..
            } => *remaining_surface_pressure_kpa,
            Self::Boundary { .. }
            | Self::SteadySeepage { .. }
            | Self::HydrostaticWater { .. }
            | Self::Birth { .. }
            | Self::MaterialSet { .. }
            | Self::ResetState { .. }
            | Self::ResetDisplay { .. } => previous,
            Self::InstallBeams { .. } | Self::RemoveBeams { .. } => previous,
            Self::InstallAnchors { .. } | Self::RemoveAnchors { .. } | Self::JackAnchors { .. } => {
                previous
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StageCase {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_system: Option<crate::joints::Configuration>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_system: Option<crate::anchors::Configuration>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub beam_system: Option<crate::beams::Configuration>,
    pub schema: String,
    pub settings: Settings,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub boundary: Option<Conditions>,
    pub monitors: BTreeMap<String, f64>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub monitor_sets: BTreeMap<String, String>,
    /// Local cuts retain surcharge only on the original ground surface.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub pressure_on_original_surface_only: bool,
    pub stages: Vec<Stage>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub initially_inactive_bodies: Vec<String>,
    /// Explicit verification-only injection, identified as such in every result.
    #[serde(default)]
    pub verification_reject_once_stages: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Committed {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mixed_mcc: Option<crate::mixed_mcc::State>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_water: Option<crate::joint_water::Field>,
    /// ONE_WAY pore pressure at every immutable integration-point index.
    /// Inactive points retain the pressure at removal. None uses hydrostatics.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pore_pressure_kpa: Option<Vec<f64>>,
    /// Prescribed hydrostatic level, or exterior water level when a separate
    /// ONE_WAY pore field is present. None preserves the fixed-water schema.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hydrostatic_water_level_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joints: Option<crate::joints::State>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchors: Option<crate::anchors::State>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub beams: Option<crate::beams::State>,
    /// Display epoch only. Constitutive displacement and geometry remain unchanged.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub display_reference_displacement_m: Vec<[f64; 3]>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub element_material_ids: Vec<u32>,
    pub displacement_m: Vec<[f64; 3]>,
    pub material_states: Vec<soil::State>,
    pub active_elements: Vec<bool>,
    pub applied_nodal_load_kn: Vec<f64>,
    pub excavation_release_load_kn: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub geostatic_release_load_kn: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transition_release_load_kn: Vec<f64>,
    pub surface_pressure_kpa: f64,
    #[serde(default)]
    pub boundary_factor: f64,
    pub stage_index: usize,
    pub stage_progress: f64,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Attempt {
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub rejected_trial_water_changed: bool,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contact_predictors: Vec<ContactPredictor>,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub rejected_trial_joint_slip_changed: bool,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub rejected_trial_structure_changed: bool,
    pub stage: String,
    pub progress_from: f64,
    pub progress_to: f64,
    pub accepted: bool,
    pub reason: String,
    pub newton_iterations: Option<usize>,
    pub relative_force_residual: Option<f64>,
    pub committed_state_unchanged_on_rejection: bool,
    pub rejected_trial_material_changed: bool,
    #[serde(default)]
    pub rejected_trial_plastic_changed: bool,
    pub rejected_trial_active_set_changed: bool,
    pub rejected_trial_load_changed: bool,
    #[serde(default)]
    pub linear_solves: Vec<LinearReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContactPredictor {
    #[serde(default)]
    pub elastic_activation: bool,
    pub projected_load: bool,
    pub elastic_neutral_correction_norm_m: f64,
    pub removed_neutral_correction_norm_m: f64,
    pub iteration: usize,
    pub zero_increment_free_dofs: Vec<usize>,
    pub activation_predictor: bool,
    pub original_relative_linear_residual: f64,
    pub temporary_linear_solve: LinearReport,
    pub full_relative_force_before: f64,
    pub full_relative_force_after: Option<f64>,
    pub accepted_alpha: Option<f64>,
    pub repartitioned: bool,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Monitor {
    pub active_node_count: usize,
    pub mean_displacement_z_m: Option<f64>,
    pub stage_increment_z_m: Option<f64>,
    pub minimum_displacement_z_m: Option<f64>,
    pub maximum_displacement_z_m: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_support: Option<crate::rigidity::Report>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_results: Option<crate::joints::Report>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_results: Option<crate::anchors::Report>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub beam_results: Option<crate::beams::Report>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub global_moment_balance_error_kn_m: Option<[f64; 3]>,
    pub name: String,
    pub state: Committed,
    pub relative_force_residual: f64,
    pub applied_force_kn: [f64; 3],
    pub support_reaction_kn: [f64; 3],
    pub force_balance_error_kn: [f64; 3],
    pub maximum_normalized_yield: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maximum_absolute_engineering_strain_from_k0: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maximum_absolute_engineering_strain_from_seed: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maximum_absolute_engineering_strain_from_initial_equilibrium: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maximum_absolute_material_engineering_strain_from_seed: Option<f64>,
    pub maximum_constraint_error_m: f64,
    pub monitors: BTreeMap<String, Monitor>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Run {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub flow_audits: Vec<FlowAudit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub joint_integration_points: Vec<crate::joints::PointInfo>,
    pub schema: String,
    pub provenance: crate::checkpoint::Identity,
    pub complete: bool,
    pub failure_reason: Option<String>,
    pub case: StageCase,
    pub ground: Ground,
    pub frames: Vec<Frame>,
    pub attempts: Vec<Attempt>,
    pub excavation_audits: Vec<ExcavationAudit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transition_audits: Vec<TransitionAudit>,
    pub last_committed: Committed,
    #[serde(default)]
    pub integration_points: Vec<PointInfo>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub continuation: Option<Continuation>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub geostatic: Option<geostatic::Report>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reference_displacement_m: Vec<[f64; 3]>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FlowAudit {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub joint_flow: Option<crate::joint_water::Flow>,
    pub stage: String,
    pub stage_index: usize,
    pub solution: crate::darcy::Report,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TransitionAudit {
    pub stage: String,
    pub selected_elements: Vec<usize>,
    /// Constitutive state after explicit replacement, before any equilibrium increment.
    pub prepared_state: Committed,
    pub initial_holding_load_kn: Vec<f64>,
    pub zero_release_relative_force_residual: f64,
}

/// Complete stage origin and load path, needed to resume partial excavation.
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Continuation {
    origin: Committed,
    active: Vec<bool>,
    end_load: Vec<f64>,
    holding_load: Vec<f64>,
    stage_index: usize,
    progress: f64,
    increment: f64,
    stage_attempts: usize,
    injected_stages: BTreeSet<String>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExcavationAudit {
    pub stage: String,
    pub removed_element_count: usize,
    pub removed_volume_m3: f64,
    /// Holding forces at lambda=0. The stage removes these continuously to zero.
    pub initial_holding_load_kn: Vec<f64>,
    pub zero_release_relative_force_residual: f64,
}

struct Element {
    volume: f64,
    points: Vec<VolumePoint>,
    point_start: usize,
    dofs: Vec<usize>,
    material: usize,
}

/// Immutable material work. The solver's observer, structures and global
/// matrices never cross a thread boundary. Results retain the original GP order.
struct MaterialWork {
    point: usize,
    material: usize,
    strain: [f64; 6],
    partition: Option<usize>,
}
type MaterialResponse = (soil::Update, Option<crate::elastic::Tensor>);

fn material_work_chunk(
    models: &[soil::Model],
    previous: Option<&[soil::Model]>,
    states: &[soil::State],
    work: &[MaterialWork],
    tangent: bool,
    perturbation: f64,
    maximum_substep: f64,
) -> Result<Vec<MaterialResponse>> {
    work.iter()
        .map(|item| {
            let material = models
                .get(item.material)
                .ok_or_else(|| Error::Input("material work model missing".into()))?;
            let old = states
                .get(item.point)
                .ok_or_else(|| Error::Input("material work history missing".into()))?;
            if let Some(previous) = previous {
                let from = previous
                    .get(item.material)
                    .ok_or_else(|| Error::Input("previous MC material missing".into()))?;
                let (
                    soil::Model::MohrCoulomb(to),
                    soil::Model::MohrCoulomb(from),
                    soil::State::MohrCoulomb(old),
                ) = (material, from, old)
                else {
                    return Err(Error::Input(
                        "material transition requires MC history and models".into(),
                    ));
                };
                let steps = item.partition.unwrap_or_else(|| {
                    (item.strain.iter().map(|x| x.abs()).fold(0., f64::max) / maximum_substep)
                        .ceil()
                        .max(1.) as usize
                });
                let (update, tensor) = if tangent {
                    let (u, d) = to.stress_tangent_transition_partition(
                        from,
                        old,
                        item.strain,
                        steps,
                        perturbation,
                    )?;
                    (u, Some(d))
                } else {
                    (
                        to.integrate_transition_partition(from, old, item.strain, steps)?,
                        None,
                    )
                };
                return Ok((
                    soil::Update {
                        state: soil::State::MohrCoulomb(update.state),
                        substeps: update.substeps,
                        plastic_substeps: update.plastic_substeps,
                    },
                    tensor,
                ));
            }
            if tangent {
                let (update, tensor) = material.stress_tangent_partition(
                    old,
                    item.strain,
                    item.partition.unwrap_or(1),
                    perturbation,
                )?;
                Ok((update, Some(tensor)))
            } else {
                let update = if let Some(partition) = item.partition {
                    material.integrate_partition(old, item.strain, partition)?
                } else {
                    material.integrate(old, item.strain, maximum_substep)?
                };
                Ok((update, None))
            }
        })
        .collect()
}

fn material_work_batch(
    models: &[soil::Model],
    previous: Option<&[soil::Model]>,
    states: &[soil::State],
    work: &[MaterialWork],
    mode: (bool, f64, f64),
    workers: usize,
) -> Result<Vec<MaterialResponse>> {
    if !(1..=4).contains(&workers) {
        return Err(Error::Input("material worker count must be1..4".into()));
    }
    let (tangent, perturbation, maximum_substep) = mode;
    if workers == 1 || work.len() < 64 {
        return material_work_chunk(
            models,
            previous,
            states,
            work,
            tangent,
            perturbation,
            maximum_substep,
        );
    }
    let chunk_size = work.len().div_ceil(workers);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        let mut failure = None;
        for (index, chunk) in work.chunks(chunk_size).enumerate() {
            let builder = std::thread::Builder::new().name(format!("soil-material-{index}"));
            match builder.spawn_scoped(scope, move || {
                material_work_chunk(
                    models,
                    previous,
                    states,
                    chunk,
                    tangent,
                    perturbation,
                    maximum_substep,
                )
            }) {
                Ok(handle) => handles.push(handle),
                Err(error) => {
                    failure = Some(Error::Control(format!(
                        "material worker creation failed: {error}"
                    )));
                    break;
                }
            }
        }
        let mut ordered = Vec::with_capacity(work.len());
        // Join every owned worker, including after an error. The first failed
        // GP in input order wins; no thread can commit any material history.
        for handle in handles {
            match handle.join() {
                Ok(Ok(values)) if failure.is_none() => ordered.extend(values),
                Ok(Err(error)) if failure.is_none() => failure = Some(error),
                Err(_) if failure.is_none() => {
                    failure = Some(Error::Control(
                        "material worker failed without committing history".into(),
                    ))
                }
                _ => {}
            }
        }
        if let Some(error) = failure {
            Err(error)
        } else {
            Ok(ordered)
        }
    })
}

struct Model<'a> {
    material_models: Vec<soil::Model>,
    material_transition_from: Option<Vec<soil::Model>>,
    observer: Option<&'a dyn crate::stage_control::Observer>,
    wet_joints: Option<crate::joint_water::Model>,
    soil_points: Vec<bool>,
    joints: Option<crate::joints::Model>,
    anchors: Option<crate::anchors::Model>,
    beams: Option<crate::beams::Model>,
    mesh: &'a Mesh,
    ground: &'a Ground,
    case: &'a StageCase,
    elements: Vec<Element>,
    bounds: [[f64; 3]; 2],
    faces: BTreeMap<[usize; 3], Vec<usize>>,
    point_count: usize,
    reference_displacement_m: Vec<[f64; 3]>,
}
struct Dofs {
    map: Vec<Option<usize>>,
    rows: Vec<Vec<(usize, f64)>>,
    transform: Option<Transform>,
    count: usize,
    active_nodes: Vec<bool>,
    sparse_pattern: Option<Arc<CsrPattern>>,
}
struct Evaluation {
    joints: Option<crate::joints::Evaluation>,
    anchors: Option<crate::anchors::Evaluation>,
    beams: Option<crate::beams::Evaluation>,
    states: Vec<soil::State>,
    internal: Vec<f64>,
    tangent: Option<Tangent>,
}
struct Trial {
    state: Committed,
    iterations: usize,
    residual: f64,
    linear_solves: Vec<LinearReport>,
    contact_predictors: Vec<ContactPredictor>,
}

#[derive(Clone)]
enum Tangent {
    Dense(DMatrix<f64>),
    Sparse(CsrMatrix),
}
impl Tangent {
    fn action(&self, x: &DVector<f64>) -> Result<DVector<f64>> {
        Ok(match self {
            Self::Dense(a) => a * x,
            Self::Sparse(a) => DVector::from_vec(a.matvec(x.as_slice())?),
        })
    }
    fn scale(&self) -> f64 {
        match self {
            Self::Dense(a) => (0..a.nrows())
                .map(|i| a.row(i).iter().map(|x| x.abs()).sum::<f64>())
                .fold(0., f64::max),
            Self::Sparse(a) => a
                .row_ptr()
                .windows(2)
                .map(|r| a.values()[r[0]..r[1]].iter().map(|x| x.abs()).sum::<f64>())
                .fold(0., f64::max),
        }
    }
    fn gauged(&self, indices: &[usize]) -> Result<Self> {
        let diagonal = self.scale().max(1.);
        Ok(match self {
            Self::Sparse(a) => Self::Sparse(a.zero_increment_gauges(indices, diagonal)?),
            Self::Dense(a) => {
                let mut b = a.clone();
                for &i in indices {
                    b.row_mut(i).fill(0.);
                    b.column_mut(i).fill(0.);
                    b[(i, i)] = diagonal;
                }
                Self::Dense(b)
            }
        })
    }
    fn add(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        match self {
            Self::Dense(a) => a[(i, j)] += value,
            Self::Sparse(a) => a.add(i, j, value)?,
        }
        Ok(())
    }
    fn solve(
        &self,
        rhs: &DVector<f64>,
        settings: &Settings,
    ) -> Result<(DVector<f64>, LinearReport)> {
        match self {
            Self::Sparse(a) => {
                let LinearSolver::SparseGmresIlu0 { options } = &settings.linear_solver else {
                    return Err(Error::Input("sparse tangent/solver mismatch".into()));
                };
                let (x, report) = linalg::gmres_ilu0(
                    a,
                    rhs.as_slice(),
                    settings.relative_linear_tolerance,
                    options,
                )?;
                Ok((DVector::from_vec(x), report))
            }
            Self::Dense(a) => {
                let started = Instant::now();
                let x = a.clone().lu().solve(rhs).ok_or_else(|| {
                    Error::Solve("singular non-symmetric tangent; inspect active supports".into())
                })?;
                let residual = (a * &x - rhs).norm() / rhs.norm().max(1e-30);
                if !residual.is_finite()
                    || residual > settings.relative_linear_tolerance
                    || !x.iter().all(|v| v.is_finite())
                {
                    return Err(Error::Solve(format!("true linear residual {residual:e}")));
                }
                let report = LinearReport {
                    method: "dense_lu".into(),
                    dimension: a.nrows(),
                    stored_nonzeros: a.len(),
                    iterations: 1,
                    restarts: 0,
                    true_relative_residual: residual,
                    matrix_storage_bytes: a.len() * size_of::<f64>(),
                    workspace_storage_bytes: (a.len() + 2 * a.nrows()) * size_of::<f64>(),
                    elapsed_seconds: started.elapsed().as_secs_f64(),
                };
                Ok((x, report))
            }
        }
    }
}
#[derive(Clone, Copy)]
struct EvaluationMode<'p> {
    structures: bool,
    anchor_control: Option<&'p crate::anchors::Control>,
    factor: Option<f64>,
    tangent: bool,
    partitions: Option<&'p [usize]>,
}
impl EvaluationMode<'_> {
    fn plain() -> Self {
        Self {
            structures: true,
            anchor_control: None,
            factor: None,
            tangent: false,
            partitions: None,
        }
    }
}
#[derive(Clone)]
struct Target {
    joint_water: Option<crate::joint_water::Field>,
    pore_pressure_kpa: Option<Vec<f64>>,
    water_level: Option<f64>,
    anchor_control: Option<crate::anchors::Control>,
    active: Vec<bool>,
    load: Vec<f64>,
    release: Vec<f64>,
    pressure: f64,
    boundary_factor: f64,
    stage: usize,
    progress: f64,
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-8
}
fn same_beam_installation(
    a: &Option<crate::beams::State>,
    b: &Option<crate::beams::State>,
) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.active_members == b.active_members
                && a.reference_global_kinematics == b.reference_global_kinematics
                && a.active_plates == b.active_plates
                && a.plate_reference_global_kinematics == b.plate_reference_global_kinematics
        }
        _ => false,
    }
}

fn joint_slip_changed(a: &Option<crate::joints::State>, b: &Option<crate::joints::State>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.points.iter().zip(&b.points).any(|(a, b)| {
            a.plastic_slip_m != b.plastic_slip_m
                || a.accumulated_plastic_slip_m != b.accumulated_plastic_slip_m
        }),
        _ => false,
    }
}
fn same_anchor_basis(a: &Option<crate::anchors::State>, b: &Option<crate::anchors::State>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.anchors.len() == b.anchors.len()
                && a.anchors.iter().zip(&b.anchors).all(|(a, b)| {
                    a.active == b.active
                        && a.reference_extensions_m == b.reference_extensions_m
                        && a.force_control_kn == b.force_control_kn
                        && a.last_lock_off_force_kn == b.last_lock_off_force_kn
                })
        }
        _ => false,
    }
}

impl StageCase {
    fn is_mixed_mcc(&self) -> bool {
        self.schema == "soil-fem-mixed-mcc-stages/1" || self.is_projected_mcc()
    }
    fn is_projected_mcc(&self) -> bool {
        self.schema == "soil-fem-projected-mcc-stages/1"
    }
    fn has_seepage(&self) -> bool {
        self.stages
            .iter()
            .any(|s| matches!(s, Stage::SteadySeepage { .. }))
    }
    fn has_staged_water(&self) -> bool {
        self.stages.iter().any(|s| {
            matches!(
                s,
                Stage::HydrostaticWater { .. } | Stage::SteadySeepage { .. }
            )
        })
    }
    fn has_lifecycle(&self) -> bool {
        self.joint_system.is_some()
            || self.anchor_system.is_some()
            || self.beam_system.is_some()
            || !self.initially_inactive_bodies.is_empty()
            || self
                .stages
                .iter()
                .any(|s| s.is_transition() || matches!(s, Stage::ResetDisplay { .. }))
    }
    pub fn validate(&self) -> Result<()> {
        let s = &self.settings;
        if let Some(c) = &self.joint_system
            && self.has_staged_water()
            && !c.joints.iter().any(|d| d.water.is_some())
        {
            return Err(Error::Input(
                "staged water with JOINT needs an explicit hydraulic interface policy".into(),
            ));
        }
        if self.joint_system.is_some()
            && (matches!(s.kinematics, Kinematics::Oedometer)
                || !self.pressure_on_original_surface_only
                || self.stages.iter().any(|s| matches!(s, Stage::Birth { .. })))
        {
            return Err(Error::Input("JOINT requires independent displacements and original-surface surcharge; BIRTH reactivation requires an explicit interface history policy".into()));
        }
        if self.anchor_system.is_some() && matches!(s.kinematics, Kinematics::Oedometer) {
            return Err(Error::Input(
                "anchors require explicit or side-roller kinematics".into(),
            ));
        }
        s.linear_solver.validate()?;
        if self.beam_system.is_some() && matches!(s.kinematics, Kinematics::Oedometer) {
            return Err(Error::Input(
                "beams require explicit or side-roller kinematics".into(),
            ));
        }
        if (matches!(s.kinematics, Kinematics::Explicit) && self.boundary.is_none())
            || (matches!(s.kinematics, Kinematics::Oedometer) && self.boundary.is_some())
        {
            return Err(Error::Input("explicit kinematics requires boundary conditions; extra conditions with tied-plane oedometer are unsupported".into()));
        }
        if (self.schema != "soil-fem-stages/1" && !self.is_mixed_mcc())
            || !s.relative_force_tolerance.is_finite()
            || !(1e-11..=1e-5).contains(&s.relative_force_tolerance)
            || !s.relative_linear_tolerance.is_finite()
            || !(1e-12..=1e-6).contains(&s.relative_linear_tolerance)
            || s.maximum_newton_iterations == 0
            || s.maximum_newton_iterations > 60
            || !s.initial_increment.is_finite()
            || !(0.0001..=1.).contains(&s.initial_increment)
            || !s.minimum_increment.is_finite()
            || s.minimum_increment <= 0.
            || s.minimum_increment > s.initial_increment
            || !s.material_maximum_substep.is_finite()
            || !(1e-7..=1e-3).contains(&s.material_maximum_substep)
            || !s.tangent_perturbation.is_finite()
            || !(1e-9..=1e-5).contains(&s.tangent_perturbation)
            || s.maximum_free_dofs == 0
            || s.maximum_free_dofs > s.linear_solver.dimension_limit()
            || self.stages.len() > 100
            || self.monitors.values().any(|v| !v.is_finite())
        {
            return Err(Error::Input("invalid staged MCC settings".into()));
        }
        let mut names = BTreeSet::new();
        if self.initially_inactive_bodies.iter().any(String::is_empty)
            || self
                .initially_inactive_bodies
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != self.initially_inactive_bodies.len()
        {
            return Err(Error::Input(
                "initially inactive bodies must be unique and nonempty".into(),
            ));
        }
        for stage in &self.stages {
            if let Stage::SteadySeepage {
                exterior_water_level_m,
                ..
            } = stage
                && !exterior_water_level_m.is_finite()
            {
                return Err(Error::Input(
                    "seepage needs a finite, separately specified exterior water level".into(),
                ));
            }
            if let Stage::HydrostaticWater { target_level_m, .. } = stage
                && !target_level_m.is_finite()
            {
                return Err(Error::Input(
                    "water stage requires a finite level in metres".into(),
                ));
            }
            if let Stage::InstallAnchors { anchors, .. }
            | Stage::RemoveAnchors { anchors, .. }
            | Stage::JackAnchors { anchors, .. } = stage
                && (self.anchor_system.is_none()
                    || anchors.is_empty()
                    || anchors.iter().collect::<BTreeSet<_>>().len() != anchors.len())
            {
                return Err(Error::Input(
                    "anchor operation requires a system and unique IDs".into(),
                ));
            }
            if let Stage::JackAnchors {
                target_force_kn, ..
            } = stage
                && (!target_force_kn.is_finite() || *target_force_kn < 0.)
            {
                return Err(Error::Input(
                    "anchor jack force must be finite and nonnegative".into(),
                ));
            }
            if let Stage::InstallBeams { members, .. } | Stage::RemoveBeams { members, .. } = stage
                && (self.beam_system.is_none()
                    || members.is_empty()
                    || members.iter().collect::<BTreeSet<_>>().len() != members.len())
            {
                return Err(Error::Input(
                    "beam operation requires a system and unique member IDs".into(),
                ));
            }
            if stage.name().is_empty()
                || !names.insert(stage.name())
                || !stage.pressure(0.).is_finite()
                || stage.pressure(0.) < 0.
                || !stage.step(s).is_finite()
                || stage.step(s) < s.minimum_increment
                || stage.step(s) > 1.
            {
                return Err(Error::Input(
                    "invalid/duplicate stage name or pressure/increment".into(),
                ));
            }
            if let Stage::Boundary { target_factor, .. } = stage
                && (!target_factor.is_finite() || self.boundary.is_none())
            {
                return Err(Error::Input(
                    "boundary stage requires conditions and a finite target factor".into(),
                ));
            }
            if let Stage::Excavate { bodies, .. }
            | Stage::Birth { bodies, .. }
            | Stage::MaterialSet { bodies, .. }
            | Stage::ResetState { bodies, .. } = stage
                && (bodies.is_empty()
                    || bodies.iter().collect::<BTreeSet<_>>().len() != bodies.len())
            {
                return Err(Error::Input(
                    "construction operation requires unique body names".into(),
                ));
            }
        }
        if self
            .verification_reject_once_stages
            .iter()
            .any(|n| !names.contains(n.as_str()))
        {
            return Err(Error::Input(
                "verification injection names unknown stage".into(),
            ));
        }
        Ok(())
    }
}

impl<'a> Model<'a> {
    fn progress(
        &self,
        kind: crate::stage_control::Kind,
        stage: usize,
        progress: f64,
    ) -> crate::stage_control::Progress {
        crate::stage_control::Progress {
            schema: "soil-fem-stage-progress/1".into(),
            kind,
            stage_index: stage,
            stage: if stage == 0 {
                "geostatic_equilibration".into()
            } else {
                self.case
                    .stages
                    .get(stage - 1)
                    .map_or("external_increment", Stage::name)
                    .into()
            },
            progress,
            newton_iteration: None,
            relative_force_residual: None,
            relative_volume_residual: None,
            relative_pressure_residual: None,
            maximum_material_partitions: None,
            linear_solve: None,
            reason: None,
        }
    }
    fn notify(&self, event: impl FnOnce() -> crate::stage_control::Progress) -> Result<()> {
        if let Some(observer) = self.observer {
            observer
                .on_progress(&event())
                .map_err(|e| Error::Control(e.to_string()))?;
        }
        Ok(())
    }
    fn stop_requested(&self) -> Result<bool> {
        self.observer.map_or(Ok(false), |o| {
            o.stop_requested()
                .map_err(|e| Error::Control(e.to_string()))
        })
    }
    fn notify_attempt(&self, stage: usize, attempt: &Attempt) -> Result<()> {
        self.notify(|| {
            let mut event = self.progress(
                if attempt.accepted {
                    crate::stage_control::Kind::Accepted
                } else {
                    crate::stage_control::Kind::Rejected
                },
                stage,
                attempt.progress_to,
            );
            event.newton_iteration = attempt.newton_iterations;
            event.relative_force_residual = attempt.relative_force_residual;
            event.reason = Some(attempt.reason.clone());
            event
        })
    }
    fn anchor_control(
        &self,
        origin: &Committed,
        stage: &Stage,
        progress: f64,
    ) -> Result<Option<crate::anchors::Control>> {
        if let Stage::JackAnchors {
            anchors,
            target_force_kn,
            ..
        } = stage
        {
            let model = self
                .anchors
                .as_ref()
                .ok_or_else(|| Error::Input("missing anchor configuration".into()))?;
            let state = origin
                .anchors
                .as_ref()
                .ok_or_else(|| Error::Input("missing anchor state".into()))?;
            Ok(Some(model.control(
                state,
                anchors,
                *target_force_kn,
                progress,
            )?))
        } else {
            Ok(None)
        }
    }
    fn material_index(&self, ids: &[u32], ei: usize) -> Result<usize> {
        if ids.is_empty() {
            return Ok(self.elements[ei].material);
        }
        if ids.len() != self.elements.len() {
            return Err(Error::Input(
                "invalid per-element material assignment".into(),
            ));
        }
        self.ground
            .materials
            .iter()
            .position(|m| m.id == ids[ei])
            .ok_or_else(|| Error::Input("unknown current element material".into()))
    }
    fn selected_bodies(&self, bodies: &[String]) -> Result<Vec<usize>> {
        if bodies.is_empty() || bodies.iter().collect::<BTreeSet<_>>().len() != bodies.len() {
            return Err(Error::Input(
                "body selection must be nonempty and unique".into(),
            ));
        }
        let mut selected = Vec::new();
        for body in bodies {
            let region = self
                .mesh
                .regions
                .iter()
                .find(|r| &r.body == body)
                .ok_or_else(|| Error::Input(format!("unknown material body {body}")))?;
            let indices: Vec<_> = self
                .mesh
                .region_ids
                .iter()
                .enumerate()
                .filter(|(_, r)| **r == region.id)
                .map(|(i, _)| i)
                .collect();
            if indices.is_empty() {
                return Err(Error::Input("body has no elements".into()));
            }
            selected.extend(indices);
        }
        selected.sort_unstable();
        Ok(selected)
    }
    fn interface_connectivity(&self) -> Vec<Vec<usize>> {
        let mut result = self
            .anchors
            .as_ref()
            .map_or_else(Vec::new, |m| m.connectivity());
        if let Some(joints) = &self.joints {
            result.extend(joints.connectivity());
        }
        result
    }
    fn joint_layout(&self) -> Vec<crate::joints::PointInfo> {
        self.joints.as_ref().map_or_else(Vec::new, |m| m.layout())
    }
    /// Potential connectivity only; actual rigid support is checked separately.
    fn explicitly_supported_elements(&self, active: &[bool]) -> Vec<usize> {
        let constrained: BTreeSet<_> = self
            .case
            .boundary
            .as_ref()
            .into_iter()
            .flat_map(|bc| &bc.constraints)
            .filter_map(|c| self.mesh.node_sets.get(&c.node_set))
            .flatten()
            .copied()
            .collect();
        (0..active.len())
            .filter(|&e| {
                active[e]
                    && self
                        .mesh
                        .cell_nodes(e)
                        .iter()
                        .any(|n| constrained.contains(n))
            })
            .collect()
    }
    fn ensure_supported(&self, active: &[bool]) -> Result<()> {
        let mut adjacency = vec![Vec::new(); active.len()];
        let mut stack = self.explicitly_supported_elements(active);
        if let Some(joints) = &self.joints {
            for &[a, b] in &joints.owner_pairs {
                if active[a] && active[b] {
                    adjacency[a].push(b);
                    adjacency[b].push(a);
                }
            }
        }
        for (face, owners) in &self.faces {
            if owners.len() == 2 && active[owners[0]] && active[owners[1]] {
                adjacency[owners[0]].push(owners[1]);
                adjacency[owners[1]].push(owners[0]);
            }
            if face
                .iter()
                .all(|&n| near(self.mesh.nodes[n][2], self.bounds[0][2]))
            {
                stack.extend(owners.iter().copied().filter(|&i| active[i]));
            }
        }
        let mut visited = BTreeSet::new();
        while let Some(ei) = stack.pop() {
            if visited.insert(ei) {
                stack.extend(&adjacency[ei]);
            }
        }
        if !active.iter().any(|x| *x)
            || active
                .iter()
                .enumerate()
                .any(|(i, a)| *a && !visited.contains(&i))
        {
            return Err(Error::Input(
                "active construction soil must connect by faces to the supported base".into(),
            ));
        }
        Ok(())
    }
    fn prepare_transition(
        &self,
        origin: &Committed,
        stage: &Stage,
    ) -> Result<(Committed, Vec<usize>)> {
        self.prepare_transition_mode(origin, stage, false)
    }
    fn prepare_transition_mode(
        &self,
        origin: &Committed,
        stage: &Stage,
        external_fluid: bool,
    ) -> Result<(Committed, Vec<usize>)> {
        if matches!(stage, Stage::Birth { .. })
            && origin.pore_pressure_kpa.is_some()
            && !external_fluid
        {
            return Err(Error::Input(
                "BIRTH after seepage needs an explicit initial hydraulic state for the new soil"
                    .into(),
            ));
        }
        if let Stage::InstallAnchors { anchors, .. } | Stage::RemoveAnchors { anchors, .. } = stage
        {
            let model = self
                .anchors
                .as_ref()
                .ok_or_else(|| Error::Input("missing anchor system".into()))?;
            let state = origin
                .anchors
                .as_ref()
                .ok_or_else(|| Error::Input("missing anchor history".into()))?;
            let mut prepared = origin.clone();
            prepared.anchors = Some(model.change(
                state,
                &origin.displacement_m,
                anchors,
                matches!(stage, Stage::InstallAnchors { .. }),
            )?);
            return Ok((prepared, vec![]));
        }
        if let Stage::InstallBeams { members, .. } | Stage::RemoveBeams { members, .. } = stage {
            let model = self
                .beams
                .as_ref()
                .ok_or_else(|| Error::Input("missing beam system".into()))?;
            let old = origin
                .beams
                .as_ref()
                .ok_or_else(|| Error::Input("missing beam state".into()))?;
            let mut prepared = origin.clone();
            prepared.beams = Some(model.change_members(
                old,
                &origin.displacement_m,
                members,
                matches!(stage, Stage::InstallBeams { .. }),
                matches!(
                    stage,
                    Stage::InstallBeams {
                        placement: BirthPlacement::CadNewNodes,
                        ..
                    }
                ),
            )?);
            return Ok((prepared, vec![]));
        }
        let (bodies, replacements, material_id, birth) = match stage {
            Stage::ResetState { bodies, states, .. } => (bodies, Some(states), None, false),
            Stage::Birth { bodies, states, .. } => (bodies, Some(states), None, true),
            Stage::MaterialSet {
                bodies,
                material_id,
                history,
                ..
            } => (
                bodies,
                match history {
                    MaterialHistory::Preserve => None,
                    MaterialHistory::Replace { states } => Some(states),
                },
                Some(*material_id),
                false,
            ),
            _ => return Ok((origin.clone(), Vec::new())),
        };
        if let Some(states) = replacements
            && states.keys().collect::<BTreeSet<_>>() != bodies.iter().collect::<BTreeSet<_>>()
        {
            return Err(Error::Input(
                "replacement states must match selected body names exactly".into(),
            ));
        }
        let selected = self.selected_bodies(bodies)?;
        let mut prepared = origin.clone();
        if matches!(
            stage,
            Stage::Birth {
                placement: BirthPlacement::CadNewNodes,
                ..
            }
        ) {
            let mut existing = vec![false; self.mesh.nodes.len()];
            for (ei, a) in origin.active_elements.iter().enumerate() {
                if *a {
                    for n in self.mesh.cell_nodes(ei) {
                        existing[n] = true;
                    }
                }
            }
            for &ei in &selected {
                for n in self.mesh.cell_nodes(ei) {
                    if !existing[n] {
                        prepared.displacement_m[n] = [0.; 3];
                    }
                }
            }
        }
        if prepared.element_material_ids.is_empty() {
            prepared.element_material_ids = self.mesh.material_ids.clone();
        }
        for &ei in &selected {
            if birth && origin.active_elements[ei] {
                return Err(Error::Input("BIRTH requires a wholly inactive body".into()));
            }
            if let Some(id) = material_id {
                prepared.element_material_ids[ei] = id;
            }
            let mi = self.material_index(&prepared.element_material_ids, ei)?;
            let model = &self.material_models[mi];
            if self.wet_joints.as_ref().is_some_and(|w| !w.is_soil(ei))
                && !matches!(model, soil::Model::MohrCoulomb(_))
            {
                return Err(Error::Input("impermeable total-stress structure cannot be rematerialised as an MCC/soil hardening model".into()));
            }
            let region = self
                .mesh
                .regions
                .iter()
                .find(|r| r.id == self.mesh.region_ids[ei])
                .unwrap();
            for (qi, q) in self.elements[ei].points.iter().enumerate() {
                let point = self.elements[ei].point_start + qi;
                if let Some(states) = replacements {
                    prepared.material_states[point] = states[&region.body].clone();
                }
                let state = &prepared.material_states[point];
                if birth {
                    if !matches!(model, soil::Model::MohrCoulomb(_))
                        || !matches!(state, soil::State::MohrCoulomb(_))
                    {
                        return Err(Error::Input("stress-free BIRTH requires a zero-stress-admissible model; MCC has zero bulk stiffness at p=0".into()));
                    }
                    if self.pore(origin.hydrostatic_water_level_m, q.position[2]) != 0.
                        || state
                            .stress()
                            .iter()
                            .chain(state.plastic_strain())
                            .any(|v| *v != 0.)
                    {
                        return Err(Error::Input(
                            "stress-free BIRTH requires dry, zero-stress, fresh material states"
                                .into(),
                        ));
                    }
                    if let soil::State::MohrCoulomb(s) = state
                        && (s.equivalent_plastic_strain != 0.
                            || s.plastic_work_kpa != 0.
                            || !matches!(s.last_return, crate::mc::Return::Elastic))
                    {
                        return Err(Error::Input(
                            "BIRTH cannot inherit prior plastic history".into(),
                        ));
                    }
                }
                model.validate_state(state)?;
            }
            if birth {
                prepared.active_elements[ei] = true;
            }
        }
        self.ensure_supported(&prepared.active_elements)?;
        if birth && !prepared.display_reference_displacement_m.is_empty() {
            let mut existing = vec![false; self.mesh.nodes.len()];
            for (ei, a) in origin.active_elements.iter().enumerate() {
                if *a {
                    for n in self.mesh.cell_nodes(ei) {
                        existing[n] = true;
                    }
                }
            }
            for &ei in &selected {
                for n in self.mesh.cell_nodes(ei) {
                    if !existing[n] {
                        prepared.display_reference_displacement_m[n] = prepared.displacement_m[n];
                    }
                }
            }
        }
        Ok((prepared, selected))
    }
    fn new(mesh: &'a Mesh, ground: &'a Ground, case: &'a StageCase) -> Result<Self> {
        if case.settings.mcc_integration.is_some()
            && !ground.materials.iter().any(|m| m.model.is_mcc())
        {
            return Err(Error::Input(
                "MCC integration override requires an MCC material".into(),
            ));
        }
        let material_models = ground
            .materials
            .iter()
            .map(|material| {
                let mut model = material.model.clone();
                if let (soil::Model::Mcc(parameters), Some(integration)) =
                    (&mut model, case.settings.mcc_integration)
                {
                    parameters.integration = integration;
                }
                model
            })
            .collect();
        if case.is_mixed_mcc() {
            mixed::validate_definition(mesh, ground, case)?;
        }
        if case.joint_system.is_some()
            && (!ground.is_general()
                || (ground.initialization.water_level_m.is_some()
                    && !case
                        .joint_system
                        .as_ref()
                        .is_some_and(|c| c.joints.iter().any(|d| d.water.is_some()))))
        {
            return Err(Error::Input(
                "dry JOINT requires general geostatic equilibration and no prescribed water level"
                    .into(),
            ));
        }
        if let Some(boundary) = &case.boundary {
            boundary.validate_definition(mesh)?;
        }
        let mut elements = Vec::new();
        let mut point_count = 0;
        for i in 0..mesh.tetrahedra.len() {
            let points = mesh.integration(i)?;
            let volume = points.iter().map(|q| q.weight_m3).sum();
            let point_start = point_count;
            point_count += points.len();
            let material = ground
                .materials
                .iter()
                .position(|m| m.id == mesh.material_ids[i])
                .ok_or_else(|| Error::Input("unknown element material".into()))?;
            elements.push(Element {
                volume,
                points,
                point_start,
                dofs: mesh
                    .cell_nodes(i)
                    .into_iter()
                    .flat_map(|n| [3 * n, 3 * n + 1, 3 * n + 2])
                    .collect(),
                material,
            });
        }
        let bounds = std::array::from_fn(|side| {
            std::array::from_fn(|j| {
                mesh.nodes.iter().map(|p| p[j]).fold(
                    if side == 0 {
                        f64::INFINITY
                    } else {
                        f64::NEG_INFINITY
                    },
                    |a, b| if side == 0 { a.min(b) } else { a.max(b) },
                )
            })
        });
        let mut faces = BTreeMap::<[usize; 3], Vec<usize>>::new();
        for (ei, tet) in mesh.tetrahedra.iter().enumerate() {
            for mut face in [
                [tet[0], tet[1], tet[2]],
                [tet[0], tet[1], tet[3]],
                [tet[0], tet[2], tet[3]],
                [tet[1], tet[2], tet[3]],
            ] {
                face.sort_unstable();
                faces.entry(face).or_default().push(ei);
            }
        }
        let joints = case
            .joint_system
            .clone()
            .map(|c| crate::joints::Model::new(mesh, c))
            .transpose()?;
        let wet_joints = joints
            .as_ref()
            .filter(|j| {
                !j.configuration.structure_bodies.is_empty()
                    || j.configuration.joints.iter().any(|d| d.water.is_some())
            })
            .map(|j| crate::joint_water::Model::new(mesh, j, ground))
            .transpose()?;
        let mut soil_points = vec![true; point_count];
        if let Some(wet) = &wet_joints {
            for (i, e) in elements.iter().enumerate() {
                for p in &mut soil_points[e.point_start..e.point_start + e.points.len()] {
                    *p = wet.is_soil(i);
                }
            }
        }
        Ok(Self {
            material_models,
            material_transition_from: None,
            wet_joints,
            observer: None,
            soil_points,
            joints,
            anchors: case
                .anchor_system
                .clone()
                .map(|c| crate::anchors::Model::new(mesh, c))
                .transpose()?,
            beams: case
                .beam_system
                .clone()
                .map(|c| crate::beams::Model::new(mesh, c))
                .transpose()?,
            mesh,
            ground,
            case,
            elements,
            bounds,
            faces,
            point_count,
            reference_displacement_m: vec![],
        })
    }

    fn dofs(&self, active: &[bool], settings: &Settings, factor: f64) -> Result<Dofs> {
        let mut active_nodes = vec![false; self.mesh.nodes.len()];
        for (i, &enabled) in active.iter().enumerate() {
            if enabled {
                for n in self.mesh.cell_nodes(i) {
                    active_nodes[n] = true;
                }
            }
        }
        let mut map = vec![None; self.mesh.nodes.len() * 3];
        let mut count = 0;
        if let Some(boundary) = &self.case.boundary {
            let mut constraints = boundary.constraints.clone();
            for c in &mut constraints {
                c.value_m *= factor;
                for g in &mut c.gradient {
                    *g *= factor;
                }
            }
            let mut extra = Vec::new();
            if matches!(settings.kinematics, Kinematics::SideRollers) {
                for (node, p) in self.mesh.nodes.iter().enumerate() {
                    if active_nodes[node] {
                        for a in 0..3 {
                            let fixed = if a == 2 {
                                near(p[2], self.bounds[0][2])
                            } else {
                                near(p[a], self.bounds[0][a]) || near(p[a], self.bounds[1][a])
                            };
                            if fixed {
                                let mut normal = Vector3::zeros();
                                normal[a] = 1.;
                                extra.push((node, normal, 0.));
                            }
                        }
                    }
                }
            }
            let transform = if self.joints.is_some() {
                Transform::with_interfaces(self.mesh, &constraints, active.to_vec(), extra)?
            } else {
                Transform::with_extra(self.mesh, &constraints, active.to_vec(), extra)?
            };
            count = transform.count;
            if count > settings.maximum_free_dofs {
                return Err(Error::Input("explicit stage DOF limit exceeded".into()));
            }
            let rows = transform.rows.clone();
            let mut beam_rows = Vec::new();
            if let Some(b) = &self.beams {
                for &n in b.attachments.iter().flatten() {
                    for a in 0..3 {
                        beam_rows.extend(rows[3 * n + a].iter().map(|(r, _)| *r));
                    }
                }
            }
            let anchor_rows: Vec<Vec<usize>> = self
                .interface_connectivity()
                .into_iter()
                .map(|c| {
                    c.iter()
                        .flat_map(|&g| rows[g].iter().map(|(r, _)| *r))
                        .collect()
                })
                .collect();
            let sparse_pattern = if count > 0 && !settings.linear_solver.is_dense() {
                Some(CsrPattern::from_connectivity(
                    count,
                    self.elements
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| active[*i])
                        .map(|(_, e)| {
                            e.dofs
                                .iter()
                                .flat_map(|&g| rows[g].iter().map(|&(r, _)| r))
                                .collect()
                        })
                        .chain(linalg::dense_block_connectivity(beam_rows))
                        .chain(anchor_rows),
                )?)
            } else {
                None
            };
            return Ok(Dofs {
                map,
                rows,
                count,
                active_nodes,
                sparse_pattern,
                transform: Some(transform),
            });
        }
        let mut planes = BTreeMap::<u64, usize>::new();
        for (i, p) in self.mesh.nodes.iter().enumerate() {
            if active_nodes[i] {
                for j in 0..3 {
                    let fixed = if j == 2 {
                        near(p[2], self.bounds[0][2])
                    } else {
                        match settings.kinematics {
                            Kinematics::Oedometer => true,
                            Kinematics::Explicit => false,
                            Kinematics::SideRollers => {
                                near(p[j], self.bounds[0][j]) || near(p[j], self.bounds[1][j])
                            }
                        }
                    };
                    if !fixed {
                        if matches!(settings.kinematics, Kinematics::Oedometer) {
                            // Height bands are shared by material nodes of a horizontal plane.
                            let key = (p[2] * 1e8).round().to_bits();
                            let r = *planes.entry(key).or_insert_with(|| {
                                let r = count;
                                count += 1;
                                r
                            });
                            map[3 * i + j] = Some(r);
                        } else {
                            map[3 * i + j] = Some(count);
                            count += 1;
                        }
                    }
                }
            }
        }
        if count == 0 || count > settings.maximum_free_dofs {
            return Err(Error::Input(format!(
                "active free DOFs {count} exceed configured solver limit {}",
                settings.maximum_free_dofs
            )));
        }
        let mut beam_rows = Vec::new();
        if let Some(b) = &self.beams {
            for &n in b.attachments.iter().flatten() {
                for a in 0..3 {
                    if let Some(r) = map[3 * n + a] {
                        beam_rows.push(r);
                    }
                }
            }
        }
        let anchor_rows: Vec<Vec<usize>> = self
            .interface_connectivity()
            .into_iter()
            .map(|c| c.iter().filter_map(|&g| map[g]).collect())
            .collect();
        let sparse_pattern = if settings.linear_solver.is_dense() {
            None
        } else {
            Some(CsrPattern::from_connectivity(
                count,
                self.elements
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| active[*i])
                    .map(|(_, e)| e.dofs.iter().filter_map(|&g| map[g]).collect())
                    .chain(linalg::dense_block_connectivity(beam_rows))
                    .chain(anchor_rows),
            )?)
        };
        Ok(Dofs {
            rows: map
                .iter()
                .map(|r| r.map(|r| vec![(r, 1.)]).unwrap_or_default())
                .collect(),
            map,
            count,
            active_nodes,
            sparse_pattern,
            transform: None,
        })
    }

    fn layout(&self) -> Vec<PointInfo> {
        self.elements
            .iter()
            .enumerate()
            .flat_map(|(ei, e)| {
                e.points.iter().enumerate().map(move |(qi, q)| PointInfo {
                    element_index: ei,
                    local_index: qi,
                    natural_coordinates: q.natural,
                    position_m: q.position,
                    weight_m3: q.weight_m3,
                    shape_functions: q.shape.clone(),
                })
            })
            .collect()
    }

    fn pore(&self, level: Option<f64>, z: f64) -> f64 {
        level
            .or(self.ground.initialization.water_level_m)
            .map_or(0., |h| {
                (h - z).max(0.) * self.ground.initialization.water_unit_weight_kn_m3
            })
    }
    fn point_pore(&self, state: &Committed, point: usize, z: f64) -> f64 {
        if !self.soil_points[point] {
            return 0.;
        }
        state.pore_pressure_kpa.as_ref().map_or_else(
            || self.pore(state.hydrostatic_water_level_m, z),
            |p| p[point],
        )
    }
    fn stage_pore(
        &self,
        origin: &Committed,
        stage: &Stage,
        active: &[bool],
        flow_target: Option<&[f64]>,
        progress: f64,
    ) -> Result<Option<Vec<f64>>> {
        if flow_target.is_none()
            && !(matches!(stage, Stage::HydrostaticWater { .. })
                && origin.pore_pressure_kpa.is_some())
        {
            return Ok(origin.pore_pressure_kpa.clone());
        }
        let mut value = Vec::with_capacity(self.point_count);
        for (ei, e) in self.elements.iter().enumerate() {
            for (qi, q) in e.points.iter().enumerate() {
                let point = e.point_start + qi;
                let start = self.point_pore(origin, point, q.position[2]);
                let end = if !self.soil_points[point] {
                    0.
                } else if let Some(target) = flow_target {
                    target[point]
                } else {
                    self.pore(self.stage_water(origin, stage, 1.), q.position[2])
                };
                let p = if active[ei] {
                    if progress == 1. {
                        end
                    } else {
                        start + (end - start) * progress
                    }
                } else {
                    start
                };
                if !p.is_finite() {
                    return Err(Error::Solve("non-finite ONE_WAY pore pressure".into()));
                }
                value.push(p);
            }
        }
        Ok(Some(value))
    }
    fn stage_water(&self, origin: &Committed, stage: &Stage, progress: f64) -> Option<f64> {
        match stage {
            Stage::HydrostaticWater { target_level_m, .. }
            | Stage::SteadySeepage {
                exterior_water_level_m: target_level_m,
                ..
            } => {
                let start = origin
                    .hydrostatic_water_level_m
                    .or(self.ground.initialization.water_level_m)
                    .unwrap_or(self.bounds[0][2]);
                Some(start + (target_level_m - start) * progress)
            }
            _ => origin.hydrostatic_water_level_m,
        }
    }
    fn stage_joint_water(
        &self,
        origin: &Committed,
        stage: &Stage,
        active: &[bool],
        flow: Option<&crate::joint_water::Field>,
        progress: f64,
    ) -> Result<Option<crate::joint_water::Field>> {
        let Some(wet) = &self.wet_joints else {
            return Ok(None);
        };
        let old = origin
            .joint_water
            .as_ref()
            .ok_or_else(|| Error::Input("missing committed JOINT water field".into()))?;
        if let Some(target) = flow {
            return Ok(Some(wet.interpolate(old, target, active, progress)?));
        }
        if matches!(stage, Stage::HydrostaticWater { .. }) {
            if origin.pore_pressure_kpa.is_some() {
                let target = wet.hydrostatic(active, self.stage_water(origin, stage, 1.))?;
                return Ok(Some(wet.interpolate(old, &target, active, progress)?));
            }
            return Ok(Some(
                wet.hydrostatic(active, self.stage_water(origin, stage, progress))?,
            ));
        }
        Ok(Some(wet.retain(old, active)?))
    }
    fn solve_flow(
        &self,
        flow: &crate::darcy::Case,
        origin: &Committed,
        stage: &Stage,
        active: &[bool],
        material_ids: &[u32],
    ) -> Result<(crate::darcy::Report, Option<crate::joint_water::Flow>)> {
        if let Some(wet) = &self.wet_joints {
            wet.solve_steady(
                flow,
                active,
                Some(material_ids),
                self.stage_water(origin, stage, 1.),
            )
            .map(|f| (f.solution.clone(), Some(f)))
        } else {
            crate::darcy::solve(self.mesh, flow, active, Some(material_ids)).map(|r| (r, None))
        }
    }
    fn reference_z(&self, state: &Committed, node: usize) -> f64 {
        state
            .display_reference_displacement_m
            .get(node)
            .or_else(|| self.reference_displacement_m.get(node))
            .map_or(0., |u| u[2])
    }
    fn reset_display(&self, origin: &Committed, stage: usize) -> Result<Trial> {
        let frame = self.frame("display_reset_check".into(), origin, origin, self.case)?;
        if frame.relative_force_residual > self.case.settings.relative_force_tolerance
            || frame.maximum_constraint_error_m > 1e-10
        {
            return Err(Error::Input(
                "display RESET requires an equilibrated committed state".into(),
            ));
        }
        let mut state = origin.clone();
        state.display_reference_displacement_m = state.displacement_m.clone();
        state.stage_index = stage;
        state.stage_progress = 1.;
        Ok(Trial {
            state,
            iterations: 0,
            residual: frame.relative_force_residual,
            linear_solves: vec![],
            contact_predictors: vec![],
        })
    }
    fn maximum_strain(&self, u: &[[f64; 3]], active: &[bool], relative: bool) -> f64 {
        let mut maximum: f64 = 0.;
        for (ei, e) in self.elements.iter().enumerate() {
            if active[ei] {
                let ue = DVector::from_fn(e.dofs.len(), |i, _| {
                    let g = e.dofs[i];
                    u[g / 3][g % 3]
                        - if relative {
                            self.reference_displacement_m
                                .get(g / 3)
                                .map_or(0., |r| r[g % 3])
                        } else {
                            0.
                        }
                });
                for q in &e.points {
                    maximum = maximum.max((&q.b * &ue).iter().map(|v| v.abs()).fold(0., f64::max));
                }
            }
        }
        maximum
    }
    fn maximum_physical_strain_from_seed(&self, state: &Committed) -> Result<f64> {
        let kinematic = self.maximum_strain(&state.displacement_m, &state.active_elements, false);
        match &state.mixed_mcc {
            Some(history) => Ok(kinematic.max(mixed::maximum_material_strain(
                self,
                &state.displacement_m,
                &state.active_elements,
                &history.cumulative_volume_strain_compression_positive,
            )?)),
            None => Ok(kinematic),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn physical_load(
        &self,
        active: &[bool],
        pressure: f64,
        original_surface: bool,
        factor: f64,
        material_ids: &[u32],
        beam_state: Option<&crate::beams::State>,
        water_level: Option<f64>,
        joint_water: Option<&crate::joint_water::Field>,
    ) -> Result<Vec<f64>> {
        let mut load = vec![0.; self.mesh.nodes.len() * 3];
        match (&self.beams, beam_state) {
            (Some(model), Some(state)) => {
                for (g, f) in model.attached_load(state, factor)? {
                    load[g] += f;
                }
            }
            (None, None) => {}
            _ => return Err(Error::Input("beam configuration/state mismatch".into())),
        }
        for (i, e) in self.elements.iter().enumerate() {
            if active[i] {
                let nodes = self.mesh.cell_nodes(i);
                let unit_weight =
                    self.ground.materials[self.material_index(material_ids, i)?].unit_weight_kn_m3;
                for q in &e.points {
                    for (j, &node) in nodes.iter().enumerate() {
                        load[3 * node + 2] -= q.weight_m3 * unit_weight * q.shape[j];
                    }
                }
            }
        }
        let mut hydrostatic_faces = BTreeSet::new();
        if let Some(cfg) = &self.ground.initialization.geostatic {
            for name in &cfg.hydrostatic_face_sets {
                for face in &self.mesh.face_sets[name] {
                    let mut face = *face;
                    face.sort_unstable();
                    hydrostatic_faces.insert(face);
                }
            }
        }
        for (face, owners) in &self.faces {
            let active_owners: Vec<_> = owners.iter().copied().filter(|&i| active[i]).collect();
            let opposite = self.joints.as_ref().and_then(|j| j.counterpart(face));
            if opposite.is_some_and(|i| active[i]) {
                continue;
            }
            if active_owners.len() == 1 {
                let nodes = self.mesh.face_nodes(active_owners[0], face)?;
                let points: Vec<_> = nodes.iter().map(|&i| self.mesh.nodes[i]).collect();
                // Hydrostatic tractions act on every new cut face, including vertical walls.
                // Outer side/base faces are supported model boundaries, not exposed water faces.
                let cut_face = owners.len() == 2 || opposite.is_some();
                let level = water_level.or(self.ground.initialization.water_level_m);
                let partial = level
                    .map(|h| crate::hydrostatic::depth_bounds(&points, h))
                    .transpose()?
                    .is_some_and(|b| b[0] < 0. && b[1] > 0.);
                let mut wet_face = false;
                for q in elements::face_points(&points)? {
                    let normal = q.area_vector;
                    let norm = normal.iter().map(|n| n * n).sum::<f64>().sqrt();
                    let horizontal_up =
                        normal[2] > 0. && normal[0].abs() + normal[1].abs() < 1e-8 * norm;
                    let surcharge = horizontal_up
                        && (!original_surface || near(q.position[2], self.bounds[1][2]));
                    let water_face = if self.ground.is_general() {
                        cut_face || hydrostatic_faces.contains(face)
                    } else {
                        cut_face || horizontal_up
                    };
                    wet_face |= water_face;
                    let water = if water_face && !partial {
                        self.pore(water_level, q.position[2])
                    } else {
                        0.
                    };
                    for (j, &node) in nodes.iter().enumerate() {
                        for axis in 0..3 {
                            load[3 * node + axis] -= q.shape[j]
                                * normal[axis]
                                * (water + if surcharge { pressure } else { 0. });
                        }
                    }
                }
                if partial && wet_face {
                    let integrated = crate::hydrostatic::face_load(
                        &points,
                        level.ok_or_else(|| Error::Input("missing partial-water level".into()))?,
                        self.ground.initialization.water_unit_weight_kn_m3,
                        crate::hydrostatic::Integration::default(),
                    )?;
                    for (j, &node) in nodes.iter().enumerate() {
                        for a in 0..3 {
                            load[3 * node + a] += integrated.nodal_force_kn[j][a];
                        }
                    }
                }
            }
        }
        if let Some(cfg) = &self.ground.initialization.geostatic {
            for (f, g) in load
                .iter_mut()
                .zip(cfg.base_loads.on_active(self.mesh, active)?)
            {
                *f += g;
            }
        }
        if let Some(wet) = &self.wet_joints {
            let field = if let Some(field) = joint_water {
                wet.retain(field, active)?
            } else {
                wet.hydrostatic(
                    active,
                    water_level.or(self.ground.initialization.water_level_m),
                )?
            };
            for (v, water) in load.iter_mut().zip(field.nodal_load_kn) {
                *v += water;
            }
        }
        if factor != 0.
            && let Some(boundary) = &self.case.boundary
        {
            let extra = boundary.loads(self.mesh, active)?;
            for (f, g) in load.iter_mut().zip(extra) {
                *f += factor * g;
            }
        }
        if !load.iter().all(|f| f.is_finite()) {
            return Err(Error::Input(
                "non-finite reconstructed physical load".into(),
            ));
        }
        Ok(load)
    }

    fn evaluate(
        &self,
        committed: &Committed,
        u: &[[f64; 3]],
        active: &[bool],
        dofs: &Dofs,
        settings: &Settings,
        mode: EvaluationMode<'_>,
    ) -> Result<Evaluation> {
        let projected_pressure = if self.case.is_projected_mcc() {
            if mode.tangent || u != committed.displacement_m.as_slice() {
                return Err(Error::Input(
                    "projected MCC increments require the explicit three-field backend".into(),
                ));
            }
            Some(mixed::work_pressure(self, committed, active)?)
        } else {
            None
        };
        let mut states = committed.material_states.clone();
        let mut internal = vec![0.; self.mesh.nodes.len() * 3];
        let mut global = if mode.tangent {
            Some(if let Some(pattern) = &dofs.sparse_pattern {
                Tangent::Sparse(CsrMatrix::zeros(pattern.clone()))
            } else {
                Tangent::Dense(DMatrix::zeros(dofs.count, dofs.count))
            })
        } else {
            None
        };
        // Bound transient responses to about4096 material points per batch.
        // Geometry and global accumulation remain in their original serial order.
        let max_points = self
            .elements
            .iter()
            .map(|e| e.points.len())
            .max()
            .unwrap_or(1)
            .max(1);
        let batch_elements = (4096 / max_points).max(1);
        for first in (0..self.elements.len()).step_by(batch_elements) {
            let last = (first + batch_elements).min(self.elements.len());
            let mut work = Vec::new();
            for (offset, e) in self.elements[first..last].iter().enumerate() {
                let ei = first + offset;
                if active[ei] {
                    let material = self.material_index(&committed.element_material_ids, ei)?;
                    let du = DVector::from_fn(e.dofs.len(), |i, _| {
                        let g = e.dofs[i];
                        u[g / 3][g % 3] - committed.displacement_m[g / 3][g % 3]
                    });
                    for (qi, q) in e.points.iter().enumerate() {
                        let point = e.point_start + qi;
                        let de = -&q.b * &du;
                        work.push(MaterialWork {
                            point,
                            material,
                            strain: std::array::from_fn(|i| de[i]),
                            partition: mode.partitions.map(|p| p[point]),
                        });
                    }
                }
            }
            let updates = material_work_batch(
                &self.material_models,
                self.material_transition_from.as_deref(),
                &committed.material_states,
                &work,
                (
                    mode.tangent,
                    settings.tangent_perturbation,
                    settings.material_maximum_substep,
                ),
                4,
            )?;
            let mut updates = updates.into_iter();
            for (offset, e) in self.elements[first..last].iter().enumerate() {
                let ei = first + offset;
                if active[ei] {
                    for (qi, q) in e.points.iter().enumerate() {
                        let point = e.point_start + qi;
                        let (update, d) = updates.next().ok_or_else(|| {
                            Error::Control("material response count differs".into())
                        })?;
                        states[point] = update.state;
                        let mut total = *states[point].stress();
                        if let Some(pressure) = &projected_pressure {
                            let mean = total.iter().take(3).sum::<f64>() / 3.;
                            for value in total.iter_mut().take(3) {
                                *value += pressure[point] - mean;
                            }
                        }
                        for s in total.iter_mut().take(3) {
                            *s += self.point_pore(committed, point, q.position[2]);
                        }
                        let force = -q.weight_m3
                            * q.b.transpose()
                            * SVector::<f64, 6>::from_row_slice(&total);
                        for i in 0..e.dofs.len() {
                            internal[e.dofs[i]] += force[i];
                        }
                        if let (Some(matrix), Some(d)) = (&mut global, d) {
                            let ke = q.weight_m3 * q.b.transpose() * d * &q.b;
                            for i in 0..e.dofs.len() {
                                for &(gi, ci) in &dofs.rows[e.dofs[i]] {
                                    for j in 0..e.dofs.len() {
                                        for &(gj, cj) in &dofs.rows[e.dofs[j]] {
                                            matrix.add(gi, gj, ci * cj * ke[(i, j)])?;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if updates.next().is_some() {
                return Err(Error::Control(
                    "material response count exceeds work".into(),
                ));
            }
        }
        let anchors = if let Some(model) = self.anchors.as_ref().filter(|_| mode.structures) {
            let state = committed
                .anchors
                .as_ref()
                .ok_or_else(|| Error::Input("missing anchor history".into()))?;
            let e = model.evaluate(state, u, active, mode.anchor_control)?;
            for (g, &f) in e.internal_kn.iter().enumerate() {
                internal[g] += f;
            }
            if let Some(matrix) = &mut global {
                for &(g, h, k) in &e.tangent {
                    for &(a, ca) in &dofs.rows[g] {
                        for &(b, cb) in &dofs.rows[h] {
                            matrix.add(a, b, ca * cb * k)?;
                        }
                    }
                }
            }
            Some(e)
        } else {
            None
        };
        let joints = if let Some(model) = self.joints.as_ref().filter(|_| mode.structures) {
            let state = committed
                .joints
                .as_ref()
                .ok_or_else(|| Error::Input("missing JOINT history".into()))?;
            let mut e = model.evaluate(
                state,
                u,
                active,
                mode.partitions.map(|p| &p[self.point_count..]),
                mode.tangent,
            )?;
            if let Some(wet) = &self.wet_joints {
                let field = wet.retain(
                    committed.joint_water.as_ref().ok_or_else(|| {
                        Error::Input("missing committed JOINT pressure field".into())
                    })?,
                    active,
                )?;
                wet.annotate(&mut e.report, &field)?;
            }
            for (g, &f) in e.internal_kn.iter().enumerate() {
                internal[g] += f;
            }
            if let Some(matrix) = &mut global {
                for &(g, h, k) in &e.tangent {
                    for &(a, ca) in &dofs.rows[g] {
                        for &(b, cb) in &dofs.rows[h] {
                            matrix.add(a, b, ca * cb * k)?;
                        }
                    }
                }
            }
            Some(e)
        } else {
            None
        };
        Ok(Evaluation {
            joints,
            anchors,
            beams: {
                if let Some(model) = self.beams.as_ref().filter(|_| mode.structures) {
                    let state = committed
                        .beams
                        .as_ref()
                        .ok_or_else(|| Error::Input("missing beam history".into()))?;
                    let beam = model.evaluate(
                        state,
                        u,
                        &dofs.active_nodes,
                        mode.factor.unwrap_or(committed.boundary_factor),
                    )?;
                    for (i, &g) in beam.interface_dofs.iter().enumerate() {
                        internal[g] += beam.interface_internal_kn[i];
                    }
                    if let Some(matrix) = &mut global {
                        for (i, &g) in beam.interface_dofs.iter().enumerate() {
                            for (j, &h) in beam.interface_dofs.iter().enumerate() {
                                let k = beam.interface_tangent_kn_m[(i, j)];
                                if k != 0. {
                                    for &(a, ca) in &dofs.rows[g] {
                                        for &(b, cb) in &dofs.rows[h] {
                                            matrix.add(a, b, ca * cb * k)?;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    Some(beam)
                } else {
                    None
                }
            },
            states,
            internal,
            tangent: global,
        })
    }

    fn reduce(&self, v: &[f64], dofs: &Dofs) -> DVector<f64> {
        let mut out = DVector::zeros(dofs.count);
        for (i, value) in v.iter().enumerate() {
            for &(r, c) in &dofs.rows[i] {
                out[r] += c * value;
            }
        }
        out
    }

    fn joint_support_space<'g>(
        &self,
        geometry: &'g crate::rigidity::Geometry,
        active: &[bool],
        dofs: &Dofs,
        eval: &Evaluation,
    ) -> Result<crate::rigidity::Space<'g>> {
        let mut space = geometry.constrained(&dofs.rows)?;
        if let (Some(model), Some(joints)) = (&self.joints, &eval.joints) {
            // Contact forces smaller than equilibrium uncertainty cannot
            // certify support (a neutral contact retains rounding pressure).
            let stress_scale = self
                .elements
                .iter()
                .enumerate()
                .filter(|(i, _)| active[*i])
                .flat_map(|(_, e)| &eval.states[e.point_start..e.point_start + e.points.len()])
                .map(|s| {
                    s.stress()
                        .iter()
                        .enumerate()
                        .map(|(i, x)| x * x * if i < 3 { 1. } else { 2. })
                        .sum::<f64>()
                        .sqrt()
                })
                .fold(1., f64::max);
            let pressure_tolerance =
                (self.case.settings.relative_force_tolerance * stress_scale).max(1e-10);
            model.add_rigid_support(&joints.state, active, pressure_tolerance, &mut space)?;
        }
        if let (Some(model), Some(anchors)) = (&self.anchors, &eval.anchors) {
            model.add_rigid_support(&anchors.state, &mut space)?;
        }
        if let Some(beams) = &eval.beams {
            let k = &beams.interface_tangent_kn_m;
            let scale = k.iter().map(|v| v.abs()).fold(0., f64::max);
            if scale > 0. {
                for i in 0..k.nrows() {
                    let row: Vec<_> = beams
                        .interface_dofs
                        .iter()
                        .enumerate()
                        .map(|(j, &g)| (g, k[(i, j)] / scale))
                        .collect();
                    space.physical_row(&row)?;
                }
            }
        }
        Ok(space)
    }
    fn predictor_partitions(
        &self,
        old: &Committed,
        u: &[[f64; 3]],
        active: &[bool],
        settings: &Settings,
    ) -> Result<Vec<usize>> {
        let mut required =
            vec![1; self.point_count + self.joints.as_ref().map_or(0, |j| j.point_count())];
        for (i, e) in self.elements.iter().enumerate() {
            if active[i] {
                let du = DVector::from_fn(e.dofs.len(), |j, _| {
                    let g = e.dofs[j];
                    u[g / 3][g % 3] - old.displacement_m[g / 3][g % 3]
                });
                for (qi, q) in e.points.iter().enumerate() {
                    required[e.point_start + qi] =
                        ((-&q.b * &du).iter().map(|v| v.abs()).sum::<f64>()
                            / settings.material_maximum_substep)
                            .ceil()
                            .max(1.) as usize;
                }
            }
        }
        if let Some(joints) = &self.joints {
            required[self.point_count..].copy_from_slice(
                &joints.partitions(
                    old.joints
                        .as_ref()
                        .ok_or_else(|| Error::Input("missing predictor history".into()))?,
                    u,
                    active,
                )?,
            );
        }
        Ok(required)
    }
    fn trial(&self, old: &Committed, target_data: Target, settings: &Settings) -> Result<Trial> {
        if self.case.is_mixed_mcc() {
            return mixed::trial(self, old, target_data, settings);
        }
        if self.joints.is_none() {
            return self.trial_with_predictor(old, target_data, settings, true, &mut false);
        }
        let mut used = false;
        let primary =
            self.trial_with_predictor(old, target_data.clone(), settings, true, &mut used);
        if used
            && matches!(
                &primary,
                Err(Error::Solve(_) | Error::Linear(_) | Error::Material(_))
            )
        {
            // An alternative gauge load projection is a numerical restart of
            // this trial. No iterate or partition from the failed trial is
            // reused; both strategies start at the same committed history.
            self.trial_with_predictor(old, target_data, settings, false, &mut false)
        } else {
            primary
        }
    }
    fn trial_with_predictor(
        &self,
        old: &Committed,
        target_data: Target,
        settings: &Settings,
        projected_load: bool,
        used: &mut bool,
    ) -> Result<Trial> {
        let Target {
            joint_water,
            pore_pressure_kpa,
            water_level,
            anchor_control,
            active,
            load: target,
            release,
            pressure,
            boundary_factor,
            stage,
            progress,
        } = target_data;
        let hydraulic_old;
        let old = if water_level != old.hydrostatic_water_level_m
            || pore_pressure_kpa != old.pore_pressure_kpa
            || joint_water != old.joint_water
        {
            hydraulic_old = Committed {
                joint_water: joint_water.clone(),
                pore_pressure_kpa: pore_pressure_kpa.clone(),
                hydrostatic_water_level_m: water_level,
                ..old.clone()
            };
            &hydraulic_old
        } else {
            old
        };
        let dofs = self.dofs(&active, settings, boundary_factor)?;
        let rigid_geometry = self
            .joints
            .as_ref()
            .map(|_| crate::rigidity::Geometry::new(self.mesh, &active))
            .transpose()?;
        let contact_trust = if let Some(geometry) = &rigid_geometry {
            !geometry.constrained(&dofs.rows)?.full()
        } else {
            false
        };
        // Use the same physical-load norm as the frame and independent audit.
        // Retaining a large previous load here can certify an under-resolved
        // endpoint after unloading or excavation to a much smaller load.
        let scale = norm(&target).max(1.);
        let mut u = old.displacement_m.clone();
        if let Some(transform) = &dofs.transform {
            let flat: Vec<_> = u.iter().flatten().copied().collect();
            let free = transform.expand(&transform.reduce(&flat));
            for g in 0..flat.len() {
                if dofs.active_nodes[g / 3] {
                    u[g / 3][g % 3] = transform.offset[g] + free[g];
                }
            }
        }
        let constrained_start = u.clone();
        let mut partitions =
            vec![1; self.point_count + self.joints.as_ref().map_or(0, |j| j.point_count())];
        if let Some(joints) = &self.joints {
            let needed = joints.partitions(
                old.joints
                    .as_ref()
                    .ok_or_else(|| Error::Input("missing JOINT history".into()))?,
                &u,
                &active,
            )?;
            partitions[self.point_count..].copy_from_slice(&needed);
        }
        if dofs.transform.is_some() {
            for (ei, e) in self.elements.iter().enumerate() {
                if active[ei] {
                    let du = DVector::from_iterator(
                        e.dofs.len(),
                        e.dofs
                            .iter()
                            .map(|&g| u[g / 3][g % 3] - old.displacement_m[g / 3][g % 3]),
                    );
                    for (qi, q) in e.points.iter().enumerate() {
                        let de = -&q.b * &du;
                        let needed = (de.iter().map(|x| x.abs()).sum::<f64>()
                            / settings.material_maximum_substep)
                            .ceil()
                            .max(1.) as usize;
                        if needed > 50000 {
                            return Err(Error::Material(
                                "prescribed increment exceeds MCC partition budget".into(),
                            ));
                        }
                        partitions[e.point_start + qi] = needed;
                    }
                }
            }
        }
        let mut linear_solves = Vec::new();
        let mut contact_predictors = Vec::new();
        for iteration in 0..=settings.maximum_newton_iterations {
            let mode = EvaluationMode {
                structures: true,
                anchor_control: anchor_control.as_ref(),
                tangent: false,
                partitions: Some(&partitions),
                factor: Some(boundary_factor),
            };
            let eval = self.evaluate(old, &u, &active, &dofs, settings, mode)?;
            let residual: Vec<_> = target
                .iter()
                .zip(&eval.internal)
                .map(|(a, b)| a - b)
                .collect();
            let reduced = self.reduce(&residual, &dofs);
            let relative = reduced.norm() / scale;
            self.notify(|| {
                let mut event = self.progress(crate::stage_control::Kind::Newton, stage, progress);
                event.newton_iteration = Some(iteration);
                event.relative_force_residual = Some(relative);
                event.maximum_material_partitions = partitions.iter().max().copied();
                event
            })?;
            if relative <= settings.relative_force_tolerance {
                if let Some(geometry) = &rigid_geometry {
                    self.joint_support_space(geometry, &active, &dofs, &eval)?
                        .require_supported()?;
                }
                if let Some(t) = &dofs.transform
                    && t.maximum_constraint_error(&u.iter().flatten().copied().collect::<Vec<_>>())
                        > 1e-10
                {
                    return Err(Error::Solve(
                        "prescribed displacement constraint residual exceeds tolerance".into(),
                    ));
                }
                return Ok(Trial {
                    state: Committed {
                        mixed_mcc: None,
                        joint_water,
                        pore_pressure_kpa,
                        hydrostatic_water_level_m: water_level,
                        joints: eval.joints.map(|j| j.state),
                        anchors: eval.anchors.map(|a| a.state),
                        beams: eval.beams.map(|b| b.state),
                        display_reference_displacement_m: old
                            .display_reference_displacement_m
                            .clone(),
                        element_material_ids: old.element_material_ids.clone(),
                        displacement_m: u,
                        material_states: eval.states,
                        active_elements: active,
                        applied_nodal_load_kn: target,
                        excavation_release_load_kn: if (stage == 0 && self.ground.is_general())
                            || (stage > 0 && self.case.stages[stage - 1].is_transition())
                        {
                            vec![0.; release.len()]
                        } else {
                            release.clone()
                        },
                        transition_release_load_kn: if stage > 0
                            && self.case.stages[stage - 1].is_transition()
                        {
                            release.clone()
                        } else {
                            vec![]
                        },
                        geostatic_release_load_kn: if stage == 0 && self.ground.is_general() {
                            release.clone()
                        } else {
                            vec![]
                        },
                        surface_pressure_kpa: pressure,
                        boundary_factor,
                        stage_index: stage,
                        stage_progress: progress,
                    },
                    iterations: iteration,
                    residual: relative,
                    linear_solves,
                    contact_predictors,
                });
            }
            if iteration == settings.maximum_newton_iterations {
                return Err(Error::Solve(format!(
                    "Newton iteration limit; material force residual={relative:e}"
                )));
            }
            let mut tangent_evaluation = self.evaluate(
                old,
                &u,
                &active,
                &dofs,
                settings,
                EvaluationMode {
                    structures: true,
                    anchor_control: anchor_control.as_ref(),
                    tangent: true,
                    partitions: Some(&partitions),
                    factor: Some(boundary_factor),
                },
            )?;
            let matrix = tangent_evaluation
                .tangent
                .take()
                .ok_or_else(|| Error::Solve("missing global tangent".into()))?;
            let gauges = if let Some(geometry) = &rigid_geometry {
                let space = self.joint_support_space(geometry, &active, &dofs, &eval)?;
                if space.full() {
                    None
                } else {
                    let modes =
                        geometry.reduced_modes(&space.null_modes(), &dofs.rows, dofs.count)?;
                    Some(crate::rigidity::gauge_coordinates(
                        &modes,
                        |x| matrix.action(x),
                        matrix.scale(),
                    )?)
                }
            } else {
                None
            };
            let cold = if iteration == 0 && gauges.is_some() {
                self.joints
                    .as_ref()
                    .map(|j| {
                        j.cold_activation_tangent(
                            old.joints.as_ref().ok_or_else(|| {
                                Error::Solve("missing cold contact history".into())
                            })?,
                            &active,
                        )
                    })
                    .transpose()?
                    .flatten()
            } else {
                None
            };
            let (correction, predictor_index) = if let Some(elastic) = cold {
                *used = true;
                let mut temporary = matrix.clone();
                let actual = tangent_evaluation
                    .joints
                    .as_ref()
                    .ok_or_else(|| Error::Solve("missing cold contact tangent".into()))?;
                for (sign, terms) in [(-1., actual.tangent.as_slice()), (1., elastic.as_slice())] {
                    for &(g, h, value) in terms {
                        for &(a, ca) in &dofs.rows[g] {
                            for &(b, cb) in &dofs.rows[h] {
                                temporary.add(a, b, sign * ca * cb * value)?;
                            }
                        }
                    }
                }
                let (x, report) = temporary.solve(&reduced, settings)?;
                let original = (matrix.action(&x)? - &reduced).norm() / reduced.norm().max(1e-30);
                let index = contact_predictors.len();
                contact_predictors.push(ContactPredictor {
                    elastic_activation: true,
                    projected_load: false,
                    elastic_neutral_correction_norm_m: x.norm(),
                    removed_neutral_correction_norm_m: 0.,
                    iteration,
                    zero_increment_free_dofs: vec![],
                    activation_predictor: original > settings.relative_linear_tolerance,
                    original_relative_linear_residual: original,
                    temporary_linear_solve: report,
                    full_relative_force_before: relative,
                    full_relative_force_after: None,
                    accepted_alpha: None,
                    repartitioned: false,
                });
                (x, Some(index))
            } else if gauges.as_ref().is_none_or(|g| g.coordinates.is_empty()) {
                let (x, report) = matrix.solve(&reduced, settings)?;
                linear_solves.push(report);
                (x, None)
            } else {
                *used = true;
                let gauges = gauges
                    .as_ref()
                    .ok_or_else(|| Error::Solve("missing contact gauges".into()))?;
                let coordinates = &gauges.coordinates;
                // First remove the unresolved rigid load from the range
                // equation. Dropping just one nodal equation would turn
                // that load into a spurious point reaction and strain.
                // The elastic neutral predictor below handles the original
                // rigid load, followed by the full physical force check.
                let mut rhs = if projected_load {
                    &reduced - &gauges.modes * (gauges.modes.transpose() * &reduced)
                } else {
                    reduced.clone()
                };
                for &i in coordinates {
                    rhs[i] = 0.;
                }
                let (mut x, temporary) = matrix.gauged(coordinates)?.solve(&rhs, settings)?;
                // Choose the minimum-norm correction, removing arbitrary
                // neutral motion induced by the temporary coordinate pins.
                let neutral = &gauges.modes * (gauges.modes.transpose() * &x);
                let removed_neutral_correction_norm_m = neutral.norm();
                x -= neutral;
                let elastic = self
                    .joints
                    .as_ref()
                    .ok_or_else(|| Error::Solve("missing predictor interface model".into()))?
                    .neutral_elastic_predictor(
                        (
                            old.joints.as_ref().ok_or_else(|| {
                                Error::Solve("missing old predictor history".into())
                            })?,
                            &eval
                                .joints
                                .as_ref()
                                .ok_or_else(|| {
                                    Error::Solve("missing current predictor history".into())
                                })?
                                .state,
                        ),
                        &active,
                        &dofs.rows,
                        &gauges.modes,
                        &x,
                        gauges.modes.transpose() * &reduced,
                    )?;
                let elastic_neutral_correction_norm_m = elastic.norm();
                x += &gauges.modes * elastic;
                let original = (matrix.action(&x)? - &reduced).norm() / reduced.norm().max(1e-30);
                if !original.is_finite() {
                    return Err(Error::Solve(
                        "nonfinite original contact predictor residual".into(),
                    ));
                }
                let activation = original > settings.relative_linear_tolerance;
                if !activation {
                    let mut report = temporary.clone();
                    report.method.push_str("_zero_increment_gauge");
                    report.true_relative_residual = original;
                    linear_solves.push(report);
                }
                let index = contact_predictors.len();
                contact_predictors.push(ContactPredictor {
                    elastic_activation: false,
                    projected_load,
                    elastic_neutral_correction_norm_m,
                    removed_neutral_correction_norm_m,
                    iteration,
                    zero_increment_free_dofs: coordinates.clone(),
                    activation_predictor: activation,
                    original_relative_linear_residual: original,
                    temporary_linear_solve: temporary,
                    full_relative_force_before: relative,
                    full_relative_force_after: None,
                    accepted_alpha: None,
                    repartitioned: false,
                });
                (x, Some(index))
            };
            self.notify(|| {
                let mut event =
                    self.progress(crate::stage_control::Kind::LinearSolved, stage, progress);
                event.newton_iteration = Some(iteration);
                if let Some(index) = predictor_index {
                    event.linear_solve =
                        Some(contact_predictors[index].temporary_linear_solve.clone());
                    event.reason = Some("auxiliary_contact_predictor".into());
                } else {
                    event.linear_solve = linear_solves.last().cloned();
                }
                event
            })?;
            let mut alpha = 1.;
            let mut accepted = false;
            let mut repartition = false;
            let mut best_candidate_relative = f64::INFINITY;
            let mut trust_rejections = 0;
            let current_required = if contact_trust {
                self.predictor_partitions(old, &u, &active, settings)?
            } else {
                vec![]
            };
            for _ in 0..16 {
                let mut candidate = u.clone();
                for (g, terms) in dofs.rows.iter().enumerate() {
                    candidate[g / 3][g % 3] +=
                        alpha * terms.iter().map(|&(r, c)| c * correction[r]).sum::<f64>();
                }
                if contact_trust {
                    let required =
                        match self.predictor_partitions(old, &candidate, &active, settings) {
                            Ok(required) => required,
                            Err(Error::Material(_)) => {
                                alpha *= 0.5;
                                continue;
                            }
                            Err(error) => return Err(error),
                        };
                    // Backtrack a wild Newton trial before it permanently
                    // expands the constitutive path partition. Current and
                    // accepted trials still use the same fixed partition.
                    if required
                        .iter()
                        .zip(&current_required)
                        .any(|(&n, &current)| n > current.saturating_mul(2).max(4))
                    {
                        trust_rejections += 1;
                        alpha *= 0.5;
                        continue;
                    }
                }
                for (ei, e) in self.elements.iter().enumerate() {
                    if active[ei] {
                        let du = DVector::from_fn(e.dofs.len(), |j, _| {
                            let g = e.dofs[j];
                            candidate[g / 3][g % 3] - old.displacement_m[g / 3][g % 3]
                        });
                        for (qi, q) in e.points.iter().enumerate() {
                            let point = e.point_start + qi;
                            let strain = -&q.b * &du;
                            let needed = (strain.iter().map(|x| x.abs()).sum::<f64>()
                                / settings.material_maximum_substep)
                                .ceil()
                                .max(1.) as usize;
                            if needed > partitions[point] {
                                if needed > 50_000 {
                                    return Err(Error::Material(
                                        "global increment exceeds MCC partition budget".into(),
                                    ));
                                }
                                partitions[point] = needed * 2;
                                repartition = true;
                            }
                        }
                    }
                }
                if let Some(joints) = &self.joints {
                    let needed =
                        joints.partitions(old.joints.as_ref().unwrap(), &candidate, &active)?;
                    for (i, count) in needed.into_iter().enumerate() {
                        let old_count = &mut partitions[self.point_count + i];
                        if count > *old_count {
                            *old_count = (count * 2).min(50_000);
                            repartition = true;
                        }
                    }
                }
                if repartition {
                    if let Some(i) = predictor_index {
                        contact_predictors[i].repartitioned = true;
                    }
                    u = constrained_start.clone();
                    break;
                }
                if let Ok(trial) = self.evaluate(
                    old,
                    &candidate,
                    &active,
                    &dofs,
                    settings,
                    EvaluationMode {
                        structures: true,
                        anchor_control: anchor_control.as_ref(),
                        tangent: false,
                        partitions: Some(&partitions),
                        factor: Some(boundary_factor),
                    },
                ) {
                    let r: Vec<_> = target
                        .iter()
                        .zip(&trial.internal)
                        .map(|(a, b)| a - b)
                        .collect();
                    best_candidate_relative =
                        best_candidate_relative.min(self.reduce(&r, &dofs).norm() / scale);
                    if self.reduce(&r, &dofs).norm() < reduced.norm() * (1. - 1e-4 * alpha) {
                        if let Some(i) = predictor_index {
                            contact_predictors[i].full_relative_force_after =
                                Some(self.reduce(&r, &dofs).norm() / scale);
                            contact_predictors[i].accepted_alpha = Some(alpha);
                        }
                        u = candidate;
                        accepted = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if repartition {
                continue;
            }
            if !accepted {
                return Err(Error::Solve(format!(
                    "global line search failed; committed state retained; iteration={iteration}, relative_force={relative:e}, best_candidate={best_candidate_relative:e}, correction_norm_m={:e}, trust_rejections={trust_rejections}, projected_load={projected_load}, predictor={:?}",
                    correction.norm(),
                    predictor_index.map(|i| {
                        let p = &contact_predictors[i];
                        (
                            p.iteration,
                            p.original_relative_linear_residual,
                            p.zero_increment_free_dofs.len(),
                            p.elastic_activation,
                        )
                    })
                )));
            }
        }
        Err(Error::Solve("Newton did not complete".into()))
    }

    fn frame(
        &self,
        name: String,
        state: &Committed,
        origin: &Committed,
        case: &StageCase,
    ) -> Result<Frame> {
        let dofs = self.dofs(
            &state.active_elements,
            &case.settings,
            state.boundary_factor,
        )?;
        let eval = self.evaluate(
            state,
            &state.displacement_m,
            &state.active_elements,
            &dofs,
            &case.settings,
            EvaluationMode::plain(),
        )?;
        let residual: Vec<_> = state
            .applied_nodal_load_kn
            .iter()
            .zip(&eval.internal)
            .map(|(a, b)| a - b)
            .collect();
        let joint_support = if self.joints.is_some() {
            Some(
                self.joint_support_space(
                    &crate::rigidity::Geometry::new(self.mesh, &state.active_elements)?,
                    &state.active_elements,
                    &dofs,
                    &eval,
                )?
                .require_supported()?,
            )
        } else {
            None
        };
        let mut applied = [0.; 3];
        let mut reactions = [0.; 3];
        let mut soil_support = vec![[0.; 3]; self.mesh.nodes.len()];
        let transformed_reactions = dofs.transform.as_ref().map(|t| {
            let r: Vec<_> = eval
                .internal
                .iter()
                .zip(&state.applied_nodal_load_kn)
                .map(|(i, f)| i - f)
                .collect();
            t.reactions(&r)
        });
        for (i, &f) in state.applied_nodal_load_kn.iter().enumerate() {
            if dofs.active_nodes[i / 3] {
                applied[i % 3] += f;
                if let Some(r) = &transformed_reactions {
                    reactions[i % 3] += r[i];
                    soil_support[i / 3][i % 3] = r[i];
                } else if dofs.map[i].is_none() {
                    reactions[i % 3] += eval.internal[i] - f;
                    soil_support[i / 3][i % 3] = eval.internal[i] - f;
                }
            }
        }
        let mut monitors = BTreeMap::new();
        let global_moment_balance_error_kn_m =
            if eval.beams.is_some() || eval.anchors.is_some() || eval.joints.is_some() {
                let mut moment = Vector3::zeros();
                for (n, p) in self.mesh.nodes.iter().enumerate() {
                    if dofs.active_nodes[n] {
                        let f = Vector3::new(
                            state.applied_nodal_load_kn[3 * n],
                            state.applied_nodal_load_kn[3 * n + 1],
                            state.applied_nodal_load_kn[3 * n + 2],
                        ) + Vector3::from(soil_support[n]);
                        moment += Vector3::from(*p).cross(&f);
                    }
                }
                if let Some(beam) = &eval.beams {
                    for a in 0..3 {
                        applied[a] += beam.report.external_force_kn[a]
                            - beam.report.interface_external_force_kn[a];
                        reactions[a] += beam.report.support_force_kn[a];
                        moment[a] += beam.report.external_moment_origin_kn_m[a]
                            - beam.report.interface_external_moment_origin_kn_m[a]
                            + beam.report.support_moment_origin_kn_m[a];
                    }
                }
                Some(moment.into())
            } else {
                None
            };
        let mut monitor_nodes = BTreeMap::new();
        for (name, &z) in &case.monitors {
            let nodes: Vec<_> = self
                .mesh
                .nodes
                .iter()
                .enumerate()
                .filter(|(i, p)| dofs.active_nodes[*i] && near(p[2], z))
                .map(|(i, _)| i)
                .collect();
            monitor_nodes.insert(name.clone(), nodes);
        }
        for (name, set) in &case.monitor_sets {
            let nodes = self
                .mesh
                .node_sets
                .get(set)
                .ok_or_else(|| Error::Input(format!("unknown monitor node set {set}")))?;
            if monitor_nodes.contains_key(name) || nodes.is_empty() {
                return Err(Error::Input(
                    "duplicate monitor name or empty monitor set".into(),
                ));
            }
            monitor_nodes.insert(
                name.clone(),
                nodes
                    .iter()
                    .copied()
                    .filter(|&i| dofs.active_nodes[i])
                    .collect(),
            );
        }
        for (name, nodes) in monitor_nodes {
            let mean = if nodes.is_empty() {
                None
            } else {
                Some(
                    nodes
                        .iter()
                        .map(|&i| state.displacement_m[i][2] - self.reference_z(state, i))
                        .sum::<f64>()
                        / nodes.len() as f64,
                )
            };
            let increment = if nodes.is_empty() {
                None
            } else {
                Some(
                    nodes
                        .iter()
                        .map(|&i| state.displacement_m[i][2] - origin.displacement_m[i][2])
                        .sum::<f64>()
                        / nodes.len() as f64,
                )
            };
            monitors.insert(
                name.clone(),
                Monitor {
                    active_node_count: nodes.len(),
                    mean_displacement_z_m: mean,
                    stage_increment_z_m: increment,
                    minimum_displacement_z_m: nodes
                        .iter()
                        .map(|&i| state.displacement_m[i][2] - self.reference_z(state, i))
                        .reduce(f64::min),
                    maximum_displacement_z_m: nodes
                        .iter()
                        .map(|&i| state.displacement_m[i][2] - self.reference_z(state, i))
                        .reduce(f64::max),
                },
            );
        }
        let max_yield = self
            .elements
            .iter()
            .enumerate()
            .filter(|(i, _)| state.active_elements[*i])
            .flat_map(|(ei, e)| {
                e.points.iter().enumerate().map(move |(qi, _)| {
                    self.material_index(&state.element_material_ids, ei)
                        .and_then(|mi| {
                            self.material_models[mi]
                                .normalized_yield(&state.material_states[e.point_start + qi])
                        })
                })
            })
            .try_fold(f64::NEG_INFINITY, |maximum, value| {
                value.map(|v| maximum.max(v))
            })?;
        let mut max_strain: f64 = 0.;
        for (ei, e) in self.elements.iter().enumerate() {
            if state.active_elements[ei] {
                let ue = DVector::from_fn(e.dofs.len(), |j, _| {
                    let g = e.dofs[j];
                    state.displacement_m[g / 3][g % 3]
                });
                for q in &e.points {
                    max_strain =
                        max_strain.max((&q.b * &ue).iter().map(|v| v.abs()).fold(0., f64::max));
                }
            }
        }
        Ok(Frame {
            joint_support,
            joint_results: eval.joints.map(|j| j.report),
            anchor_results: eval.anchors.map(|a| a.report),
            beam_results: eval.beams.map(|b| b.report),
            global_moment_balance_error_kn_m,
            name,
            state: state.clone(),
            relative_force_residual: self.reduce(&residual, &dofs).norm()
                / norm(&state.applied_nodal_load_kn).max(1.),
            applied_force_kn: applied,
            support_reaction_kn: reactions,
            force_balance_error_kn: std::array::from_fn(|i| applied[i] + reactions[i]),
            maximum_normalized_yield: max_yield,
            maximum_absolute_engineering_strain_from_k0: (!self.ground.is_general())
                .then_some(max_strain),
            maximum_absolute_material_engineering_strain_from_seed: state
                .mixed_mcc
                .as_ref()
                .map(|history| {
                    mixed::maximum_material_strain(
                        self,
                        &state.displacement_m,
                        &state.active_elements,
                        &history.cumulative_volume_strain_compression_positive,
                    )
                })
                .transpose()?,
            maximum_absolute_engineering_strain_from_seed: self
                .ground
                .is_general()
                .then_some(max_strain),
            maximum_absolute_engineering_strain_from_initial_equilibrium: self
                .ground
                .is_general()
                .then(|| self.maximum_strain(&state.displacement_m, &state.active_elements, true)),
            maximum_constraint_error_m: dofs.transform.as_ref().map_or(0., |t| {
                t.maximum_constraint_error(
                    &state
                        .displacement_m
                        .iter()
                        .flatten()
                        .copied()
                        .collect::<Vec<_>>(),
                )
            }),
            monitors,
        })
    }

    fn excavation(&self, old: &Committed, bodies: &[String]) -> Result<Vec<bool>> {
        let mut active = old.active_elements.clone();
        let mut removed = 0;
        for body in bodies {
            let indices = if let Some(set) = self.mesh.excavation_sets.get(body) {
                set.clone()
            } else {
                let region = self
                    .mesh
                    .regions
                    .iter()
                    .find(|r| &r.body == body)
                    .ok_or_else(|| Error::Input(format!("unknown excavation body {body}")))?;
                self.mesh
                    .region_ids
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| **r == region.id)
                    .map(|(i, _)| i)
                    .collect()
            };
            if indices.iter().any(|&i| !old.active_elements[i]) {
                return Err(Error::Input(format!(
                    "excavation body {body} overlaps previously removed soil"
                )));
            }
            let mut found = false;
            for i in indices {
                if active[i] {
                    active[i] = false;
                    found = true;
                    removed += 1;
                }
            }
            if !found {
                return Err(Error::Input(format!("body {body} already inactive")));
            }
        }
        if removed == 0 || !active.iter().any(|x| *x) {
            return Err(Error::Input("excavation must leave supported soil".into()));
        }
        let mut adjacency = vec![vec![]; active.len()];
        let mut accessible = vec![];
        let mut supported = self.explicitly_supported_elements(&active);
        if let Some(joints) = &self.joints {
            for &[a, b] in &joints.owner_pairs {
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
        }
        for (face, owners) in &self.faces {
            let opposite = self.joints.as_ref().and_then(|j| j.counterpart(face));
            if owners.len() == 2 {
                adjacency[owners[0]].push(owners[1]);
                adjacency[owners[1]].push(owners[0]);
            }
            for &ei in owners {
                if active[ei]
                    && face
                        .iter()
                        .all(|&n| near(self.mesh.nodes[n][2], self.bounds[0][2]))
                {
                    supported.push(ei);
                }
                if old.active_elements[ei] && !active[ei] {
                    let was_cut = owners.iter().any(|&i| !old.active_elements[i])
                        || opposite.is_some_and(|i| !old.active_elements[i]);
                    let was_top = if opposite.is_some() {
                        false
                    } else if self.ground.is_general() && owners.len() == 1 {
                        let nodes = self.mesh.face_nodes(ei, face)?;
                        let points: Vec<_> = nodes.iter().map(|&n| self.mesh.nodes[n]).collect();
                        elements::face_points(&points)?.iter().all(|q| {
                            let norm = q.area_vector.iter().map(|x| x * x).sum::<f64>().sqrt();
                            q.area_vector[2] > 1e-8 * norm
                        })
                    } else {
                        owners.len() == 1
                            && face
                                .iter()
                                .all(|&n| near(self.mesh.nodes[n][2], self.bounds[1][2]))
                    };
                    if was_cut || was_top {
                        accessible.push(ei);
                    }
                }
            }
        }
        for (mut stack, retained) in [(supported, true), (accessible, false)] {
            let mut visited = BTreeSet::new();
            while let Some(i) = stack.pop() {
                if !visited.insert(i) {
                    continue;
                }
                stack.extend(
                    adjacency[i]
                        .iter()
                        .copied()
                        .filter(|&j| old.active_elements[j] && active[j] == retained),
                );
            }
            if (0..active.len())
                .any(|i| old.active_elements[i] && active[i] == retained && !visited.contains(&i))
            {
                return Err(Error::Input(
                    if retained {
                        "excavation leaves soil disconnected from the supported base"
                    } else {
                        "excavation contains a buried cavity; remove soil from an exposed surface"
                    }
                    .into(),
                ));
            }
        }
        Ok(active)
    }
}

/// Explicit material epochs after a completed, unreduced construction run.
/// This experimental kernel performs equilibrium trials, not stability or Fs classification.
/// Each independent strength factor must construct its own controller from the same origin.
pub struct McParameterMechanics<'a> {
    observer: Option<&'a dyn crate::stage_control::Observer>,
    mesh: &'a Mesh,
    origin: &'a Run,
    parameters: Vec<crate::mc::Parameters>,
    committed: Committed,
    origin_sha256: String,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct McParameterStep {
    pub schema: String,
    pub origin_sha256: String,
    /// Material table order is exactly the original Ground material table order.
    pub parameters: Vec<crate::mc::Parameters>,
    pub frame: Frame,
    pub iterations: usize,
    pub linear_solves: Vec<LinearReport>,
}

/// Dedicated material-epoch checkpoint. It references the immutable original
/// completed Run by full hash, rather than rewriting a normal Ground/Checkpoint.
/// The original Run must be retained alongside this checkpoint for resumption.
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct McEpochCheckpoint {
    pub schema: String,
    pub payload_sha256: String,
    pub payload: McEpochPayload,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct McEpochPayload {
    pub identity: crate::checkpoint::Identity,
    pub origin_run_sha256: String,
    pub origin_state_sha256: String,
    pub parameters: Vec<crate::mc::Parameters>,
    pub committed: Committed,
}

pub(crate) fn validate_mc_epoch_parameters(
    from: &[crate::mc::Parameters],
    to: &[crate::mc::Parameters],
) -> Result<()> {
    if from.len() != to.len() {
        return Err(Error::Input(
            "MC transition material table length differs".into(),
        ));
    }
    for (from, to) in from.iter().zip(to) {
        to.validate()?;
        if from.young_kpa != to.young_kpa || from.poisson != to.poisson || from.flow != to.flow {
            return Err(Error::Input(
                "MC continuation holds elastic constants and flow family fixed".into(),
            ));
        }
    }
    Ok(())
}

impl<'a> McParameterMechanics<'a> {
    pub fn new(mesh: &'a Mesh, origin: &'a Run) -> Result<Self> {
        mesh.validate()?;
        origin.ground.validate()?;
        origin.case.validate()?;
        if (!origin.reference_displacement_m.is_empty()
            && origin.reference_displacement_m.len() != mesh.nodes.len())
            || !origin
                .reference_displacement_m
                .iter()
                .flatten()
                .all(|v| v.is_finite())
        {
            return Err(Error::Input(
                "invalid MC origin displacement reference".into(),
            ));
        }
        if !origin.complete
            || origin.failure_reason.is_some()
            || origin.continuation.is_some()
            || origin.case.is_mixed_mcc()
            || origin.case.joint_system.is_some()
            || origin.case.beam_system.is_some()
            || origin.case.anchor_system.is_some()
            || origin.case.has_seepage()
            || origin.case.has_staged_water()
            || origin.ground.initialization.water_level_m.is_some()
            || origin.provenance
                != crate::checkpoint::Identity::current(mesh, &origin.ground, &origin.case)?
        {
            return Err(Error::Input("MC continuation requires a complete dry, all-MC construction origin from this exact solver; structures and flow are not enabled".into()));
        }
        let parameters = origin
            .ground
            .materials
            .iter()
            .map(|m| match &m.model {
                soil::Model::MohrCoulomb(p) => Ok(p.clone()),
                _ => Err(Error::Input(
                    "MC continuation never converts MCC/IDH materials".into(),
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        let model = Model::new(mesh, &origin.ground, &origin.case)?;
        validate_committed(&model, &origin.last_committed, &origin.case)?;
        Ok(Self {
            observer: None,
            mesh,
            origin,
            parameters,
            committed: origin.last_committed.clone(),
            origin_sha256: crate::checkpoint::hash_json(&origin.last_committed)?,
        })
    }
    pub fn committed(&self) -> &Committed {
        &self.committed
    }
    pub fn parameters(&self) -> &[crate::mc::Parameters] {
        &self.parameters
    }
    pub fn set_observer(&mut self, observer: &'a dyn crate::stage_control::Observer) {
        self.observer = Some(observer);
    }

    /// Only private accepted state is serialised; no failed Newton iterate is reachable.
    pub fn checkpoint(&self) -> Result<McEpochCheckpoint> {
        let payload = McEpochPayload {
            identity: self.origin.provenance.clone(),
            origin_run_sha256: crate::checkpoint::hash_json(self.origin)?,
            origin_state_sha256: self.origin_sha256.clone(),
            parameters: self.parameters.clone(),
            committed: self.committed.clone(),
        };
        Ok(McEpochCheckpoint {
            schema: "soil-fem-mc-epoch-checkpoint/1".into(),
            payload_sha256: crate::checkpoint::hash_json(&payload)?,
            payload,
        })
    }

    /// Validate the original run, exact executable, seal, frozen construction
    /// data, current material admissibility and current global equilibrium.
    pub fn resume(mesh: &'a Mesh, origin: &'a Run, checkpoint: &McEpochCheckpoint) -> Result<Self> {
        let mut controller = Self::new(mesh, origin)?;
        let payload = &checkpoint.payload;
        if checkpoint.schema != "soil-fem-mc-epoch-checkpoint/1"
            || checkpoint.payload_sha256 != crate::checkpoint::hash_json(payload)?
            || payload.identity != origin.provenance
            || payload.origin_run_sha256 != crate::checkpoint::hash_json(origin)?
            || payload.origin_state_sha256 != controller.origin_sha256
        {
            return Err(Error::Input(
                "MC epoch checkpoint seal/origin/executable mismatch".into(),
            ));
        }
        validate_mc_epoch_parameters(&controller.parameters, &payload.parameters)?;
        let state = &payload.committed;
        if state.displacement_m.len() != mesh.nodes.len()
            || state.material_states.len() != origin.last_committed.material_states.len()
            || !state.displacement_m.iter().flatten().all(|v| v.is_finite())
        {
            return Err(Error::Input(
                "invalid MC epoch checkpoint state dimensions/nonfinite displacement".into(),
            ));
        }
        let permitted = Committed {
            displacement_m: state.displacement_m.clone(),
            material_states: state.material_states.clone(),
            ..origin.last_committed.clone()
        };
        if crate::checkpoint::hash_json(&permitted)? != crate::checkpoint::hash_json(state)? {
            return Err(Error::Input(
                "MC epoch checkpoint changed frozen construction data".into(),
            ));
        }
        let mut model = Model::new(mesh, &origin.ground, &origin.case)?;
        model.reference_displacement_m = origin.reference_displacement_m.clone();
        model.material_models = payload
            .parameters
            .iter()
            .cloned()
            .map(soil::Model::MohrCoulomb)
            .collect();
        for (ei, element) in model.elements.iter().enumerate() {
            let material = model.material_index(&state.element_material_ids, ei)?;
            for offset in 0..element.points.len() {
                let point = element.point_start + offset;
                let now = &state.material_states[point];
                let seed = &origin.last_committed.material_states[point];
                if !state.active_elements[ei] {
                    // Removed material retains its original dormant epoch; this
                    // controller never performs later reactivation/construction.
                    if crate::checkpoint::hash_json(now)? != crate::checkpoint::hash_json(seed)? {
                        return Err(Error::Input(
                            "MC epoch checkpoint changed dormant material history".into(),
                        ));
                    }
                    continue;
                }
                model.material_models[material].validate_state(now)?;
                let (soil::State::MohrCoulomb(now), soil::State::MohrCoulomb(seed)) = (now, seed)
                else {
                    return Err(Error::Input(
                        "MC epoch checkpoint history family differs".into(),
                    ));
                };
                if now.equivalent_plastic_strain + 1e-12 < seed.equivalent_plastic_strain {
                    return Err(Error::Input(
                        "MC epoch checkpoint lost cumulative plastic history".into(),
                    ));
                }
            }
        }
        let audit = model.frame(
            "mc_epoch_resume_audit".into(),
            state,
            &origin.last_committed,
            &origin.case,
        )?;
        if !audit.relative_force_residual.is_finite()
            || audit.relative_force_residual > origin.case.settings.relative_force_tolerance
            || !audit.maximum_normalized_yield.is_finite()
            || audit.maximum_normalized_yield > 1e-8
            || !audit.maximum_constraint_error_m.is_finite()
            || audit.maximum_constraint_error_m > 1e-10
        {
            return Err(Error::Input(
                "MC epoch checkpoint failed current-material equilibrium/yield/constraint audit"
                    .into(),
            ));
        }
        controller.parameters = payload.parameters.clone();
        controller.committed = payload.committed.clone();
        Ok(controller)
    }

    /// Immutable Newton trials use the accepted old material and history together.
    /// Only a converged, audited trial replaces BOTH; every error retains both.
    pub fn try_parameters(
        &mut self,
        parameters: Vec<crate::mc::Parameters>,
    ) -> Result<McParameterStep> {
        validate_mc_epoch_parameters(&self.parameters, &parameters)?;
        let mut model = Model::new(self.mesh, &self.origin.ground, &self.origin.case)?;
        model.observer = self.observer;
        model.reference_displacement_m = self.origin.reference_displacement_m.clone();
        model.material_models = parameters
            .iter()
            .cloned()
            .map(soil::Model::MohrCoulomb)
            .collect();
        model.material_transition_from = Some(
            self.parameters
                .iter()
                .cloned()
                .map(soil::Model::MohrCoulomb)
                .collect(),
        );
        let old = &self.committed;
        let trial = model.trial(
            old,
            Target {
                joint_water: old.joint_water.clone(),
                pore_pressure_kpa: old.pore_pressure_kpa.clone(),
                water_level: old.hydrostatic_water_level_m,
                anchor_control: None,
                active: old.active_elements.clone(),
                load: old.applied_nodal_load_kn.clone(),
                release: old.excavation_release_load_kn.clone(),
                pressure: old.surface_pressure_kpa,
                boundary_factor: old.boundary_factor,
                stage: old.stage_index,
                progress: old.stage_progress,
            },
            &self.origin.case.settings,
        )?;
        // Strength continuation is NOT another construction/load stage.
        // All construction bookkeeping and external data stay byte-for-byte frozen.
        let accepted = Committed {
            displacement_m: trial.state.displacement_m,
            material_states: trial.state.material_states,
            ..old.clone()
        };
        model.material_transition_from = None;
        let frame = model.frame(
            "mc_parameter_epoch".into(),
            &accepted,
            &self.origin.last_committed,
            &self.origin.case,
        )?;
        if frame.relative_force_residual > self.origin.case.settings.relative_force_tolerance
            || !frame.relative_force_residual.is_finite()
            || frame.maximum_normalized_yield > 1e-8
        {
            return Err(Error::Solve(
                "MC parameter epoch failed accepted-state force/yield audit".into(),
            ));
        }
        let result = McParameterStep {
            schema: "soil-fem-mc-parameter-step/1".into(),
            origin_sha256: self.origin_sha256.clone(),
            parameters: parameters.clone(),
            frame,
            iterations: trial.iterations,
            linear_solves: trial.linear_solves,
        };
        self.committed = result.frame.state.clone();
        self.parameters = parameters;
        Ok(result)
    }
}

/// Internal force-equilibrium kernel for a physical-time ONE_WAY controller.
/// Its RHS is explicitly F_physical - I_fluid. Material integration remains
/// the existing MCC/MC/IDH Newton. No hydraulic pressure is added a second time.
pub(crate) struct ExternalMechanics<'a> {
    model: Model<'a>,
}
pub(crate) struct ExternalTrial {
    pub frame: Frame,
    pub iterations: usize,
    pub linear_solves: Vec<LinearReport>,
}
impl<'a> ExternalMechanics<'a> {
    pub fn new(mesh: &'a Mesh, ground: &'a Ground, case: &'a StageCase) -> Result<Self> {
        mesh.validate()?;
        ground.validate()?;
        case.validate()?;
        if !case.stages.is_empty()
            || case.joint_system.is_some()
            || case.beam_system.is_some()
            || case.anchor_system.is_some()
            || !case.initially_inactive_bodies.is_empty()
            || !case.verification_reject_once_stages.is_empty()
            || ground.initialization.water_level_m.is_some()
            || ground
                .initialization
                .geostatic
                .as_ref()
                .is_some_and(|g| !g.hydrostatic_face_sets.is_empty() || g.verification_reject_once)
        {
            return Err(Error::Input("ONE_WAY requires dry-reference weights, explicit fluid/exterior loads and fixed active soil; construction/structures need their coupled policies".into()));
        }
        Ok(Self {
            model: Model::new(mesh, ground, case)?,
        })
    }
    pub fn seed(&self) -> Result<Committed> {
        let m = &self.model;
        let states = if m.ground.is_general() {
            geostatic::states(m.mesh, m.ground)?
        } else {
            ground::initialize(m.mesh, m.ground.clone())?
                .point_states
                .into_iter()
                .map(|p| p.constitutive_state.into_state())
                .collect()
        };
        let active = vec![true; m.elements.len()];
        let factor = m
            .ground
            .initialization
            .geostatic
            .as_ref()
            .map_or(0., |g| g.initial_boundary_factor);
        let load = m.physical_load(
            &active,
            m.ground.initialization.surface_surcharge_kpa,
            m.case.pressure_on_original_surface_only,
            factor,
            &[],
            None,
            None,
            None,
        )?;
        Ok(Committed {
            mixed_mcc: None,
            joint_water: None,
            pore_pressure_kpa: Some(vec![0.; m.point_count]),
            hydrostatic_water_level_m: None,
            joints: None,
            anchors: None,
            beams: None,
            display_reference_displacement_m: vec![],
            element_material_ids: vec![],
            displacement_m: vec![[0.; 3]; m.mesh.nodes.len()],
            material_states: states,
            active_elements: active,
            applied_nodal_load_kn: load,
            excavation_release_load_kn: vec![0.; m.mesh.nodes.len() * 3],
            geostatic_release_load_kn: vec![],
            transition_release_load_kn: vec![],
            surface_pressure_kpa: m.ground.initialization.surface_surcharge_kpa,
            boundary_factor: factor,
            stage_index: 0,
            stage_progress: 0.,
        })
    }
    pub fn layout(&self) -> Vec<PointInfo> {
        self.model.layout()
    }
    pub fn seed_on_domain(&self, active: &[bool]) -> Result<Committed> {
        if active.len() != self.model.elements.len() || !active.iter().any(|a| *a) {
            return Err(Error::Input("invalid initial mechanical domain".into()));
        }
        self.model.ensure_supported(active)?;
        let mut seed = self.seed()?;
        seed.active_elements = active.to_vec();
        seed.applied_nodal_load_kn = self.dry_load(&seed, active)?;
        Ok(seed)
    }
    pub fn body_selection(&self, bodies: &[String]) -> Result<Vec<usize>> {
        self.model.selected_bodies(bodies)
    }
    pub fn prepare_external_transition(
        &self,
        old: &Committed,
        stage: &Stage,
    ) -> Result<(Committed, Vec<usize>)> {
        if old.hydrostatic_water_level_m.is_some()
            || old
                .pore_pressure_kpa
                .as_ref()
                .is_none_or(|p| p.len() != self.model.point_count || p.iter().any(|p| *p != 0.))
        {
            return Err(Error::Input(
                "external transition requires separate fluid stress, not kernel pore pressure"
                    .into(),
            ));
        }
        if !matches!(stage, Stage::Birth { .. } | Stage::MaterialSet { .. }) {
            return Err(Error::Input("unsupported external transition".into()));
        }
        self.model.prepare_transition_mode(old, stage, true)
    }
    pub fn transition_holding(&self, prepared: &Committed, rhs: &[f64]) -> Result<(Vec<f64>, f64)> {
        let internal = self.internal(prepared)?;
        let dofs = self.model.dofs(
            &prepared.active_elements,
            &self.model.case.settings,
            prepared.boundary_factor,
        )?;
        let mut holding = internal
            .iter()
            .zip(rhs)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        for (g, h) in holding.iter_mut().enumerate() {
            if !dofs.active_nodes[g / 3] {
                *h = 0.;
            }
        }
        let residual = rhs
            .iter()
            .zip(&holding)
            .zip(internal)
            .map(|((f, h), i)| f + h - i)
            .collect::<Vec<_>>();
        Ok((
            holding,
            self.model.reduce(&residual, &dofs).norm() / norm(rhs).max(1.),
        ))
    }
    pub fn transition_trial(
        &self,
        old: &Committed,
        rhs: Vec<f64>,
        holding: Vec<f64>,
        progress: f64,
        name: String,
    ) -> Result<ExternalTrial> {
        let mut trial = self.excavation_trial(
            old,
            old.active_elements.clone(),
            rhs,
            holding.clone(),
            progress,
            name,
        )?;
        trial.frame.state.excavation_release_load_kn.fill(0.);
        trial.frame.state.transition_release_load_kn = holding;
        Ok(trial)
    }
    pub fn set_reference(&mut self, value: Vec<[f64; 3]>) {
        self.model.reference_displacement_m = value;
    }
    pub fn internal(&self, state: &Committed) -> Result<Vec<f64>> {
        let m = &self.model;
        let dofs = m.dofs(
            &state.active_elements,
            &m.case.settings,
            state.boundary_factor,
        )?;
        Ok(m.evaluate(
            state,
            &state.displacement_m,
            &state.active_elements,
            &dofs,
            &m.case.settings,
            EvaluationMode::plain(),
        )?
        .internal)
    }
    pub fn maximum_strain(&self, state: &Committed) -> f64 {
        self.model
            .maximum_strain(&state.displacement_m, &state.active_elements, false)
    }
    pub fn excavation_active(&self, old: &Committed, bodies: &[String]) -> Result<Vec<bool>> {
        self.model.excavation(old, bodies)
    }
    pub fn dry_load(&self, state: &Committed, active: &[bool]) -> Result<Vec<f64>> {
        self.model.physical_load(
            active,
            state.surface_pressure_kpa,
            self.model.case.pressure_on_original_surface_only,
            state.boundary_factor,
            &state.element_material_ids,
            None,
            None,
            None,
        )
    }
    pub fn cut_holding(
        &self,
        old: &Committed,
        active: &[bool],
        rhs: &[f64],
    ) -> Result<(Vec<f64>, f64)> {
        let m = &self.model;
        let dofs = m.dofs(active, &m.case.settings, old.boundary_factor)?;
        let removed = old
            .active_elements
            .iter()
            .zip(active)
            .map(|(a, b)| *a && !b)
            .collect::<Vec<_>>();
        let removed_force = m
            .evaluate(
                old,
                &old.displacement_m,
                &removed,
                &dofs,
                &m.case.settings,
                EvaluationMode::plain(),
            )?
            .internal;
        let mut holding = old
            .applied_nodal_load_kn
            .iter()
            .zip(rhs)
            .zip(removed_force)
            .map(|((a, b), r)| a - b - r)
            .collect::<Vec<_>>();
        for (g, h) in holding.iter_mut().enumerate() {
            if !dofs.active_nodes[g / 3] {
                *h = 0.;
            }
        }
        let retained = m
            .evaluate(
                old,
                &old.displacement_m,
                active,
                &dofs,
                &m.case.settings,
                EvaluationMode::plain(),
            )?
            .internal;
        let residual = rhs
            .iter()
            .zip(&holding)
            .zip(retained)
            .map(|((f, h), i)| f + h - i)
            .collect::<Vec<_>>();
        Ok((
            holding,
            m.reduce(&residual, &dofs).norm() / norm(rhs).max(1.),
        ))
    }
    pub fn excavation_trial(
        &self,
        old: &Committed,
        active: Vec<bool>,
        rhs: Vec<f64>,
        holding: Vec<f64>,
        progress: f64,
        name: String,
    ) -> Result<ExternalTrial> {
        let load = rhs.iter().zip(&holding).map(|(f, h)| f + h).collect();
        let mut trial = self.model.trial(
            old,
            Target {
                joint_water: None,
                pore_pressure_kpa: old.pore_pressure_kpa.clone(),
                water_level: None,
                anchor_control: None,
                active,
                load,
                release: holding.clone(),
                pressure: old.surface_pressure_kpa,
                boundary_factor: old.boundary_factor,
                stage: 0,
                progress,
            },
            &self.model.case.settings,
        )?;
        // This is construction, even when the seed used general geostatic
        // initialization. Do not label its holding force as initial gravity.
        trial.state.excavation_release_load_kn = holding;
        trial.state.geostatic_release_load_kn = vec![];
        let frame = self.model.frame(name, &trial.state, old, self.model.case)?;
        Ok(ExternalTrial {
            frame,
            iterations: trial.iterations,
            linear_solves: trial.linear_solves,
        })
    }
    pub fn trial(&self, old: &Committed, rhs: Vec<f64>, name: String) -> Result<ExternalTrial> {
        if rhs.len() != self.model.mesh.nodes.len() * 3 || !rhs.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("invalid external mechanical RHS".into()));
        }
        let trial = self.model.trial(
            old,
            Target {
                joint_water: None,
                pore_pressure_kpa: old.pore_pressure_kpa.clone(),
                water_level: None,
                anchor_control: None,
                active: old.active_elements.clone(),
                load: rhs,
                release: vec![0.; self.model.mesh.nodes.len() * 3],
                pressure: old.surface_pressure_kpa,
                boundary_factor: old.boundary_factor,
                stage: 0,
                progress: 0.,
            },
            &self.model.case.settings,
        )?;
        let frame = self.model.frame(name, &trial.state, old, self.model.case)?;
        Ok(ExternalTrial {
            frame,
            iterations: trial.iterations,
            linear_solves: trial.linear_solves,
        })
    }
}

pub fn solve(mesh: &Mesh, ground: Ground, case: StageCase) -> Result<Run> {
    solve_controlled(mesh, ground, case, None, None, None)
}

pub fn solve_observed(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    observer: &dyn crate::stage_control::Observer,
) -> Result<Run> {
    solve_controlled(mesh, ground, case, None, None, Some(observer))
}

/// Stop only at an accepted increment boundary and retain the full stage path.
pub fn solve_until(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    accepted_increments: usize,
) -> Result<Run> {
    if accepted_increments == 0 {
        return Err(Error::Input("checkpoint interval must be positive".into()));
    }
    solve_controlled(mesh, ground, case, None, Some(accepted_increments), None)
}

pub(crate) fn resume(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    previous: Run,
    limit: Option<usize>,
) -> Result<Run> {
    solve_controlled(mesh, ground, case, Some(previous), limit, None)
}

pub(crate) fn resume_observed(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    previous: Run,
    observer: &dyn crate::stage_control::Observer,
) -> Result<Run> {
    solve_controlled(mesh, ground, case, Some(previous), None, Some(observer))
}

pub(crate) fn solve_until_observed(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    accepted_increments: usize,
    observer: &dyn crate::stage_control::Observer,
) -> Result<Run> {
    if accepted_increments == 0 {
        return Err(Error::Input("checkpoint interval must be positive".into()));
    }
    solve_controlled(
        mesh,
        ground,
        case,
        None,
        Some(accepted_increments),
        Some(observer),
    )
}

pub(crate) fn resume_limited_observed(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    previous: Run,
    limit: Option<usize>,
    observer: &dyn crate::stage_control::Observer,
) -> Result<Run> {
    if limit == Some(0) {
        return Err(Error::Input("checkpoint interval must be positive".into()));
    }
    solve_controlled(mesh, ground, case, Some(previous), limit, Some(observer))
}

fn equilibrate(
    model: &Model<'_>,
    mut seed: Committed,
    previous: Option<geostatic::Report>,
    limit: Option<usize>,
) -> Result<(geostatic::Report, usize)> {
    let cfg = model
        .ground
        .initialization
        .geostatic
        .as_ref()
        .ok_or_else(|| Error::Input("missing geostatic configuration".into()))?;
    if cfg.initial_boundary_factor != 0. && model.case.boundary.is_none() {
        return Err(Error::Input(
            "geostatic boundary factor requires boundary conditions".into(),
        ));
    }
    let mut settings = model.case.settings.clone();
    settings.relative_force_tolerance = settings
        .relative_force_tolerance
        .min(model.ground.initialization.relative_equilibrium_tolerance);
    let dofs = model.dofs(&seed.active_elements, &settings, 0.)?;
    let internal = model
        .evaluate(
            &seed,
            &seed.displacement_m,
            &seed.active_elements,
            &dofs,
            &settings,
            EvaluationMode::plain(),
        )?
        .internal;
    let physical = model.physical_load(
        &seed.active_elements,
        seed.surface_pressure_kpa,
        model.case.pressure_on_original_surface_only,
        0.,
        &seed.element_material_ids,
        seed.beams.as_ref(),
        seed.hydrostatic_water_level_m,
        seed.joint_water.as_ref(),
    )?;
    let holding: Vec<_> = internal.iter().zip(&physical).map(|(a, b)| a - b).collect();
    seed.applied_nodal_load_kn = internal;
    seed.geostatic_release_load_kn = holding.clone();
    let mut report = if let Some(mut saved) = previous {
        if serde_json::to_vec(&saved.seed)? != serde_json::to_vec(&seed)?
            || saved.initial_holding_load_kn != holding
            || !saved.progress.is_finite()
            || !(0. ..=1.).contains(&saved.progress)
            || !saved.next_increment.is_finite()
            || !(settings.minimum_increment..=settings.initial_increment)
                .contains(&saved.next_increment)
            || saved.last_committed.stage_index != 0
            || saved.last_committed.stage_progress != saved.progress
            || saved.complete != (saved.progress == 1.)
            || saved.attempts.len() > 2000
            || (saved.verification_rejection_done && !cfg.verification_reject_once)
        {
            return Err(Error::Input(
                "invalid geostatic checkpoint origin/progress/continuation".into(),
            ));
        }
        let state = &saved.last_committed;
        let factor = saved.progress * cfg.initial_boundary_factor;
        let expected = model.physical_load(
            &state.active_elements,
            state.surface_pressure_kpa,
            model.case.pressure_on_original_surface_only,
            factor,
            &state.element_material_ids,
            state.beams.as_ref(),
            state.hydrostatic_water_level_m,
            state.joint_water.as_ref(),
        )?;
        let release: Vec<_> = holding.iter().map(|h| (1. - saved.progress) * h).collect();
        let expected: Vec<_> = expected.iter().zip(&release).map(|(f, h)| f + h).collect();
        if state.active_elements != seed.active_elements
            || !same_beam_installation(&state.beams, &seed.beams)
            || !same_anchor_basis(&state.anchors, &seed.anchors)
            || !state.display_reference_displacement_m.is_empty()
            || state.element_material_ids != seed.element_material_ids
            || !state.transition_release_load_kn.is_empty()
            || state.surface_pressure_kpa != seed.surface_pressure_kpa
            || state.boundary_factor != factor
            || state.geostatic_release_load_kn != release
            || !state.excavation_release_load_kn.iter().all(|v| *v == 0.)
            || state.applied_nodal_load_kn != expected
        {
            return Err(Error::Input(
                "invalid geostatic checkpoint holding/physical load path".into(),
            ));
        }
        validate_committed(model, state, model.case)?;
        let strain = model.maximum_physical_strain_from_seed(state)?;
        if (strain - saved.maximum_absolute_engineering_strain_from_seed).abs() > 1e-10
            || strain > cfg.maximum_initial_engineering_strain
        {
            return Err(Error::Input(
                "geostatic checkpoint strain/reference mismatch".into(),
            ));
        }
        let mut cursor = 0.;
        let mut increment = settings.initial_increment;
        let mut injected = false;
        for attempt in &saved.attempts {
            let target = if cursor + increment >= 1. - 1e-12 {
                1.
            } else {
                cursor + increment
            };
            if attempt.stage != "geostatic_equilibration"
                || attempt.progress_from != cursor
                || attempt.progress_to != target
            {
                return Err(Error::Input(
                    "geostatic checkpoint increment path mismatch".into(),
                ));
            }
            if attempt.accepted {
                cursor = target;
            } else {
                if !attempt.committed_state_unchanged_on_rejection {
                    return Err(Error::Input(
                        "geostatic checkpoint rejected trial changed history".into(),
                    ));
                }
                if attempt.reason == "verification_injected_rejection_after_convergence" {
                    if injected || !cfg.verification_reject_once {
                        return Err(Error::Input(
                            "invalid geostatic verification rejection".into(),
                        ));
                    }
                    injected = true;
                }
                increment *= 0.5;
            }
        }
        if cursor != saved.progress
            || increment != saved.next_increment
            || injected != saved.verification_rejection_done
        {
            return Err(Error::Input(
                "geostatic checkpoint path endpoint mismatch".into(),
            ));
        }
        saved.failure_reason = None;
        saved
    } else {
        geostatic::Report {
            complete: false,
            progress: 0.,
            next_increment: settings.initial_increment,
            verification_rejection_done: false,
            seed: seed.clone(),
            last_committed: seed,
            initial_holding_load_kn: holding.clone(),
            attempts: vec![],
            maximum_absolute_engineering_strain_from_seed: 0.,
            final_relative_physical_force_residual: None,
            failure_reason: None,
        }
    };
    let mut accepted = 0;
    while report.progress < 1. - 1e-12 {
        if report.attempts.len() >= 2000 {
            report.failure_reason = Some("geostatic attempt budget exhausted".into());
            break;
        }
        let sum = report.progress + report.next_increment;
        let next = if sum >= 1. - 1e-12 { 1. } else { sum };
        let active = report.last_committed.active_elements.clone();
        let factor = next * cfg.initial_boundary_factor;
        let physical = model.physical_load(
            &active,
            report.last_committed.surface_pressure_kpa,
            model.case.pressure_on_original_surface_only,
            factor,
            &report.last_committed.element_material_ids,
            report.last_committed.beams.as_ref(),
            report.last_committed.hydrostatic_water_level_m,
            report.last_committed.joint_water.as_ref(),
        )?;
        let release: Vec<_> = holding.iter().map(|h| (1. - next) * h).collect();
        let load: Vec<_> = physical.iter().zip(&release).map(|(f, h)| f + h).collect();
        let before = serde_json::to_vec(&report.last_committed)?;
        let load_changed = load != report.last_committed.applied_nodal_load_kn;
        model.notify(|| model.progress(crate::stage_control::Kind::AttemptStarted, 0, next))?;
        let trial = model.trial(
            &report.last_committed,
            Target {
                joint_water: report.last_committed.joint_water.clone(),
                pore_pressure_kpa: report.last_committed.pore_pressure_kpa.clone(),
                water_level: report.last_committed.hydrostatic_water_level_m,
                anchor_control: None,
                active,
                load,
                release,
                pressure: report.last_committed.surface_pressure_kpa,
                boundary_factor: factor,
                stage: 0,
                progress: next,
            },
            &settings,
        );
        let trial = match trial {
            Err(error @ Error::Control(_)) => return Err(error),
            other => other,
        };
        let changed_structure = trial.as_ref().is_ok_and(|t| {
            t.state.beams != report.last_committed.beams
                || t.state.anchors != report.last_committed.anchors
                || t.state.joints != report.last_committed.joints
        });
        let changed_joint_slip = trial
            .as_ref()
            .is_ok_and(|t| joint_slip_changed(&t.state.joints, &report.last_committed.joints));
        let injected =
            cfg.verification_reject_once && !report.verification_rejection_done && trial.is_ok();
        let (iterations, residual, changed, plastic, strain, linear) = if let Ok(t) = &trial {
            (
                Some(t.iterations),
                Some(t.residual),
                serde_json::to_vec(&t.state.material_states)?
                    != serde_json::to_vec(&report.last_committed.material_states)?,
                t.state
                    .material_states
                    .iter()
                    .zip(&report.last_committed.material_states)
                    .any(|(a, b)| a.plastic_strain() != b.plastic_strain()),
                model.maximum_physical_strain_from_seed(&t.state)?,
                t.linear_solves.clone(),
            )
        } else {
            (None, None, false, false, 0., vec![])
        };
        let strain_limit = strain > cfg.maximum_initial_engineering_strain;
        let contact_predictors = trial
            .as_ref()
            .map_or_else(|_| vec![], |t| t.contact_predictors.clone());
        match trial {
            Ok(t) if !injected && !strain_limit => {
                report.attempts.push(Attempt {
                    rejected_trial_water_changed: false,
                    contact_predictors,
                    stage: "geostatic_equilibration".into(),
                    progress_from: report.progress,
                    progress_to: next,
                    accepted: true,
                    reason: "converged_material_force".into(),
                    newton_iterations: iterations,
                    relative_force_residual: residual,
                    committed_state_unchanged_on_rejection: false,
                    rejected_trial_joint_slip_changed: false,
                    rejected_trial_structure_changed: false,
                    rejected_trial_material_changed: false,
                    rejected_trial_plastic_changed: false,
                    rejected_trial_active_set_changed: false,
                    rejected_trial_load_changed: false,
                    linear_solves: linear,
                });
                report.last_committed = t.state;
                report.progress = next;
                report.maximum_absolute_engineering_strain_from_seed = strain;
                accepted += 1;
                model.notify_attempt(0, report.attempts.last().unwrap())?;
                if limit.is_some_and(|n| accepted >= n) || model.stop_requested()? {
                    report.failure_reason = Some("checkpoint_requested".into());
                    break;
                }
            }
            result => {
                let reason = if strain_limit {
                    format!(
                        "geostatic initial strain {strain:e} exceeds configured limit {}",
                        cfg.maximum_initial_engineering_strain
                    )
                } else if injected {
                    report.verification_rejection_done = true;
                    "verification_injected_rejection_after_convergence".into()
                } else {
                    match result {
                        Err(e) => e.to_string(),
                        Ok(_) => "rejected".into(),
                    }
                };
                let unchanged = before == serde_json::to_vec(&report.last_committed)?;
                report.attempts.push(Attempt {
                    rejected_trial_water_changed: false,
                    contact_predictors,
                    stage: "geostatic_equilibration".into(),
                    progress_from: report.progress,
                    progress_to: next,
                    accepted: false,
                    reason: reason.clone(),
                    newton_iterations: iterations,
                    relative_force_residual: residual,
                    committed_state_unchanged_on_rejection: unchanged,
                    rejected_trial_joint_slip_changed: changed_joint_slip,
                    rejected_trial_structure_changed: changed_structure,
                    rejected_trial_material_changed: changed,
                    rejected_trial_plastic_changed: plastic,
                    rejected_trial_active_set_changed: false,
                    rejected_trial_load_changed: load_changed,
                    linear_solves: linear,
                });
                if !unchanged {
                    return Err(Error::Solve(
                        "geostatic rejection changed committed history".into(),
                    ));
                }
                model.notify_attempt(0, report.attempts.last().unwrap())?;
                if strain_limit {
                    report.failure_reason = Some(reason);
                    break;
                }
                report.next_increment *= 0.5;
                if report.next_increment < settings.minimum_increment {
                    report.failure_reason = Some(format!("geostatic cutback exhausted: {reason}"));
                    break;
                }
            }
        }
    }
    report.complete = report.progress == 1.;
    if report.complete {
        let physical = model.physical_load(
            &report.last_committed.active_elements,
            report.last_committed.surface_pressure_kpa,
            model.case.pressure_on_original_surface_only,
            report.last_committed.boundary_factor,
            &report.last_committed.element_material_ids,
            report.last_committed.beams.as_ref(),
            report.last_committed.hydrostatic_water_level_m,
            report.last_committed.joint_water.as_ref(),
        )?;
        let mut state = report.last_committed.clone();
        state.applied_nodal_load_kn = physical;
        let frame = model.frame(
            "geostatic_physical_equilibrium".into(),
            &state,
            &state,
            model.case,
        )?;
        report.final_relative_physical_force_residual = Some(frame.relative_force_residual);
        if frame.relative_force_residual > settings.relative_force_tolerance {
            report.complete = false;
            report.failure_reason = Some(
                "geostatic endpoint is not in physical equilibrium; no initial reference accepted"
                    .into(),
            );
        }
    }
    Ok((report, accepted))
}

fn result_schema(mesh: &Mesh, ground: &Ground, case: &StageCase) -> &'static str {
    if case.is_projected_mcc() {
        return if case.beam_system.is_some() {
            "soil-fem-projected-mcc-stage-result/2"
        } else {
            "soil-fem-projected-mcc-stage-result/1"
        };
    }
    if case.is_mixed_mcc() {
        return if case.beam_system.is_some() {
            "soil-fem-mixed-mcc-stage-result/2"
        } else {
            "soil-fem-mixed-mcc-stage-result/1"
        };
    }
    if case
        .beam_system
        .as_ref()
        .is_some_and(|s| !s.plates.is_empty())
    {
        return "soil-fem-stage-result/14";
    }
    if case
        .joint_system
        .as_ref()
        .is_some_and(|c| c.joints.iter().any(|d| d.water.is_some()))
    {
        "soil-fem-stage-result/13"
    } else if case.has_seepage() {
        "soil-fem-stage-result/12"
    } else if case.has_staged_water() {
        "soil-fem-stage-result/11"
    } else if ground
        .materials
        .iter()
        .any(|m| matches!(m.model, soil::Model::InvariantDoubleHardening(_)))
    {
        "soil-fem-stage-result/10"
    } else if case.joint_system.is_some() {
        "soil-fem-stage-result/9"
    } else if case.anchor_system.is_some() {
        "soil-fem-stage-result/8"
    } else if case.beam_system.is_some() {
        "soil-fem-stage-result/7"
    } else if case
        .stages
        .iter()
        .any(|s| matches!(s, Stage::ResetDisplay { .. } | Stage::ResetState { .. }))
    {
        "soil-fem-stage-result/6"
    } else if case.has_lifecycle() {
        "soil-fem-stage-result/5"
    } else if ground.is_general() {
        "soil-fem-stage-result/4"
    } else if ground.materials.iter().any(|m| !m.model.is_mcc()) {
        "soil-fem-stage-result/3"
    } else if mesh.tetrahedra10.is_empty() {
        "soil-fem-stage-result/1"
    } else {
        "soil-fem-stage-result/2"
    }
}

fn solve_controlled(
    mesh: &Mesh,
    ground: Ground,
    case: StageCase,
    mut previous: Option<Run>,
    stop_after: Option<usize>,
    observer: Option<&dyn crate::stage_control::Observer>,
) -> Result<Run> {
    case.validate()?;
    let states = if ground.is_general() {
        geostatic::states(mesh, &ground)?
    } else {
        ground::initialize(mesh, ground.clone())?
            .point_states
            .into_iter()
            .map(|p| p.constitutive_state.into_state())
            .collect()
    };
    let provenance = crate::checkpoint::Identity::current(mesh, &ground, &case)?;
    if let Some(p) = &previous {
        if p.schema != result_schema(mesh, &ground, &case)
            || p.complete
            || p.provenance != provenance
            || p.failure_reason.as_deref() != Some("checkpoint_requested")
            || serde_json::to_vec(&p.case)? != serde_json::to_vec(&case)?
            || serde_json::to_vec(&p.ground)? != serde_json::to_vec(&ground)?
            || p.geostatic.is_some() != ground.is_general()
        {
            return Err(Error::Input(
                "checkpoint stage result/inputs mismatch".into(),
            ));
        }
        if let Some(g) = &p.geostatic {
            if !g.complete
                && (p.continuation.is_some()
                    || !p.frames.is_empty()
                    || serde_json::to_vec(&p.last_committed)?
                        != serde_json::to_vec(&g.last_committed)?
                    || !p.reference_displacement_m.is_empty())
            {
                return Err(Error::Input(
                    "geostatic checkpoint phase/state mismatch".into(),
                ));
            }
            if g.complete && p.reference_displacement_m != g.last_committed.displacement_m {
                return Err(Error::Input(
                    "initial reference displacement mismatch".into(),
                ));
            }
        }
    }
    let mut model = Model::new(mesh, &ground, &case)?;
    model.observer = observer;
    if let Some(p) = &previous
        && serde_json::to_vec(&p.joint_integration_points)?
            != serde_json::to_vec(&model.joint_layout())?
    {
        return Err(Error::Input(
            "checkpoint JOINT integration-point layout mismatch".into(),
        ));
    }
    // A restored Run is already an accepted state. Do not replay its geostatic
    // preamble after a stop has arrived: that would discard construction state.
    if previous.is_some() && model.stop_requested()? {
        return previous.ok_or_else(|| Error::Control("accepted resume state disappeared".into()));
    }
    let resuming_construction = previous.as_ref().is_some_and(|p| p.continuation.is_some());
    let mut active = vec![true; mesh.tetrahedra.len()];
    if !case.initially_inactive_bodies.is_empty() {
        if !ground.is_general() {
            return Err(Error::Input(
                "initially inactive soil requires general equilibration".into(),
            ));
        }
        for i in model.selected_bodies(&case.initially_inactive_bodies)? {
            active[i] = false;
        }
        model.ensure_supported(&active)?;
    }
    let element_material_ids = if case.has_lifecycle() {
        mesh.material_ids.clone()
    } else {
        vec![]
    };
    if case.stages.iter().any(|s| matches!(s, Stage::Excavate { bodies, .. } if bodies.iter().any(|b| mesh.excavation_sets.contains_key(b)))) && !case.pressure_on_original_surface_only {
        return Err(Error::Input("local excavation requires pressure_on_original_surface_only=true; floor surcharge is zero".into()));
    }
    let initial_beams = model.beams.as_ref().map(|m| m.initial_state());
    let initial_water = if case.has_staged_water() {
        Some(
            ground
                .initialization
                .water_level_m
                .unwrap_or(model.bounds[0][2]),
        )
    } else {
        None
    };
    let initial_joint_water = model
        .wet_joints
        .as_ref()
        .map(|w| {
            w.hydrostatic(
                &active,
                initial_water.or(ground.initialization.water_level_m),
            )
        })
        .transpose()?;
    let load = model.physical_load(
        &active,
        ground.initialization.surface_surcharge_kpa,
        case.pressure_on_original_surface_only,
        0.,
        &element_material_ids,
        initial_beams.as_ref(),
        initial_water,
        initial_joint_water.as_ref(),
    )?;
    let mut committed = Committed {
        mixed_mcc: None,
        joint_water: initial_joint_water,
        pore_pressure_kpa: None,
        hydrostatic_water_level_m: initial_water,
        joints: model
            .joints
            .as_ref()
            .map(|j| j.initial_state())
            .transpose()?,
        anchors: model.anchors.as_ref().map(|a| a.initial_state()),
        beams: initial_beams,
        display_reference_displacement_m: vec![],
        element_material_ids,
        displacement_m: vec![[0.; 3]; mesh.nodes.len()],
        material_states: states,
        active_elements: active,
        applied_nodal_load_kn: load,
        excavation_release_load_kn: vec![0.; mesh.nodes.len() * 3],
        geostatic_release_load_kn: vec![],
        transition_release_load_kn: vec![],
        surface_pressure_kpa: ground.initialization.surface_surcharge_kpa,
        boundary_factor: 0.,
        stage_index: 0,
        stage_progress: 0.,
    };
    if case.is_mixed_mcc() {
        mixed::initialize(&model, &mut committed)?;
    }
    if let Some(model_beams) = &model.beams {
        let dofs = model.dofs(&committed.active_elements, &case.settings, 0.)?;
        committed.beams = Some(
            model_beams
                .evaluate(
                    committed.beams.as_ref().unwrap(),
                    &committed.displacement_m,
                    &dofs.active_nodes,
                    0.,
                )?
                .state,
        );
    }
    let mut geostatic_report = None;
    let mut accepted_initialization = 0;
    if ground.is_general() {
        let was_initializing = previous
            .as_ref()
            .and_then(|p| p.geostatic.as_ref())
            .is_some_and(|g| !g.complete);
        let saved = previous.as_mut().and_then(|p| p.geostatic.take());
        let (report, count) = equilibrate(&model, committed, saved, stop_after)?;
        accepted_initialization = count;
        committed = report.last_committed.clone();
        if report.complete {
            model.reference_displacement_m = committed.displacement_m.clone();
            committed.stage_progress = 0.;
            committed.geostatic_release_load_kn.clear();
        }
        if !report.complete
            || stop_after.is_some_and(|n| count >= n)
            || (!resuming_construction && model.stop_requested()?)
        {
            let frames = if report.complete {
                vec![model.frame("initial_equilibrated".into(), &committed, &committed, &case)?]
            } else {
                vec![]
            };
            let joint_integration_points = model.joint_layout();
            let layout = model.layout();
            let reference_displacement_m = model.reference_displacement_m.clone();
            let failure_reason = report.failure_reason.clone().or_else(|| {
                Some(
                    if report.complete {
                        "checkpoint_requested"
                    } else {
                        "geostatic initialization incomplete"
                    }
                    .into(),
                )
            });
            return Ok(Run {
                flow_audits: vec![],
                joint_integration_points,
                schema: result_schema(mesh, &ground, &case).into(),
                provenance,
                complete: false,
                failure_reason,
                case,
                ground,
                frames,
                attempts: vec![],
                excavation_audits: vec![],
                transition_audits: vec![],
                last_committed: committed,
                integration_points: layout,
                continuation: None,
                geostatic: Some(report),
                reference_displacement_m,
            });
        }
        geostatic_report = Some(report);
        if was_initializing {
            previous = None;
        }
    }
    let mut frames = vec![model.frame(
        if ground.is_general() {
            "initial_equilibrated".into()
        } else {
            "initial_k0".into()
        },
        &committed,
        &committed,
        &case,
    )?];
    if previous
        .as_ref()
        .is_some_and(|p| ground.is_general() && p.continuation.is_none())
    {
        let p = previous
            .as_ref()
            .ok_or_else(|| Error::Input("missing checkpoint".into()))?;
        if p.frames.len() != 1
            || serde_json::to_vec(&p.frames)? != serde_json::to_vec(&frames)?
            || serde_json::to_vec(&p.last_committed)? != serde_json::to_vec(&committed)?
            || !p.attempts.is_empty()
            || !p.excavation_audits.is_empty()
            || !p.flow_audits.is_empty()
        {
            return Err(Error::Input(
                "invalid checkpoint at initialization/construction boundary".into(),
            ));
        }
        previous = None;
    }
    if frames[0].relative_force_residual > case.settings.relative_force_tolerance {
        return Err(Error::Input(
            "initial K0 state is not in equilibrium for selected kinematics".into(),
        ));
    }
    let mut geometry_check = committed.clone();
    for stage in case.stages.iter().filter(|_| !case.has_lifecycle()) {
        geometry_check.surface_pressure_kpa = stage.pressure(geometry_check.surface_pressure_kpa);
        if let Stage::Boundary { target_factor, .. } = stage {
            geometry_check.boundary_factor = *target_factor;
        }
        geometry_check.hydrostatic_water_level_m = model.stage_water(&geometry_check, stage, 1.);
        if let Stage::Excavate { bodies, .. } = stage {
            geometry_check.active_elements = model.excavation(&geometry_check, bodies)?;
        }
        model.dofs(
            &geometry_check.active_elements,
            &case.settings,
            geometry_check.boundary_factor,
        )?;
        model.physical_load(
            &geometry_check.active_elements,
            geometry_check.surface_pressure_kpa,
            case.pressure_on_original_surface_only,
            geometry_check.boundary_factor,
            &geometry_check.element_material_ids,
            geometry_check.beams.as_ref(),
            geometry_check.hydrostatic_water_level_m,
            geometry_check.joint_water.as_ref(),
        )?;
    }
    let mut attempts = Vec::new();
    let mut flow_audits = Vec::new();
    let mut excavation_audits = Vec::new();
    let mut transition_audits = Vec::new();
    let mut failure = None;
    let mut injected = BTreeSet::new();
    let mut continuation = None;
    if let Some(mut previous) = previous {
        if previous.schema != result_schema(mesh, &ground, &case)
            || previous.complete
            || previous.provenance != provenance
            || previous.failure_reason.as_deref() != Some("checkpoint_requested")
            || serde_json::to_vec(&previous.case)? != serde_json::to_vec(&case)?
            || serde_json::to_vec(&previous.ground)? != serde_json::to_vec(&ground)?
        {
            return Err(Error::Input(
                "checkpoint stage result/inputs mismatch".into(),
            ));
        }
        if serde_json::to_vec(&previous.integration_points)? != serde_json::to_vec(&model.layout())?
        {
            return Err(Error::Input(
                "checkpoint integration-point layout mismatch".into(),
            ));
        }
        let runtime = previous
            .continuation
            .take()
            .ok_or_else(|| Error::Input("checkpoint missing stage continuation".into()))?;
        if runtime.stage_index >= case.stages.len()
            || !runtime.progress.is_finite()
            || !(0. ..=1.).contains(&runtime.progress)
            || !runtime.increment.is_finite()
            || runtime.increment < case.settings.minimum_increment
            || runtime.increment > 1.
            || runtime.stage_attempts > 2000
            || runtime
                .injected_stages
                .iter()
                .any(|name| !case.verification_reject_once_stages.contains(name))
            || previous.frames.len() != runtime.stage_index + 1
            || previous.last_committed.stage_index != runtime.stage_index + 1
            || previous.last_committed.stage_progress != runtime.progress
            || previous.frames.is_empty()
            || serde_json::to_vec(&previous.frames[0].state)? != serde_json::to_vec(&committed)?
            || serde_json::to_vec(&previous.frames[runtime.stage_index].state)?
                != serde_json::to_vec(&runtime.origin)?
        {
            return Err(Error::Input(
                "invalid checkpoint stage/progress/origin".into(),
            ));
        }
        validate_committed(&model, &previous.last_committed, &case)?;
        mixed::validate_domain_origin(&model, &runtime.origin, &previous.last_committed)?;
        if case.has_seepage() {
            validate_committed(&model, &runtime.origin, &case)?;
        }
        committed = previous.last_committed;
        frames = previous.frames;
        attempts = previous.attempts;
        excavation_audits = previous.excavation_audits;
        transition_audits = previous.transition_audits;
        flow_audits = previous.flow_audits;
        injected = runtime.injected_stages.clone();
        continuation = Some(runtime);
    }
    let start_stage = continuation.as_ref().map_or(0, |c| c.stage_index);
    if continuation.is_some() {
        let expected: Vec<_> = case
            .stages
            .iter()
            .enumerate()
            .take(start_stage + 1)
            .filter_map(|(i, s)| matches!(s, Stage::SteadySeepage { .. }).then_some(i))
            .collect();
        if flow_audits
            .iter()
            .map(|a| a.stage_index)
            .collect::<Vec<_>>()
            != expected
        {
            return Err(Error::Input(
                "checkpoint hydraulic audit sequence mismatch".into(),
            ));
        }
        for audit in flow_audits.iter().filter(|a| a.stage_index < start_stage) {
            let Stage::SteadySeepage { flow, .. } = &case.stages[audit.stage_index] else {
                return Err(Error::Input("checkpoint hydraulic stage mismatch".into()));
            };
            let start = &frames[audit.stage_index].state;
            let (check, joint_check) = model.solve_flow(
                flow,
                start,
                &case.stages[audit.stage_index],
                &start.active_elements,
                &start.element_material_ids,
            )?;
            if audit.stage != case.stages[audit.stage_index].name()
                || !crate::darcy::same_solution(&check, &audit.solution)?
                || !crate::joint_water::same_flow(&joint_check, &audit.joint_flow)?
            {
                return Err(Error::Input(
                    "checkpoint completed hydraulic solution mismatch".into(),
                ));
            }
        }
    }
    let mut accepted_this_run = accepted_initialization;
    'stages: for (stage_index, stage) in case.stages.iter().enumerate().skip(start_stage) {
        let resumed = continuation.take();
        let origin = resumed
            .as_ref()
            .map_or_else(|| committed.clone(), |c| c.origin.clone());
        let (prepared, selected) = model.prepare_transition(&origin, stage)?;
        let active = match stage {
            Stage::Pressure { .. }
            | Stage::SteadySeepage { .. }
            | Stage::HydrostaticWater { .. }
            | Stage::Boundary { .. }
            | Stage::ResetDisplay { .. }
            | Stage::JackAnchors { .. } => origin.active_elements.clone(),
            Stage::Birth { .. }
            | Stage::MaterialSet { .. }
            | Stage::ResetState { .. }
            | Stage::InstallBeams { .. }
            | Stage::RemoveBeams { .. }
            | Stage::InstallAnchors { .. }
            | Stage::RemoveAnchors { .. } => prepared.active_elements.clone(),
            Stage::Excavate { bodies, .. } => model.excavation(&origin, bodies)?,
        };
        let flow_target = if let Stage::SteadySeepage { flow, .. } = stage {
            let solved = model.solve_flow(
                flow,
                &origin,
                stage,
                &active,
                &prepared.element_material_ids,
            );
            let (solution, joint_flow) = match solved {
                Ok(value) => value,
                Err(Error::Input(message)) => return Err(Error::Input(message)),
                Err(error) => {
                    failure = Some(format!("{}: hydraulic solve failed: {error}", stage.name()));
                    break;
                }
            };
            let audit = FlowAudit {
                joint_flow,
                stage: stage.name().into(),
                stage_index,
                solution,
            };
            if resumed.is_some() {
                let saved = flow_audits.last().ok_or_else(|| {
                    Error::Input("checkpoint missing current flow solution".into())
                })?;
                if saved.stage != audit.stage
                    || saved.stage_index != stage_index
                    || !crate::darcy::same_solution(&saved.solution, &audit.solution)?
                    || !crate::joint_water::same_flow(&saved.joint_flow, &audit.joint_flow)?
                {
                    return Err(Error::Input(
                        "checkpoint hydraulic target reconstruction mismatch".into(),
                    ));
                }
            } else {
                flow_audits.push(audit);
            }
            let solution = &flow_audits
                .last()
                .ok_or_else(|| Error::Input("missing flow solution".into()))?
                .solution;
            if !solution.complete {
                failure = Some(format!(
                    "{}: {}",
                    stage.name(),
                    solution
                        .failure_reason
                        .as_deref()
                        .unwrap_or("incomplete steady flow")
                ));
                break;
            }
            let mut target = vec![0.; model.point_count];
            for p in &solution.points {
                target[model.elements[p.element_index].point_start + p.local_index] =
                    p.pressure_head_m * ground.initialization.water_unit_weight_kn_m3;
            }
            Some(target)
        } else {
            None
        };
        let end_factor = if let Stage::Boundary { target_factor, .. } = stage {
            *target_factor
        } else {
            origin.boundary_factor
        };
        let end_pressure = stage.pressure(origin.surface_pressure_kpa);
        let joint_flow_target = if matches!(stage, Stage::SteadySeepage { .. }) {
            flow_audits
                .last()
                .and_then(|a| a.joint_flow.as_ref())
                .map(|f| {
                    model
                        .wet_joints
                        .as_ref()
                        .ok_or_else(|| Error::Input("unexpected JOINT flow".into()))?
                        .from_heads(
                            &active,
                            &f.original_field,
                            model.stage_water(&origin, stage, 1.),
                        )
                })
                .transpose()?
        } else {
            None
        };
        let end_joint_water =
            model.stage_joint_water(&origin, stage, &active, joint_flow_target.as_ref(), 1.)?;
        let endload = if matches!(stage, Stage::ResetDisplay { .. }) {
            // Display RESET owns no physical load update. The last pressure
            // increment may differ from a fresh assembly by a few ulps; retain
            // its exact committed load, including in checkpoint reconstruction.
            origin.applied_nodal_load_kn.clone()
        } else {
            model.physical_load(
                &active,
                end_pressure,
                case.pressure_on_original_surface_only,
                end_factor,
                &prepared.element_material_ids,
                prepared.beams.as_ref(),
                model.stage_water(&origin, stage, 1.),
                end_joint_water.as_ref(),
            )?
        };
        let mut release = if matches!(stage, Stage::Excavate { .. }) && case.is_projected_mcc() {
            // Projection is domain-dependent: full minus removed work forces
            // need not equal the restricted-domain work force. Reconstruct
            // the retained force without changing any GP material history.
            let dofs = model.dofs(&active, &case.settings, end_factor)?;
            model
                .evaluate(
                    &origin,
                    &origin.displacement_m,
                    &active,
                    &dofs,
                    &case.settings,
                    EvaluationMode::plain(),
                )?
                .internal
                .iter()
                .zip(&endload)
                .map(|(i, f)| i - f)
                .collect::<Vec<_>>()
        } else if matches!(stage, Stage::Excavate { .. }) {
            let removed: Vec<_> = origin
                .active_elements
                .iter()
                .zip(&active)
                .map(|(a, b)| *a && !*b)
                .collect();
            let dofs = model.dofs(&active, &case.settings, end_factor)?;
            let removed_force = model
                .evaluate(
                    &origin,
                    &origin.displacement_m,
                    &removed,
                    &dofs,
                    &case.settings,
                    EvaluationMode {
                        structures: false,
                        ..EvaluationMode::plain()
                    },
                )?
                .internal;
            origin
                .applied_nodal_load_kn
                .iter()
                .zip(&endload)
                .zip(removed_force)
                .map(|((a, b), f)| a - b - f)
                .collect::<Vec<_>>()
        } else if stage.is_transition() {
            let dofs = model.dofs(&active, &case.settings, end_factor)?;
            let retained = model.evaluate(
                &prepared,
                &prepared.displacement_m,
                &active,
                &dofs,
                &case.settings,
                EvaluationMode::plain(),
            )?;
            retained
                .internal
                .iter()
                .zip(&endload)
                .map(|(a, b)| a - b)
                .collect()
        } else {
            vec![0.; endload.len()]
        };
        if matches!(stage, Stage::Excavate { .. })
            && let Some(joints) = &model.joints
        {
            // Split interfaces have different soil DOFs on their two sides.
            // The disappearing contact force is not part of removed-soil
            // internal force and must be removed explicitly from the old
            // equilibrium identity before constructing the holding load.
            let history = origin
                .joints
                .as_ref()
                .ok_or_else(|| Error::Input("missing JOINT history".into()))?;
            let before = joints.evaluate(
                history,
                &origin.displacement_m,
                &origin.active_elements,
                None,
                false,
            )?;
            let retained =
                joints.evaluate(history, &origin.displacement_m, &active, None, false)?;
            for ((h, a), b) in release
                .iter_mut()
                .zip(before.internal_kn)
                .zip(retained.internal_kn)
            {
                *h -= a - b;
            }
        }
        let active_dofs = model.dofs(&active, &case.settings, end_factor)?;
        for (i, r) in release.iter_mut().enumerate() {
            if !active_dofs.active_nodes[i / 3] {
                *r = 0.;
            }
        }
        if matches!(stage, Stage::Excavate { .. }) {
            let retained = model.evaluate(
                &origin,
                &origin.displacement_m,
                &active,
                &active_dofs,
                &case.settings,
                EvaluationMode::plain(),
            )?;
            let residual: Vec<_> = endload
                .iter()
                .zip(&release)
                .zip(&retained.internal)
                .map(|((f, h), r)| f + h - r)
                .collect();
            let selected: Vec<_> = (0..active.len())
                .filter(|&i| origin.active_elements[i] && !active[i])
                .collect();
            let audit = ExcavationAudit {
                stage: stage.name().into(),
                removed_element_count: selected.len(),
                removed_volume_m3: selected.iter().map(|&i| model.elements[i].volume).sum(),
                initial_holding_load_kn: release.clone(),
                zero_release_relative_force_residual: model.reduce(&residual, &active_dofs).norm()
                    / norm(&endload).max(1.),
            };
            if resumed.is_some() {
                let previous = excavation_audits
                    .last()
                    .ok_or_else(|| Error::Input("checkpoint missing excavation audit".into()))?;
                if serde_json::to_vec(previous)? != serde_json::to_vec(&audit)? {
                    return Err(Error::Input(
                        "checkpoint excavation holding-force audit mismatch".into(),
                    ));
                }
            } else {
                excavation_audits.push(audit);
            }
        }
        if stage.is_transition() {
            let retained = model.evaluate(
                &prepared,
                &prepared.displacement_m,
                &active,
                &active_dofs,
                &case.settings,
                EvaluationMode::plain(),
            )?;
            let residual: Vec<_> = endload
                .iter()
                .zip(&release)
                .zip(&retained.internal)
                .map(|((f, h), i)| f + h - i)
                .collect();
            let audit = TransitionAudit {
                stage: stage.name().into(),
                selected_elements: selected,
                prepared_state: prepared.clone(),
                initial_holding_load_kn: release.clone(),
                zero_release_relative_force_residual: model.reduce(&residual, &active_dofs).norm()
                    / norm(&endload).max(1.),
            };
            if let Some(old) = resumed.as_ref().and_then(|_| transition_audits.last()) {
                if serde_json::to_vec(old)? != serde_json::to_vec(&audit)? {
                    return Err(Error::Input(
                        "checkpoint construction transition mismatch".into(),
                    ));
                }
            } else if resumed.is_some() {
                return Err(Error::Input("checkpoint missing transition audit".into()));
            } else {
                transition_audits.push(audit);
            }
        }
        let mut progress = resumed.as_ref().map_or(0., |c| c.progress);
        let mut increment = resumed
            .as_ref()
            .map_or_else(|| stage.step(&case.settings), |c| c.increment);
        let mut count = resumed.as_ref().map_or(0, |c| c.stage_attempts);
        if let Some(runtime) = &resumed {
            if runtime.active != active
                || runtime.end_load != endload
                || runtime.holding_load != release
            {
                return Err(Error::Input(
                    "checkpoint geometry/load reconstruction mismatch".into(),
                ));
            }
            let expected_load: Vec<_> = if matches!(stage, Stage::ResetDisplay { .. }) {
                origin.applied_nodal_load_kn.clone()
            } else if matches!(
                stage,
                Stage::HydrostaticWater { .. } | Stage::SteadySeepage { .. }
            ) {
                model.physical_load(
                    &active,
                    origin.surface_pressure_kpa,
                    case.pressure_on_original_surface_only,
                    origin.boundary_factor,
                    &prepared.element_material_ids,
                    prepared.beams.as_ref(),
                    model.stage_water(&origin, stage, progress),
                    model
                        .stage_joint_water(
                            &origin,
                            stage,
                            &active,
                            joint_flow_target.as_ref(),
                            progress,
                        )?
                        .as_ref(),
                )?
            } else if !matches!(stage, Stage::Excavate { .. }) && !stage.is_transition() {
                origin
                    .applied_nodal_load_kn
                    .iter()
                    .zip(&endload)
                    .map(|(a, b)| a + (b - a) * progress)
                    .collect()
            } else {
                endload
                    .iter()
                    .zip(&release)
                    .map(|(a, h)| a + (1. - progress) * h)
                    .collect()
            };
            let expected_release: Vec<_> = release.iter().map(|h| (1. - progress) * h).collect();
            let expected_pressure =
                if !matches!(stage, Stage::Excavate { .. }) && !stage.is_transition() {
                    origin.surface_pressure_kpa
                        + (end_pressure - origin.surface_pressure_kpa) * progress
                } else {
                    end_pressure
                };
            let expected_anchor_basis =
                if let Some(control) = model.anchor_control(&origin, stage, progress)? {
                    let m = model.anchors.as_ref().unwrap();
                    Some(
                        m.evaluate(
                            prepared.anchors.as_ref().unwrap(),
                            &committed.displacement_m,
                            &active,
                            Some(&control),
                        )?
                        .state,
                    )
                } else {
                    prepared.anchors.clone()
                };
            if committed.active_elements != active
                || committed.pore_pressure_kpa
                    != model.stage_pore(
                        &origin,
                        stage,
                        &active,
                        flow_target.as_deref(),
                        progress,
                    )?
                || committed.hydrostatic_water_level_m
                    != model.stage_water(&origin, stage, progress)
                || committed.joint_water
                    != model.stage_joint_water(
                        &origin,
                        stage,
                        &active,
                        joint_flow_target.as_ref(),
                        progress,
                    )?
                || !same_anchor_basis(&committed.anchors, &expected_anchor_basis)
                || !same_beam_installation(&committed.beams, &prepared.beams)
                || committed.display_reference_displacement_m
                    != if matches!(stage, Stage::ResetDisplay { .. }) {
                        origin.displacement_m.clone()
                    } else {
                        prepared.display_reference_displacement_m.clone()
                    }
                || committed.surface_pressure_kpa != expected_pressure
                || committed.boundary_factor
                    != origin.boundary_factor + (end_factor - origin.boundary_factor) * progress
                || committed.applied_nodal_load_kn != expected_load
                || (stage.is_transition()
                    && (committed.element_material_ids != prepared.element_material_ids
                        || committed.transition_release_load_kn != expected_release
                        || committed
                            .excavation_release_load_kn
                            .iter()
                            .any(|x| *x != 0.)))
                || (!stage.is_transition()
                    && (committed.excavation_release_load_kn != expected_release
                        || (if matches!(stage, Stage::ResetDisplay { .. }) {
                            committed.transition_release_load_kn
                                != origin.transition_release_load_kn
                        } else {
                            !committed.transition_release_load_kn.is_empty()
                        })))
            {
                return Err(Error::Input(
                    "checkpoint current active set/load path mismatch".into(),
                ));
            }
        }
        while progress < 1. - 1e-12 {
            count += 1;
            if count > 2000 {
                failure = Some("stage attempt budget exceeded".into());
                break;
            }
            let next = if matches!(stage, Stage::ResetDisplay { .. }) {
                1.
            } else {
                (progress + increment).min(1.)
            };
            let target: Vec<_> = if matches!(stage, Stage::ResetDisplay { .. }) {
                origin.applied_nodal_load_kn.clone()
            } else if matches!(
                stage,
                Stage::HydrostaticWater { .. } | Stage::SteadySeepage { .. }
            ) {
                model.physical_load(
                    &active,
                    origin.surface_pressure_kpa,
                    case.pressure_on_original_surface_only,
                    origin.boundary_factor,
                    &prepared.element_material_ids,
                    prepared.beams.as_ref(),
                    model.stage_water(&origin, stage, next),
                    model
                        .stage_joint_water(
                            &origin,
                            stage,
                            &active,
                            joint_flow_target.as_ref(),
                            next,
                        )?
                        .as_ref(),
                )?
            } else if !matches!(stage, Stage::Excavate { .. }) && !stage.is_transition() {
                origin
                    .applied_nodal_load_kn
                    .iter()
                    .zip(&endload)
                    .map(|(a, b)| a + (b - a) * next)
                    .collect()
            } else {
                endload
                    .iter()
                    .zip(&release)
                    .map(|(f, r)| f + (1. - next) * r)
                    .collect()
            };
            let releasing = release.iter().map(|r| (1. - next) * r).collect::<Vec<_>>();
            let pressure = if !matches!(stage, Stage::Excavate { .. }) && !stage.is_transition() {
                origin.surface_pressure_kpa + (end_pressure - origin.surface_pressure_kpa) * next
            } else {
                end_pressure
            };
            let before = serde_json::to_vec(&committed)?;
            model.notify(|| {
                model.progress(
                    crate::stage_control::Kind::AttemptStarted,
                    stage_index + 1,
                    next,
                )
            })?;
            let trial = if matches!(stage, Stage::ResetDisplay { .. }) {
                if target != origin.applied_nodal_load_kn {
                    return Err(Error::Input(
                        "display RESET must not change physical loads".into(),
                    ));
                }
                model.reset_display(&origin, stage_index + 1)
            } else {
                model.trial(
                    if stage.is_transition() && progress == 0. {
                        &prepared
                    } else {
                        &committed
                    },
                    Target {
                        joint_water: model.stage_joint_water(
                            &origin,
                            stage,
                            &active,
                            joint_flow_target.as_ref(),
                            next,
                        )?,
                        pore_pressure_kpa: model.stage_pore(
                            &origin,
                            stage,
                            &active,
                            flow_target.as_deref(),
                            next,
                        )?,
                        water_level: model.stage_water(&origin, stage, next),
                        anchor_control: model.anchor_control(&origin, stage, next)?,
                        active: active.clone(),
                        load: target,
                        release: releasing,
                        pressure,
                        boundary_factor: origin.boundary_factor
                            + (end_factor - origin.boundary_factor) * next,
                        stage: stage_index + 1,
                        progress: next,
                    },
                    &case.settings,
                )
            };
            let trial = match trial {
                Err(error @ Error::Control(_)) => return Err(error),
                other => other,
            };
            let changed_structure = trial.as_ref().is_ok_and(|t| {
                t.state.beams != committed.beams
                    || t.state.anchors != committed.anchors
                    || t.state.joints != committed.joints
            });
            let changed_water = trial.as_ref().is_ok_and(|t| {
                t.state.hydrostatic_water_level_m != committed.hydrostatic_water_level_m
                    || t.state.pore_pressure_kpa != committed.pore_pressure_kpa
                    || t.state.joint_water != committed.joint_water
            });
            let changed_joint_slip = trial
                .as_ref()
                .is_ok_and(|t| joint_slip_changed(&t.state.joints, &committed.joints));
            let force_reject = case
                .verification_reject_once_stages
                .iter()
                .any(|n| n == stage.name())
                && !injected.contains(stage.name())
                && trial.is_ok();
            let (
                iterations,
                residual,
                changed_material,
                changed_plastic,
                changed_active,
                changed_load,
            ) = if let Ok(t) = &trial {
                (
                    Some(t.iterations),
                    Some(t.residual),
                    t.state.element_material_ids != committed.element_material_ids
                        || serde_json::to_vec(&t.state.material_states)?
                            != serde_json::to_vec(&committed.material_states)?,
                    t.state
                        .material_states
                        .iter()
                        .zip(&committed.material_states)
                        .any(|(a, b)| a.plastic_strain() != b.plastic_strain()),
                    t.state.active_elements != committed.active_elements,
                    t.state.applied_nodal_load_kn != committed.applied_nodal_load_kn,
                )
            } else {
                (None, None, false, false, false, false)
            };
            let from = progress;
            let linear_solves = trial
                .as_ref()
                .map_or_else(|_| Vec::new(), |t| t.linear_solves.clone());
            let contact_predictors = trial
                .as_ref()
                .map_or_else(|_| vec![], |t| t.contact_predictors.clone());
            match trial {
                Ok(t) if !force_reject => {
                    committed = t.state;
                    progress = next;
                    attempts.push(Attempt {
                        contact_predictors,
                        stage: stage.name().into(),
                        progress_from: from,
                        progress_to: next,
                        accepted: true,
                        rejected_trial_water_changed: false,
                        reason: "converged_material_force".into(),
                        newton_iterations: iterations,
                        relative_force_residual: residual,
                        committed_state_unchanged_on_rejection: false,
                        rejected_trial_joint_slip_changed: false,
                        rejected_trial_structure_changed: false,
                        rejected_trial_material_changed: false,
                        rejected_trial_plastic_changed: false,
                        rejected_trial_active_set_changed: false,
                        rejected_trial_load_changed: false,
                        linear_solves,
                    });
                    accepted_this_run += 1;
                    model.notify_attempt(stage_index + 1, attempts.last().unwrap())?;
                    if stop_after.is_some_and(|limit| accepted_this_run >= limit)
                        || model.stop_requested()?
                    {
                        continuation = Some(Continuation {
                            origin,
                            active,
                            end_load: endload,
                            holding_load: release,
                            stage_index,
                            progress,
                            increment,
                            stage_attempts: count,
                            injected_stages: injected,
                        });
                        failure = Some("checkpoint_requested".into());
                        break 'stages;
                    }
                }
                result => {
                    let reason = if force_reject {
                        injected.insert(stage.name().to_string());
                        "verification_injected_rejection_after_convergence".into()
                    } else {
                        match result {
                            Err(e) => e.to_string(),
                            Ok(_) => "rejected".into(),
                        }
                    };
                    let unchanged = before == serde_json::to_vec(&committed)?;
                    attempts.push(Attempt {
                        contact_predictors,
                        stage: stage.name().into(),
                        progress_from: progress,
                        progress_to: next,
                        accepted: false,
                        rejected_trial_water_changed: changed_water,
                        reason: reason.clone(),
                        newton_iterations: iterations,
                        relative_force_residual: residual,
                        committed_state_unchanged_on_rejection: unchanged,
                        rejected_trial_joint_slip_changed: changed_joint_slip,
                        rejected_trial_structure_changed: changed_structure,
                        rejected_trial_material_changed: changed_material,
                        rejected_trial_plastic_changed: changed_plastic,
                        rejected_trial_active_set_changed: changed_active,
                        rejected_trial_load_changed: changed_load,
                        linear_solves,
                    });
                    if !unchanged {
                        return Err(Error::Solve(
                            "committed history changed during rejected trial".into(),
                        ));
                    }
                    model.notify_attempt(stage_index + 1, attempts.last().unwrap())?;
                    increment *= 0.5;
                    if increment < case.settings.minimum_increment {
                        failure = Some(format!("{}: cutback exhausted: {reason}", stage.name()));
                        break;
                    }
                }
            }
        }
        if failure.is_some() {
            break;
        }
        frames.push(model.frame(
            stage.name().into(),
            &committed,
            if matches!(stage, Stage::Birth { .. }) {
                &prepared
            } else {
                &origin
            },
            &case,
        )?);
    }
    let joint_integration_points = model.joint_layout();
    let integration_points = model.layout();
    let reference_displacement_m = model.reference_displacement_m.clone();
    Ok(Run {
        flow_audits,
        joint_integration_points,
        provenance,
        schema: result_schema(mesh, &ground, &case).into(),
        complete: failure.is_none(),
        failure_reason: failure,
        case,
        ground,
        frames,
        attempts,
        excavation_audits,
        transition_audits,
        last_committed: committed,
        continuation,
        integration_points,
        geostatic: geostatic_report,
        reference_displacement_m,
    })
}

fn validate_committed(model: &Model<'_>, state: &Committed, case: &StageCase) -> Result<()> {
    if case.is_mixed_mcc() {
        mixed::validate_state(model, state)?;
    } else if state.mixed_mcc.is_some() {
        return Err(Error::Input(
            "mixed state supplied to displacement formulation".into(),
        ));
    }
    match (&model.wet_joints, &state.joint_water) {
        (None, None) => {}
        (Some(w), Some(field)) => {
            if w.retain(field, &state.active_elements)? != *field {
                return Err(Error::Input(
                    "checkpoint JOINT water retains an inactive pair load".into(),
                ));
            }
            if state.pore_pressure_kpa.is_none()
                && w.hydrostatic(
                    &state.active_elements,
                    state
                        .hydrostatic_water_level_m
                        .or(model.ground.initialization.water_level_m),
                )? != *field
            {
                return Err(Error::Input(
                    "checkpoint hydrostatic JOINT field differs from its water level".into(),
                ));
            }
        }
        _ => {
            return Err(Error::Input(
                "checkpoint hydraulic JOINT policy/state mismatch".into(),
            ));
        }
    }
    let nodes = model.mesh.nodes.len();
    let elements = model.elements.len();
    if state.pore_pressure_kpa.as_ref().is_some_and(|p| {
        !case.has_seepage()
            || p.len() != model.point_count
            || !p.iter().all(|v| v.is_finite())
            || p.iter()
                .zip(&model.soil_points)
                .any(|(v, soil)| !soil && *v != 0.)
    }) {
        return Err(Error::Input(
            "invalid checkpoint ONE_WAY pore pressures".into(),
        ));
    }
    if state
        .hydrostatic_water_level_m
        .is_some_and(|h| !h.is_finite())
        || state.hydrostatic_water_level_m.is_some() != case.has_staged_water()
        || state.displacement_m.len() != nodes
        || (!state.display_reference_displacement_m.is_empty()
            && state.display_reference_displacement_m.len() != nodes)
        || (!state.element_material_ids.is_empty() && state.element_material_ids.len() != elements)
        || (case.has_lifecycle() && state.element_material_ids.is_empty())
        || state.material_states.len() != model.point_count
        || state.active_elements.len() != elements
        || state.applied_nodal_load_kn.len() != 3 * nodes
        || state.excavation_release_load_kn.len() != 3 * nodes
        || (!state.geostatic_release_load_kn.is_empty()
            && state.geostatic_release_load_kn.len() != 3 * nodes)
        || (!state.transition_release_load_kn.is_empty()
            && state.transition_release_load_kn.len() != 3 * nodes)
        || !state
            .displacement_m
            .iter()
            .flatten()
            .chain(state.display_reference_displacement_m.iter().flatten())
            .chain(&state.applied_nodal_load_kn)
            .chain(&state.excavation_release_load_kn)
            .chain(&state.geostatic_release_load_kn)
            .chain(&state.transition_release_load_kn)
            .all(|v| v.is_finite())
    {
        return Err(Error::Input(
            "invalid checkpoint state dimensions/nonfinite data".into(),
        ));
    }
    for (ei, e) in model.elements.iter().enumerate() {
        for qi in 0..e.points.len() {
            model.ground.materials[model.material_index(&state.element_material_ids, ei)?]
                .model
                .validate_state(&state.material_states[e.point_start + qi])?;
        }
    }
    let dofs = model.dofs(
        &state.active_elements,
        &case.settings,
        state.boundary_factor,
    )?;
    match (&model.joints, &state.joints) {
        (None, None) => {}
        (Some(j), Some(s)) => j.validate_state(s, &state.displacement_m, &state.active_elements)?,
        _ => {
            return Err(Error::Input(
                "checkpoint JOINT configuration/history mismatch".into(),
            ));
        }
    }
    match (&model.anchors, &state.anchors) {
        (None, None) => {}
        (Some(model), Some(anchors)) => {
            model.validate_state(anchors)?;
            let evaluated =
                model.evaluate(anchors, &state.displacement_m, &state.active_elements, None)?;
            for (old, new) in anchors.anchors.iter().zip(&evaluated.state.anchors) {
                if old.segment_taut != new.segment_taut
                    || old
                        .segment_forces_kn
                        .iter()
                        .zip(&new.segment_forces_kn)
                        .any(|(a, b)| (a - b).abs() > 1e-9 * (1. + a.abs()))
                {
                    return Err(Error::Input(
                        "checkpoint anchor force/taut cache is inconsistent".into(),
                    ));
                }
            }
        }
        _ => {
            return Err(Error::Input(
                "checkpoint anchor configuration/state mismatch".into(),
            ));
        }
    }
    match (&model.beams, &state.beams) {
        (None, None) => {}
        (Some(model), Some(beams)) => {
            model.validate_state(beams)?;
            let evaluated = model.evaluate(
                beams,
                &state.displacement_m,
                &dofs.active_nodes,
                state.boundary_factor,
            )?;
            if beams
                .node_displacement_m
                .iter()
                .flatten()
                .zip(evaluated.state.node_displacement_m.iter().flatten())
                .chain(
                    beams
                        .node_rotation_rad
                        .iter()
                        .flatten()
                        .zip(evaluated.state.node_rotation_rad.iter().flatten()),
                )
                .any(|(a, b)| (a - b).abs() > 1e-10 * (1. + a.abs()))
            {
                return Err(Error::Input(
                    "checkpoint beam interior kinematics are inconsistent".into(),
                ));
            }
        }
        _ => {
            return Err(Error::Input(
                "checkpoint beam configuration/state mismatch".into(),
            ));
        }
    }
    let mut tied = BTreeMap::new();
    if let Some(t) = &dofs.transform {
        if t.maximum_constraint_error(
            &state
                .displacement_m
                .iter()
                .flatten()
                .copied()
                .collect::<Vec<_>>(),
        ) > 1e-10
        {
            return Err(Error::Input(
                "checkpoint explicit displacement mismatch".into(),
            ));
        }
    } else {
        for (g, &mapped) in dofs.map.iter().enumerate() {
            let value = state.displacement_m[g / 3][g % 3];
            if dofs.active_nodes[g / 3] {
                if let Some(r) = mapped {
                    if let Some(old) = tied.insert(r, value)
                        && (old - value).abs() > 1e-12
                    {
                        return Err(Error::Input("checkpoint tied displacement mismatch".into()));
                    }
                } else if value.abs() > 1e-12 {
                    return Err(Error::Input(
                        "checkpoint support displacement mismatch".into(),
                    ));
                }
            }
        }
    }
    let frame = model.frame("checkpoint_validation".into(), state, state, case)?;
    if frame.relative_force_residual > case.settings.relative_force_tolerance
        || frame.maximum_constraint_error_m > 1e-10
    {
        return Err(Error::Input(
            "checkpoint is not an equilibrated committed state".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod predictor_tests {
    use super::*;
    fn material_fixture() -> Result<(Vec<soil::Model>, Vec<soil::State>, Vec<MaterialWork>)> {
        let mut parameters = crate::mcc::Parameters {
            critical_state_m: 1.2,
            lambda: 0.2,
            kappa: 0.04,
            poisson: 0.3,
            reference_pressure_kpa: 1.,
            specific_volume_ncl_at_reference: 2.7,
            integration: crate::mcc::Integration::Midpoint,
        };
        let midpoint = parameters.clone();
        parameters.integration = crate::mcc::Integration::BackwardEuler;
        let models = vec![
            soil::Model::Mcc(midpoint),
            soil::Model::Mcc(parameters.clone()),
        ];
        let state = crate::mcc::State {
            effective_stress_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
            void_ratio: parameters.reference_void_ratio(100., 120.)?,
            preconsolidation_pressure_kpa: 120.,
            plastic_strain_compression_positive: [0.; 6],
        };
        let states = (0..128).map(|_| soil::State::Mcc(state.clone())).collect();
        let work = (0..128)
            .map(|i| MaterialWork {
                point: i,
                material: i % 2,
                partition: Some(32),
                strain: match i % 4 {
                    0 => [1e-5, -2e-6, 3e-6, 4e-6, 0., 0.],
                    1 => [0.003, 0.003, 0.003, 0., 0., 0.],
                    2 => [0.0001, -0.00005, -0.00005, 0.001, 0., 0.],
                    _ => [0.00003, 0., 0., 0., 0., 0.],
                },
            })
            .collect();
        Ok((models, states, work))
    }

    #[test]
    fn parallel_material_states_tangents_partitions_and_order_match_serial_without_history_mutation()
    -> Result<()> {
        let (models, states, work) = material_fixture()?;
        let before = serde_json::to_value(&states)?;
        for tangent in [false, true] {
            let serial =
                material_work_batch(&models, None, &states, &work, (tangent, 1e-7, 5e-5), 1)?;
            let parallel =
                material_work_batch(&models, None, &states, &work, (tangent, 1e-7, 5e-5), 4)?;
            assert_eq!(serial.len(), 128);
            assert_eq!(parallel.len(), serial.len());
            for ((a, ta), (b, tb)) in serial.iter().zip(&parallel) {
                assert_eq!(
                    serde_json::to_value(&a.state)?,
                    serde_json::to_value(&b.state)?
                );
                assert_eq!(
                    (a.substeps, a.plastic_substeps),
                    (b.substeps, b.plastic_substeps)
                );
                assert_eq!(ta, tb);
            }
            assert_eq!(serde_json::to_value(&states)?, before);
        }
        Ok(())
    }

    #[test]
    fn parallel_material_errors_choose_first_input_failure_and_leave_all_committed_histories_unchanged()
    -> Result<()> {
        let (models, states, mut work) = material_fixture()?;
        let before = serde_json::to_value(&states)?;
        work[10].partition = Some(0);
        work[100].material = usize::MAX;
        let message = |result: Result<Vec<MaterialResponse>>| match result {
            Err(error) => error.to_string(),
            Ok(_) => panic!("invalid material work accepted"),
        };
        let serial = message(material_work_batch(
            &models,
            None,
            &states,
            &work,
            (true, 1e-7, 5e-5),
            1,
        ));
        let parallel = message(material_work_batch(
            &models,
            None,
            &states,
            &work,
            (true, 1e-7, 5e-5),
            4,
        ));
        assert_eq!(serial, parallel);
        assert!(serial.contains("partition"));
        assert_eq!(serde_json::to_value(&states)?, before);
        assert!(material_work_batch(&models, None, &states, &work, (true, 1e-7, 5e-5), 0).is_err());
        assert!(material_work_batch(&models, None, &states, &work, (true, 1e-7, 5e-5), 5).is_err());
        Ok(())
    }

    #[test]
    fn temporary_sparse_and_dense_gauges_do_not_certify_the_original_equations() -> Result<()> {
        for sparse in [false, true] {
            let mut settings: Settings = serde_json::from_value(
                serde_json::from_str::<serde_json::Value>(include_str!(
                    "../examples/stages-free-contact.json"
                ))?["settings"]
                    .clone(),
            )?;
            if !sparse {
                settings.linear_solver = LinearSolver::DenseLu;
            }
            let mut a = if sparse {
                Tangent::Sparse(CsrMatrix::zeros(CsrPattern::from_connectivity(
                    2,
                    [vec![0, 1]],
                )?))
            } else {
                Tangent::Dense(DMatrix::zeros(2, 2))
            };
            a.add(0, 1, -3.)?;
            a.add(1, 1, 10.)?;
            let gauges = crate::rigidity::gauge_coordinates(
                &DMatrix::identity(2, 2),
                |v| a.action(v),
                a.scale(),
            )?;
            assert_eq!(gauges.coordinates, vec![0]);
            let rhs = DVector::from_vec(vec![0., 10.]);
            let (x, report) = a.gauged(&gauges.coordinates)?.solve(&rhs, &settings)?;
            assert!(report.true_relative_residual < 1e-12);
            assert!(((a.action(&x)? - &rhs).norm() / rhs.norm() - 0.3).abs() < 1e-12);
            assert!(a.action(&DVector::from_vec(vec![1., 0.]))?.norm() == 0.);
            if let Tangent::Sparse(original) = &a {
                assert_eq!(original.pattern().nonzeros(), 4);
                assert_eq!(original.values()[0], 0.);
            }
        }
        Ok(())
    }
}
