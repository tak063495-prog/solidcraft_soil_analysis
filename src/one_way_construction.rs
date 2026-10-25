//! Physical-time ONE_WAY with instantaneous CAD DEATH and separate mechanical
//! stress-release continuation. No water time is consumed by release progress.
use crate::{
    Error, Result, boundary, checkpoint, elements,
    ground::Ground,
    hydraulic_transition as water,
    mesh::Mesh,
    one_way::{self, Fluid, MechanicalAudit, State, Transfer},
    richards, soil, stages,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CutWaterPolicy {
    ExplicitPeriodBoundaries,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CutTractionPolicy {
    ExplicitPeriodReservoirs,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Excavation {
    pub name: String,
    pub before_period: String,
    pub bodies: Vec<String>,
    pub water_boundary_policy: CutWaterPolicy,
    pub cut_traction_policy: CutTractionPolicy,
    #[serde(default)]
    pub initial_increment: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub schema: String,
    pub one_way: one_way::Case,
    pub excavations: Vec<Excavation>,
    #[serde(default)]
    pub verification_reject_once_excavations: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub initially_inactive_bodies: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transitions: Vec<Transition>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub verification_reject_once_transitions: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Action {
    Birth {
        bodies: Vec<String>,
        states: BTreeMap<String, soil::State>,
        placement: stages::BirthPlacement,
        initial_head: richards::AffineHead,
    },
    MaterialSet {
        bodies: Vec<String>,
        material_id: u32,
        history: stages::MaterialHistory,
        water: water::WaterPolicy,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Transition {
    pub name: String,
    pub before_period: String,
    pub action: Action,
    pub water_boundary_policy: CutWaterPolicy,
    pub surface_traction_policy: CutTractionPolicy,
    #[serde(default)]
    pub initial_increment: Option<f64>,
}
impl Transition {
    fn stage(&self) -> stages::Stage {
        match &self.action {
            Action::Birth {
                bodies,
                states,
                placement,
                ..
            } => stages::Stage::Birth {
                name: self.name.clone(),
                bodies: bodies.clone(),
                states: states.clone(),
                placement: placement.clone(),
                initial_increment: self.initial_increment,
            },
            Action::MaterialSet {
                bodies,
                material_id,
                history,
                ..
            } => stages::Stage::MaterialSet {
                name: self.name.clone(),
                bodies: bodies.clone(),
                material_id: *material_id,
                history: history.clone(),
                initial_increment: self.initial_increment,
            },
        }
    }
    fn water(&self) -> water::WaterPolicy {
        match &self.action {
            Action::Birth { initial_head, .. } => water::WaterPolicy::InitializeHead {
                value: initial_head.clone(),
            },
            Action::MaterialSet { water, .. } => water.clone(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArchivedPoint {
    pub integration_point_index: usize,
    pub material_id: u32,
    pub state: soil::State,
    pub active_before: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TransitionRecord {
    pub name: String,
    pub time_s: f64,
    pub origin_pair_sha256: String,
    pub prepared_mechanical: stages::Committed,
    pub old_material_points: Vec<ArchivedPoint>,
    pub water: water::WaterChange,
    pub initial_holding_load_kn: Vec<f64>,
    pub zero_release_relative_force_residual: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PendingTransition {
    pub event_index: usize,
    pub progress: f64,
    pub increment: f64,
    pub record: TransitionRecord,
    pub fluid: Fluid,
    pub porosity_reference: Vec<f64>,
    pub archived: bool,
}
fn zero(v: &f64) -> bool {
    *v == 0.
}
fn zero_usize(v: &usize) -> bool {
    *v == 0
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExcavationRecord {
    pub name: String,
    pub time_s: f64,
    pub origin_pair_sha256: String,
    pub origin_hydraulic: richards::State,
    pub removal: richards::DomainRemoval,
    pub initial_holding_load_kn: Vec<f64>,
    pub zero_release_relative_force_residual: f64,
    pub removed_mechanical_points: Vec<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub removed_material_points: Vec<ArchivedPoint>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    pub name: String,
    pub phase: String,
    pub progress: f64,
    pub state: State,
    pub hydraulic_faces: Vec<richards::Face>,
    pub hydraulic_cell_element_indices: Vec<usize>,
    /// None for a boundary jump at fixed time: it is not a flow equilibrium.
    pub hydraulic_balance: Option<richards::Balance>,
    pub mechanical_audit: MechanicalAudit,
    pub holding_load_kn: Vec<f64>,
    pub cumulative_removed_water_m3: f64,
    pub cumulative_net_water_inflow_m3: f64,
    pub signed_global_water_budget_error_m3: f64,
    #[serde(default, skip_serializing_if = "zero")]
    pub cumulative_added_water_m3: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub hydraulic_material_ids: Vec<u32>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Attempt {
    pub name: String,
    pub phase: String,
    pub time_s: f64,
    pub dt_s: Option<f64>,
    pub progress_from: f64,
    pub progress_to: f64,
    pub accepted: bool,
    pub reason: Option<String>,
    pub candidate_changed: bool,
    pub plastic_candidate_changed: bool,
    pub committed_pair_sha256_before: String,
    pub committed_pair_sha256_after: String,
    pub hydraulic_iterations: Vec<richards::Iteration>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pending {
    pub event_index: usize,
    pub progress: f64,
    pub increment: f64,
    pub active: Vec<bool>,
    pub fluid: Fluid,
    pub record: ExcavationRecord,
    pub archived: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Runtime {
    pub period_index: usize,
    pub event_index: usize,
    pub pending: Option<Pending>,
    pub next_dt_s: f64,
    pub injected_excavations: BTreeSet<String>,
    pub injected_times: BTreeSet<usize>,
    pub removed_water_m3: f64,
    pub net_water_inflow_m3: f64,
    #[serde(default, skip_serializing_if = "zero_usize")]
    pub transition_index: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pending_transition: Option<PendingTransition>,
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub injected_transitions: BTreeSet<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub porosity_reference: Vec<f64>,
    #[serde(default, skip_serializing_if = "zero")]
    pub added_water_m3: f64,
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
    pub initialization: one_way::Report,
    pub initial_water_m3: f64,
    pub frames: Vec<Frame>,
    pub attempts: Vec<Attempt>,
    pub excavations: Vec<ExcavationRecord>,
    pub last_committed: State,
    pub runtime: Option<Runtime>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transitions: Vec<TransitionRecord>,
}

/// Existing named sets keep their original CAD membership, intersected with
/// the currently exposed faces. No coordinate or centroid selection is used.
pub fn domain_mesh(mesh: &Mesh, active: &[bool]) -> Result<Mesh> {
    if active.len() != mesh.tetrahedra.len() {
        return Err(Error::Input("invalid construction mask".into()));
    }
    let owners = boundary::surface_owners(mesh);
    let exposed = owners
        .iter()
        .filter(|(_, o)| o.iter().filter(|&&e| active[e]).count() == 1)
        .map(|(k, _)| *k)
        .collect::<BTreeSet<_>>();
    let mut result = mesh.clone();
    for set in result.face_sets.values_mut() {
        set.retain(|f| {
            let mut key = *f;
            key.sort_unstable();
            exposed.contains(&key)
        });
    }
    result.face_sets.retain(|_, set| !set.is_empty());
    result.validate()?;
    Ok(result)
}

fn periods(case: &Case) -> Result<&[richards::Period]> {
    one_way::validate(&case.one_way)?;
    if case.schema != "soil-fem-one-way-construction/1"
        || (case.excavations.is_empty()
            && case.transitions.is_empty()
            && case.initially_inactive_bodies.is_empty())
    {
        return Err(Error::Input(
            "construction requires a nonempty CAD excavation schedule".into(),
        ));
    }
    let richards::Analysis::Transient { periods } = &case.one_way.flow.analysis else {
        return Err(Error::Input("construction requires transient flow".into()));
    };
    if periods.is_empty()
        || case.one_way.exterior_periods.len() != periods.len()
        || case
            .one_way
            .exterior_periods
            .iter()
            .zip(periods)
            .any(|(a, b)| a.name != b.name)
    {
        return Err(Error::Input(
            "construction exterior/flow history mismatch".into(),
        ));
    }
    let map = periods
        .iter()
        .enumerate()
        .map(|(i, p)| (p.name.as_str(), i))
        .collect::<BTreeMap<_, _>>();
    let mut previous = None;
    let mut names = BTreeSet::new();
    for event in &case.excavations {
        if event.initial_increment.is_some_and(|v| {
            !v.is_finite() || v < case.one_way.mechanics.settings.minimum_increment || v > 1.
        }) {
            return Err(Error::Input("invalid excavation release increment".into()));
        }
        let i = *map
            .get(event.before_period.as_str())
            .ok_or_else(|| Error::Input("unknown excavation period".into()))?;
        if event.name.is_empty()
            || !names.insert(event.name.clone())
            || previous.is_some_and(|j| i <= j)
            || event.bodies.is_empty()
            || event.bodies.iter().any(String::is_empty)
            || event.bodies.iter().collect::<BTreeSet<_>>().len() != event.bodies.len()
        {
            return Err(Error::Input(
                "invalid or unordered excavation schedule".into(),
            ));
        }
        previous = Some(i);
    }
    let mut occupied = case
        .excavations
        .iter()
        .map(|e| e.before_period.clone())
        .collect::<BTreeSet<_>>();
    let mut transition_names = BTreeSet::new();
    let mut prior = None;
    for event in &case.transitions {
        let i = *map
            .get(event.before_period.as_str())
            .ok_or_else(|| Error::Input("unknown transition period".into()))?;
        if event.name.is_empty()
            || !names.insert(event.name.clone())
            || !transition_names.insert(event.name.clone())
            || !occupied.insert(event.before_period.clone())
            || prior.is_some_and(|p| i <= p)
            || event.initial_increment.is_some_and(|v| {
                !v.is_finite() || v < case.one_way.mechanics.settings.minimum_increment || v > 1.
            })
        {
            return Err(Error::Input(
                "invalid, overlapping or unordered lifecycle schedule".into(),
            ));
        }
        prior = Some(i);
        match event.water() {
            water::WaterPolicy::InitializeHead { value } => value.validate()?,
            water::WaterPolicy::PreserveWater { reference_head } => reference_head.validate()?,
            _ => {}
        }
    }
    if case
        .verification_reject_once_transitions
        .iter()
        .any(|n| !transition_names.contains(n))
        || case
            .verification_reject_once_transitions
            .iter()
            .collect::<BTreeSet<_>>()
            .len()
            != case.verification_reject_once_transitions.len()
    {
        return Err(Error::Input("invalid lifecycle rejection targets".into()));
    }
    if case
        .verification_reject_once_excavations
        .iter()
        .any(|n| !names.contains(n))
        || case
            .verification_reject_once_excavations
            .iter()
            .collect::<BTreeSet<_>>()
            .len()
            != case.verification_reject_once_excavations.len()
    {
        return Err(Error::Input(
            "invalid construction rejection targets".into(),
        ));
    }
    Ok(periods)
}
fn material_ids(mesh: &Mesh, state: &stages::Committed) -> Vec<u32> {
    if state.element_material_ids.is_empty() {
        mesh.material_ids.clone()
    } else {
        state.element_material_ids.clone()
    }
}
fn initial_active(
    kernel: &stages::ExternalMechanics<'_>,
    case: &Case,
    n: usize,
) -> Result<Vec<bool>> {
    let mut active = vec![true; n];
    if !case.initially_inactive_bodies.is_empty() {
        for e in kernel.body_selection(&case.initially_inactive_bodies)? {
            active[e] = false;
        }
    }
    Ok(active)
}
fn new_exposed(mesh: &Mesh, before: &[bool], after: &[bool]) -> BTreeSet<[usize; 3]> {
    boundary::surface_owners(mesh)
        .into_iter()
        .filter(|(_, o)| {
            o.iter().filter(|&&e| after[e]).count() == 1 && o.iter().all(|&e| !before[e])
        })
        .map(|(k, _)| k)
        .collect()
}
fn prepare_pending(
    kernel: &stages::ExternalMechanics<'_>,
    mesh: &Mesh,
    ground: &Ground,
    case: &Case,
    old: &State,
    event: &Transition,
    porosity: &[f64],
) -> Result<PendingTransition> {
    let (prepared, selected) =
        kernel.prepare_external_transition(&old.mechanical, &event.stage())?;
    let old_ids = material_ids(mesh, &old.mechanical);
    let ids = material_ids(mesh, &prepared);
    let change = water::change(
        mesh,
        &case.one_way.flow,
        water::Domain {
            active: &old.mechanical.active_elements,
            material_ids: &old_ids,
        },
        water::Domain {
            active: &prepared.active_elements,
            material_ids: &ids,
        },
        &old.hydraulic,
        &selected,
        &event.water(),
    )?;
    let layout = kernel.layout();
    let chosen = selected.into_iter().collect::<BTreeSet<_>>();
    if porosity.len() != layout.len() {
        return Err(Error::Input("missing lifecycle porosity origins".into()));
    }
    let mut reference = porosity.to_vec();
    let mut archived = vec![];
    for (i, p) in layout.iter().enumerate() {
        if chosen.contains(&p.element_index) {
            let e = prepared.material_states[i].void_ratio();
            reference[i] = e / (1. + e);
            archived.push(ArchivedPoint {
                integration_point_index: i,
                material_id: old_ids[p.element_index],
                state: old.mechanical.material_states[i].clone(),
                active_before: old.mechanical.active_elements[p.element_index],
            });
        }
    }
    let projected = domain_mesh(mesh, &prepared.active_elements)?;
    let mut transfer = Transfer::new(&projected, ground, &case.one_way, &prepared, layout)?;
    transfer.set_porosity_reference(&reference)?;
    let reservoirs = &case
        .one_way
        .exterior_periods
        .iter()
        .find(|p| p.name == event.before_period)
        .ok_or_else(|| Error::Input("missing lifecycle reservoirs".into()))?
        .reservoirs;
    let fluid = transfer.field_with_base(
        &change.state,
        reservoirs,
        &kernel.dry_load(&prepared, &prepared.active_elements)?,
    )?;
    let (holding, residual) =
        kernel.transition_holding(&prepared, &fluid.equivalent_mechanical_rhs_kn)?;
    if residual > case.one_way.mechanics.settings.relative_force_tolerance * 10. {
        return Err(Error::Solve(
            "lifecycle zero-release force identity failed".into(),
        ));
    }
    Ok(PendingTransition {
        event_index: 0,
        progress: 0.,
        increment: event
            .initial_increment
            .unwrap_or(case.one_way.mechanics.settings.initial_increment),
        record: TransitionRecord {
            name: event.name.clone(),
            time_s: old.hydraulic.time_s,
            origin_pair_sha256: checkpoint::hash_json(old)?,
            prepared_mechanical: prepared,
            old_material_points: archived,
            water: change,
            initial_holding_load_kn: holding,
            zero_release_relative_force_residual: residual,
        },
        fluid,
        porosity_reference: reference,
        archived: false,
    })
}

fn validate_domain(
    mesh: &Mesh,
    case: &Case,
    active: &[bool],
    period: usize,
    projected: &Mesh,
) -> Result<()> {
    let p = &periods(case)?[period];
    let explicit = p
        .boundaries
        .iter()
        .flat_map(|b| projected.face_sets.get(&b.face_set).into_iter().flatten())
        .map(|f| {
            let mut key = *f;
            key.sort_unstable();
            key
        })
        .collect::<BTreeSet<_>>();
    for (key, owners) in boundary::surface_owners(mesh) {
        if matches!(
            case.one_way.mechanics.settings.kinematics,
            stages::Kinematics::Oedometer
        ) && owners.len() == 2
            && owners.iter().filter(|&&e| active[e]).count() == 1
            && key
                .iter()
                .any(|&n| (mesh.nodes[n][2] - mesh.nodes[key[0]][2]).abs() > 1e-8)
        {
            return Err(Error::Input("oedometer domain changes require horizontal cuts; use 3D supports for local excavation".into()));
        }
        if owners.len() == 2
            && owners.iter().filter(|&&e| active[e]).count() == 1
            && !explicit.contains(&key)
        {
            return Err(Error::Input(
                "new cut faces require explicit water boundaries in every subsequent period".into(),
            ));
        }
    }
    Ok(())
}

pub fn solve(mesh: &Mesh, ground: Ground, case: Case) -> Result<Report> {
    advance(mesh, ground, case, None, None)
}
pub fn solve_until(mesh: &Mesh, ground: Ground, case: Case, count: usize) -> Result<Report> {
    if count == 0 {
        return Err(Error::Input(
            "construction checkpoint count must be positive".into(),
        ));
    }
    advance(mesh, ground, case, None, Some(count))
}
pub fn normalized(report: &Report) -> Result<serde_json::Value> {
    fn strip(value: &mut serde_json::Value) {
        match value {
            serde_json::Value::Object(map) => {
                for (k, v) in map {
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
    let mut value = serde_json::to_value(report)?;
    strip(&mut value);
    Ok(value)
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
            "construction resume requires an accepted intermediate state".into(),
        ));
    }
    let count = saved.attempts.iter().filter(|a| a.accepted).count();
    let canonical = solve_until(mesh, ground.clone(), case.clone(), count)?;
    if normalized(&canonical)? != normalized(&saved)? {
        return Err(Error::Input(
            "construction checkpoint input, executable, domain, archive or history changed".into(),
        ));
    }
    drop(canonical);
    advance(mesh, ground, case, Some(saved), additional)
}

fn advance(
    mesh: &Mesh,
    ground: Ground,
    case: Case,
    saved: Option<Report>,
    limit: Option<usize>,
) -> Result<Report> {
    if limit == Some(0) {
        return Err(Error::Input(
            "construction checkpoint count must be positive".into(),
        ));
    }
    let ps = periods(&case)?;
    let mut kernel = stages::ExternalMechanics::new(mesh, &ground, &case.one_way.mechanics)?;
    let starting_active = initial_active(&kernel, &case, mesh.tetrahedra.len())?;
    let seed = if case.initially_inactive_bodies.is_empty() {
        kernel.seed()?
    } else {
        kernel.seed_on_domain(&starting_active)?
    };
    let lifecycle = !case.transitions.is_empty() || !case.initially_inactive_bodies.is_empty();
    // Validate ALL future domains, boundaries, exterior loads and named CAD
    // selections before initialization or any accepted removal is produced.
    let mut dummy = seed.clone();
    let mut domains: Vec<(Arc<Mesh>, richards::Stepper)> = Vec::new();
    let mut previous_active = vec![];
    let mut previous_ids = vec![];
    let mut required_exposed = BTreeSet::new();
    for (i, p) in ps.iter().enumerate() {
        let before = dummy.active_elements.clone();
        if let Some(event) = case.excavations.iter().find(|e| e.before_period == p.name) {
            dummy.active_elements = kernel.excavation_active(&dummy, &event.bodies)?;
        }
        if let Some(event) = case.transitions.iter().find(|e| e.before_period == p.name) {
            let (prepared, selected) =
                kernel.prepare_external_transition(&dummy, &event.stage())?;
            if let Action::MaterialSet { water, .. } = &event.action {
                let inactive = matches!(water, water::WaterPolicy::Inactive);
                if selected
                    .iter()
                    .any(|&e| dummy.active_elements[e] == inactive)
                {
                    return Err(Error::Input(
                        "MATSET water policy does not match selected active/inactive soil".into(),
                    ));
                }
            }
            dummy = prepared;
            required_exposed.extend(new_exposed(mesh, &before, &dummy.active_elements));
        }
        let ids = material_ids(mesh, &dummy);
        let (projected, flow) = if let Some((prior_mesh, prior_flow)) = domains
            .last()
            .filter(|_| previous_active == dummy.active_elements && previous_ids == ids)
        {
            (
                Arc::clone(prior_mesh),
                prior_flow.with_period_boundaries(prior_mesh, ps, i)?,
            )
        } else {
            let projected = Arc::new(domain_mesh(mesh, &dummy.active_elements)?);
            let flow = if lifecycle {
                richards::Stepper::deferred_domain(
                    &projected,
                    &case.one_way.flow,
                    &dummy.active_elements,
                    Some(&ids),
                    i,
                )?
            } else {
                richards::Stepper::for_domain(
                    &projected,
                    &case.one_way.flow,
                    &dummy.active_elements,
                    i,
                )?
            };
            (projected, flow)
        };
        validate_domain(mesh, &case, &dummy.active_elements, i, &projected)?;
        let explicit = ps[i]
            .boundaries
            .iter()
            .flat_map(|b| projected.face_sets.get(&b.face_set).into_iter().flatten())
            .map(|f| {
                let mut k = *f;
                k.sort_unstable();
                k
            })
            .collect::<BTreeSet<_>>();
        let exposed = boundary::surface_owners(mesh)
            .into_iter()
            .filter(|(_, o)| o.iter().filter(|&&e| dummy.active_elements[e]).count() == 1)
            .map(|(k, _)| k)
            .collect::<BTreeSet<_>>();
        if required_exposed
            .intersection(&exposed)
            .any(|k| !explicit.contains(k))
        {
            return Err(Error::Input(
                "new BIRTH surfaces require explicit water conditions in every exposed period"
                    .into(),
            ));
        }
        let t = Transfer::new(&projected, &ground, &case.one_way, &dummy, kernel.layout())?;
        if lifecycle {
            t.exterior_on_domain(
                &case.one_way.exterior_periods[i].reservoirs,
                &dummy.active_elements,
            )?;
        } else {
            t.field_with_base(
                &flow.initial()?.state,
                &case.one_way.exterior_periods[i].reservoirs,
                &kernel.dry_load(&dummy, &dummy.active_elements)?,
            )?;
        }
        domains.push((projected, flow));
        previous_active = dummy.active_elements.clone();
        previous_ids = ids;
    }
    let mut report = if let Some(saved) = saved {
        saved
    } else {
        let mut initial_case = case.one_way.clone();
        initial_case.flow.analysis = richards::Analysis::Transient {
            periods: vec![ps[0].clone()],
        };
        // Initial conditions precede every excavation, including a time-zero cut.
        // Future cut boundaries must not be imposed on this full soil domain.
        let at_zero = case
            .excavations
            .first()
            .is_some_and(|e| e.before_period == ps[0].name)
            || case
                .transitions
                .first()
                .is_some_and(|e| e.before_period == ps[0].name);
        if at_zero || !case.initially_inactive_bodies.is_empty() {
            if let richards::Analysis::Transient { periods } = &mut initial_case.flow.analysis {
                periods[0].boundaries.clear();
            }
            initial_case.exterior_periods[0].reservoirs =
                initial_case.initial_exterior_water.clone();
        }
        initial_case.exterior_periods.truncate(1);
        initial_case.verification_reject_once_at_times_s.clear();
        let initial_mesh = domain_mesh(mesh, &seed.active_elements)?;
        let initialization = if case.initially_inactive_bodies.is_empty() {
            one_way::initialize(&initial_mesh, ground.clone(), initial_case)?
        } else {
            one_way::initialize_on_domain(
                &initial_mesh,
                ground.clone(),
                initial_case,
                &starting_active,
            )?
        };
        if initialization.failure_reason.is_some()
            || initialization.initial_reference_displacement_m.is_empty()
        {
            return Err(Error::Solve(format!(
                "construction initial equilibrium failed: {:?}",
                initialization.failure_reason
            )));
        }
        let initial = initialization.last_committed.clone();
        let water = initial.fluid.water_equivalent_mass_volume_m3;
        Report {
            schema: "soil-fem-one-way-construction-result/1".into(),
            complete: false,
            failure_reason: None,
            mesh_sha256: checkpoint::hash_json(mesh)?,
            ground_sha256: checkpoint::hash_json(&ground)?,
            case_sha256: checkpoint::hash_json(&case)?,
            solver_binary_sha256: checkpoint::hash_file(&std::env::current_exe()?)?,
            ground: ground.clone(),
            case: case.clone(),
            mechanical_integration_points: kernel.layout(),
            initialization,
            initial_water_m3: water,
            frames: vec![],
            attempts: vec![],
            excavations: vec![],
            transitions: vec![],
            last_committed: initial,
            runtime: Some(Runtime {
                period_index: 0,
                event_index: 0,
                pending: None,
                next_dt_s: case.one_way.flow.settings.initial_time_step_s,
                injected_excavations: BTreeSet::new(),
                injected_times: BTreeSet::new(),
                removed_water_m3: 0.,
                net_water_inflow_m3: 0.,
                transition_index: 0,
                pending_transition: None,
                injected_transitions: BTreeSet::new(),
                added_water_m3: 0.,
                porosity_reference: if lifecycle {
                    seed.material_states
                        .iter()
                        .map(|s| s.void_ratio() / (1. + s.void_ratio()))
                        .collect()
                } else {
                    vec![]
                },
            }),
        }
    };
    kernel.set_reference(
        report
            .initialization
            .initial_reference_displacement_m
            .clone(),
    );
    let mut rt = report
        .runtime
        .take()
        .ok_or_else(|| Error::Input("missing construction runtime".into()))?;
    let mut accepted = 0;
    while rt.period_index < ps.len() {
        if report.attempts.len() >= case.one_way.flow.settings.maximum_steps {
            report.failure_reason = Some("construction attempt budget exhausted".into());
            break;
        }
        let old = &report.last_committed;
        if rt.pending.is_none()
            && rt.pending_transition.is_none()
            && case
                .excavations
                .get(rt.event_index)
                .is_some_and(|e| e.before_period == ps[rt.period_index].name)
        {
            let event = &case.excavations[rt.event_index];
            let active = kernel.excavation_active(&old.mechanical, &event.bodies)?;
            let ids = material_ids(mesh, &old.mechanical);
            let removal = richards::remove_domain_materials(
                mesh,
                &case.one_way.flow,
                &old.mechanical.active_elements,
                &active,
                Some(&ids),
                &old.hydraulic,
            )?;
            let projected = &domains[rt.period_index].0;
            let mut transfer =
                Transfer::new(projected, &ground, &case.one_way, &seed, kernel.layout())?;
            transfer.set_material_ids(&ids)?;
            if !rt.porosity_reference.is_empty() {
                transfer.set_porosity_reference(&rt.porosity_reference)?;
            }
            let fluid = transfer.field_with_base(
                &removal.state,
                &case.one_way.exterior_periods[rt.period_index].reservoirs,
                &kernel.dry_load(&old.mechanical, &active)?,
            )?;
            let (holding, residual) = kernel.cut_holding(
                &old.mechanical,
                &active,
                &fluid.equivalent_mechanical_rhs_kn,
            )?;
            if residual > case.one_way.mechanics.settings.relative_force_tolerance * 10. {
                return Err(Error::Solve(
                    "construction zero-release force identity failed".into(),
                ));
            }
            let removed: Vec<usize> = report
                .mechanical_integration_points
                .iter()
                .enumerate()
                .filter(|(_, p)| {
                    old.mechanical.active_elements[p.element_index] && !active[p.element_index]
                })
                .map(|(i, _)| i)
                .collect();
            let removed_material_points = if lifecycle {
                removed
                    .iter()
                    .map(|&i| ArchivedPoint {
                        integration_point_index: i,
                        material_id: ids[report.mechanical_integration_points[i].element_index],
                        state: old.mechanical.material_states[i].clone(),
                        active_before: true,
                    })
                    .collect()
            } else {
                vec![]
            };
            rt.pending = Some(Pending {
                event_index: rt.event_index,
                progress: 0.,
                increment: event
                    .initial_increment
                    .unwrap_or(case.one_way.mechanics.settings.initial_increment),
                active,
                fluid,
                record: ExcavationRecord {
                    name: event.name.clone(),
                    time_s: old.hydraulic.time_s,
                    origin_pair_sha256: checkpoint::hash_json(old)?,
                    origin_hydraulic: old.hydraulic.clone(),
                    removal,
                    initial_holding_load_kn: holding,
                    zero_release_relative_force_residual: residual,
                    removed_mechanical_points: removed,
                    removed_material_points,
                },
                archived: false,
            });
        }
        if rt.pending.is_none()
            && rt.pending_transition.is_none()
            && case
                .transitions
                .get(rt.transition_index)
                .is_some_and(|e| e.before_period == ps[rt.period_index].name)
        {
            let event = &case.transitions[rt.transition_index];
            match prepare_pending(
                &kernel,
                mesh,
                &ground,
                &case,
                old,
                event,
                &rt.porosity_reference,
            ) {
                Ok(mut p) => {
                    p.event_index = rt.transition_index;
                    rt.pending_transition = Some(p);
                }
                Err(error) => {
                    let hash = checkpoint::hash_json(old)?;
                    let reason = error.to_string();
                    report.attempts.push(Attempt {
                        name: event.name.clone(),
                        phase: "transition_prepare".into(),
                        time_s: old.hydraulic.time_s,
                        dt_s: None,
                        progress_from: 0.,
                        progress_to: 0.,
                        accepted: false,
                        reason: Some(reason.clone()),
                        candidate_changed: false,
                        plastic_candidate_changed: false,
                        committed_pair_sha256_before: hash.clone(),
                        committed_pair_sha256_after: hash,
                        hydraulic_iterations: vec![],
                    });
                    report.failure_reason = Some(reason);
                    break;
                }
            }
        }
        let cutting = rt.pending.is_some();
        let transitioning = rt.pending_transition.is_some();
        let (name, from, to, dt) = if let Some(p) = &rt.pending {
            (
                p.record.name.clone(),
                p.progress,
                (p.progress + p.increment).min(1.),
                None,
            )
        } else if let Some(p) = &rt.pending_transition {
            (
                p.record.name.clone(),
                p.progress,
                (p.progress + p.increment).min(1.),
                None,
            )
        } else {
            (
                ps[rt.period_index].name.clone(),
                0.,
                1.,
                Some(
                    rt.next_dt_s
                        .min(ps[rt.period_index].end_time_s - old.hydraulic.time_s),
                ),
            )
        };
        let hash = checkpoint::hash_json(old)?;
        let mut log = vec![];
        let mut candidate = None;
        let mut changed = false;
        let mut plastic = false;
        let result = (|| -> Result<()> {
            let active = if let Some(p) = &rt.pending_transition {
                &p.record.prepared_mechanical.active_elements
            } else {
                rt.pending
                    .as_ref()
                    .map_or(&old.mechanical.active_elements, |p| &p.active)
            };
            let (projected, flow) = &domains[rt.period_index];
            let (hydraulic, balance, fluid, holding, trial) = if let Some(p) = &rt.pending {
                let holding = p
                    .record
                    .initial_holding_load_kn
                    .iter()
                    .map(|h| if to == 1. { 0. } else { (1. - to) * h })
                    .collect::<Vec<_>>();
                let trial = kernel.excavation_trial(
                    &old.mechanical,
                    p.active.clone(),
                    p.fluid.equivalent_mechanical_rhs_kn.clone(),
                    holding.clone(),
                    to,
                    name.clone(),
                )?;
                (
                    p.record.removal.state.clone(),
                    None,
                    p.fluid.clone(),
                    holding,
                    trial,
                )
            } else if let Some(p) = &rt.pending_transition {
                let holding = p
                    .record
                    .initial_holding_load_kn
                    .iter()
                    .map(|h| if to == 1. { 0. } else { (1. - to) * h })
                    .collect::<Vec<_>>();
                let origin = if p.progress == 0. {
                    &p.record.prepared_mechanical
                } else {
                    &old.mechanical
                };
                let trial = kernel.transition_trial(
                    origin,
                    p.fluid.equivalent_mechanical_rhs_kn.clone(),
                    holding.clone(),
                    to,
                    name.clone(),
                )?;
                (
                    p.record.water.state.clone(),
                    None,
                    p.fluid.clone(),
                    holding,
                    trial,
                )
            } else {
                let (hyd, iterations) = flow.trial(
                    &old.hydraulic,
                    rt.period_index,
                    dt.ok_or_else(|| Error::Input("missing physical dt".into()))?,
                );
                log = iterations;
                let hyd = hyd?;
                let mut transfer =
                    Transfer::new(projected, &ground, &case.one_way, &seed, kernel.layout())?;
                transfer.set_material_ids(&material_ids(mesh, &old.mechanical))?;
                if !rt.porosity_reference.is_empty() {
                    transfer.set_porosity_reference(&rt.porosity_reference)?;
                }
                let fluid = transfer.field_with_base(
                    &hyd.state,
                    &case.one_way.exterior_periods[rt.period_index].reservoirs,
                    &kernel.dry_load(&old.mechanical, active)?,
                )?;
                let trial = kernel.trial(
                    &old.mechanical,
                    fluid.equivalent_mechanical_rhs_kn.clone(),
                    name.clone(),
                )?;
                (
                    hyd.state,
                    Some(hyd.balance),
                    fluid,
                    vec![0.; mesh.nodes.len() * 3],
                    trial,
                )
            };
            plastic = old
                .mechanical
                .material_states
                .iter()
                .zip(&trial.frame.state.material_states)
                .any(|(a, b)| a.plastic_strain() != b.plastic_strain());
            let mut audit_fluid = fluid.clone();
            for (f, h) in audit_fluid.physical_nodal_load_kn.iter_mut().zip(&holding) {
                *f += h;
            }
            let audit = one_way::audit(&trial, &audit_fluid);
            let state = State {
                hydraulic,
                mechanical: trial.frame.state,
                fluid,
            };
            changed = checkpoint::hash_json(&state)? != hash;
            let removed = rt.removed_water_m3
                + rt.pending
                    .as_ref()
                    .filter(|p| !p.archived)
                    .map_or(0., |p| p.record.removal.water_removed_m3)
                + rt.pending_transition
                    .as_ref()
                    .filter(|p| !p.archived)
                    .map_or(0., |p| p.record.water.water_removed_m3);
            let added = rt.added_water_m3
                + rt.pending_transition
                    .as_ref()
                    .filter(|p| !p.archived)
                    .map_or(0., |p| p.record.water.water_added_m3);
            let inflow = rt.net_water_inflow_m3
                + balance
                    .as_ref()
                    .zip(dt)
                    .map_or(0., |(b, dt)| dt * (b.source_m3_s - b.boundary_outflow_m3_s));
            let water_error = state.fluid.water_equivalent_mass_volume_m3 + removed
                - added
                - report.initial_water_m3
                - inflow;
            candidate = Some(Frame {
                name: name.clone(),
                phase: if transitioning {
                    "material_transition"
                } else if cutting {
                    "excavation_release"
                } else {
                    "transient"
                }
                .into(),
                progress: to,
                hydraulic_faces: flow.faces(),
                hydraulic_cell_element_indices: active
                    .iter()
                    .enumerate()
                    .filter_map(|(e, a)| a.then_some(e))
                    .collect(),
                hydraulic_balance: balance,
                mechanical_audit: audit,
                holding_load_kn: holding,
                cumulative_removed_water_m3: removed,
                cumulative_net_water_inflow_m3: inflow,
                signed_global_water_budget_error_m3: water_error,
                cumulative_added_water_m3: added,
                hydraulic_material_ids: if lifecycle {
                    material_ids(mesh, &state.mechanical)
                } else {
                    vec![]
                },
                state,
            });
            if cutting
                && case.verification_reject_once_excavations.contains(&name)
                && !rt.injected_excavations.contains(&name)
            {
                rt.injected_excavations.insert(name.clone());
                return Err(Error::Solve(
                    "verification rejection after changed domain/fluid/mechanical candidate".into(),
                ));
            }
            if transitioning
                && case.verification_reject_once_transitions.contains(&name)
                && !rt.injected_transitions.contains(&name)
            {
                rt.injected_transitions.insert(name.clone());
                return Err(Error::Solve("verification rejection after lifecycle candidate changed; restoring water, material, domain and ledger".into()));
            }
            if !cutting && !transitioning {
                for (i, t) in case
                    .one_way
                    .verification_reject_once_at_times_s
                    .iter()
                    .enumerate()
                {
                    if !rt.injected_times.contains(&i)
                        && candidate
                            .as_ref()
                            .is_some_and(|f| f.state.hydraulic.time_s == *t)
                    {
                        rt.injected_times.insert(i);
                        return Err(Error::Solve(
                            "verification rejection after physical coupled candidate".into(),
                        ));
                    }
                }
            }
            Ok(())
        })();
        let success = result.is_ok();
        let reason = result.as_ref().err().map(ToString::to_string);
        let invalid = matches!(result, Err(Error::Input(_)));
        let start_time = old.hydraulic.time_s;
        if success {
            let frame =
                candidate.ok_or_else(|| Error::Solve("missing construction candidate".into()))?;
            rt.removed_water_m3 = frame.cumulative_removed_water_m3;
            rt.added_water_m3 = frame.cumulative_added_water_m3;
            rt.net_water_inflow_m3 = frame.cumulative_net_water_inflow_m3;
            report.last_committed = frame.state.clone();
            report.frames.push(frame);
            accepted += 1;
            if let Some(p) = &mut rt.pending {
                if !p.archived {
                    report.excavations.push(p.record.clone());
                    p.archived = true;
                }
                p.progress = to;
                p.increment =
                    (p.increment * 1.5).min(case.one_way.mechanics.settings.initial_increment);
                if to == 1. {
                    rt.pending = None;
                    rt.event_index += 1;
                }
            } else if let Some(p) = &mut rt.pending_transition {
                if !p.archived {
                    report.transitions.push(p.record.clone());
                    p.archived = true;
                    rt.porosity_reference = p.porosity_reference.clone();
                }
                p.progress = to;
                p.increment =
                    (p.increment * 1.5).min(case.one_way.mechanics.settings.initial_increment);
                if to == 1. {
                    rt.pending_transition = None;
                    rt.transition_index += 1;
                }
            } else {
                if report.last_committed.hydraulic.time_s >= ps[rt.period_index].end_time_s {
                    rt.period_index += 1;
                }
                rt.next_dt_s =
                    (rt.next_dt_s * 1.5).min(case.one_way.flow.settings.maximum_time_step_s);
            }
        }
        let after = checkpoint::hash_json(&report.last_committed)?;
        if !success && after != hash {
            return Err(Error::Solve(
                "rejected construction changed committed triple".into(),
            ));
        }
        report.attempts.push(Attempt {
            name,
            phase: if transitioning {
                "material_transition"
            } else if cutting {
                "excavation_release"
            } else {
                "transient"
            }
            .into(),
            time_s: start_time,
            dt_s: dt,
            progress_from: from,
            progress_to: to,
            accepted: success,
            reason: reason.clone(),
            candidate_changed: changed,
            plastic_candidate_changed: plastic,
            committed_pair_sha256_before: hash,
            committed_pair_sha256_after: after,
            hydraulic_iterations: log,
        });
        if success {
            if limit == Some(accepted) {
                break;
            }
        } else if !invalid
            && rt.pending_transition.as_ref().is_some_and(|p| {
                p.increment / 2. >= case.one_way.mechanics.settings.minimum_increment
            })
        {
            if let Some(p) = &mut rt.pending_transition {
                p.increment /= 2.;
            }
        } else if !invalid
            && rt.pending.as_ref().is_some_and(|p| {
                p.increment / 2. >= case.one_way.mechanics.settings.minimum_increment
            })
        {
            if let Some(p) = &mut rt.pending {
                p.increment /= 2.;
            }
        } else if !invalid
            && !cutting
            && !transitioning
            && dt.is_some_and(|d| d / 2. >= case.one_way.flow.settings.minimum_time_step_s)
        {
            rt.next_dt_s = dt.unwrap_or(0.) / 2.;
        } else {
            report.failure_reason = reason;
            break;
        }
    }
    report.complete = report.failure_reason.is_none()
        && rt.period_index == ps.len()
        && rt.pending.is_none()
        && rt.pending_transition.is_none()
        && rt.event_index == case.excavations.len()
        && rt.transition_index == case.transitions.len();
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
                "construction checkpoint needs an accepted intermediate state".into(),
            ));
        }
        Ok(Self {
            schema: "soil-fem-one-way-construction-checkpoint/1".into(),
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
        if self.schema != "soil-fem-one-way-construction-checkpoint/1"
            || self.run_sha256 != checkpoint::hash_json(&self.run)?
        {
            return Err(Error::Input(
                "construction checkpoint payload hash changed".into(),
            ));
        }
        resume(mesh, ground, case, self.run, additional)
    }
}
