//! Hybrid RT0/P0 mixed Richards equation, conservative backward Euler.
//! Total head is cellwise constant; conductivity/storage use psi=H-z at
//! quadrature points. Atmospheric boundaries have explicit NO-PONDING policy.
//! This hydraulic module does not yet modify mechanical effective stresses.
use crate::{
    Error, Result, checkpoint, darcy, elements,
    linalg::{CsrMatrix, CsrPattern, LinearReport, LinearSolver},
    mesh::Mesh,
    mixed_flow,
    retention::VanGenuchten,
};
use nalgebra::{Matrix3, Matrix4, Vector4};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub id: u32,
    pub saturated_conductivity_m_s: [[f64; 3]; 3],
    pub retention: VanGenuchten,
    #[serde(default)]
    pub volume_source_per_s: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineHead {
    pub head_m: f64,
    #[serde(default)]
    pub gradient: [f64; 3],
    #[serde(default)]
    pub reference_point_m: [f64; 3],
}
impl AffineHead {
    pub(crate) fn validate(&self) -> Result<()> {
        if !self
            .gradient
            .iter()
            .chain(&self.reference_point_m)
            .chain(std::iter::once(&self.head_m))
            .all(|v| v.is_finite())
        {
            return Err(Error::Input("non-finite hydraulic head definition".into()));
        }
        Ok(())
    }
    fn relative(&self, p: [f64; 3], datum: f64) -> f64 {
        self.head_m - datum
            + (0..3)
                .map(|a| self.gradient[a] * (p[a] - self.reference_point_m[a]))
                .sum::<f64>()
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Condition {
    Head {
        value: AffineHead,
    },
    Pressure {
        pressure_head_m: f64,
    },
    OutwardFlux {
        flux_m_s: f64,
    },
    /// Vertical rainfall [m/s per horizontal area], immediately routed away
    /// when infiltration is limited. No ponding/storage is assumed.
    AtmosphericNoPonding {
        rainfall_vertical_m_s: f64,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Boundary {
    pub face_set: String,
    pub condition: Condition,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Period {
    pub name: String,
    pub end_time_s: f64,
    pub boundaries: Vec<Boundary>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Analysis {
    Steady { boundaries: Vec<Boundary> },
    Transient { periods: Vec<Period> },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub schema: String,
    pub units: String,
    pub materials: Vec<Material>,
    /// Subtracted from all total-head unknowns; not a hydraulic boundary.
    #[serde(default)]
    pub head_datum_m: f64,
    pub initial_head: AffineHead,
    pub analysis: Analysis,
    #[serde(default)]
    pub settings: Settings,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Settings {
    pub relative_flow_tolerance: f64,
    pub absolute_flow_tolerance_m3_s: f64,
    pub head_tolerance_m: f64,
    pub maximum_iterations: usize,
    pub initial_time_step_s: f64,
    pub minimum_time_step_s: f64,
    pub maximum_time_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_free_faces: usize,
    pub linear_solver: LinearSolver,
    pub store_integration_points: bool,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            relative_flow_tolerance: 1e-8,
            absolute_flow_tolerance_m3_s: 1e-16,
            head_tolerance_m: 1e-8,
            maximum_iterations: 100,
            initial_time_step_s: 10.,
            minimum_time_step_s: 1e-4,
            maximum_time_step_s: 100.,
            maximum_steps: 10000,
            maximum_free_faces: 300000,
            linear_solver: LinearSolver::SparseGmresIlu0 {
                options: Default::default(),
            },
            store_integration_points: true,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Face {
    pub corners: [usize; 3],
    /// [position in Report.cell_element_indices, local face opposite vertex].
    pub owners: Vec<[usize; 2]>,
    pub area_m2: f64,
    pub projected_horizontal_area_m2: f64,
    #[serde(default)]
    pub projected_horizontal_area_error_bound_m2: f64,
    pub reference_elevation_m: f64,
    pub reference_position_m: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Iteration {
    pub iteration: usize,
    pub scaled_residual: f64,
    pub line_search_fraction: f64,
    pub linear: Option<LinearReport>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Attempt {
    pub period: String,
    pub start_time_s: f64,
    pub time_step_s: Option<f64>,
    pub accepted: bool,
    pub iterations: Vec<Iteration>,
    pub failure_reason: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub time_s: f64,
    pub relative_cell_head_m: Vec<Option<f64>>,
    pub relative_face_head_m: Vec<f64>,
    pub cell_outflow_m3_s: Vec<Option<[f64; 4]>>,
    pub cell_storage_offset_m3: Vec<Option<f64>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Point {
    pub element_index: usize,
    pub local_index: usize,
    pub position_m: [f64; 3],
    pub weight_m3: f64,
    pub total_head_m: f64,
    pub pressure_head_m: f64,
    pub retention: crate::retention::State,
    pub darcy_flux_m_s: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Balance {
    pub maximum_absolute_cell_mass_residual_m3_s: f64,
    pub maximum_absolute_face_flux_residual_m3_s: f64,
    pub maximum_scaled_residual: f64,
    pub source_m3_s: f64,
    pub boundary_outflow_m3_s: f64,
    pub storage_rate_m3_s: f64,
    pub absolute_global_mass_residual_m3_s: f64,
    pub precipitation_m3_s: f64,
    pub surface_routed_water_m3_s: f64,
    pub atmospheric_faces: Vec<AtmosphericFace>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtmosphericFace {
    pub face_index: usize,
    pub projected_pressure_head_m: f64,
    pub outward_flux_m3_s: f64,
    pub precipitation_m3_s: f64,
    pub surface_routed_water_m3_s: f64,
    pub wet: bool,
    /// Positive outward seepage, distinct from a zero-pressure rain inlet.
    pub seeping: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    pub name: String,
    pub state: State,
    pub balance: Balance,
    pub points: Vec<Point>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub complete: bool,
    pub failure_reason: Option<String>,
    pub mesh_sha256: String,
    pub case_sha256: String,
    pub solver_binary_sha256: String,
    pub case: Case,
    pub active_elements: Vec<bool>,
    pub material_ids: Vec<u32>,
    pub cell_element_indices: Vec<usize>,
    pub faces: Vec<Face>,
    pub frames: Vec<Frame>,
    pub attempts: Vec<Attempt>,
    pub last_committed: State,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub runtime: Option<Runtime>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Runtime {
    pub period_index: usize,
    pub next_time_step_s: f64,
}

struct Model {
    cells: Vec<mixed_flow::Cell>,
    materials: Vec<Material>,
    inverse_k: Vec<Matrix3<f64>>,
    faces: Vec<Face>,
    cell_faces: Vec<[usize; 4]>,
    active: Vec<bool>,
    ids: Vec<u32>,
    datum: f64,
}
#[derive(Clone)]
enum Bc {
    Internal,
    Flux(f64),
    Head(f64),
    Atmosphere(f64),
}
struct Local {
    c: Matrix4<f64>,
    mass_derivative: Matrix4<f64>,
    b: Vector4<f64>,
    storage: f64,
    capacity: f64,
    source: f64,
}
struct Evaluated {
    locals: Vec<Local>,
    fluxes: Vec<Vector4<f64>>,
    balance: Balance,
}

impl Model {
    fn new(mesh: &Mesh, case: &Case, active: &[bool], ids: Option<&[u32]>) -> Result<Self> {
        mesh.validate()?;
        case.initial_head.validate()?;
        case.settings.linear_solver.validate()?;
        let s = &case.settings;
        if case.schema != "soil-fem-richards/1"
            || case.units != "m-s"
            || !case.head_datum_m.is_finite()
            || active.len() != mesh.tetrahedra.len()
            || !active.iter().any(|v| *v)
            || ![
                s.relative_flow_tolerance,
                s.absolute_flow_tolerance_m3_s,
                s.head_tolerance_m,
                s.initial_time_step_s,
                s.minimum_time_step_s,
                s.maximum_time_step_s,
            ]
            .iter()
            .all(|v| v.is_finite() && *v > 0.)
            || !(1e-12..=1e-4).contains(&s.relative_flow_tolerance)
            || s.head_tolerance_m > 1e-3
            || s.minimum_time_step_s > s.initial_time_step_s
            || s.initial_time_step_s > s.maximum_time_step_s
            || s.maximum_iterations == 0
            || s.maximum_iterations > 1000
            || s.maximum_steps == 0
            || s.maximum_steps > 100000
            || s.maximum_free_faces == 0
            || s.maximum_free_faces > s.linear_solver.dimension_limit()
        {
            return Err(Error::Input(
                "invalid Richards schema, units, domain or settings".into(),
            ));
        }
        let ids = ids.unwrap_or(&mesh.material_ids).to_vec();
        if ids.len() != active.len() {
            return Err(Error::Input("Richards material map length mismatch".into()));
        }
        let mut map = BTreeMap::new();
        for m in &case.materials {
            m.retention.validate()?;
            let k = darcy::tensor(&darcy::Material {
                id: m.id,
                hydraulic_conductivity_m_s: m.saturated_conductivity_m_s,
                volume_source_per_s: m.volume_source_per_s,
            })?;
            let inverse = k
                .try_inverse()
                .ok_or_else(|| Error::Input("singular Richards conductivity".into()))?;
            if map.insert(m.id, (m.clone(), inverse)).is_some() {
                return Err(Error::Input("duplicate Richards material id".into()));
            }
        }
        let mut cells = Vec::new();
        let mut materials = Vec::new();
        let mut inverse_k = Vec::new();
        let mut owners: BTreeMap<[usize; 3], Vec<[usize; 2]>> = BTreeMap::new();
        for (e, &live) in active.iter().enumerate() {
            if live {
                let (mat, inv) = map
                    .get(&ids[e])
                    .ok_or_else(|| Error::Input("missing Richards material".into()))?;
                let i = cells.len();
                cells.push(mixed_flow::Cell::new(mesh, e)?);
                materials.push(mat.clone());
                inverse_k.push(*inv);
                for j in 0..4 {
                    let mut key = elements::FACES[j].map(|v| mesh.tetrahedra[e][v]);
                    key.sort_unstable();
                    owners.entry(key).or_default().push([i, j]);
                }
            }
        }
        let mut faces = Vec::new();
        let mut cell_faces = vec![[0; 4]; cells.len()];
        for (key, o) in owners {
            if o.len() > 2 {
                return Err(Error::Input("nonmanifold hydraulic face".into()));
            }
            let [i, j] = o[0];
            let index = faces.len();
            for &[ci, cj] in &o {
                cell_faces[ci][cj] = index;
            }
            let coords: Vec<_> = mesh
                .face_nodes(cells[i].element_index, &key)?
                .iter()
                .map(|&n| mesh.nodes[n])
                .collect();
            let horizontal = crate::surface_projection::upward_area(&coords)?;
            let (surface_area, _) = crate::surface_projection::surface_area(&coords)?;
            let range = if coords.len() == 3 { 0..3 } else { 3..6 };
            let mean =
                std::array::from_fn(|a| range.clone().map(|j| coords[j][a]).sum::<f64>() / 3.);
            faces.push(Face {
                corners: key,
                owners: o,
                area_m2: surface_area,
                projected_horizontal_area_m2: horizontal.area_m2,
                projected_horizontal_area_error_bound_m2: horizontal.domain_error_bound_m2,
                reference_elevation_m: cells[i].face_elevation_m[j],
                reference_position_m: mean,
            });
        }
        Ok(Self {
            cells,
            materials,
            inverse_k,
            faces,
            cell_faces,
            active: active.to_vec(),
            ids,
            datum: case.head_datum_m,
        })
    }
    fn boundaries(&self, mesh: &Mesh, definitions: &[Boundary]) -> Result<Vec<Bc>> {
        let face_map: BTreeMap<_, _> = self
            .faces
            .iter()
            .enumerate()
            .map(|(i, f)| (f.corners, i))
            .collect();
        let mut bc: Vec<_> = self
            .faces
            .iter()
            .map(|f| {
                if f.owners.len() == 2 {
                    Bc::Internal
                } else {
                    Bc::Flux(0.)
                }
            })
            .collect();
        let mut used = BTreeSet::new();
        for b in definitions {
            let set = mesh
                .face_sets
                .get(&b.face_set)
                .filter(|v| !v.is_empty())
                .ok_or_else(|| {
                    Error::Input(format!("unknown/empty Richards face set {}", b.face_set))
                })?;
            match &b.condition {
                Condition::Head { value } => value.validate()?,
                Condition::Pressure { pressure_head_m } if !pressure_head_m.is_finite() => {
                    return Err(Error::Input("non-finite pressure boundary".into()));
                }
                Condition::OutwardFlux { flux_m_s } if !flux_m_s.is_finite() => {
                    return Err(Error::Input("non-finite flux boundary".into()));
                }
                Condition::AtmosphericNoPonding {
                    rainfall_vertical_m_s,
                } if !rainfall_vertical_m_s.is_finite() || *rainfall_vertical_m_s < 0. => {
                    return Err(Error::Input(
                        "rainfall must be finite and nonnegative".into(),
                    ));
                }
                _ => {}
            }
            for mut key in set.clone() {
                key.sort_unstable();
                let &i = face_map.get(&key).ok_or_else(|| {
                    Error::Input("inactive/unknown Richards boundary face".into())
                })?;
                let f = &self.faces[i];
                if f.owners.len() != 1 || !used.insert(i) {
                    return Err(Error::Input(
                        "overlapping/internal Richards boundary".into(),
                    ));
                }
                bc[i] = match &b.condition {
                    Condition::Pressure { pressure_head_m } => {
                        Bc::Head(f.reference_elevation_m - self.datum + pressure_head_m)
                    }
                    Condition::Head { value } => {
                        Bc::Head(value.relative(f.reference_position_m, self.datum))
                    }
                    Condition::OutwardFlux { flux_m_s } => Bc::Flux(flux_m_s * f.area_m2),
                    Condition::AtmosphericNoPonding {
                        rainfall_vertical_m_s,
                    } => Bc::Atmosphere(-rainfall_vertical_m_s * f.projected_horizontal_area_m2),
                };
                let value = match bc[i] {
                    Bc::Head(v) | Bc::Flux(v) | Bc::Atmosphere(v) => v,
                    Bc::Internal => 0.,
                };
                if !value.is_finite() {
                    return Err(Error::Input("Richards boundary value/flow overflow".into()));
                }
            }
        }
        Ok(bc)
    }
    fn local(&self, i: usize, h: f64) -> Result<Local> {
        if !h.is_finite() || !(self.datum + h).is_finite() {
            return Err(Error::Solve("hydraulic head overflow".into()));
        }
        let cell = &self.cells[i];
        let law = &self.materials[i].retention;
        let mut inverses = Vec::new();
        let mut derivatives = Vec::new();
        let mut storage = 0.;
        let mut capacity = 0.;
        for q in &cell.points {
            let state = law.evaluate((self.datum - q.position_m[2]) + h)?;
            if state.relative_conductivity == 0. {
                return Err(Error::Solve("relative conductivity underflow: refine/change the physical initial state; no numerical floor was applied".into()));
            }
            inverses.push(self.inverse_k[i] / state.relative_conductivity);
            derivatives.push(
                -self.inverse_k[i]
                    * (state.d_log_relative_conductivity_per_m / state.relative_conductivity),
            );
            storage += q.weight_m3 * state.storage_offset;
            capacity += q.weight_m3 * state.capacity_per_m;
        }
        let mass = cell.mass(&inverses)?;
        let mass_derivative = cell.weighted_mass(&derivatives)?;
        let c = mass
            .cholesky()
            .ok_or_else(|| Error::Solve("RT0 mass Cholesky failed".into()))?
            .inverse();
        let b = c * Vector4::repeat(1.);
        let a = b.sum();
        let source = cell.volume_m3 * self.materials[i].volume_source_per_s;
        if !c
            .iter()
            .chain(b.iter())
            .chain([a, storage, capacity, source].iter())
            .all(|v| v.is_finite())
            || a <= 0.
            || capacity < 0.
        {
            return Err(Error::Solve("invalid Richards local condensation".into()));
        }
        Ok(Local {
            c,
            mass_derivative,
            b,
            storage,
            capacity,
            source,
        })
    }
    fn evaluate(
        &self,
        state: &State,
        old: &State,
        dt: Option<f64>,
        bc: &[Bc],
        s: &Settings,
    ) -> Result<Evaluated> {
        let mut locals = Vec::new();
        let mut fluxes = Vec::new();
        let mut mass_error: f64 = 0.;
        let mut face_error: f64 = 0.;
        let mut merit: f64 = 0.;
        let mut source = 0.;
        let mut storage_rate = 0.;
        for (i, cell) in self.cells.iter().enumerate() {
            let e = cell.element_index;
            let h = state.relative_cell_head_m[e]
                .ok_or_else(|| Error::Input("undefined hydraulic cell head".into()))?;
            let l = self.local(i, h)?;
            let difference =
                Vector4::from_fn(|j, _| h - state.relative_face_head_m[self.cell_faces[i][j]]);
            let flux = l.c * difference;
            let rate = dt.map_or(0., |dt| {
                (l.storage - old.cell_storage_offset_m3[e].unwrap_or(0.)) / dt
            });
            let r = flux.sum() + rate - l.source;
            let scale = flux.abs().sum() + rate.abs() + l.source.abs();
            mass_error = mass_error.max(r.abs());
            merit = merit.max(
                r.abs() / (s.absolute_flow_tolerance_m3_s + s.relative_flow_tolerance * scale),
            );
            source += l.source;
            storage_rate += rate;
            locals.push(l);
            fluxes.push(flux);
        }
        let mut boundary_outflow = 0.;
        let mut precipitation = 0.;
        let mut routed = 0.;
        let mut atmospheric = Vec::new();
        for (fi, f) in self.faces.iter().enumerate() {
            let q: f64 = f.owners.iter().map(|&[i, j]| fluxes[i][j]).sum();
            let scale: f64 = f.owners.iter().map(|&[i, j]| fluxes[i][j].abs()).sum();
            if f.owners.len() == 1 {
                boundary_outflow += q;
            }
            let r = match bc[fi] {
                Bc::Internal => q,
                Bc::Flux(target) => q - target,
                Bc::Head(target) => {
                    merit = merit
                        .max((state.relative_face_head_m[fi] - target).abs() / s.head_tolerance_m);
                    0.
                }
                Bc::Atmosphere(rain) => {
                    let [i, j] = f.owners[0];
                    let t = locals[i].c[(j, j)];
                    let psi =
                        (self.datum - f.reference_elevation_m) + state.relative_face_head_m[fi];
                    let excess = q - rain;
                    merit = merit.max((-psi).min(excess / t).abs() / s.head_tolerance_m);
                    let flow_violation = (-excess).max(0.);
                    merit = merit.max(
                        flow_violation
                            / (s.absolute_flow_tolerance_m3_s
                                + s.relative_flow_tolerance * (scale + rain.abs())),
                    );
                    merit = merit.max(psi.max(0.) / s.head_tolerance_m);
                    precipitation -= rain;
                    routed += excess;
                    atmospheric.push(AtmosphericFace {
                        face_index: fi,
                        projected_pressure_head_m: psi,
                        outward_flux_m3_s: q,
                        precipitation_m3_s: -rain,
                        surface_routed_water_m3_s: excess,
                        wet: psi.abs() <= s.head_tolerance_m,
                        seeping: psi.abs() <= s.head_tolerance_m
                            && q > s.absolute_flow_tolerance_m3_s
                                + s.relative_flow_tolerance * (q.abs() + rain.abs()),
                    });
                    0.
                }
            };
            face_error = face_error.max(r.abs());
            let target = match bc[fi] {
                Bc::Flux(q) => q.abs(),
                _ => 0.,
            };
            merit = merit.max(
                r.abs()
                    / (s.absolute_flow_tolerance_m3_s
                        + s.relative_flow_tolerance * (scale + target)),
            );
        }
        let global = (boundary_outflow + storage_rate - source).abs();
        if ![
            merit,
            mass_error,
            face_error,
            source,
            storage_rate,
            boundary_outflow,
            precipitation,
            routed,
            global,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(Error::Solve("Richards residual overflow".into()));
        }
        Ok(Evaluated {
            locals,
            fluxes,
            balance: Balance {
                maximum_absolute_cell_mass_residual_m3_s: mass_error,
                maximum_absolute_face_flux_residual_m3_s: face_error,
                maximum_scaled_residual: merit,
                source_m3_s: source,
                boundary_outflow_m3_s: boundary_outflow,
                storage_rate_m3_s: storage_rate,
                absolute_global_mass_residual_m3_s: global,
                precipitation_m3_s: precipitation,
                surface_routed_water_m3_s: routed,
                atmospheric_faces: atmospheric,
            },
        })
    }
    fn initial(&self, case: &Case) -> Result<State> {
        let mut state = State {
            time_s: 0.,
            relative_cell_head_m: vec![None; self.active.len()],
            relative_face_head_m: vec![0.; self.faces.len()],
            cell_outflow_m3_s: vec![None; self.active.len()],
            cell_storage_offset_m3: vec![None; self.active.len()],
        };
        for (i, c) in self.cells.iter().enumerate() {
            let anchor = c.points[0].position_m;
            let correction: f64 = c
                .points
                .iter()
                .map(|q| {
                    q.weight_m3
                        * (0..3)
                            .map(|a| case.initial_head.gradient[a] * (q.position_m[a] - anchor[a]))
                            .sum::<f64>()
                })
                .sum::<f64>()
                / c.volume_m3;
            let h = case.initial_head.relative(anchor, self.datum) + correction;
            state.relative_cell_head_m[c.element_index] = Some(h);
            state.cell_storage_offset_m3[c.element_index] = Some(self.local(i, h)?.storage);
        }
        // Projection of the specified initial affine field on each reference face.
        for (fi, f) in self.faces.iter().enumerate() {
            state.relative_face_head_m[fi] = case
                .initial_head
                .relative(f.reference_position_m, self.datum);
        }
        for (i, c) in self.cells.iter().enumerate() {
            let h = state.relative_cell_head_m[c.element_index].unwrap_or(0.);
            let q = self.local(i, h)?.c
                * Vector4::from_fn(|j, _| h - state.relative_face_head_m[self.cell_faces[i][j]]);
            state.cell_outflow_m3_s[c.element_index] = Some(q.into());
        }
        Ok(state)
    }
    fn newton_step(
        &self,
        state: &State,
        old: &State,
        dt: Option<f64>,
        bc: &[Bc],
        eval: &Evaluated,
        s: &Settings,
    ) -> Result<(State, Option<LinearReport>)> {
        let mut fixed = vec![None; self.faces.len()];
        let mut targets = vec![0.; self.faces.len()];
        for (fi, f) in self.faces.iter().enumerate() {
            match bc[fi] {
                Bc::Head(h) => fixed[fi] = Some(h),
                Bc::Flux(q) => targets[fi] = q,
                Bc::Atmosphere(rain) => {
                    let [i, j] = f.owners[0];
                    let q = eval.fluxes[i][j];
                    let t = eval.locals[i].c[(j, j)];
                    let psi =
                        (self.datum - f.reference_elevation_m) + state.relative_face_head_m[fi];
                    if psi + (q - rain) / t >= 0. {
                        fixed[fi] = Some(f.reference_elevation_m - self.datum);
                    } else {
                        targets[fi] = rain;
                    }
                }
                Bc::Internal => {}
            }
        }
        // Solve increments from the ORIGINAL local mass/face residuals.
        // Absolute-head lifting can conceal small physical flows behind the
        // much larger storage*H/dt right-hand side.
        for (fi, h) in fixed.iter_mut().enumerate() {
            if let Some(h) = h {
                *h -= state.relative_face_head_m[fi];
            }
        }
        let mut free = vec![None; self.faces.len()];
        let mut count = 0;
        for (i, value) in fixed.iter().enumerate() {
            if value.is_none() {
                free[i] = Some(count);
                count += 1;
            }
        }
        if count > s.maximum_free_faces {
            return Err(Error::Input("Richards free face budget exceeded".into()));
        }
        let mut condensed = Vec::new();
        for (i, l) in eval.locals.iter().enumerate() {
            let e = self.cells[i].element_index;
            let capacity = dt.map_or(0., |d| l.capacity / d);
            let rate = dt.map_or(0., |d| {
                (l.storage - old.cell_storage_offset_m3[e].unwrap_or(0.)) / d
            });
            let v = l.b - l.c * (l.mass_derivative * eval.fluxes[i]);
            let a = v.sum() + capacity;
            if !a.is_finite() || a.abs() <= 1e-14 * (v.abs().sum() + capacity) {
                return Err(Error::Solve(
                    "singular hydraulic cell Newton condensation".into(),
                ));
            }
            let rhs = -(eval.fluxes[i].sum() + rate - l.source);
            condensed.push((l.c - v * l.b.transpose() / a, v * (rhs / a), a, rhs));
        }
        let mut candidate = state.clone();
        let mut linear = None;
        if count > 0 {
            let pattern = CsrPattern::from_connectivity(
                count,
                self.cell_faces
                    .iter()
                    .map(|faces| faces.iter().filter_map(|&f| free[f]).collect()),
            )?;
            let mut matrix = CsrMatrix::zeros(pattern);
            let mut rhs = vec![0.; count];
            for (i, (d, load, _, _)) in condensed.iter().enumerate() {
                for j in 0..4 {
                    let fj = self.cell_faces[i][j];
                    if let Some(row) = free[fj] {
                        rhs[row] += load[j];
                        for k in 0..4 {
                            let fk = self.cell_faces[i][k];
                            if let Some(col) = free[fk] {
                                matrix.add(row, col, d[(j, k)])?;
                            } else {
                                rhs[row] -= d[(j, k)] * fixed[fk].unwrap_or(0.);
                            }
                        }
                    }
                }
            }
            for (fi, row) in free.iter().enumerate() {
                if let Some(row) = row {
                    rhs[*row] += self.faces[fi]
                        .owners
                        .iter()
                        .map(|&[i, j]| eval.fluxes[i][j])
                        .sum::<f64>()
                        - targets[fi];
                }
            }
            let settings = darcy::Settings {
                relative_linear_tolerance: 1e-11,
                pressure_head_tolerance_m: 0.,
                maximum_free_heads: s.maximum_free_faces,
                linear_solver: s.linear_solver.clone(),
            };
            let (solution, report) = darcy::linear_solve(&matrix, &rhs, &settings)?;
            linear = Some(report);
            for (fi, row) in free.iter().enumerate() {
                if let Some(row) = row {
                    candidate.relative_face_head_m[fi] += solution[*row];
                }
            }
        }
        for (fi, h) in fixed.into_iter().enumerate() {
            if let Some(h) = h {
                candidate.relative_face_head_m[fi] += h;
            }
        }
        for (i, (_, _, a, rhs)) in condensed.iter().enumerate() {
            let lambda = Vector4::from_fn(|j, _| {
                let f = self.cell_faces[i][j];
                candidate.relative_face_head_m[f] - state.relative_face_head_m[f]
            });
            candidate.relative_cell_head_m[self.cells[i].element_index] = Some(
                state.relative_cell_head_m[self.cells[i].element_index].unwrap_or(0.)
                    + (eval.locals[i].b.dot(&lambda) + rhs) / a,
            );
        }
        Ok((candidate, linear))
    }
    // Globalization uses a head-valued L2 residual with transmissivities
    // frozen at the current Newton state. The local physical tolerances in
    // evaluate remain the acceptance gate. Recomputing flow-relative scales
    // at each trial can make a Newton descent direction look worse solely
    // because a small trial flux shrinks its denominator.
    fn search_norm(
        &self,
        state: &State,
        old: &State,
        dt: Option<f64>,
        bc: &[Bc],
        value: &Evaluated,
        scales: &Evaluated,
    ) -> f64 {
        let mut norm: f64 = 0.;
        for (i, cell) in self.cells.iter().enumerate() {
            let e = cell.element_index;
            let rate = dt.map_or(0., |d| {
                (value.locals[i].storage - old.cell_storage_offset_m3[e].unwrap_or(0.)) / d
            });
            let residual = value.fluxes[i].sum() + rate - value.locals[i].source;
            let scale = scales.locals[i].b.sum() + dt.map_or(0., |d| scales.locals[i].capacity / d);
            norm = norm.hypot(residual / scale);
        }
        for (fi, face) in self.faces.iter().enumerate() {
            let q: f64 = face.owners.iter().map(|&[i, j]| value.fluxes[i][j]).sum();
            let t: f64 = face
                .owners
                .iter()
                .map(|&[i, j]| scales.locals[i].c[(j, j)])
                .sum();
            let r = match bc[fi] {
                Bc::Internal => q / t,
                Bc::Flux(target) => (q - target) / t,
                Bc::Head(target) => state.relative_face_head_m[fi] - target,
                Bc::Atmosphere(rain) => {
                    let psi =
                        self.datum - face.reference_elevation_m + state.relative_face_head_m[fi];
                    (-psi).min((q - rain) / t)
                }
            };
            norm = norm.hypot(r);
        }
        norm
    }
    fn trial(
        &self,
        old: &State,
        dt: Option<f64>,
        bc: &[Bc],
        s: &Settings,
        log: &mut Vec<Iteration>,
    ) -> Result<State> {
        let mut state = old.clone();
        for (fi, b) in bc.iter().enumerate() {
            if let Bc::Head(h) = b {
                state.relative_face_head_m[fi] = *h;
            }
        }
        let mut eval = self.evaluate(&state, old, dt, bc, s)?;
        for iteration in 0..=s.maximum_iterations {
            if eval.balance.maximum_scaled_residual <= 1. {
                self.pressure_reference(&state, &eval, bc, s, dt.is_some())?;
                for (i, c) in self.cells.iter().enumerate() {
                    state.cell_outflow_m3_s[c.element_index] = Some(eval.fluxes[i].into());
                    state.cell_storage_offset_m3[c.element_index] = Some(eval.locals[i].storage);
                }
                if let Some(dt) = dt {
                    state.time_s = old.time_s + dt;
                }
                log.push(Iteration {
                    iteration,
                    scaled_residual: eval.balance.maximum_scaled_residual,
                    line_search_fraction: 0.,
                    linear: None,
                });
                return Ok(state);
            }
            if iteration == s.maximum_iterations {
                break;
            }
            let (candidate, linear) = self.newton_step(&state, old, dt, bc, &eval, s)?;
            let old_norm = self.search_norm(&state, old, dt, bc, &eval, &eval);
            let mut accepted = None;
            for half in 0..14 {
                let fraction = 2f64.powi(-half);
                let mut test = state.clone();
                for c in &self.cells {
                    let e = c.element_index;
                    let h = state.relative_cell_head_m[e].unwrap_or(0.);
                    test.relative_cell_head_m[e] =
                        Some(h + fraction * (candidate.relative_cell_head_m[e].unwrap_or(0.) - h));
                }
                for (i, h) in test.relative_face_head_m.iter_mut().enumerate() {
                    *h += fraction * (candidate.relative_face_head_m[i] - *h);
                }
                if let Ok(next) = self.evaluate(&test, old, dt, bc, s)
                    && (next.balance.maximum_scaled_residual <= 1.
                        || self.search_norm(&test, old, dt, bc, &next, &eval)
                            < old_norm * (1. - 1e-4 * fraction))
                {
                    accepted = Some((test, next, fraction));
                    break;
                }
            }
            let Some((next, next_eval, fraction)) = accepted else {
                log.push(Iteration {
                    iteration,
                    scaled_residual: eval.balance.maximum_scaled_residual,
                    line_search_fraction: 0.,
                    linear,
                });
                return Err(Error::Solve(format!(
                    "Richards line search failed at residual {}; cell mass {} m3/s, face flux {} m3/s",
                    eval.balance.maximum_scaled_residual,
                    eval.balance.maximum_absolute_cell_mass_residual_m3_s,
                    eval.balance.maximum_absolute_face_flux_residual_m3_s
                )));
            };
            log.push(Iteration {
                iteration,
                scaled_residual: next_eval.balance.maximum_scaled_residual,
                line_search_fraction: fraction,
                linear,
            });
            state = next;
            eval = next_eval;
        }
        Err(Error::Solve(format!(
            "Richards iteration limit; scaled residual {}",
            eval.balance.maximum_scaled_residual
        )))
    }
    fn pressure_reference(
        &self,
        state: &State,
        eval: &Evaluated,
        bc: &[Bc],
        s: &Settings,
        transient: bool,
    ) -> Result<()> {
        let mut seen = vec![false; self.cells.len()];
        for start in 0..self.cells.len() {
            if !seen[start] {
                let mut stack = vec![start];
                seen[start] = true;
                let mut anchored = false;
                while let Some(i) = stack.pop() {
                    anchored |= transient && eval.locals[i].capacity > 0.;
                    for &fi in &self.cell_faces[i] {
                        anchored |= matches!(bc[fi], Bc::Head(_));
                        if let Bc::Atmosphere(rain) = bc[fi] {
                            let [ci, j] = self.faces[fi].owners[0];
                            let q = eval.fluxes[ci][j];
                            let psi = self.datum - self.faces[fi].reference_elevation_m
                                + state.relative_face_head_m[fi];
                            anchored |= psi.abs() <= s.head_tolerance_m
                                && q - rain
                                    > s.absolute_flow_tolerance_m3_s
                                        + s.relative_flow_tolerance * (q.abs() + rain.abs());
                        }
                        for &[j, _] in &self.faces[fi].owners {
                            if !seen[j] {
                                seen[j] = true;
                                stack.push(j);
                            }
                        }
                    }
                }
                if !anchored {
                    return Err(Error::Solve(if transient {
                        "unresolved pressure reference in zero-storage transient hydraulic component: prescribe head or physical storage/strict active atmosphere".into()
                    } else {
                        "non-unique steady hydraulic component: no prescribed head or strictly active atmospheric reference".into()
                    }));
                }
            }
        }
        Ok(())
    }
    fn frame(
        &self,
        name: String,
        state: State,
        old: &State,
        dt: Option<f64>,
        bc: &[Bc],
        s: &Settings,
    ) -> Result<Frame> {
        let eval = self.evaluate(&state, old, dt, bc, s)?;
        let mut points = Vec::new();
        for (i, c) in self
            .cells
            .iter()
            .enumerate()
            .filter(|_| s.store_integration_points)
        {
            let h = state.relative_cell_head_m[c.element_index].unwrap_or(0.);
            for (j, q) in c.points.iter().enumerate() {
                let psi = (self.datum - q.position_m[2]) + h;
                let flux = (0..4)
                    .map(|k| q.basis_per_m2[k] * eval.fluxes[i][k])
                    .fold(nalgebra::Vector3::zeros(), |a, b| a + b);
                points.push(Point {
                    element_index: c.element_index,
                    local_index: j,
                    position_m: q.position_m,
                    weight_m3: q.weight_m3,
                    total_head_m: self.datum + h,
                    pressure_head_m: psi,
                    retention: self.materials[i].retention.evaluate(psi)?,
                    darcy_flux_m_s: flux.into(),
                });
            }
        }
        Ok(Frame {
            name,
            state,
            balance: eval.balance,
            points,
        })
    }
}

pub fn solve(
    mesh: &Mesh,
    case: Case,
    active: &[bool],
    material_ids: Option<&[u32]>,
) -> Result<Report> {
    advance(mesh, case, active, material_ids, None, None)
}
/// Crate-internal physical-time trial service. The caller owns commit/cutback
/// of both hydraulic and mechanical states. It cannot silently advance a
/// hydraulic substep after mechanics has rejected the same physical interval.
type BoundaryPeriods = Vec<(String, Option<f64>, Vec<Bc>)>;
pub(crate) struct Stepper {
    model: Arc<Model>,
    boundaries: BoundaryPeriods,
    settings: Settings,
    initial: Option<Arc<Frame>>,
}
impl Stepper {
    pub fn new(mesh: &Mesh, case: &Case) -> Result<Self> {
        Self::initial_on_domain(mesh, case, &vec![true; mesh.tetrahedra.len()])
    }
    pub fn initial_on_domain(mesh: &Mesh, case: &Case, active: &[bool]) -> Result<Self> {
        let model = Model::new(mesh, case, active, None)?;
        let boundaries = validated_boundaries(&model, mesh, case)?;
        if !matches!(case.analysis, Analysis::Transient { .. }) {
            return Err(Error::Input(
                "ONE_WAY physical time requires transient Richards periods".into(),
            ));
        }
        let state = model.initial(case)?;
        let initial = model.frame(
            "initial".into(),
            state.clone(),
            &state,
            None,
            &boundaries[0].2,
            &case.settings,
        )?;
        Ok(Self {
            model: Arc::new(model),
            boundaries,
            settings: case.settings.clone(),
            initial: Some(Arc::new(initial)),
        })
    }
    pub fn initial(&self) -> Result<Frame> {
        self.initial
            .as_ref()
            .map(|f| f.as_ref().clone())
            .ok_or_else(|| {
                Error::Input("deferred hydraulic domain has no imposed initial field".into())
            })
    }
    /// Validate only the period applicable to this construction domain. Earlier
    /// and later periods may prescribe faces that are not exposed here.
    pub fn for_domain(mesh: &Mesh, case: &Case, active: &[bool], period: usize) -> Result<Self> {
        Self::for_domain_materials(mesh, case, active, None, period)
    }
    pub fn for_domain_materials(
        mesh: &Mesh,
        case: &Case,
        active: &[bool],
        ids: Option<&[u32]>,
        period: usize,
    ) -> Result<Self> {
        let Analysis::Transient { periods } = &case.analysis else {
            return Err(Error::Input("construction needs transient flow".into()));
        };
        let model = Model::new(mesh, case, active, ids)?;
        let boundaries = one_period_boundaries(&model, mesh, periods, period)?;
        let state = model.initial(case)?;
        let initial = model.frame(
            "domain_initial_guess".into(),
            state.clone(),
            &state,
            None,
            &boundaries[period].2,
            &case.settings,
        )?;
        Ok(Self {
            model: Arc::new(model),
            boundaries,
            settings: case.settings.clone(),
            initial: Some(Arc::new(initial)),
        })
    }
    pub fn deferred_domain(
        mesh: &Mesh,
        case: &Case,
        active: &[bool],
        ids: Option<&[u32]>,
        period: usize,
    ) -> Result<Self> {
        let Analysis::Transient { periods } = &case.analysis else {
            return Err(Error::Input("construction needs transient flow".into()));
        };
        let model = Model::new(mesh, case, active, ids)?;
        let boundaries = one_period_boundaries(&model, mesh, periods, period)?;
        Ok(Self {
            model: Arc::new(model),
            boundaries,
            settings: case.settings.clone(),
            initial: None,
        })
    }
    /// Immutable domain and its specified initial field are shared across
    /// periods. Only this period's boundary operator is prepared and stored.
    pub fn with_period_boundaries(
        &self,
        mesh: &Mesh,
        periods: &[Period],
        period: usize,
    ) -> Result<Self> {
        Ok(Self {
            model: Arc::clone(&self.model),
            boundaries: one_period_boundaries(&self.model, mesh, periods, period)?,
            settings: self.settings.clone(),
            initial: self.initial.clone(),
        })
    }
    pub fn faces(&self) -> Vec<Face> {
        self.model.faces.clone()
    }
    pub fn trial(&self, old: &State, period: usize, dt: f64) -> (Result<Frame>, Vec<Iteration>) {
        let mut log = Vec::new();
        let result = (|| {
            let (name, end, bc) = self
                .boundaries
                .get(period)
                .ok_or_else(|| Error::Input("invalid hydraulic period".into()))?;
            if bc.len() != self.model.faces.len() {
                return Err(Error::Input("unprepared hydraulic domain period".into()));
            }
            let start = if period == 0 {
                0.
            } else {
                self.boundaries[period - 1].1.unwrap_or(0.)
            };
            if !dt.is_finite()
                || dt <= 0.
                || old.time_s < start
                || old.time_s + dt == old.time_s
                || end.is_none_or(|end| old.time_s + dt > end)
            {
                return Err(Error::Input(
                    "hydraulic trial crosses physical period boundary".into(),
                ));
            }
            let next = self
                .model
                .trial(old, Some(dt), bc, &self.settings, &mut log)?;
            self.model
                .frame(name.clone(), next, old, Some(dt), bc, &self.settings)
        })();
        (result, log)
    }
}

fn one_period_boundaries(
    model: &Model,
    mesh: &Mesh,
    periods: &[Period],
    period: usize,
) -> Result<BoundaryPeriods> {
    if period >= periods.len() {
        return Err(Error::Input("invalid construction period".into()));
    }
    let mut end = 0.;
    let mut names = BTreeSet::new();
    let mut out = vec![];
    for (i, p) in periods.iter().enumerate() {
        if p.name.is_empty()
            || !names.insert(p.name.clone())
            || !p.end_time_s.is_finite()
            || p.end_time_s <= end
        {
            return Err(Error::Input(
                "Richards period names/times must be unique and increasing".into(),
            ));
        }
        end = p.end_time_s;
        out.push((
            p.name.clone(),
            Some(end),
            if i == period {
                model.boundaries(mesh, &p.boundaries)?
            } else {
                vec![]
            },
        ));
    }
    Ok(out)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RemovedCell {
    pub element_index: usize,
    pub relative_head_m: f64,
    pub outflow_m3_s: [f64; 4],
    pub storage_offset_m3: f64,
    pub volume_m3: f64,
    pub water_equivalent_volume_m3: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DomainRemoval {
    pub state: State,
    pub old_faces: Vec<Face>,
    pub new_faces: Vec<Face>,
    /// New face index -> original face index. Corner keys, never coordinates.
    pub new_to_old_face: Vec<usize>,
    pub newly_exposed_corners: Vec<[usize; 3]>,
    pub removed_cells: Vec<RemovedCell>,
    pub water_before_m3: f64,
    pub water_remaining_m3: f64,
    pub water_removed_m3: f64,
    pub absolute_mass_transfer_error_m3: f64,
}

/// Instantaneous DEATH transfer, with no flow solve or hydraulic equilibration.
/// Remaining primary H/W and local fluxes are copied exactly. New boundary
/// traces inherit the old interior trace; their new BC applies to the NEXT dt.
pub fn remove_domain(
    mesh: &Mesh,
    case: &Case,
    before: &[bool],
    after: &[bool],
    old: &State,
) -> Result<DomainRemoval> {
    remove_domain_materials(mesh, case, before, after, None, old)
}
pub(crate) fn remove_domain_materials(
    mesh: &Mesh,
    case: &Case,
    before: &[bool],
    after: &[bool],
    ids: Option<&[u32]>,
    old: &State,
) -> Result<DomainRemoval> {
    if before.len() != after.len()
        || before.len() != mesh.tetrahedra.len()
        || before.iter().zip(after).any(|(a, b)| !a && *b)
        || !before.iter().zip(after).any(|(a, b)| *a && !b)
    {
        return Err(Error::Input(
            "hydraulic DEATH requires a strict remaining subdomain".into(),
        ));
    }
    let original = Model::new(mesh, case, before, ids)?;
    let remaining = Model::new(mesh, case, after, ids)?;
    let n = before.len();
    if !old.time_s.is_finite()
        || old.time_s < 0.
        || old.relative_cell_head_m.len() != n
        || old.cell_outflow_m3_s.len() != n
        || old.cell_storage_offset_m3.len() != n
        || old.relative_face_head_m.len() != original.faces.len()
        || !old.relative_face_head_m.iter().all(|v| v.is_finite())
    {
        return Err(Error::Input("invalid hydraulic origin for DEATH".into()));
    }
    let mut state = old.clone();
    let mut removed_cells = vec![];
    let mut water_before = 0.;
    let mut water_remaining = 0.;
    let mut water_removed = 0.;
    for (e, &live) in before.iter().enumerate() {
        if !live
            && (old.relative_cell_head_m[e].is_some()
                || old.cell_storage_offset_m3[e].is_some()
                || old.cell_outflow_m3_s[e].is_some())
        {
            return Err(Error::Input(
                "inactive hydraulic origin contains live cell data".into(),
            ));
        }
    }
    for (i, c) in original.cells.iter().enumerate() {
        let e = c.element_index;
        let h = old.relative_cell_head_m[e]
            .filter(|v| v.is_finite())
            .ok_or_else(|| Error::Input("missing hydraulic origin head".into()))?;
        let w = old.cell_storage_offset_m3[e]
            .filter(|v| v.is_finite())
            .ok_or_else(|| Error::Input("missing hydraulic origin storage".into()))?;
        let q = old.cell_outflow_m3_s[e]
            .filter(|v| v.iter().all(|v| v.is_finite()))
            .ok_or_else(|| Error::Input("missing hydraulic origin flux".into()))?;
        let expected = original.local(i, h)?.storage;
        if (w - expected).abs() > 1e-11 * c.volume_m3.max(w.abs()).max(1e-30) {
            return Err(Error::Input("hydraulic origin H/storage mismatch".into()));
        }
        let water = original.materials[i].retention.saturated_water_content * c.volume_m3 + w;
        if !water.is_finite() || water < 0. {
            return Err(Error::Input("invalid hydraulic origin water volume".into()));
        }
        water_before += water;
        if after[e] {
            water_remaining += water;
        } else {
            water_removed += water;
            removed_cells.push(RemovedCell {
                element_index: e,
                relative_head_m: h,
                outflow_m3_s: q,
                storage_offset_m3: w,
                volume_m3: c.volume_m3,
                water_equivalent_volume_m3: water,
            });
            state.relative_cell_head_m[e] = None;
            state.cell_storage_offset_m3[e] = None;
            state.cell_outflow_m3_s[e] = None;
        }
    }
    let old_map: BTreeMap<_, _> = original
        .faces
        .iter()
        .enumerate()
        .map(|(i, f)| (f.corners, i))
        .collect();
    let mut mapping = vec![];
    let mut newly_exposed = vec![];
    for f in &remaining.faces {
        let i = *old_map
            .get(&f.corners)
            .ok_or_else(|| Error::Input("DEATH created an unknown face".into()))?;
        mapping.push(i);
        if f.owners.len() == 1 && original.faces[i].owners.len() == 2 {
            newly_exposed.push(f.corners);
        }
    }
    state.relative_face_head_m = mapping
        .iter()
        .map(|&i| old.relative_face_head_m[i])
        .collect();
    let error = (water_before - water_remaining - water_removed).abs();
    if error > 1e-12 * water_before.max(1e-30) {
        return Err(Error::Solve("DEATH water mass transfer failed".into()));
    }
    Ok(DomainRemoval {
        state,
        old_faces: original.faces,
        new_faces: remaining.faces,
        new_to_old_face: mapping,
        newly_exposed_corners: newly_exposed,
        removed_cells,
        water_before_m3: water_before,
        water_remaining_m3: water_remaining,
        water_removed_m3: water_removed,
        absolute_mass_transfer_error_m3: error,
    })
}

fn storage_only(model: &Model, i: usize, h: f64) -> Result<f64> {
    let c = &model.cells[i];
    let law = &model.materials[i].retention;
    if !h.is_finite() || !(h + model.datum).is_finite() {
        return Err(Error::Input("non-finite transition head".into()));
    }
    c.points.iter().try_fold(0., |sum, q| {
        Ok(sum
            + q.weight_m3
                * law
                    .evaluate((model.datum - q.position_m[2]) + h)?
                    .storage_offset)
    })
}
fn projected_head(cell: &mixed_flow::Cell, value: &AffineHead, datum: f64) -> Result<f64> {
    value.validate()?;
    let anchor = cell.points[0].position_m;
    Ok(value.relative(anchor, datum)
        + cell
            .points
            .iter()
            .map(|q| {
                q.weight_m3
                    * (0..3)
                        .map(|a| value.gradient[a] * (q.position_m[a] - anchor[a]))
                        .sum::<f64>()
            })
            .sum::<f64>()
            / cell.volume_m3)
}
fn storage_head(model: &Model, i: usize, target: f64, reference: &AffineHead) -> Result<f64> {
    let cell = &model.cells[i];
    let law = &model.materials[i].retention;
    let v = cell.volume_m3;
    let tolerance = 5e-14 * v + 16. * f64::EPSILON * target.abs();
    let anchor = projected_head(cell, reference, model.datum)?;
    if !target.is_finite()
        || target <= (law.residual_water_content - law.saturated_water_content) * v
    {
        return Err(Error::Input(
            "water-preserving MATSET has no finite retention head".into(),
        ));
    }
    if law.saturated_specific_storage_per_m == 0. && target >= -tolerance {
        if target > tolerance || (storage_only(model, i, anchor)? - target).abs() > tolerance {
            return Err(Error::Input(
                "saturated zero-storage MATSET needs a compatible explicit pressure reference"
                    .into(),
            ));
        }
        return Ok(anchor);
    }
    let mut lo = anchor;
    let mut hi = anchor;
    let mut step = 1.;
    for _ in 0..100 {
        if storage_only(model, i, lo)? <= target && storage_only(model, i, hi)? >= target {
            break;
        }
        if storage_only(model, i, lo)? > target {
            lo -= step;
        }
        if storage_only(model, i, hi)? < target {
            hi += step;
        }
        step *= 2.;
    }
    if storage_only(model, i, lo)? > target || storage_only(model, i, hi)? < target {
        return Err(Error::Input(
            "water-preserving MATSET head could not be bracketed".into(),
        ));
    }
    for _ in 0..200 {
        let h = lo + (hi - lo) / 2.;
        let w = storage_only(model, i, h)?;
        if (w - target).abs() <= tolerance {
            return Ok(h);
        }
        if w < target {
            lo = h;
        } else {
            hi = h;
        }
    }
    Err(Error::Solve(
        "water-preserving MATSET inversion exhausted".into(),
    ))
}
pub(crate) fn transition_domain(
    mesh: &Mesh,
    case: &Case,
    before: crate::hydraulic_transition::Domain<'_>,
    after: crate::hydraulic_transition::Domain<'_>,
    old: &State,
    selected: &[usize],
    policy: &crate::hydraulic_transition::WaterPolicy,
) -> Result<crate::hydraulic_transition::WaterChange> {
    use crate::hydraulic_transition::{ChangedCell, WaterChange, WaterPolicy};
    let n = mesh.tetrahedra.len();
    let chosen = selected.iter().copied().collect::<BTreeSet<_>>();
    if before.active.len() != n
        || after.active.len() != n
        || before.material_ids.len() != n
        || after.material_ids.len() != n
        || chosen.is_empty()
        || chosen.len() != selected.len()
        || chosen.iter().any(|e| *e >= n)
        || (0..n).any(|e| {
            (before.active[e] && !after.active[e])
                || (chosen.contains(&e)
                    && !after.active[e]
                    && !matches!(policy, WaterPolicy::Inactive))
                || (!chosen.contains(&e)
                    && (before.active[e] != after.active[e]
                        || before.material_ids[e] != after.material_ids[e]))
        })
    {
        return Err(Error::Input(
            "invalid BIRTH/MATSET hydraulic domain or selection".into(),
        ));
    }
    if matches!(policy, WaterPolicy::Inactive)
        && chosen.iter().any(|&e| before.active[e] || after.active[e])
    {
        return Err(Error::Input(
            "inactive water policy requires wholly inactive selected soil".into(),
        ));
    }
    if (0..n).any(|e| !before.active[e] && after.active[e])
        && !matches!(policy, WaterPolicy::InitializeHead { .. })
    {
        return Err(Error::Input("BIRTH requires explicit initial head".into()));
    }
    match policy {
        WaterPolicy::InitializeHead { value } => value.validate()?,
        WaterPolicy::PreserveWater { reference_head } => reference_head.validate()?,
        _ => {}
    }
    let original = Model::new(mesh, case, before.active, Some(before.material_ids))?;
    let new = Model::new(mesh, case, after.active, Some(after.material_ids))?;
    if !old.time_s.is_finite()
        || old.time_s < 0.
        || old.relative_cell_head_m.len() != n
        || old.cell_storage_offset_m3.len() != n
        || old.cell_outflow_m3_s.len() != n
        || old.relative_face_head_m.len() != original.faces.len()
        || old.relative_face_head_m.iter().any(|v| !v.is_finite())
    {
        return Err(Error::Input("invalid hydraulic transition origin".into()));
    }
    let mut snapshots = BTreeMap::new();
    let mut water_before = 0.;
    for e in 0..n {
        if !before.active[e]
            && (old.relative_cell_head_m[e].is_some()
                || old.cell_storage_offset_m3[e].is_some()
                || old.cell_outflow_m3_s[e].is_some())
        {
            return Err(Error::Input(
                "inactive transition origin has live water data".into(),
            ));
        }
    }
    for (i, c) in original.cells.iter().enumerate() {
        let e = c.element_index;
        let h = old.relative_cell_head_m[e]
            .filter(|h| h.is_finite())
            .ok_or_else(|| Error::Input("missing transition origin head".into()))?;
        let w = old.cell_storage_offset_m3[e]
            .filter(|w| w.is_finite())
            .ok_or_else(|| Error::Input("missing transition origin storage".into()))?;
        let q = old.cell_outflow_m3_s[e]
            .filter(|q| q.iter().all(|v| v.is_finite()))
            .ok_or_else(|| Error::Input("missing transition origin flux".into()))?;
        if (w - storage_only(&original, i, h)?).abs() > 1e-11 * c.volume_m3.max(w.abs()) {
            return Err(Error::Input("transition origin H/storage mismatch".into()));
        }
        let water = original.materials[i].retention.saturated_water_content * c.volume_m3 + w;
        water_before += water;
        snapshots.insert(
            e,
            RemovedCell {
                element_index: e,
                relative_head_m: h,
                outflow_m3_s: q,
                storage_offset_m3: w,
                volume_m3: c.volume_m3,
                water_equivalent_volume_m3: water,
            },
        );
    }
    let map = original
        .faces
        .iter()
        .enumerate()
        .map(|(i, f)| (f.corners, i))
        .collect::<BTreeMap<_, _>>();
    let mapping = new
        .faces
        .iter()
        .map(|f| map.get(&f.corners).copied())
        .collect::<Vec<_>>();
    let mut state = old.clone();
    state.relative_face_head_m = new
        .faces
        .iter()
        .zip(&mapping)
        .map(|(f, i)| {
            if let Some(i) = i {
                Ok(old.relative_face_head_m[*i])
            } else if let WaterPolicy::InitializeHead { value } = policy {
                Ok(value.relative(f.reference_position_m, new.datum))
            } else {
                Err(Error::Input(
                    "new hydraulic face requires initial head".into(),
                ))
            }
        })
        .collect::<Result<Vec<_>>>()?;
    let mut changed = vec![];
    let mut added = 0.;
    let mut removed = 0.;
    let mut water_after = 0.;
    for (i, c) in new.cells.iter().enumerate() {
        let e = c.element_index;
        if chosen.contains(&e) {
            let prior = snapshots.get(&e);
            let old_head = prior.map(|p| p.relative_head_m);
            let h = match policy {
                WaterPolicy::Inactive => {
                    return Err(Error::Input(
                        "inactive hydraulic transition selected live cell".into(),
                    ));
                }
                WaterPolicy::InitializeHead { value } => projected_head(c, value, new.datum)?,
                WaterPolicy::PreserveHead => {
                    old_head.ok_or_else(|| Error::Input("missing MATSET head".into()))?
                }
                WaterPolicy::PreserveWater { reference_head } => {
                    let p = prior.ok_or_else(|| Error::Input("missing MATSET water".into()))?;
                    let old_law = &case
                        .materials
                        .iter()
                        .find(|m| m.id == before.material_ids[e])
                        .ok_or_else(|| Error::Input("missing old hydraulic material".into()))?
                        .retention;
                    let target = p.storage_offset_m3
                        + (old_law.saturated_water_content
                            - new.materials[i].retention.saturated_water_content)
                            * c.volume_m3;
                    storage_head(&new, i, target, reference_head)?
                }
            };
            let local = new.local(i, h)?;
            let q: [f64; 4] = (local.c
                * Vector4::from_fn(|j, _| h - state.relative_face_head_m[new.cell_faces[i][j]]))
            .into();
            let water =
                new.materials[i].retention.saturated_water_content * c.volume_m3 + local.storage;
            if !water.is_finite() || water < 0. {
                return Err(Error::Input("invalid changed water mass".into()));
            }
            let delta = if let Some(p) = prior {
                water - p.water_equivalent_volume_m3
            } else {
                water
            };
            if !matches!(policy, WaterPolicy::PreserveWater { .. }) {
                added += delta.max(0.);
                removed += (-delta).max(0.);
            }
            state.relative_cell_head_m[e] = Some(h);
            state.cell_storage_offset_m3[e] = Some(local.storage);
            state.cell_outflow_m3_s[e] = Some(q);
            changed.push(ChangedCell {
                element_index: e,
                old_material_id: prior.map(|_| before.material_ids[e]),
                new_material_id: after.material_ids[e],
                before: prior.cloned(),
                after: RemovedCell {
                    element_index: e,
                    relative_head_m: h,
                    outflow_m3_s: q,
                    storage_offset_m3: local.storage,
                    volume_m3: c.volume_m3,
                    water_equivalent_volume_m3: water,
                },
            });
            water_after += water;
        } else {
            water_after += snapshots[&e].water_equivalent_volume_m3;
        }
    }
    let error = (water_after - water_before - added + removed).abs();
    if error > 1e-11 * water_before.max(water_after).max(1e-30) {
        return Err(Error::Solve(
            "BIRTH/MATSET water mass transfer failed".into(),
        ));
    }
    Ok(WaterChange {
        state,
        old_faces: original.faces,
        new_faces: new.faces,
        new_to_old_face: mapping,
        changed_cells: changed,
        water_before_m3: water_before,
        water_after_m3: water_after,
        water_added_m3: added,
        water_removed_m3: removed,
        absolute_mass_transfer_error_m3: error,
    })
}

fn validated_boundaries(model: &Model, mesh: &Mesh, case: &Case) -> Result<BoundaryPeriods> {
    match &case.analysis {
        Analysis::Steady { boundaries } => Ok(vec![(
            "steady".into(),
            None,
            model.boundaries(mesh, boundaries)?,
        )]),
        Analysis::Transient { periods } => {
            if periods.is_empty() {
                return Err(Error::Input("Richards transient requires periods".into()));
            }
            let mut end = 0.;
            let mut names = BTreeSet::new();
            let mut list = Vec::new();
            for p in periods {
                if p.name.is_empty()
                    || !names.insert(p.name.clone())
                    || !p.end_time_s.is_finite()
                    || p.end_time_s <= end
                {
                    return Err(Error::Input(
                        "Richards period names/times must be unique and increasing".into(),
                    ));
                }
                end = p.end_time_s;
                list.push((
                    p.name.clone(),
                    Some(end),
                    model.boundaries(mesh, &p.boundaries)?,
                ));
            }
            Ok(list)
        }
    }
}
pub fn solve_until(
    mesh: &Mesh,
    case: Case,
    active: &[bool],
    material_ids: Option<&[u32]>,
    accepted_steps: usize,
) -> Result<Report> {
    advance(mesh, case, active, material_ids, None, Some(accepted_steps))
}
/// Exact input/binary/history replay is an intentionally conservative resume
/// validator. Its cost grows with the saved prefix; it is not a fast restart.
pub fn resume(
    mesh: &Mesh,
    case: Case,
    active: &[bool],
    material_ids: Option<&[u32]>,
    saved: Report,
    accepted_steps: Option<usize>,
) -> Result<Report> {
    if saved.complete || saved.failure_reason.is_some() || saved.runtime.is_none() {
        return Err(Error::Input(
            "Richards resume requires an accepted intermediate checkpoint".into(),
        ));
    }
    let count = saved.attempts.iter().filter(|a| a.accepted).count();
    if count == 0 {
        return Err(Error::Input(
            "Richards checkpoint has no accepted step".into(),
        ));
    }
    let canonical = solve_until(mesh, case.clone(), active, material_ids, count)?;
    if normalized(&canonical)? != normalized(&saved)? {
        return Err(Error::Input(
            "Richards checkpoint inputs, executable, runtime, frames or hydraulic history changed"
                .into(),
        ));
    }
    advance(
        mesh,
        case,
        active,
        material_ids,
        Some(saved),
        accepted_steps,
    )
}
fn normalized(run: &Report) -> Result<Vec<u8>> {
    let mut copy = run.clone();
    for attempt in &mut copy.attempts {
        for iteration in &mut attempt.iterations {
            if let Some(linear) = &mut iteration.linear {
                linear.elapsed_seconds = 0.;
            }
        }
    }
    Ok(serde_json::to_vec(&copy)?)
}
fn advance(
    mesh: &Mesh,
    case: Case,
    active: &[bool],
    material_ids: Option<&[u32]>,
    saved: Option<Report>,
    limit: Option<usize>,
) -> Result<Report> {
    if limit == Some(0) {
        return Err(Error::Input(
            "Richards checkpoint count must be positive".into(),
        ));
    }
    let model = Model::new(mesh, &case, active, material_ids)?;
    let boundaries = validated_boundaries(&model, mesh, &case)?;
    let (mut state, mut frames, mut attempts, mut index, mut step) = if let Some(saved) = saved {
        let runtime = saved
            .runtime
            .ok_or_else(|| Error::Input("missing Richards runtime".into()))?;
        (
            saved.last_committed,
            saved.frames,
            saved.attempts,
            runtime.period_index,
            runtime.next_time_step_s,
        )
    } else {
        let state = model.initial(&case)?;
        let frame = model.frame(
            "initial".into(),
            state.clone(),
            &state,
            None,
            &boundaries[0].2,
            &case.settings,
        )?;
        (
            state,
            vec![frame],
            Vec::new(),
            0,
            case.settings.initial_time_step_s,
        )
    };
    let mut failure = None;
    let mut accepted = 0;
    while index < boundaries.len() {
        let (name, end, bc) = &boundaries[index];
        if attempts.len() >= case.settings.maximum_steps {
            failure = Some("Richards step/attempt budget exceeded".into());
            break;
        }
        let dt = end.map(|end| step.min(end - state.time_s));
        if dt.is_some_and(|dt| dt <= 0. || state.time_s + dt == state.time_s) {
            failure = Some("Richards physical time cannot advance".into());
            break;
        }
        let mut iterations = Vec::new();
        let trial = model.trial(&state, dt, bc, &case.settings, &mut iterations);
        let success = trial.is_ok();
        let reason = trial.as_ref().err().map(ToString::to_string);
        let invalid = matches!(&trial, Err(Error::Input(_)));
        attempts.push(Attempt {
            period: name.clone(),
            start_time_s: state.time_s,
            time_step_s: dt,
            accepted: success,
            iterations,
            failure_reason: reason.clone(),
        });
        match trial {
            Ok(next) => {
                frames.push(model.frame(
                    name.clone(),
                    next.clone(),
                    &state,
                    dt,
                    bc,
                    &case.settings,
                )?);
                state = next;
                accepted += 1;
                if end.is_none() || end.is_some_and(|end| state.time_s >= end) {
                    index += 1;
                } else {
                    step = (step * 1.5).min(case.settings.maximum_time_step_s);
                }
                if limit == Some(accepted) {
                    break;
                }
            }
            Err(_) => {
                let smaller = dt.map(|dt| dt / 2.);
                if !invalid && smaller.is_some_and(|v| v >= case.settings.minimum_time_step_s) {
                    step = smaller.unwrap_or(step);
                } else {
                    failure = reason;
                    break;
                }
            }
        }
    }
    let complete = failure.is_none() && index == boundaries.len();
    let runtime = if !complete && failure.is_none() {
        Some(Runtime {
            period_index: index,
            next_time_step_s: step,
        })
    } else {
        None
    };
    Ok(Report {
        schema: "soil-fem-richards-result/1".into(),
        complete,
        failure_reason: failure,
        mesh_sha256: checkpoint::hash_json(mesh)?,
        case_sha256: checkpoint::hash_json(&case)?,
        solver_binary_sha256: checkpoint::hash_file(&std::env::current_exe()?)?,
        case,
        active_elements: model.active.clone(),
        material_ids: model.ids.clone(),
        cell_element_indices: model.cells.iter().map(|c| c.element_index).collect(),
        faces: model.faces.clone(),
        frames,
        attempts,
        last_committed: state,
        runtime,
    })
}
