//! Structural shell path. Accepted install/load/remove events follow the dry
//! offset composite convention. Each installation owns its reference geometry.
use super::{
    Field, Kernel, Q36, Tying,
    geometry::{Director, Face, Geometry, V, geometries},
};
use crate::{Error, Result, checkpoint, offset_composite::path::Step, plates::Section};
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FaceLoad {
    pub face: String,
    pub z_m: f64,
    pub traction_kpa: [f64; 3],
    pub pressure_kpa: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub nodes: Vec<[f64; 3]>,
    pub faces: Vec<Face>,
    pub directors: Vec<Director>,
    pub section: Section,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub soil: Option<super::soil::Attachment>,
    pub steps: Vec<Step>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub face_loads: BTreeMap<usize, Vec<FaceLoad>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElementState {
    pub face: Face,
    pub soil_owner: Option<crate::offset_composite::Face>,
    pub exposed_side: Option<i8>,
    pub initial_directors: [[f64; 3]; 6],
    pub positions_m: [[f64; 3]; 6],
    pub directors: [[f64; 3]; 6],
    pub explicit: [bool; 6],
    pub director_changes: [[f64; 3]; 6],
    pub fields: Vec<Field>,
    pub resultants: Vec<super::Resultant>,
    pub physical_energy_kn_m: f64,
    pub drilling_energy_kn_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Generation {
    pub generation: u64,
    pub installation: Vec<f64>,
    pub elements: Vec<ElementState>,
    pub plate_cumulative_load_work_kn_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    pub accepted_index: usize,
    pub event: String,
    pub holder: Option<String>,
    pub configuration: Vec<f64>,
    pub active: Option<Generation>,
    pub archived: Vec<Generation>,
    pub soil_stress_kpa: Vec<Vec<[f64; 6]>>,
    pub soil_translations_m: Vec<[f64; 3]>,
    pub applied: Vec<f64>,
    pub plate_load: Vec<f64>,
    pub reactions: Vec<f64>,
    pub cumulative_applied_work_kn_m: f64,
    pub cumulative_reaction_work_kn_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Run {
    pub schema: String,
    pub units: String,
    pub section: Section,
    pub rotation_convention: String,
    pub input_sha256: String,
    pub solver_binary_sha256: String,
    pub frames: Vec<Frame>,
}
fn v(x: &[f64]) -> DVector<f64> {
    DVector::from_column_slice(x)
}
fn list(x: &DVector<f64>) -> Vec<f64> {
    x.as_slice().to_vec()
}
fn indices(face: &Face) -> [usize; 36] {
    std::array::from_fn(|j| face.nodes[j / 6] * 6 + j % 6)
}
fn gather(x: &DVector<f64>, rows: &[usize; 36]) -> Q36 {
    Q36::from_fn(|i, _| x[rows[i]])
}
struct Model {
    elements: Vec<Kernel>,
    geometry: Vec<Geometry>,
    stiffness: DMatrix<f64>,
}
impl Model {
    fn new(case: &Case, g: &[Geometry], n: usize) -> Result<Self> {
        let mut stiffness = DMatrix::zeros(n, n);
        let mut elements = vec![];
        for (geo, f) in g.iter().zip(&case.faces) {
            let k = Kernel::new(geo, &case.section, Tying::default())?;
            let idx = indices(f);
            for (i, &a) in idx.iter().enumerate() {
                for (j, &b) in idx.iter().enumerate() {
                    stiffness[(a, b)] += k.stiffness[(i, j)];
                }
            }
            elements.push(k);
        }
        Ok(Self {
            elements,
            geometry: g.to_vec(),
            stiffness,
        })
    }
    fn loads(&self, case: &Case, traction: [f64; 3]) -> DVector<f64> {
        let mut load = DVector::zeros(self.stiffness.nrows());
        for (k, f) in self.elements.iter().zip(&case.faces) {
            let element = k.traction.transpose() * V::from(traction) + k.gravity;
            for (i, r) in indices(f).into_iter().enumerate() {
                load[r] += element[i];
            }
        }
        load
    }
    fn face_actions(
        &self,
        case: &Case,
        index: usize,
        uniform: [f64; 3],
        load: &mut DVector<f64>,
    ) -> Result<()> {
        let mut owned = BTreeSet::new();
        for entry in case.face_loads.get(&index).into_iter().flatten() {
            let i = case
                .faces
                .iter()
                .position(|f| f.id == entry.face)
                .ok_or_else(|| Error::Input("unknown shell face load owner".into()))?;
            let layer = if entry.z_m == 0. {
                0
            } else {
                entry.z_m.to_bits()
            };
            if !owned.insert((i, layer)) || (layer == 0 && uniform != [0.; 3]) {
                return Err(Error::Input(
                    "duplicate shell face load on the same layer".into(),
                ));
            }
            let f = super::surface_action(
                &self.geometry[i],
                case.section.thickness_m,
                entry.z_m,
                V::from(entry.traction_kpa),
                entry.pressure_kpa,
            )?;
            for (j, row) in indices(&case.faces[i]).into_iter().enumerate() {
                load[row] += f[j];
            }
        }
        Ok(())
    }
    fn recover(&self, g: &mut Generation, q: &DVector<f64>) -> DVector<f64> {
        let delta = q - v(&g.installation);
        let mut internal = DVector::zeros(q.len());
        for (k, e) in self.elements.iter().zip(&mut g.elements) {
            let rows = indices(&e.face);
            let local = gather(&delta, &rows);
            e.fields = k.fields(&local);
            e.resultants = k.resultants(&local, e.face.offset_m);
            e.drilling_energy_kn_m = 0.5 * local.dot(&(k.drill * local));
            e.physical_energy_kn_m = 0.5 * local.dot(&((k.stiffness - k.drill) * local));
            e.director_changes = std::array::from_fn(|i| {
                let theta = V::from_fn(|j, _| local[6 * i + 3 + j]);
                theta.cross(&V::from(e.directors[i])).into()
            });
            let force = k.internal(&local);
            for (i, r) in rows.into_iter().enumerate() {
                internal[r] += force[i];
            }
        }
        internal
    }
}

fn install(case: &Case, old: &Frame) -> Result<(Model, Generation)> {
    let x: Vec<_> = case
        .nodes
        .iter()
        .enumerate()
        .map(|(i, x)| std::array::from_fn(|a| x[a] + old.configuration[6 * i + a]))
        .collect();
    let initial = geometries(&case.nodes, &case.faces, &case.directors)?;
    let explicit: Vec<_> = case
        .directors
        .iter()
        .map(|director| {
            let mut updated = director.clone();
            let j = 6 * updated.node + 3;
            let spin = V::from_column_slice(&old.configuration[j..j + 3]);
            let d = V::from(updated.direction);
            updated.direction = (d + spin.cross(&d)).normalize().into();
            updated
        })
        .collect();
    let g = geometries(&x, &case.faces, &explicit)?;
    let model = Model::new(case, &g, old.configuration.len())?;
    let elements = case
        .faces
        .iter()
        .zip(g)
        .zip(initial)
        .map(|((f, g), initial)| ElementState {
            initial_directors: initial.directors.map(Into::into),
            face: f.clone(),
            soil_owner: case.soil.as_ref().map(|input| input.owners[&f.id].clone()),
            exposed_side: case.soil.as_ref().map(|input| input.exposed_side[&f.id]),
            positions_m: g.x.map(Into::into),
            directors: g.directors.map(Into::into),
            explicit: f.nodes.map(|n| {
                case.directors
                    .iter()
                    .any(|d| d.group == f.group && d.node == n)
            }),
            director_changes: [[0.; 3]; 6],
            fields: vec![],
            resultants: vec![],
            physical_energy_kn_m: 0.,
            drilling_energy_kn_m: 0.,
        })
        .collect();
    Ok((
        model,
        Generation {
            generation: 1 + old.archived.len() as u64,
            installation: old.configuration.clone(),
            elements,
            plate_cumulative_load_work_kn_m: 0.,
        },
    ))
}
fn solve_supported(k: DMatrix<f64>, rhs: &DVector<f64>) -> Result<DVector<f64>> {
    let scale = k.diagonal().map(|v| 1. / v.sqrt());
    let scaled = DMatrix::from_fn(k.nrows(), k.ncols(), |i, j| scale[i] * k[(i, j)] * scale[j]);
    let factor = scaled
        .cholesky()
        .ok_or_else(|| Error::Solve("singular supported curved shell".into()))?;
    let correction = |r: &DVector<f64>| factor.solve(&scale.component_mul(r)).component_mul(&scale);
    let mut x = correction(rhs);
    let limit = 1e-10 * rhs.norm().max(1.);
    for _ in 0..4 {
        // Evaluate the residual of the original, unscaled system with
        // product and summation compensation for thin supported shells.
        let residual = DVector::from_fn(rhs.len(), |i, _| {
            let mut sum = rhs[i];
            let mut low = 0.;
            for j in 0..x.len() {
                let a = -k[(i, j)];
                let product = a * x[j];
                let next = sum + product;
                let v = next - sum;
                low += (sum - (next - v)) + (product - v) + a.mul_add(x[j], -product);
                sum = next;
            }
            sum + low
        });
        let norm = residual.norm();
        if norm <= limit {
            return Ok(x);
        }
        if !norm.is_finite() {
            return Err(Error::Solve(
                "nonfinite curved shell linear residual".into(),
            ));
        }
        let next = &x + correction(&residual);
        if next == x {
            return Err(Error::Solve(format!(
                "curved shell linear residual {norm:e} exceeds {limit:e}"
            )));
        }
        x = next;
    }
    Err(Error::Solve(
        "curved shell linear refinement exhausted at original tolerance".into(),
    ))
}

fn step(
    case: &Case,
    soil: &super::soil::Domain,
    model: &mut Option<Model>,
    old: &Frame,
    step: &Step,
) -> Result<Frame> {
    let mut next = old.clone();
    let mut q = v(&old.configuration);
    let n = q.len();
    let mut applied = v(&old.applied);
    let mut plate = v(&old.plate_load);
    let mut free = &[][..];
    next.holder = None;
    match step {
        Step::Install { contacts, holder } => {
            if model.is_some() || holder.trim().is_empty() || !contacts.is_empty() {
                return Err(Error::Input("shell installation requires absent shell, holding support and no contact pairs".into()));
            }
            let (m, g) = install(case, old)?;
            plate = m.loads(case, [0.; 3]);
            applied += &plate;
            *model = Some(m);
            next.active = Some(g);
            next.event = "install".into();
            next.holder = Some(holder.clone());
        }
        Step::Remove { holder } => {
            if holder.trim().is_empty() {
                return Err(Error::Input("shell removal needs a holding support".into()));
            }
            next.archived.push(
                next.active
                    .take()
                    .ok_or_else(|| Error::Input("no shell to remove".into()))?,
            );
            applied -= &plate;
            plate.fill(0.);
            *model = None;
            next.event = "remove".into();
            next.holder = Some(holder.clone());
        }
        Step::Load {
            prescribed_configuration,
            free_dofs,
            nodal_load,
            plate_traction_kpa,
        } => {
            if prescribed_configuration.len() != n
                || nodal_load.len() != n
                || free_dofs.iter().any(|&i| i >= n)
                || free_dofs.iter().collect::<BTreeSet<_>>().len() != free_dofs.len()
                || prescribed_configuration
                    .iter()
                    .chain(nodal_load)
                    .chain(plate_traction_kpa)
                    .any(|x| !x.is_finite())
            {
                return Err(Error::Input("invalid curved shell load layout".into()));
            }
            q = v(prescribed_configuration);
            applied = v(nodal_load);
            free = free_dofs;
            let mut rhs = applied.clone() - soil.force(&q);
            if let (Some(m), Some(g)) = (model.as_ref(), &next.active) {
                plate = m.loads(case, *plate_traction_kpa);
                m.face_actions(
                    case,
                    old.accepted_index + 1,
                    *plate_traction_kpa,
                    &mut plate,
                )?;
                applied += &plate;
                rhs += &plate - &m.stiffness * (&q - v(&g.installation));
            } else if case.face_loads.contains_key(&(old.accepted_index + 1))
                || *plate_traction_kpa != [0.; 3]
                || (soil.stiffness.is_none() && (!free.is_empty() || applied.norm() != 0.))
                || free.iter().any(|i| *i < 6 * case.nodes.len() && i % 6 >= 3)
                || (0..6 * case.nodes.len()).any(|i| i % 6 >= 3 && applied[i] != 0.)
            {
                return Err(Error::Input(
                    "uninstalled shell cannot carry load or free rotations".into(),
                ));
            }
            if !free.is_empty() {
                let sub = |k: &DMatrix<f64>| k.select_rows(free).select_columns(free);
                let mut k = soil
                    .stiffness
                    .as_ref()
                    .map_or_else(|| DMatrix::zeros(free.len(), free.len()), sub);
                if let Some(m) = model.as_ref() {
                    k += sub(&m.stiffness);
                }
                let dx = solve_supported(k, &rhs.select_rows(free))?;
                for (j, &i) in free.iter().enumerate() {
                    q[i] += dx[j];
                }
            }
            next.event = "load".into();
        }
    }
    let delta = &q - v(&old.configuration);
    let mut reaction = soil.force(&q) - &applied;
    if let (Some(m), Some(g)) = (model.as_ref(), &mut next.active) {
        reaction += m.recover(g, &q);
        g.plate_cumulative_load_work_kn_m += 0.5 * (&plate + v(&old.plate_load)).dot(&delta);
    }
    if reaction.select_rows(free).norm() > 1e-8 * applied.select_rows(free).norm().max(1.) {
        return Err(Error::Solve(
            "curved shell reassembled equilibrium residual exceeded tolerance".into(),
        ));
    }
    next.cumulative_applied_work_kn_m += 0.5 * (&applied + v(&old.applied)).dot(&delta);
    next.cumulative_reaction_work_kn_m += 0.5 * (&reaction + v(&old.reactions)).dot(&delta);
    next.soil_stress_kpa = soil.stress(&q);
    next.soil_translations_m = soil.translations(&q);
    next.configuration = list(&q);
    next.reactions = list(&reaction);
    next.applied = list(&applied);
    next.plate_load = list(&plate);
    next.accepted_index += 1;
    Ok(next)
}
pub fn run(case: &Case, until: usize, resume: Option<&Run>) -> Result<Run> {
    if until > case.steps.len() || case.faces.is_empty() || case.nodes.is_empty() {
        return Err(Error::Input(
            "empty curved shell or invalid accepted index".into(),
        ));
    }
    if case.face_loads.keys().any(|i| {
        i.checked_sub(1)
            .and_then(|j| case.steps.get(j))
            .is_none_or(|s| !matches!(s, Step::Load { .. }))
    }) {
        return Err(Error::Input("face load must belong to a load event".into()));
    }
    case.section.validate()?;
    geometries(&case.nodes, &case.faces, &case.directors)?;
    let soil = super::soil::Domain::new(case)?;
    let n = soil.size;
    let initial = Frame {
        accepted_index: 0,
        event: "initial".into(),
        soil_stress_kpa: soil.stress(&DVector::zeros(n)),
        soil_translations_m: soil.translations(&DVector::zeros(n)),
        holder: None,
        configuration: vec![0.; n],
        active: None,
        archived: vec![],
        applied: vec![0.; n],
        reactions: vec![0.; n],
        plate_load: vec![0.; n],
        cumulative_applied_work_kn_m: 0.,
        cumulative_reaction_work_kn_m: 0.,
    };
    let mut run = Run {
        schema: "curved-tri6-shell/1".into(),
        units: "m-kN-kPa".into(),
        section: case.section.clone(),
        rotation_convention: "global-small-rotation".into(),
        input_sha256: checkpoint::hash_json(case)?,
        solver_binary_sha256: checkpoint::hash_file(&std::env::current_exe()?)?,
        frames: vec![initial],
    };
    let mut model = None;
    if let Some(saved) = resume {
        if saved.frames.is_empty()
            || saved.frames.len() - 1 > until
            || saved.input_sha256 != run.input_sha256
            || saved.solver_binary_sha256 != run.solver_binary_sha256
        {
            return Err(Error::Input(
                "curved shell checkpoint identity or index mismatch".into(),
            ));
        }
        for event in &case.steps[..saved.frames.len() - 1] {
            run.frames.push(step(
                case,
                &soil,
                &mut model,
                run.frames.last().unwrap(),
                event,
            )?);
        }
        if checkpoint::hash_json(saved)? != checkpoint::hash_json(&run)? {
            return Err(Error::Input(
                "curved shell accepted path replay mismatch".into(),
            ));
        }
        run = saved.clone();
    }
    for event in &case.steps[run.frames.len() - 1..until] {
        run.frames.push(step(
            case,
            &soil,
            &mut model,
            run.frames.last().unwrap(),
            event,
        )?);
    }
    Ok(run)
}
