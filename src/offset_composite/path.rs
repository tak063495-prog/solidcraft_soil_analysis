//! Quasi-static dry Tet10 patch path. Event frames hold the installation
//! configuration. A following load frame specifies which reactions are released.
use super::{Patch, Plate};
use crate::{
    Error, Result, checkpoint, elastic, joint_law::PointState, joints, mesh::Mesh, plates,
};
use nalgebra::{DMatrix, DVector, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Soil {
    pub young_kpa: f64,
    pub poisson: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "policy", rename_all = "snake_case", deny_unknown_fields)]
pub enum Placement {
    StressFree,
    TransferForce {
        supplier: String,
        traction_kpa: Vec<[f64; 3]>,
        /// In local order A translations/rotations (36), then B translations (18).
        requested_nodal_actions: Vec<f64>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "event", rename_all = "snake_case", deny_unknown_fields)]
pub enum Step {
    Install {
        contacts: BTreeMap<String, Placement>,
        holder: String,
    },
    Remove {
        holder: String,
    },
    Load {
        /// All remaining translations and rotations have prescribed values.
        prescribed_configuration: Vec<f64>,
        free_dofs: Vec<usize>,
        nodal_load: Vec<f64>,
        /// Applied dead traction on the neutral plane, in addition to plate weight.
        plate_traction_kpa: [f64; 3],
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub plate: Plate,
    pub soil: BTreeMap<u32, Soil>,
    pub steps: Vec<Step>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContactState {
    pub id: String,
    pub region: joints::FaceRegion,
    pub soil_a: super::Face,
    pub soil_b: super::Face,
    pub arm_m: [f64; 3],
    pub local_instantaneous_work_kn_m: f64,
    pub nodal_instantaneous_work_kn_m: f64,
    pub nodal_actions: Vec<f64>,
    pub generation: u64,
    pub plate_generation: u64,
    pub soil_generations: [u64; 2],
    pub policy: joints::BirthPolicy,
    pub supplier: Option<String>,
    pub points: Vec<PointState>,
    pub endpoint_work_kn_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Generation {
    pub generation: u64,
    pub installation: Vec<f64>,
    pub plate_id: String,
    pub last_sections: Option<plates::MemberResult>,
    pub attachment_bending_kn: Vec<[f64; 3]>,
    pub plate_cumulative_load_work_kn_m: f64,
    pub contacts: Vec<ContactState>,
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
    pub applied: Vec<f64>,
    pub reactions: Vec<f64>,
    pub plate_traction_kpa: [f64; 3],
    pub instantaneous_external_kn_m: f64,
    pub cumulative_applied_work_kn_m: f64,
    pub cumulative_reaction_work_kn_m: f64,
    pub soil_stress_kpa: Vec<Vec<[f64; 6]>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Run {
    pub input_sha256: String,
    pub solver_binary_sha256: String,
    pub frames: Vec<Frame>,
}
struct Domain<'a> {
    case: &'a Case,
    patch: Patch,
    k_soil: DMatrix<f64>,
    cells: Vec<Cell>,
    plate_rows: Vec<usize>,
    pair_rows: Vec<Vec<usize>>,
}
fn gather(v: &DVector<f64>, rows: &[usize]) -> DVector<f64> {
    DVector::from_iterator(rows.len(), rows.iter().map(|&r| v[r]))
}
fn add_vector(v: &mut DVector<f64>, local: &DVector<f64>, rows: &[usize]) {
    for (&r, &f) in rows.iter().zip(local.iter()) {
        v[r] += f;
    }
}
fn add_matrix(k: &mut DMatrix<f64>, local: &DMatrix<f64>, rows: &[usize]) {
    for (i, &a) in rows.iter().enumerate() {
        for (j, &b) in rows.iter().enumerate() {
            k[(a, b)] += local[(i, j)];
        }
    }
}
pub(crate) struct Cell {
    rows: Vec<usize>,
    stress: Vec<DMatrix<f64>>,
}
pub(crate) fn assemble_soil(
    mesh: &Mesh,
    materials: &BTreeMap<u32, Soil>,
    nd: usize,
    node_rows: &[[usize; 3]],
) -> Result<(DMatrix<f64>, Vec<Cell>)> {
    let mut k_soil = DMatrix::zeros(nd, nd);
    let mut cells = Vec::new();
    for (i, id) in mesh.material_ids.iter().enumerate() {
        let soil = materials
            .get(id)
            .ok_or_else(|| Error::Input("missing dry elastic soil properties".into()))?;
        let d = elastic::constitutive(soil.young_kpa, soil.poisson)?;
        let rows: Vec<_> = mesh
            .cell_nodes(i)
            .into_iter()
            .flat_map(|n| node_rows[n])
            .collect();
        let mut stress = Vec::new();
        for gp in mesh.integration(i)? {
            let stress_operator = d * &gp.b;
            let stiffness = gp.b.transpose() * &stress_operator * gp.weight_m3;
            stress.push(DMatrix::from_column_slice(
                6,
                rows.len(),
                stress_operator.as_slice(),
            ));
            add_matrix(&mut k_soil, &stiffness, &rows);
        }
        cells.push(Cell { rows, stress });
    }
    Ok((k_soil, cells))
}

impl Cell {
    pub(crate) fn stresses(&self, q: &DVector<f64>) -> Vec<[f64; 6]> {
        let u = gather(q, &self.rows);
        self.stress
            .iter()
            .map(|b| {
                let sig = b * &u;
                std::array::from_fn(|i| sig[i])
            })
            .collect()
    }
}
struct Trial {
    force: DVector<f64>,
    tangent: DMatrix<f64>,
    active: Option<Generation>,
}
impl<'a> Domain<'a> {
    fn new(mesh: &Mesh, case: &'a Case) -> Result<Self> {
        if mesh.tetrahedra.len() > 24 {
            return Err(Error::Input(
                "dry offset composite patch supports up to 24 Tet10 cells".into(),
            ));
        }
        let patch = Patch::new(mesh, &case.plate)?;
        let soil_dofs = mesh.nodes.len() * 3;
        let nd = soil_dofs + 18;
        let node_rows: Vec<_> = (0..mesh.nodes.len())
            .map(|i| [3 * i, 3 * i + 1, 3 * i + 2])
            .collect();
        let (k_soil, cells) = assemble_soil(mesh, &case.soil, nd, &node_rows)?;
        let mut plate_rows = Vec::with_capacity(36);
        for (i, &node) in case.plate.soil_a.nodes.iter().enumerate() {
            plate_rows.extend(3 * node..3 * node + 3);
            plate_rows.extend(soil_dofs + 3 * i..soil_dofs + 3 * i + 3);
        }
        let pair_rows = case
            .plate
            .contacts
            .iter()
            .map(|c| {
                let mut rows = plate_rows.clone();
                rows.extend(c.soil_b.nodes.iter().flat_map(|n| 3 * n..3 * n + 3));
                rows
            })
            .collect();
        Ok(Self {
            case,
            patch,
            k_soil,
            cells,
            plate_rows,
            pair_rows,
        })
    }
    fn install(
        &self,
        d: &DVector<f64>,
        number: u64,
        placements: &BTreeMap<String, Placement>,
    ) -> Result<Generation> {
        if placements.len() != self.case.plate.contacts.len() {
            return Err(Error::Input(
                "each physical contact needs its own placement".into(),
            ));
        }
        let mut contacts = Vec::new();
        for (i, c) in self.case.plate.contacts.iter().enumerate() {
            let placement = placements
                .get(&c.id)
                .ok_or_else(|| Error::Input("missing contact placement".into()))?;
            let local = gather(d, &self.pair_rows[i]);
            let mut force = DVector::zeros(54);
            let mut points = Vec::new();
            let (policy, supplier) = match placement {
                Placement::StressFree => (joints::BirthPolicy::StressFree, None),
                Placement::TransferForce {
                    supplier,
                    traction_kpa,
                    requested_nodal_actions,
                } => {
                    if supplier.trim().is_empty()
                        || traction_kpa.len() != self.patch.points.len()
                        || requested_nodal_actions.len() != 54
                        || requested_nodal_actions.iter().any(|x| !x.is_finite())
                    {
                        return Err(Error::Input("initial force transfer needs a supplier, every GP and all nodal actions".into()));
                    }
                    (joints::BirthPolicy::TransferForce, Some(supplier.clone()))
                }
            };
            for (j, p) in self.patch.points.iter().enumerate() {
                let jump = &p.jump * &local;
                let traction = match placement {
                    Placement::StressFree => [0.; 3],
                    Placement::TransferForce { traction_kpa, .. } => traction_kpa[j],
                };
                let state = c.material.birth(
                    self.patch.normal.into(),
                    [jump[0], jump[1], jump[2]],
                    traction,
                )?;
                force += p.jump.transpose() * Vector3::from(traction) * p.weight;
                points.push(state);
            }
            if let Placement::TransferForce {
                requested_nodal_actions,
                ..
            } = placement
                && force
                    .iter()
                    .zip(requested_nodal_actions)
                    .any(|(f, w)| (f - w).abs() > 1e-10)
            {
                return Err(Error::Input(
                    "initial force transfer misses a translation or arm moment".into(),
                ));
            }
            contacts.push(ContactState {
                id: c.id.clone(),
                soil_a: self.case.plate.soil_a.clone(),
                soil_b: c.soil_b.clone(),
                arm_m: self.patch.arm.into(),
                local_instantaneous_work_kn_m: 0.,
                nodal_instantaneous_work_kn_m: 0.,
                nodal_actions: vec![],
                region: joints::FaceRegion {
                    face_index: i,
                    owner_a: self.case.plate.soil_a.cell,
                    owner_b: c.soil_b.cell,
                    area_m2: self.patch.points.iter().map(|p| p.weight).sum(),
                },
                generation: number,
                plate_generation: number,
                soil_generations: [1, 1],
                policy,
                supplier,
                points,
                endpoint_work_kn_m: 0.,
            });
        }
        Ok(Generation {
            generation: number,
            plate_id: self.case.plate.id.clone(),
            last_sections: None,
            attachment_bending_kn: vec![],
            plate_cumulative_load_work_kn_m: 0.,
            installation: gather(d, &self.plate_rows).as_slice().to_vec(),
            contacts,
        })
    }

    fn trial(&self, d: &DVector<f64>, active: Option<&Generation>) -> Result<Trial> {
        let mut t = Trial {
            force: &self.k_soil * d,
            tangent: self.k_soil.clone(),
            active: active.cloned(),
        };
        if let Some(generation) = &mut t.active {
            let relative =
                gather(d, &self.plate_rows) - DVector::from_column_slice(&generation.installation);
            let k = &self.patch.plate.stiffness;
            add_vector(&mut t.force, &(k * &relative), &self.plate_rows);
            add_matrix(&mut t.tangent, k, &self.plate_rows);
            let sections = self
                .patch
                .plate
                .recover(self.case.plate.id.clone(), &relative);
            generation.attachment_bending_kn = sections
                .samples
                .iter()
                .map(|s| {
                    std::array::from_fn(|i| {
                        s.bending_moment_kn[i] - self.case.plate.offset_m * s.membrane_force_kn_m[i]
                    })
                })
                .collect();
            generation.last_sections = Some(sections);
            for (i, c) in self.case.plate.contacts.iter().enumerate() {
                let old = &mut generation.contacts[i];
                let rows = &self.pair_rows[i];
                let ct = self
                    .patch
                    .assemble_contact(&c.material, &old.points, &gather(d, rows))?;
                add_vector(&mut t.force, &ct.force, rows);
                add_matrix(&mut t.tangent, &ct.tangent, rows);
                old.endpoint_work_kn_m = ct
                    .points
                    .iter()
                    .zip(&self.patch.points)
                    .map(|(s, p)| s.endpoint_work_kpa_m * p.weight)
                    .sum();
                old.local_instantaneous_work_kn_m = ct
                    .points
                    .iter()
                    .zip(&self.patch.points)
                    .map(|(s, p)| {
                        Vector3::from(s.traction_kpa).dot(&Vector3::from(s.jump_m)) * p.weight
                    })
                    .sum();
                old.nodal_instantaneous_work_kn_m = ct.force.dot(&gather(d, rows));
                old.nodal_actions = ct.force.as_slice().to_vec();
                old.points = ct.points;
            }
        }
        Ok(t)
    }

    fn solve(
        &self,
        d: &mut DVector<f64>,
        active: Option<&Generation>,
        f: &DVector<f64>,
        free: &[usize],
    ) -> Result<Trial> {
        for _ in 0..50 {
            let t = self.trial(d, active)?;
            let r = gather(&(f - &t.force), free);
            if r.norm() <= 1e-8 * gather(f, free).norm().max(1.) {
                return Ok(t);
            }
            let k = t.tangent.select_rows(free).select_columns(free);
            let dx = k
                .clone()
                .lu()
                .solve(&r)
                .ok_or_else(|| Error::Solve("singular free offset composite dofs".into()))?;
            if (&k * &dx - &r).norm() > 1e-10 * r.norm().max(1.) {
                return Err(Error::Solve(
                    "offset composite linear residual exceeded tolerance".into(),
                ));
            }
            add_vector(d, &dx, free);
        }
        Err(Error::Solve(
            "offset composite equilibrium did not converge".into(),
        ))
    }

    fn total_traction(&self, mut applied: [f64; 3]) -> [f64; 3] {
        let s = &self.case.plate.section;
        applied[2] -= s.thickness_m * s.unit_weight_kn_m3;
        applied
    }
    fn step(&self, old: &Frame, step: &Step) -> Result<Frame> {
        let mut next = old.clone();
        let mut d = DVector::from_column_slice(&old.configuration);
        let mut applied = DVector::from_column_slice(&old.applied);
        let mut free = &[][..];
        next.holder = None;
        match step {
            Step::Install { contacts, holder } => {
                if old.active.is_some() || holder.trim().is_empty() {
                    return Err(Error::Input(
                        "installation needs absent plate and named holding support".into(),
                    ));
                }
                next.active = Some(self.install(&d, 1 + old.archived.len() as u64, contacts)?);
                next.plate_traction_kpa = self.total_traction([0.; 3]);
                add_vector(
                    &mut applied,
                    &self.patch.plate.uniform_load(next.plate_traction_kpa),
                    &self.plate_rows,
                );
                next.event = "install".into();
                next.holder = Some(holder.clone());
            }
            Step::Remove { holder } => {
                if holder.trim().is_empty() {
                    return Err(Error::Input("removal needs a named holding support".into()));
                }
                next.archived.push(
                    next.active
                        .take()
                        .ok_or_else(|| Error::Input("no installed plate to remove".into()))?,
                );
                let load = self.patch.plate.uniform_load(old.plate_traction_kpa);
                add_vector(&mut applied, &(-load), &self.plate_rows);
                applied.rows_mut(d.len() - 18, 18).fill(0.);
                next.plate_traction_kpa = [0.; 3];
                next.event = "remove".into();
                next.holder = Some(holder.clone());
            }
            Step::Load {
                prescribed_configuration,
                free_dofs,
                nodal_load,
                plate_traction_kpa,
            } => {
                if prescribed_configuration.len() != d.len()
                    || nodal_load.len() != d.len()
                    || free_dofs.iter().any(|&r| r >= d.len())
                    || free_dofs.iter().collect::<BTreeSet<_>>().len() != free_dofs.len()
                    || prescribed_configuration
                        .iter()
                        .chain(nodal_load)
                        .chain(plate_traction_kpa)
                        .any(|x| !x.is_finite())
                {
                    return Err(Error::Input(
                        "invalid offset composite load/constraint layout".into(),
                    ));
                }
                if old.active.is_none()
                    && (*plate_traction_kpa != [0.; 3]
                        || free_dofs.iter().any(|&r| r >= d.len() - 18)
                        || nodal_load[d.len() - 18..].iter().any(|&x| x != 0.))
                {
                    return Err(Error::Input(
                        "uninstalled plate cannot carry load or free rotations".into(),
                    ));
                }
                d.copy_from_slice(prescribed_configuration);
                free = free_dofs;
                for &r in free {
                    d[r] = old.configuration[r];
                }
                applied.copy_from_slice(nodal_load);
                next.plate_traction_kpa = if next.active.is_some() {
                    self.total_traction(*plate_traction_kpa)
                } else {
                    [0.; 3]
                };
                add_vector(
                    &mut applied,
                    &self.patch.plate.uniform_load(next.plate_traction_kpa),
                    &self.plate_rows,
                );
                next.event = "load".into();
            }
        }
        let trial = self.solve(&mut d, next.active.as_ref(), &applied, free)?;
        let mut reactions = &trial.force - &applied;
        for &r in free {
            reactions[r] = 0.;
        }
        let delta = &d - DVector::from_column_slice(&old.configuration);
        next.cumulative_applied_work_kn_m +=
            0.5 * (&applied + DVector::from_column_slice(&old.applied)).dot(&delta);
        next.cumulative_reaction_work_kn_m +=
            0.5 * (&reactions + DVector::from_column_slice(&old.reactions)).dot(&delta);
        next.instantaneous_external_kn_m = (&applied + &reactions).dot(&d);
        next.active = trial.active;
        if let Some(g) = &mut next.active {
            let average = std::array::from_fn(|i| {
                0.5 * (old.plate_traction_kpa[i] + next.plate_traction_kpa[i])
            });
            g.plate_cumulative_load_work_kn_m += self
                .patch
                .plate
                .uniform_load(average)
                .dot(&gather(&delta, &self.plate_rows));
        }
        next.soil_stress_kpa = self.cells.iter().map(|c| c.stresses(&d)).collect();
        next.applied = applied.as_slice().to_vec();
        next.reactions = reactions.as_slice().to_vec();
        next.configuration = d.as_slice().to_vec();
        next.accepted_index += 1;
        Ok(next)
    }
}

pub fn run(mesh: &Mesh, case: &Case, until: usize, resume: Option<&Run>) -> Result<Run> {
    let model = Domain::new(mesh, case)?;
    let n = model.k_soil.nrows();
    let input_sha256 = checkpoint::hash_json(&(mesh, case))?;
    let initial = Frame {
        accepted_index: 0,
        event: "initial".into(),
        holder: None,
        configuration: vec![0.; n],
        active: None,
        archived: vec![],
        applied: vec![0.; n],
        reactions: vec![0.; n],
        plate_traction_kpa: [0.; 3],
        instantaneous_external_kn_m: 0.,
        cumulative_applied_work_kn_m: 0.,
        cumulative_reaction_work_kn_m: 0.,
        soil_stress_kpa: model
            .cells
            .iter()
            .map(|c| vec![[0.; 6]; c.stress.len()])
            .collect(),
    };
    if until > case.steps.len() {
        return Err(Error::Input(
            "offset composite checkpoint index out of range".into(),
        ));
    }
    let mut run = Run {
        input_sha256,
        solver_binary_sha256: checkpoint::hash_file(&std::env::current_exe()?)?,
        frames: vec![initial],
    };
    if let Some(saved) = resume {
        if saved.solver_binary_sha256 != run.solver_binary_sha256
            || saved.input_sha256 != run.input_sha256
            || saved.frames.is_empty()
            || saved.frames.len() - 1 > until
        {
            return Err(Error::Input(
                "offset composite checkpoint input/index mismatch".into(),
            ));
        }
        for step in &case.steps[..saved.frames.len() - 1] {
            run.frames
                .push(model.step(run.frames.last().unwrap(), step)?);
        }
        if checkpoint::hash_json(&run)? != checkpoint::hash_json(saved)? {
            return Err(Error::Input(
                "offset composite accepted path replay mismatch".into(),
            ));
        }
        run = saved.clone();
    }
    for step in &case.steps[run.frames.len() - 1..until] {
        run.frames
            .push(model.step(run.frames.last().unwrap(), step)?);
    }
    Ok(run)
}
