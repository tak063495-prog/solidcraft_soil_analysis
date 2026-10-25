//! Linear 3D frame substructure. Interior translations and rotations are solved
//! and statically condensed exactly; no rotational spring or diagonal shift.
use crate::{Error, Result, mesh::Mesh};
use nalgebra::{DMatrix, DVector, Dyn, Matrix3, SMatrix, SVector, Vector3, linalg::Cholesky};
use serde::{Deserialize, Serialize};
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
    rc::Rc,
};
type M12 = SMatrix<f64, 12, 12>;
type V12 = SVector<f64, 12>;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Node {
    pub id: String,
    pub position_m: [f64; 3],
    pub soil_node_set: Option<String>,
    #[serde(default)]
    pub match_soil_node: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Formulation {
    EulerBernoulli,
    Timoshenko {
        shear_area_y_m2: f64,
        shear_area_z_m2: f64,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Section {
    pub id: String,
    pub young_kpa: f64,
    pub poisson: f64,
    pub area_m2: f64,
    pub iy_m4: f64,
    pub iz_m4: f64,
    pub torsion_j_m4: f64,
    pub unit_weight_kn_m3: f64,
    pub formulation: Formulation,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Member {
    pub id: String,
    pub nodes: [String; 2],
    pub section_id: String,
    pub local_y: [f64; 3],
    pub initially_active: bool,
    #[serde(default)]
    pub initial_axial_force_kn: f64,
    #[serde(default)]
    pub constant_line_load_kn_m: [f64; 3],
    #[serde(default)]
    pub boundary_line_load_kn_m: [f64; 3],
    #[serde(default)]
    pub rotational_releases: [[bool; 3]; 2],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Constraint {
    pub node_id: String,
    pub direction: [f64; 3],
    pub value: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NodalLoad {
    pub node_id: String,
    #[serde(default)]
    pub constant_force_kn: [f64; 3],
    #[serde(default)]
    pub boundary_force_kn: [f64; 3],
    #[serde(default)]
    pub constant_moment_kn_m: [f64; 3],
    #[serde(default)]
    pub boundary_moment_kn_m: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    pub schema: String,
    pub nodes: Vec<Node>,
    pub sections: Vec<Section>,
    pub members: Vec<Member>,
    pub rotation_length_scale_m: f64,
    pub maximum_condensed_dofs: usize,
    #[serde(default)]
    pub rotation_constraints: Vec<Constraint>,
    #[serde(default)]
    pub translation_constraints: Vec<Constraint>,
    #[serde(default)]
    pub nodal_loads: Vec<NodalLoad>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub plate_sections: Vec<crate::plates::Section>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub plates: Vec<crate::plates::Member>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub active_members: Vec<bool>,
    pub reference_global_kinematics: Vec<[f64; 12]>,
    pub node_displacement_m: Vec<[f64; 3]>,
    pub node_rotation_rad: Vec<[f64; 3]>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub active_plates: Vec<bool>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub plate_reference_global_kinematics: Vec<[f64; 18]>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MemberResult {
    pub id: String,
    pub local_axes: [[f64; 3]; 3],
    pub length_m: f64,
    pub local_end_actions: [f64; 12],
    pub local_end_relative_rotation_rad: [[f64; 3]; 2],
    pub section_forces_at_start_mid_end: [[f64; 6]; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub active_members: usize,
    pub condensed_dofs: usize,
    pub interface_translational_dofs: usize,
    pub true_relative_interior_residual: f64,
    pub node_positions_m: Vec<[f64; 3]>,
    pub node_displacement_m: Vec<[f64; 3]>,
    pub node_rotation_rad: Vec<[f64; 3]>,
    pub node_active: Vec<bool>,
    pub external_force_kn: [f64; 3],
    pub external_moment_origin_kn_m: [f64; 3],
    pub interface_external_force_kn: [f64; 3],
    pub interface_external_moment_origin_kn_m: [f64; 3],
    pub support_force_kn: [f64; 3],
    pub support_moment_origin_kn_m: [f64; 3],
    pub members: Vec<MemberResult>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub plates: Vec<crate::plates::MemberResult>,
}
pub struct Evaluation {
    pub state: State,
    pub report: Report,
    pub interface_dofs: Vec<usize>,
    pub interface_internal_kn: DVector<f64>,
    pub interface_tangent_kn_m: DMatrix<f64>,
}

pub struct Kernel {
    pub length: f64,
    pub axes: Matrix3<f64>,
    pub transform: M12,
    pub stiffness: M12,
    local: M12,
    retained: Vec<usize>,
    released: Vec<usize>,
    release_inverse: DMatrix<f64>,
}
impl Section {
    pub fn validate(&self) -> Result<()> {
        crate::elastic::constitutive(self.young_kpa, self.poisson)?;
        if self.id.is_empty()
            || [self.area_m2, self.iy_m4, self.iz_m4, self.torsion_j_m4]
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.)
            || !self.unit_weight_kn_m3.is_finite()
            || self.unit_weight_kn_m3 < 0.
        {
            return Err(Error::Input("invalid beam section/weight".into()));
        }
        if let Formulation::Timoshenko {
            shear_area_y_m2,
            shear_area_z_m2,
        } = self.formulation
            && [shear_area_y_m2, shear_area_z_m2]
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.)
        {
            return Err(Error::Input(
                "positive effective shear areas required".into(),
            ));
        }
        Ok(())
    }
}
impl Kernel {
    pub fn new(
        a: [f64; 3],
        b: [f64; 3],
        local_y: [f64; 3],
        section: &Section,
        releases: [[bool; 3]; 2],
    ) -> Result<Self> {
        section.validate()?;
        let dx = Vector3::from(b) - Vector3::from(a);
        let length = dx.norm();
        let ey = Vector3::from(local_y);
        if !length.is_finite()
            || length < 1e-8
            || !ey.iter().all(|v| v.is_finite())
            || ey.norm() < 1e-12
        {
            return Err(Error::Input(
                "zero/invalid beam length or local axis".into(),
            ));
        }
        let ex = dx / length;
        let projected = ey - ex * ey.dot(&ex);
        if projected.norm() < 1e-8 * ey.norm() {
            return Err(Error::Input("beam local_y is parallel to its axis".into()));
        }
        let ey = projected.normalize();
        let ez = ex.cross(&ey);
        let axes = Matrix3::from_rows(&[ex.transpose(), ey.transpose(), ez.transpose()]);
        let mut transform = M12::zeros();
        for block in 0..4 {
            transform
                .fixed_view_mut::<3, 3>(3 * block, 3 * block)
                .copy_from(&axes);
        }
        let e = section.young_kpa;
        let g = e / (2. * (1. + section.poisson));
        let mut local = M12::zeros();
        for (ids, value) in [
            ([0, 6], e * section.area_m2 / length),
            ([3, 9], g * section.torsion_j_m4 / length),
        ] {
            local[(ids[0], ids[0])] += value;
            local[(ids[1], ids[1])] += value;
            local[(ids[0], ids[1])] -= value;
            local[(ids[1], ids[0])] -= value;
        }
        let (ay, az) = match section.formulation {
            Formulation::EulerBernoulli => (None, None),
            Formulation::Timoshenko {
                shear_area_y_m2,
                shear_area_z_m2,
            } => (Some(shear_area_y_m2), Some(shear_area_z_m2)),
        };
        for (ids, ei, area, sign) in [
            ([1, 5, 7, 11], e * section.iz_m4, ay, 1.),
            ([2, 4, 8, 10], e * section.iy_m4, az, -1.),
        ] {
            let phi = area.map_or(0., |a| 12. * ei / (g * a * length * length));
            let l = length;
            let values = [
                [12., 6. * l * sign, -12., 6. * l * sign],
                [
                    6. * l * sign,
                    (4. + phi) * l * l,
                    -6. * l * sign,
                    (2. - phi) * l * l,
                ],
                [-12., -6. * l * sign, 12., -6. * l * sign],
                [
                    6. * l * sign,
                    (2. - phi) * l * l,
                    -6. * l * sign,
                    (4. + phi) * l * l,
                ],
            ];
            for i in 0..4 {
                for j in 0..4 {
                    local[(ids[i], ids[j])] += values[i][j] * ei / (l.powi(3) * (1. + phi));
                }
            }
        }
        let released: Vec<_> = (0..2)
            .flat_map(|end| {
                (0..3)
                    .filter(move |&a| releases[end][a])
                    .map(move |a| 6 * end + 3 + a)
            })
            .collect();
        let retained: Vec<_> = (0..12).filter(|i| !released.contains(i)).collect();
        let inverse = if released.is_empty() {
            DMatrix::zeros(0, 0)
        } else {
            let qq = DMatrix::from_fn(released.len(), released.len(), |i, j| {
                local[(released[i], released[j])]
            });
            let inverse = qq.clone().try_inverse().ok_or_else(|| {
                Error::Input(
                    "singular beam release block (including double torsion release)".into(),
                )
            })?;
            if (&qq * &inverse - DMatrix::identity(qq.nrows(), qq.nrows())).norm() > 1e-9 {
                return Err(Error::Input("ill-conditioned beam release block".into()));
            }
            inverse
        };
        let mut condensed = local;
        if !released.is_empty() {
            condensed.fill(0.);
            for &i in &retained {
                for &j in &retained {
                    let correction = (0..released.len())
                        .map(|p| {
                            (0..released.len())
                                .map(|q| {
                                    local[(i, released[p])]
                                        * inverse[(p, q)]
                                        * local[(released[q], j)]
                                })
                                .sum::<f64>()
                        })
                        .sum::<f64>();
                    condensed[(i, j)] = local[(i, j)] - correction;
                }
            }
        }
        let stiffness = transform.transpose() * condensed * transform;
        Ok(Self {
            length,
            axes,
            transform,
            stiffness,
            local,
            retained,
            released,
            release_inverse: inverse,
        })
    }
    fn condense_load(&self, f: &V12) -> V12 {
        if self.released.is_empty() {
            return *f;
        }
        let q = DVector::from_iterator(self.released.len(), self.released.iter().map(|&i| f[i]));
        let answer = &self.release_inverse * q;
        let mut out = V12::zeros();
        for &i in &self.retained {
            out[i] = f[i]
                - (0..self.released.len())
                    .map(|p| self.local[(i, self.released[p])] * answer[p])
                    .sum::<f64>();
        }
        out
    }
    pub fn uniform_load(&self, global: [f64; 3]) -> V12 {
        let q = self.axes * Vector3::from(global);
        let l = self.length;
        let mut f = V12::zeros();
        for end in 0..2 {
            for a in 0..3 {
                f[6 * end + a] = q[a] * l / 2.;
            }
        }
        f[4] = -q[2] * l * l / 12.;
        f[5] = q[1] * l * l / 12.;
        f[10] = -f[4];
        f[11] = -f[5];
        f
    }
    pub fn recover(&self, global_relative: &V12, line: [f64; 3], prestress: f64) -> (V12, V12) {
        let load = self.uniform_load(line);
        let mut local = self.transform * global_relative;
        if !self.released.is_empty() {
            let rhs = DVector::from_iterator(
                self.released.len(),
                self.released.iter().map(|&i| {
                    load[i]
                        - self
                            .retained
                            .iter()
                            .map(|&j| self.local[(i, j)] * local[j])
                            .sum::<f64>()
                }),
            );
            let recovered = &self.release_inverse * rhs;
            for (p, &i) in self.released.iter().enumerate() {
                local[i] = recovered[p];
            }
        }
        let mut action = self.local * local - load;
        action[0] -= prestress;
        action[6] += prestress;
        (local, action)
    }
}

struct MemberData {
    nodes: [usize; 2],
    section: usize,
    kernel: Kernel,
}
struct PlateData {
    nodes: [usize; 3],
    section: usize,
    kernel: crate::plates::Kernel,
}
struct Compiled {
    active_nodes: Vec<bool>,
    k: DMatrix<f64>,
    t: DMatrix<f64>,
    q: DMatrix<f64>,
    offset: DVector<f64>,
    chol: Option<Cholesky<f64, Dyn>>,
    kqq: DMatrix<f64>,
    kc: DMatrix<f64>,
    interface_dofs: Vec<usize>,
    normals: Vec<(usize, Vector3<f64>)>,
}
pub struct Model {
    pub configuration: Configuration,
    pub positions: Vec<[f64; 3]>,
    pub attachments: Vec<Option<usize>>,
    members: Vec<MemberData>,
    plates: Vec<PlateData>,
    cache: RefCell<BTreeMap<Vec<bool>, Rc<Compiled>>>,
}

type ConstraintBasis = (Vector3<f64>, Vec<Vector3<f64>>, Vec<Vector3<f64>>);
type AssembledLoads = (DVector<f64>, DVector<f64>, Vec<[f64; 3]>);

fn basis(eqs: &[(Vector3<f64>, f64)]) -> Result<ConstraintBasis> {
    let mut normals = Vec::<Vector3<f64>>::new();
    let mut offset = Vector3::zeros();
    for (a, value) in eqs {
        let mut n = *a;
        for q in &normals {
            n -= q * n.dot(q);
        }
        let norm = n.norm();
        let r = value - a.dot(&offset);
        if norm < 1e-10 {
            if r.abs() > 1e-11 * (1. + value.abs()) {
                return Err(Error::Input("conflicting structural constraints".into()));
            }
        } else {
            let q = n / norm;
            offset += q * (r / norm);
            normals.push(q);
        }
    }
    let mut free = Vec::<Vector3<f64>>::new();
    for a in 0..3 {
        let mut v = Vector3::zeros();
        v[a] = 1.;
        for _ in 0..2 {
            for q in normals.iter().chain(&free) {
                v -= q * v.dot(q);
            }
        }
        if v.norm() > 1e-10 {
            free.push(v.normalize());
        }
    }
    if normals.len() + free.len() != 3 {
        return Err(Error::Input("structural constraint basis failed".into()));
    }
    Ok((offset, free, normals))
}

impl Model {
    pub fn new(mesh: &Mesh, c: Configuration) -> Result<Self> {
        if !["soil-fem-beams/1", "soil-fem-structures/1"].contains(&c.schema.as_str())
            || c.nodes.is_empty()
            || c.nodes.len() > 500
            || (c.members.is_empty() && c.plates.is_empty())
            || (!c.members.is_empty() && c.sections.is_empty())
            || (!c.plates.is_empty() && c.plate_sections.is_empty())
            || c.members.len() + c.plates.len() > 2000
            || !c.rotation_length_scale_m.is_finite()
            || c.rotation_length_scale_m <= 0.
            || c.maximum_condensed_dofs == 0
            || c.maximum_condensed_dofs > 3000
        {
            return Err(Error::Input("invalid beam system/scaling/budget".into()));
        }
        let node_ids: BTreeMap<_, _> = c
            .nodes
            .iter()
            .enumerate()
            .map(|(i, n)| (n.id.clone(), i))
            .collect();
        let sections: BTreeMap<_, _> = c
            .sections
            .iter()
            .enumerate()
            .map(|(i, s)| (s.id.clone(), i))
            .collect();
        if node_ids.len() != c.nodes.len()
            || sections.len() != c.sections.len()
            || c.nodes
                .iter()
                .any(|n| n.id.is_empty() || !n.position_m.iter().all(|x| x.is_finite()))
        {
            return Err(Error::Input(
                "invalid/duplicate structural identifiers".into(),
            ));
        }
        for s in &c.sections {
            s.validate()?;
        }
        let plate_sections: BTreeMap<_, _> = c
            .plate_sections
            .iter()
            .enumerate()
            .map(|(i, s)| (s.id.clone(), i))
            .collect();
        if plate_sections.len() != c.plate_sections.len() {
            return Err(Error::Input("duplicate plate section identifier".into()));
        }
        for section in &c.plate_sections {
            section.validate()?;
        }
        let mut positions = Vec::new();
        let mut attachments = Vec::new();
        for n in &c.nodes {
            if n.match_soil_node && n.soil_node_set.is_some() {
                return Err(Error::Input(
                    "choose coordinate matching or a soil node set".into(),
                ));
            }
            let attached = if n.match_soil_node {
                let found: Vec<_> = mesh
                    .nodes
                    .iter()
                    .enumerate()
                    .filter(|(_, p)| {
                        (Vector3::from(**p) - Vector3::from(n.position_m)).norm() <= 1e-8
                    })
                    .map(|(i, _)| i)
                    .collect();
                if found.len() != 1 {
                    return Err(Error::Input(
                        "beam coordinate must match exactly one soil mesh node".into(),
                    ));
                }
                Some(found[0])
            } else if let Some(set) = &n.soil_node_set {
                let nodes = mesh
                    .node_sets
                    .get(set)
                    .filter(|s| s.len() == 1)
                    .ok_or_else(|| {
                        Error::Input(
                            "beam soil_node_set must identify exactly one mesh node".into(),
                        )
                    })?;
                if (Vector3::from(mesh.nodes[nodes[0]]) - Vector3::from(n.position_m)).norm() > 1e-8
                {
                    return Err(Error::Input("beam attachment coordinate mismatch".into()));
                }
                Some(nodes[0])
            } else {
                None
            };
            positions.push(attached.map_or(n.position_m, |i| mesh.nodes[i]));
            attachments.push(attached);
        }
        let mut names = BTreeSet::new();
        let mut members = Vec::new();
        for m in &c.members {
            if m.id.is_empty()
                || !names.insert(m.id.clone())
                || !m.initial_axial_force_kn.is_finite()
                || !m
                    .constant_line_load_kn_m
                    .iter()
                    .chain(&m.boundary_line_load_kn_m)
                    .all(|x| x.is_finite())
            {
                return Err(Error::Input("invalid/duplicate beam member".into()));
            }
            let nodes = [
                *node_ids
                    .get(&m.nodes[0])
                    .ok_or_else(|| Error::Input("unknown beam node".into()))?,
                *node_ids
                    .get(&m.nodes[1])
                    .ok_or_else(|| Error::Input("unknown beam node".into()))?,
            ];
            let section = *sections
                .get(&m.section_id)
                .ok_or_else(|| Error::Input("unknown beam section".into()))?;
            let kernel = Kernel::new(
                positions[nodes[0]],
                positions[nodes[1]],
                m.local_y,
                &c.sections[section],
                m.rotational_releases,
            )?;
            members.push(MemberData {
                nodes,
                section,
                kernel,
            });
        }
        let mut plates = Vec::new();
        for plate in &c.plates {
            if plate.id.is_empty()
                || !names.insert(plate.id.clone())
                || !plate
                    .constant_traction_kpa
                    .iter()
                    .chain(&plate.boundary_traction_kpa)
                    .all(|v| v.is_finite())
            {
                return Err(Error::Input("invalid/duplicate plate member".into()));
            }
            let mut nodes = [0; 3];
            for (i, name) in plate.nodes.iter().enumerate() {
                nodes[i] = *node_ids
                    .get(name)
                    .ok_or_else(|| Error::Input("unknown plate node".into()))?;
            }
            let section = *plate_sections
                .get(&plate.section_id)
                .ok_or_else(|| Error::Input("unknown plate section".into()))?;
            let kernel = crate::plates::Kernel::new(
                nodes.map(|i| positions[i]),
                &c.plate_sections[section],
            )?;
            plates.push(PlateData {
                nodes,
                section,
                kernel,
            });
        }
        for constraint in c
            .rotation_constraints
            .iter()
            .chain(&c.translation_constraints)
        {
            let i = *node_ids
                .get(&constraint.node_id)
                .ok_or_else(|| Error::Input("unknown structural constraint node".into()))?;
            if !constraint.value.is_finite()
                || !constraint.direction.iter().all(|v| v.is_finite())
                || Vector3::from(constraint.direction).norm() < 1e-12
            {
                return Err(Error::Input("invalid structural constraint".into()));
            }
            if c.translation_constraints
                .iter()
                .any(|a| std::ptr::eq(a, constraint))
                && attachments[i].is_some()
            {
                return Err(Error::Input(
                    "attached translations must use soil boundary constraints".into(),
                ));
            }
        }
        for load in &c.nodal_loads {
            if !node_ids.contains_key(&load.node_id)
                || !load
                    .constant_force_kn
                    .iter()
                    .chain(&load.boundary_force_kn)
                    .chain(&load.constant_moment_kn_m)
                    .chain(&load.boundary_moment_kn_m)
                    .all(|v| v.is_finite())
            {
                return Err(Error::Input("invalid structural load".into()));
            }
        }
        Ok(Self {
            configuration: c,
            positions,
            attachments,
            members,
            plates,
            cache: RefCell::new(BTreeMap::new()),
        })
    }
    pub fn initial_state(&self) -> State {
        State {
            active_members: self
                .configuration
                .members
                .iter()
                .map(|m| m.initially_active)
                .collect(),
            reference_global_kinematics: vec![[0.; 12]; self.members.len()],
            node_displacement_m: vec![[0.; 3]; self.positions.len()],
            node_rotation_rad: vec![[0.; 3]; self.positions.len()],
            active_plates: self
                .configuration
                .plates
                .iter()
                .map(|p| p.initially_active)
                .collect(),
            plate_reference_global_kinematics: vec![[0.; 18]; self.plates.len()],
        }
    }
    pub fn validate_state(&self, s: &State) -> Result<()> {
        if s.active_members.len() != self.members.len()
            || s.active_plates.len() != self.plates.len()
            || s.plate_reference_global_kinematics.len() != self.plates.len()
            || s.reference_global_kinematics.len() != self.members.len()
            || s.node_displacement_m.len() != self.positions.len()
            || s.node_rotation_rad.len() != self.positions.len()
            || !s
                .reference_global_kinematics
                .iter()
                .flatten()
                .chain(s.plate_reference_global_kinematics.iter().flatten())
                .chain(s.node_displacement_m.iter().flatten())
                .chain(s.node_rotation_rad.iter().flatten())
                .all(|x| x.is_finite())
        {
            return Err(Error::Input("invalid committed beam state".into()));
        }
        Ok(())
    }
    pub fn change_members(
        &self,
        s: &State,
        soil_u: &[[f64; 3]],
        names: &[String],
        install: bool,
        cad_new_nodes: bool,
    ) -> Result<State> {
        self.validate_state(s)?;
        if names.is_empty() || names.iter().collect::<BTreeSet<_>>().len() != names.len() {
            return Err(Error::Input("unique beam member selection required".into()));
        }
        let mut next = s.clone();
        let mut old_nodes = vec![false; self.positions.len()];
        for (i, m) in self.members.iter().enumerate() {
            if s.active_members[i] {
                for n in m.nodes {
                    old_nodes[n] = true;
                }
            }
        }
        for (i, p) in self.plates.iter().enumerate() {
            if s.active_plates[i] {
                for n in p.nodes {
                    old_nodes[n] = true;
                }
            }
        }
        for name in names {
            if let Some(i) = self.configuration.plates.iter().position(|p| &p.id == name) {
                if s.active_plates[i] == install {
                    return Err(Error::Input(
                        "plate is already in requested active state".into(),
                    ));
                }
                next.active_plates[i] = install;
                if install {
                    for (end, n) in self.plates[i].nodes.iter().copied().enumerate() {
                        if let Some(soil) = self.attachments[n] {
                            next.node_displacement_m[n] = soil_u[soil];
                        } else if cad_new_nodes && !old_nodes[n] {
                            next.node_displacement_m[n] = [0.; 3];
                        }
                        if cad_new_nodes && !old_nodes[n] {
                            next.node_rotation_rad[n] = [0.; 3];
                        }
                        for a in 0..3 {
                            next.plate_reference_global_kinematics[i][6 * end + a] =
                                next.node_displacement_m[n][a];
                            next.plate_reference_global_kinematics[i][6 * end + 3 + a] =
                                next.node_rotation_rad[n][a];
                        }
                    }
                }
                continue;
            }
            let i = self
                .configuration
                .members
                .iter()
                .position(|m| &m.id == name)
                .ok_or_else(|| Error::Input("unknown beam member operation".into()))?;
            if s.active_members[i] == install {
                return Err(Error::Input(
                    "beam is already in requested active state".into(),
                ));
            }
            next.active_members[i] = install;
            if install {
                for (end, n) in self.members[i].nodes.iter().copied().enumerate() {
                    if let Some(soil) = self.attachments[n] {
                        next.node_displacement_m[n] = soil_u[soil];
                    } else if cad_new_nodes && !old_nodes[n] {
                        next.node_displacement_m[n] = [0.; 3];
                    }
                    if cad_new_nodes && !old_nodes[n] {
                        next.node_rotation_rad[n] = [0.; 3];
                    }
                    for a in 0..3 {
                        next.reference_global_kinematics[i][6 * end + a] =
                            next.node_displacement_m[n][a];
                        next.reference_global_kinematics[i][6 * end + 3 + a] =
                            next.node_rotation_rad[n][a];
                    }
                }
            }
        }
        Ok(next)
    }
    fn compiled(&self, state: &State) -> Result<Rc<Compiled>> {
        let active = &state.active_members;
        let key: Vec<bool> = active.iter().chain(&state.active_plates).copied().collect();
        if let Some(c) = self.cache.borrow().get(&key) {
            return Ok(c.clone());
        }
        let n = self.positions.len();
        let nd = 6 * n;
        let mut active_nodes = vec![false; n];
        let mut k = DMatrix::zeros(nd, nd);
        for (i, m) in self.members.iter().enumerate() {
            if active[i] {
                for node in m.nodes {
                    active_nodes[node] = true;
                }
                for a in 0..12 {
                    for b in 0..12 {
                        k[(6 * m.nodes[a / 6] + a % 6, 6 * m.nodes[b / 6] + b % 6)] +=
                            m.kernel.stiffness[(a, b)];
                    }
                }
            }
        }
        for (i, p) in self.plates.iter().enumerate() {
            if state.active_plates[i] {
                for node in p.nodes {
                    active_nodes[node] = true;
                }
                for a in 0..18 {
                    for b in 0..18 {
                        k[(6 * p.nodes[a / 6] + a % 6, 6 * p.nodes[b / 6] + b % 6)] +=
                            p.kernel.stiffness[(a, b)];
                    }
                }
            }
        }
        let mesh_nodes: BTreeSet<_> = self
            .attachments
            .iter()
            .enumerate()
            .filter(|(i, _)| active_nodes[*i])
            .filter_map(|(_, a)| *a)
            .collect();
        let interface_dofs: Vec<_> = mesh_nodes
            .iter()
            .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        let mut t = DMatrix::zeros(nd, interface_dofs.len());
        let mut offset = DVector::zeros(nd);
        let mut columns = Vec::<DVector<f64>>::new();
        let mut normals = Vec::new();
        for i in 0..n {
            if !active_nodes[i] {
                continue;
            }
            for part in 0..2 {
                if part == 0
                    && let Some(node) = self.attachments[i]
                {
                    for a in 0..3 {
                        let col = interface_dofs
                            .iter()
                            .position(|g| *g == 3 * node + a)
                            .unwrap();
                        t[(6 * i + a, col)] = 1.;
                    }
                    continue;
                }
                let constraints = if part == 0 {
                    &self.configuration.translation_constraints
                } else {
                    &self.configuration.rotation_constraints
                };
                let mut eqs = Vec::new();
                for c in constraints
                    .iter()
                    .filter(|c| c.node_id == self.configuration.nodes[i].id)
                {
                    let d = Vector3::from(c.direction);
                    eqs.push((d / d.norm(), c.value / d.norm()));
                }
                let (u0, free, fixed) = basis(&eqs)?;
                for a in 0..3 {
                    offset[6 * i + 3 * part + a] = u0[a];
                }
                for normal in fixed {
                    normals.push((6 * i + 3 * part, normal));
                }
                for b in free {
                    let mut v = DVector::zeros(nd);
                    for a in 0..3 {
                        v[6 * i + 3 * part + a] = b[a]
                            / if part == 1 {
                                self.configuration.rotation_length_scale_m
                            } else {
                                1.
                            };
                    }
                    columns.push(v);
                }
            }
        }
        if columns.len() > self.configuration.maximum_condensed_dofs {
            return Err(Error::Input("beam interior DOF budget exceeded".into()));
        }
        let q = DMatrix::from_fn(nd, columns.len(), |i, j| columns[j][i]);
        let kqq = q.transpose() * &k * &q;
        let kqt = q.transpose() * &k * &t;
        let chol = if columns.is_empty() {
            None
        } else {
            let c = kqq.clone().cholesky().ok_or_else(|| {
                Error::Solve(
                    "unconstrained beam mechanism; inspect torsion/releases/supports".into(),
                )
            })?;
            let max = (0..kqq.nrows())
                .map(|i| kqq[(i, i)].abs())
                .fold(0., f64::max);
            if c.l().diagonal().iter().any(|d| d * d <= 1e-12 * max) {
                return Err(Error::Solve(
                    "singular/ill-conditioned beam interior; no stabilization applied".into(),
                ));
            }
            Some(c)
        };
        let kc = t.transpose() * &k * &t
            - chol.as_ref().map_or_else(
                || DMatrix::zeros(t.ncols(), t.ncols()),
                |c| kqt.transpose() * c.solve(&kqt),
            );
        let value = Rc::new(Compiled {
            active_nodes,
            k,
            t,
            q,
            offset,
            chol,
            kqq,
            kc,
            interface_dofs,
            normals,
        });
        self.cache.borrow_mut().insert(key, value.clone());
        Ok(value)
    }
    fn loads(&self, s: &State, factor: f64) -> Result<AssembledLoads> {
        let c = self.compiled(s)?;
        let mut loads = DVector::zeros(6 * self.positions.len());
        let mut rest = DVector::zeros(loads.len());
        let mut lines = vec![[0.; 3]; self.members.len()];
        for (i, m) in self.members.iter().enumerate() {
            if !s.active_members[i] {
                continue;
            }
            let d = &self.configuration.members[i];
            let section = &self.configuration.sections[m.section];
            let mut line = std::array::from_fn(|a| {
                d.constant_line_load_kn_m[a] + factor * d.boundary_line_load_kn_m[a]
            });
            line[2] -= section.unit_weight_kn_m3 * section.area_m2;
            lines[i] = line;
            let local_load = m.kernel.uniform_load(line);
            let global_load = m.kernel.transform.transpose() * m.kernel.condense_load(&local_load);
            let reference = V12::from_row_slice(&s.reference_global_kinematics[i]);
            let mut prestress = V12::zeros();
            prestress[0] = -d.initial_axial_force_kn;
            prestress[6] = d.initial_axial_force_kn;
            let constant = -m.kernel.stiffness * reference
                + m.kernel.transform.transpose() * m.kernel.condense_load(&prestress);
            for a in 0..12 {
                loads[6 * m.nodes[a / 6] + a % 6] += global_load[a];
                rest[6 * m.nodes[a / 6] + a % 6] += constant[a];
            }
        }
        for (i, p) in self.plates.iter().enumerate() {
            if !s.active_plates[i] {
                continue;
            }
            let config = &self.configuration.plates[i];
            let section = &self.configuration.plate_sections[p.section];
            let mut traction = std::array::from_fn(|a| {
                config.constant_traction_kpa[a] + factor * config.boundary_traction_kpa[a]
            });
            traction[2] -= section.unit_weight_kn_m3 * section.thickness_m;
            let load = p.kernel.uniform_load(traction);
            let reference =
                crate::plates::V18::from_row_slice(&s.plate_reference_global_kinematics[i]);
            let constant = -p.kernel.stiffness * reference;
            for a in 0..18 {
                loads[6 * p.nodes[a / 6] + a % 6] += load[a];
                rest[6 * p.nodes[a / 6] + a % 6] += constant[a];
            }
        }
        for load in &self.configuration.nodal_loads {
            let n = self
                .configuration
                .nodes
                .iter()
                .position(|n| n.id == load.node_id)
                .unwrap();
            for a in 0..3 {
                let f = load.constant_force_kn[a] + factor * load.boundary_force_kn[a];
                let m = load.constant_moment_kn_m[a] + factor * load.boundary_moment_kn_m[a];
                if !c.active_nodes[n] && (f != 0. || m != 0.) {
                    return Err(Error::Input("load on inactive beam node".into()));
                }
                loads[6 * n + a] += f;
                loads[6 * n + 3 + a] += m;
            }
        }
        Ok((loads, rest, lines))
    }
    pub fn attached_load(&self, s: &State, factor: f64) -> Result<Vec<(usize, f64)>> {
        self.validate_state(s)?;
        let c = self.compiled(s)?;
        let (load, _, _) = self.loads(s, factor)?;
        let f = c.t.transpose() * load;
        Ok(c.interface_dofs
            .iter()
            .copied()
            .zip(f.iter().copied())
            .collect())
    }
    pub fn evaluate(
        &self,
        s: &State,
        soil_u: &[[f64; 3]],
        soil_active_nodes: &[bool],
        factor: f64,
    ) -> Result<Evaluation> {
        self.validate_state(s)?;
        let c = self.compiled(s)?;
        if c.interface_dofs.iter().any(|g| !soil_active_nodes[g / 3]) {
            return Err(Error::Input("active beam attachment lost its soil node; remove or explicitly detach the member before soil excavation".into()));
        }
        let (loads, rest, lines) = self.loads(s, factor)?;
        let u = DVector::from_iterator(
            c.interface_dofs.len(),
            c.interface_dofs.iter().map(|g| soil_u[g / 3][g % 3]),
        );
        let base = &c.t * &u + &c.offset;
        let rhs = c.q.transpose() * (&loads - &rest - &c.k * &base);
        let interior = c
            .chol
            .as_ref()
            .map_or_else(|| DVector::zeros(0), |chol| chol.solve(&rhs));
        let relative = if interior.is_empty() {
            0.
        } else {
            (&c.kqq * &interior - &rhs).norm() / rhs.norm().max(1.)
        };
        if !relative.is_finite() || relative > 1e-9 {
            return Err(Error::Solve(
                "beam interior true residual exceeded tolerance".into(),
            ));
        }
        let d = base + &c.q * interior;
        let internal = &c.k * &d + rest;
        let residual = &internal - &loads;
        let relative =
            (c.q.transpose() * &residual).norm() / (c.q.transpose() * &loads).norm().max(1.);
        if !relative.is_finite() || relative > 1e-9 {
            return Err(Error::Solve(
                "beam physical interior equilibrium residual exceeded tolerance".into(),
            ));
        }
        let interface_internal_kn = c.t.transpose() * &internal;
        let mut state = s.clone();
        for i in 0..self.positions.len() {
            if c.active_nodes[i] {
                for a in 0..3 {
                    state.node_displacement_m[i][a] = d[6 * i + a];
                    state.node_rotation_rad[i][a] = d[6 * i + 3 + a];
                }
            }
        }
        let mut external_f = Vector3::zeros();
        let mut external_m = Vector3::zeros();
        let mut attached_f = Vector3::zeros();
        let mut attached_m = Vector3::zeros();
        let mut support_f = Vector3::zeros();
        let mut support_m = Vector3::zeros();
        for i in 0..self.positions.len() {
            let f = Vector3::new(loads[6 * i], loads[6 * i + 1], loads[6 * i + 2]);
            let m = Vector3::new(loads[6 * i + 3], loads[6 * i + 4], loads[6 * i + 5]);
            let x = Vector3::from(self.positions[i]);
            external_f += f;
            external_m += x.cross(&f) + m;
            if self.attachments[i].is_some() {
                attached_f += f;
                attached_m += x.cross(&f);
            }
        }
        for (start, n) in &c.normals {
            let r = Vector3::new(residual[*start], residual[*start + 1], residual[*start + 2]);
            let projected = n * n.dot(&r);
            if start % 6 == 0 {
                support_f += projected;
                support_m += Vector3::from(self.positions[start / 6]).cross(&projected);
            } else {
                support_m += projected;
            }
        }
        let mut results = Vec::new();
        for (i, m) in self.members.iter().enumerate() {
            if !s.active_members[i] {
                continue;
            }
            let global = V12::from_fn(|a, _| {
                d[6 * m.nodes[a / 6] + a % 6] - s.reference_global_kinematics[i][a]
            });
            let (local, action) = m.kernel.recover(
                &global,
                lines[i],
                self.configuration.members[i].initial_axial_force_kn,
            );
            let q = m.kernel.axes * Vector3::from(lines[i]);
            let samples = std::array::from_fn(|j| {
                let x = m.kernel.length * j as f64 / 2.;
                [
                    -action[0] - q[0] * x,
                    -action[1] - q[1] * x,
                    -action[2] - q[2] * x,
                    -action[3],
                    -action[4] - x * action[2] - q[2] * x * x / 2.,
                    -action[5] + x * action[1] + q[1] * x * x / 2.,
                ]
            });
            results.push(MemberResult {
                id: self.configuration.members[i].id.clone(),
                local_axes: std::array::from_fn(|a| std::array::from_fn(|b| m.kernel.axes[(a, b)])),
                length_m: m.kernel.length,
                local_end_actions: std::array::from_fn(|a| action[a]),
                local_end_relative_rotation_rad: std::array::from_fn(|end| {
                    std::array::from_fn(|a| local[6 * end + 3 + a])
                }),
                section_forces_at_start_mid_end: samples,
            });
        }
        let plates = self
            .plates
            .iter()
            .enumerate()
            .filter(|(i, _)| s.active_plates[*i])
            .map(|(i, p)| {
                let relative = crate::plates::V18::from_fn(|a, _| {
                    d[6 * p.nodes[a / 6] + a % 6] - s.plate_reference_global_kinematics[i][a]
                });
                p.kernel
                    .recover(self.configuration.plates[i].id.clone(), &relative)
            })
            .collect();
        let report = Report {
            active_members: s.active_members.iter().filter(|a| **a).count(),
            condensed_dofs: c.q.ncols(),
            interface_translational_dofs: c.t.ncols(),
            true_relative_interior_residual: relative,
            node_positions_m: self.positions.clone(),
            node_displacement_m: state.node_displacement_m.clone(),
            node_rotation_rad: state.node_rotation_rad.clone(),
            node_active: c.active_nodes.clone(),
            external_force_kn: external_f.into(),
            external_moment_origin_kn_m: external_m.into(),
            interface_external_force_kn: attached_f.into(),
            interface_external_moment_origin_kn_m: attached_m.into(),
            support_force_kn: support_f.into(),
            support_moment_origin_kn_m: support_m.into(),
            members: results,
            plates,
        };
        Ok(Evaluation {
            state,
            report,
            interface_dofs: c.interface_dofs.clone(),
            interface_internal_kn,
            interface_tangent_kn_m: c.kc.clone(),
        })
    }
}
