//! Affine displacement constraints and work-conjugate nullspace reduction.
use crate::{Error, Result, elements, mesh::Mesh};
use nalgebra::{DMatrix, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Constraint {
    pub node_set: String,
    pub direction: [f64; 3],
    pub value_m: f64,
    #[serde(default)]
    pub gradient: [f64; 3],
    #[serde(default)]
    pub reference_point_m: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NodalLoad {
    pub node_set: String,
    /// Total force, shared equally by the selected nodes.
    pub total_force_kn: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum SurfaceLoad {
    Traction {
        face_set: String,
        traction_kpa: [f64; 3],
        #[serde(default)]
        gradient_kpa_per_m: [[f64; 3]; 3],
        #[serde(default)]
        reference_point_m: [f64; 3],
    },
    Pressure {
        face_set: String,
        pressure_kpa: f64,
        #[serde(default)]
        gradient_kpa_per_m: [f64; 3],
        #[serde(default)]
        reference_point_m: [f64; 3],
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Conditions {
    pub constraints: Vec<Constraint>,
    #[serde(default)]
    pub nodal_loads: Vec<NodalLoad>,
    #[serde(default)]
    pub surface_loads: Vec<SurfaceLoad>,
}

/// Loads that remain present during initialization and every construction
/// stage, independently of the proportional boundary-stage multiplier.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Loads {
    #[serde(default)]
    pub nodal_loads: Vec<NodalLoad>,
    #[serde(default)]
    pub surface_loads: Vec<SurfaceLoad>,
}
impl Loads {
    pub fn is_empty(&self) -> bool {
        self.nodal_loads.is_empty() && self.surface_loads.is_empty()
    }
    pub fn conditions(&self) -> Conditions {
        Conditions {
            constraints: vec![],
            nodal_loads: self.nodal_loads.clone(),
            surface_loads: self.surface_loads.clone(),
        }
    }
    pub fn on_active(&self, mesh: &Mesh, active: &[bool]) -> Result<Vec<f64>> {
        self.conditions().loads_impl(mesh, active, true)
    }
}

#[derive(Debug, Clone)]
pub struct Transform {
    /// u = offset + T q. Each physical component has at most three terms.
    pub rows: Vec<Vec<(usize, f64)>>,
    pub offset: Vec<f64>,
    pub count: usize,
    pub active_nodes: Vec<bool>,
    equations: Vec<(usize, Vector3<f64>, f64)>,
}
impl Transform {
    pub fn new(
        mesh: &Mesh,
        constraints: &[Constraint],
        active_elements: Vec<bool>,
    ) -> Result<Self> {
        Self::with_extra(mesh, constraints, active_elements, Vec::new())
    }
    pub(crate) fn with_extra(
        mesh: &Mesh,
        constraints: &[Constraint],
        active_elements: Vec<bool>,
        extra: Vec<(usize, Vector3<f64>, f64)>,
    ) -> Result<Self> {
        Self::build(mesh, constraints, active_elements, extra, false)
    }
    /// Interface-supported components are checked from their current contact
    /// state after equilibrium, rather than treating every component as free.
    pub(crate) fn with_interfaces(
        mesh: &Mesh,
        constraints: &[Constraint],
        active_elements: Vec<bool>,
        extra: Vec<(usize, Vector3<f64>, f64)>,
    ) -> Result<Self> {
        Self::build(mesh, constraints, active_elements, extra, true)
    }
    fn build(
        mesh: &Mesh,
        constraints: &[Constraint],
        active_elements: Vec<bool>,
        extra: Vec<(usize, Vector3<f64>, f64)>,
        interfaces: bool,
    ) -> Result<Self> {
        if active_elements.len() != mesh.tetrahedra.len() {
            return Err(Error::Input(
                "constraint active-node dimension mismatch".into(),
            ));
        }
        let mut active = vec![false; mesh.nodes.len()];
        for (ei, &enabled) in active_elements.iter().enumerate() {
            if enabled {
                for n in mesh.cell_nodes(ei) {
                    active[n] = true;
                }
            }
        }
        let mut by_node = vec![Vec::new(); mesh.nodes.len()];
        let mut equations = extra;
        for c in constraints {
            if !c
                .direction
                .iter()
                .chain(&c.gradient)
                .chain(&c.reference_point_m)
                .chain([&c.value_m])
                .all(|x| x.is_finite())
            {
                return Err(Error::Input("non-finite displacement constraint".into()));
            }
            let direction = Vector3::from_row_slice(&c.direction);
            let norm = direction.norm();
            if !norm.is_finite() || norm < 1e-12 {
                return Err(Error::Input(
                    "zero displacement-constraint direction".into(),
                ));
            }
            let nodes = mesh
                .node_sets
                .get(&c.node_set)
                .filter(|n| !n.is_empty())
                .ok_or_else(|| {
                    Error::Input(format!("unknown/empty constraint node set {}", c.node_set))
                })?;
            for &node in nodes {
                if !active[node] {
                    continue;
                }
                let value = (c.value_m
                    + (0..3)
                        .map(|a| c.gradient[a] * (mesh.nodes[node][a] - c.reference_point_m[a]))
                        .sum::<f64>())
                    / norm;
                equations.push((node, direction / norm, value));
            }
        }
        for &(node, direction, value) in &equations {
            if node >= mesh.nodes.len()
                || !value.is_finite()
                || !direction.iter().all(|v| v.is_finite())
            {
                return Err(Error::Input(
                    "invalid expanded displacement constraint".into(),
                ));
            }
            if active[node] {
                by_node[node].push((direction, value));
            }
        }
        let mut rows = vec![vec![]; 3 * mesh.nodes.len()];
        let mut offset = vec![0.; rows.len()];
        let mut count = 0;
        for (node, eqs) in by_node.iter().enumerate() {
            if !active[node] {
                continue;
            }
            let mut normals = Vec::<Vector3<f64>>::new();
            let mut u0 = Vector3::zeros();
            for (a, value) in eqs {
                let mut normal = *a;
                for q in &normals {
                    normal -= q * normal.dot(q);
                }
                let norm = normal.norm();
                let residual = value - a.dot(&u0);
                if norm < 1e-10 {
                    if residual.abs() > 1e-11 * (1. + value.abs()) {
                        return Err(Error::Input("conflicting prescribed displacements".into()));
                    }
                } else {
                    let q = normal / norm;
                    u0 += q * (residual / norm);
                    normals.push(q);
                }
            }
            for a in 0..3 {
                offset[3 * node + a] = u0[a];
            }
            let mut basis = Vec::<Vector3<f64>>::new();
            for a in 0..3 {
                let mut v = Vector3::zeros();
                v[a] = 1.;
                for _ in 0..2 {
                    for q in normals.iter().chain(&basis) {
                        v -= q * v.dot(q);
                    }
                }
                if v.norm() > 1e-10 {
                    basis.push(v.normalize());
                }
            }
            if basis.len() + normals.len() != 3 {
                return Err(Error::Input(
                    "constraint nullspace construction failed".into(),
                ));
            }
            for v in basis {
                for a in 0..3 {
                    if v[a].abs() > 1e-15 {
                        rows[3 * node + a].push((count, v[a]));
                    }
                }
                count += 1;
            }
        }
        let result = Self {
            rows,
            offset,
            count,
            active_nodes: active,
            equations,
        };
        if !interfaces {
            result.check_rigid_modes(mesh, &active_elements)?;
        }
        Ok(result)
    }
    fn check_rigid_modes(&self, mesh: &Mesh, active: &[bool]) -> Result<()> {
        fn root(parents: &[usize], mut i: usize) -> usize {
            while parents[i] != i {
                i = parents[i];
            }
            i
        }
        let mut parents: Vec<_> = (0..mesh.nodes.len()).collect();
        for (ei, &enabled) in active.iter().enumerate() {
            if enabled {
                let nodes = mesh.cell_nodes(ei);
                for &n in nodes.iter().skip(1) {
                    let a = root(&parents, nodes[0]);
                    let b = root(&parents, n);
                    parents[b] = a;
                }
            }
        }
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for n in 0..mesh.nodes.len() {
            if self.active_nodes[n] {
                groups.entry(root(&parents, n)).or_default().push(n);
            }
        }
        if groups.is_empty() {
            return Err(Error::Input("no active mechanics nodes".into()));
        }
        for (component, nodes) in groups {
            let equations: Vec<_> = self
                .equations
                .iter()
                .filter(|(n, _, _)| root(&parents, *n) == component)
                .collect();
            if equations.len() < 6 {
                return Err(Error::Input(
                    "mechanics component leaves unconstrained rigid-body modes".into(),
                ));
            }
            let center = nodes
                .iter()
                .map(|&n| Vector3::from_row_slice(&mesh.nodes[n]))
                .sum::<Vector3<f64>>()
                / nodes.len() as f64;
            let radius = nodes
                .iter()
                .map(|&n| (Vector3::from_row_slice(&mesh.nodes[n]) - center).norm())
                .fold(0., f64::max)
                .max(1e-30);
            let mut modes = DMatrix::zeros(equations.len(), 6);
            for (i, &&(node, direction, _)) in equations.iter().enumerate() {
                let x = (Vector3::from_row_slice(&mesh.nodes[node]) - center) / radius;
                for a in 0..3 {
                    modes[(i, a)] = direction[a];
                    let mut axis = Vector3::zeros();
                    axis[a] = 1.;
                    modes[(i, a + 3)] = direction.dot(&axis.cross(&x));
                }
            }
            let sv = modes.svd(false, false).singular_values;
            if sv.len() < 6 || sv.iter().filter(|&&v| v > 1e-10 * sv.max().max(1.)).count() != 6 {
                return Err(Error::Input(
                    "mechanics component leaves unconstrained rigid-body modes".into(),
                ));
            }
        }
        Ok(())
    }
    pub fn reduce(&self, physical: &[f64]) -> Vec<f64> {
        let mut result = vec![0.; self.count];
        for (g, terms) in self.rows.iter().enumerate() {
            for &(r, c) in terms {
                result[r] += c * physical[g];
            }
        }
        result
    }
    pub fn expand(&self, reduced: &[f64]) -> Vec<f64> {
        self.rows
            .iter()
            .map(|terms| terms.iter().map(|&(r, c)| c * reduced[r]).sum())
            .collect()
    }
    pub fn maximum_constraint_error(&self, u: &[f64]) -> f64 {
        self.equations
            .iter()
            .map(|(node, d, value)| {
                (d.dot(&Vector3::from_row_slice(&u[3 * node..3 * node + 3])) - value).abs()
            })
            .fold(0., f64::max)
    }
    pub fn reactions(&self, internal_minus_external: &[f64]) -> Vec<f64> {
        let projected = self.expand(&self.reduce(internal_minus_external));
        internal_minus_external
            .iter()
            .zip(projected)
            .enumerate()
            .map(|(g, (r, p))| if self.active_nodes[g / 3] { r - p } else { 0. })
            .collect()
    }
}

pub fn surface_owners(mesh: &Mesh) -> BTreeMap<[usize; 3], Vec<usize>> {
    let mut faces = BTreeMap::new();
    for (ei, tet) in mesh.tetrahedra.iter().enumerate() {
        for local in elements::FACES {
            let mut face = local.map(|i| tet[i]);
            face.sort_unstable();
            faces.entry(face).or_insert_with(Vec::new).push(ei);
        }
    }
    faces
}
impl Conditions {
    pub fn validate_definition(&self, mesh: &Mesh) -> Result<()> {
        for c in &self.constraints {
            let direction = Vector3::from_row_slice(&c.direction);
            if mesh.node_sets.get(&c.node_set).is_none_or(|s| s.is_empty())
                || !c
                    .direction
                    .iter()
                    .chain(&c.gradient)
                    .chain(&c.reference_point_m)
                    .chain([&c.value_m])
                    .all(|v| v.is_finite())
                || !direction.norm().is_finite()
                || direction.norm() < 1e-12
            {
                return Err(Error::Input("invalid named displacement constraint".into()));
            }
        }
        for load in &self.nodal_loads {
            if mesh
                .node_sets
                .get(&load.node_set)
                .is_none_or(|s| s.is_empty())
                || !load.total_force_kn.iter().all(|x| x.is_finite())
            {
                return Err(Error::Input("invalid named nodal load".into()));
            }
        }
        for load in &self.surface_loads {
            let (name, finite) = match load {
                SurfaceLoad::Pressure {
                    face_set,
                    pressure_kpa,
                    gradient_kpa_per_m,
                    reference_point_m,
                } => (
                    face_set,
                    gradient_kpa_per_m
                        .iter()
                        .chain(reference_point_m)
                        .chain([pressure_kpa])
                        .all(|x| x.is_finite()),
                ),
                SurfaceLoad::Traction {
                    face_set,
                    traction_kpa,
                    gradient_kpa_per_m,
                    reference_point_m,
                } => (
                    face_set,
                    traction_kpa
                        .iter()
                        .chain(gradient_kpa_per_m.iter().flatten())
                        .chain(reference_point_m)
                        .all(|x| x.is_finite()),
                ),
            };
            if !finite || mesh.face_sets.get(name).is_none_or(|s| s.is_empty()) {
                return Err(Error::Input("invalid named surface load".into()));
            }
        }
        Ok(())
    }
    pub fn loads(&self, mesh: &Mesh, active: &[bool]) -> Result<Vec<f64>> {
        self.loads_impl(mesh, active, false)
    }
    fn loads_impl(
        &self,
        mesh: &Mesh,
        active: &[bool],
        skip_removed_surfaces: bool,
    ) -> Result<Vec<f64>> {
        self.validate_definition(mesh)?;
        if active.len() != mesh.tetrahedra.len() {
            return Err(Error::Input(
                "load active-element dimension mismatch".into(),
            ));
        }
        let mut used = vec![false; mesh.nodes.len()];
        for (i, &enabled) in active.iter().enumerate() {
            if enabled {
                for n in mesh.cell_nodes(i) {
                    used[n] = true;
                }
            }
        }
        let mut f = vec![0.; mesh.nodes.len() * 3];
        for load in &self.nodal_loads {
            let nodes = mesh
                .node_sets
                .get(&load.node_set)
                .filter(|n| !n.is_empty())
                .ok_or_else(|| Error::Input("unknown/empty nodal load set".into()))?;
            if !load.total_force_kn.iter().all(|v| v.is_finite()) || nodes.iter().any(|&n| !used[n])
            {
                return Err(Error::Input(
                    "invalid nodal load or load on inactive node".into(),
                ));
            }
            for &node in nodes {
                for a in 0..3 {
                    f[3 * node + a] += load.total_force_kn[a] / nodes.len() as f64;
                }
            }
        }
        let owners = surface_owners(mesh);
        for load in &self.surface_loads {
            let name = match load {
                SurfaceLoad::Traction { face_set, .. } | SurfaceLoad::Pressure { face_set, .. } => {
                    face_set
                }
            };
            let set = mesh
                .face_sets
                .get(name)
                .filter(|f| !f.is_empty())
                .ok_or_else(|| Error::Input(format!("unknown/empty load face set {name}")))?;
            for face in set {
                let mut key = *face;
                key.sort_unstable();
                let enabled: Vec<_> = owners
                    .get(&key)
                    .ok_or_else(|| Error::Input("unknown loaded face".into()))?
                    .iter()
                    .copied()
                    .filter(|&i| active[i])
                    .collect();
                if skip_removed_surfaces && enabled.is_empty() {
                    continue;
                }
                if enabled.len() != 1 {
                    return Err(Error::Input(
                        "loaded face must be on the active exterior".into(),
                    ));
                }
                let nodes = mesh.face_nodes(enabled[0], face)?;
                let points: Vec<_> = nodes.iter().map(|&i| mesh.nodes[i]).collect();
                for q in elements::face_points(&points)? {
                    let area = Vector3::from_row_slice(&q.area_vector);
                    let mut force = Vector3::zeros();
                    match load {
                        SurfaceLoad::Pressure {
                            pressure_kpa,
                            gradient_kpa_per_m,
                            reference_point_m,
                            ..
                        } => {
                            let p = pressure_kpa
                                + (0..3)
                                    .map(|a| {
                                        gradient_kpa_per_m[a]
                                            * (q.position[a] - reference_point_m[a])
                                    })
                                    .sum::<f64>();
                            force -= area * p;
                        }
                        SurfaceLoad::Traction {
                            traction_kpa,
                            gradient_kpa_per_m,
                            reference_point_m,
                            ..
                        } => {
                            for a in 0..3 {
                                force[a] = (traction_kpa[a]
                                    + (0..3)
                                        .map(|b| {
                                            gradient_kpa_per_m[a][b]
                                                * (q.position[b] - reference_point_m[b])
                                        })
                                        .sum::<f64>())
                                    * area.norm();
                            }
                        }
                    }
                    if !force.iter().all(|v| v.is_finite()) {
                        return Err(Error::Input("non-finite surface load".into()));
                    }
                    for (j, &node) in nodes.iter().enumerate() {
                        for a in 0..3 {
                            f[3 * node + a] += q.shape[j] * force[a];
                        }
                    }
                }
            }
        }
        Ok(f)
    }
}
