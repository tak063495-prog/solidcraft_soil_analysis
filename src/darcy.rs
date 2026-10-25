//! Steady saturated hydraulic-head FEM. q=-K grad(H), div(q)=source.
//! Positive outward Neumann flux. Prescribed-head reaction outflow is -r.
//! No guessed free surface, pressure clipping, diagonal shifts or fallback.
use crate::{
    Error, Result, boundary, elements,
    linalg::{self, CsrMatrix, CsrPattern, LinearReport, LinearSolver},
    mesh::Mesh,
};
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub id: u32,
    pub hydraulic_conductivity_m_s: [[f64; 3]; 3],
    #[serde(default)]
    pub volume_source_per_s: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Head {
    pub face_set: String,
    pub head_m: f64,
    #[serde(default)]
    pub gradient: [f64; 3],
    #[serde(default)]
    pub reference_point_m: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Flux {
    Outward {
        face_set: String,
        flux_m_s: f64,
    },
    Vector {
        face_set: String,
        flux_m_s: [f64; 3],
        #[serde(default)]
        gradient_per_s: [[f64; 3]; 3],
        #[serde(default)]
        reference_point_m: [f64; 3],
    },
}
impl Flux {
    fn face_set(&self) -> &str {
        match self {
            Self::Outward { face_set, .. } | Self::Vector { face_set, .. } => face_set,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub schema: String,
    pub units: String,
    pub materials: Vec<Material>,
    pub heads: Vec<Head>,
    #[serde(default)]
    pub fluxes: Vec<Flux>,
    #[serde(default)]
    pub settings: Settings,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Settings {
    pub relative_linear_tolerance: f64,
    pub pressure_head_tolerance_m: f64,
    pub maximum_free_heads: usize,
    pub linear_solver: LinearSolver,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            relative_linear_tolerance: 1e-11,
            pressure_head_tolerance_m: 1e-9,
            maximum_free_heads: 300_000,
            linear_solver: LinearSolver::SparseGmresIlu0 {
                options: Default::default(),
            },
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct HeadField {
    /// Primary values retain small head differences against a large datum.
    pub relative_head_m: Vec<f64>,
    pub datum_head_m: Vec<f64>,
    pub head_m: Vec<f64>,
    pub defined_nodes: Vec<bool>,
}
impl HeadField {
    pub fn validate(&self, mesh: &Mesh, active: &[bool]) -> Result<()> {
        mesh.validate()?;
        if self.head_m.len() != mesh.nodes.len()
            || self.defined_nodes.len() != mesh.nodes.len()
            || self.relative_head_m.len() != mesh.nodes.len()
            || self.datum_head_m.len() != mesh.nodes.len()
            || !self
                .head_m
                .iter()
                .chain(&self.relative_head_m)
                .chain(&self.datum_head_m)
                .all(|v| v.is_finite())
            || active.len() != mesh.tetrahedra.len()
        {
            return Err(Error::Input(
                "invalid hydraulic head field dimensions or values".into(),
            ));
        }
        for (i, live) in active.iter().enumerate() {
            if *live && mesh.cell_nodes(i).iter().any(|&n| !self.defined_nodes[n]) {
                return Err(Error::Input("active soil references an undefined hydraulic head; solve flow on this geometry".into()));
            }
            let nodes = mesh.cell_nodes(i);
            if *live
                && nodes.iter().any(|&n| {
                    self.datum_head_m[n] != self.datum_head_m[nodes[0]]
                        || self.head_m[n] != self.datum_head_m[n] + self.relative_head_m[n]
                })
            {
                return Err(Error::Input(
                    "inconsistent hydraulic head datum or absolute field".into(),
                ));
            }
        }
        Ok(())
    }
    pub fn pressure_head(&self, nodes: &[usize], shape: &[f64], z: f64) -> Result<f64> {
        if nodes.is_empty()
            || nodes.len() != shape.len()
            || self.defined_nodes.len() != self.head_m.len()
            || self.relative_head_m.len() != self.head_m.len()
            || self.datum_head_m.len() != self.head_m.len()
            || !z.is_finite()
            || !shape.iter().all(|v| v.is_finite())
            || nodes.iter().any(|&n| {
                n >= self.head_m.len()
                    || !self.defined_nodes[n]
                    || !self.head_m[n].is_finite()
                    || !self.relative_head_m[n].is_finite()
                    || !self.datum_head_m[n].is_finite()
            })
        {
            return Err(Error::Input("undefined hydraulic interpolation".into()));
        }
        if nodes
            .iter()
            .any(|&n| self.datum_head_m[n] != self.datum_head_m[nodes[0]])
            || (shape.iter().sum::<f64>() - 1.).abs() > 1e-10
        {
            return Err(Error::Input(
                "interpolation crosses hydraulic datums or has invalid shape weights".into(),
            ));
        }
        let origin = self.relative_head_m[nodes[0]];
        let v = self.datum_head_m[nodes[0]] - z
            + origin
            + nodes
                .iter()
                .zip(shape)
                .map(|(&n, &w)| (self.relative_head_m[n] - origin) * w)
                .sum::<f64>();
        if !v.is_finite() {
            return Err(Error::Solve("pore head interpolation overflow".into()));
        }
        Ok(v)
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Point {
    pub element_index: usize,
    pub local_index: usize,
    pub position_m: [f64; 3],
    pub weight_m3: f64,
    pub head_m: f64,
    pub pressure_head_m: f64,
    pub head_gradient: [f64; 3],
    pub darcy_flux_m_s: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FaceFlux {
    pub corners: [usize; 3],
    pub owner: usize,
    pub boundary_kind: String,
    /// Direct integration of the discontinuous H1 element gradient trace.
    pub gradient_outflow_m3_s: f64,
    pub prescribed_outflow_m3_s: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Balance {
    pub reference_head_m: f64,
    pub relative_free_head_residual: f64,
    pub linear_solves: Vec<LinearReport>,
    pub absolute_source_m3_s: f64,
    pub absolute_neumann_flow_m3_s: f64,
    pub component: usize,
    pub elements: usize,
    pub source_m3_s: f64,
    pub neumann_outflow_m3_s: f64,
    pub dirichlet_reaction_outflow_m3_s: f64,
    pub relative_nodal_balance_error: f64,
    pub gradient_boundary_outflow_m3_s: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub refinement_solves: Vec<LinearReport>,
    pub schema: String,
    pub complete: bool,
    pub failure_reason: Option<String>,
    pub mesh_sha256: String,
    pub case_sha256: String,
    pub active_elements: Vec<bool>,
    pub material_ids: Vec<u32>,
    pub field: HeadField,
    pub points: Vec<Point>,
    pub boundary_fluxes: Vec<FaceFlux>,
    pub dirichlet_nodal_outflow_m3_s: Vec<Option<f64>>,
    pub balances: Vec<Balance>,
    pub relative_free_head_residual: f64,
    pub minimum_sampled_pressure_head_m: f64,
    pub bernstein_pressure_head_lower_bound_m: f64,
    /// Initial solve for a single region; per-region solves are in balances.
    pub linear: Option<LinearReport>,
}
fn norm(v: &[f64]) -> f64 {
    v.iter().fold(0., |a, b| a.hypot(*b))
}
fn gradients(q: &elements::VolumePoint) -> Vec<Vector3<f64>> {
    (0..q.shape.len())
        .map(|i| Vector3::new(q.b[(0, 3 * i)], q.b[(1, 3 * i + 1)], q.b[(2, 3 * i + 2)]))
        .collect()
}
pub(crate) fn tensor(m: &Material) -> Result<Matrix3<f64>> {
    let k = Matrix3::from_fn(|i, j| m.hydraulic_conductivity_m_s[i][j]);
    let scale = k.iter().map(|v| v.abs()).fold(0., f64::max);
    if !k.iter().all(|v| v.is_finite())
        || scale == 0.
        || scale > 1.
        || !m.volume_source_per_s.is_finite()
        || (k - k.transpose()).norm() > 1e-12 * scale
    {
        return Err(Error::Input(
            "hydraulic conductivity must be finite and symmetric".into(),
        ));
    }
    let eig = k.symmetric_eigen().eigenvalues;
    if eig.min() < 1e-20 || eig.max() > 1. || eig.max() / eig.min() > 1e12 {
        return Err(Error::Input(
            "hydraulic conductivity must be positive definite, 1e-20..1 m/s, condition <=1e12"
                .into(),
        ));
    }
    Ok((k + k.transpose()) * 0.5)
}
fn root(parent: &mut [usize], mut n: usize) -> usize {
    while parent[n] != n {
        parent[n] = parent[parent[n]];
        n = parent[n];
    }
    n
}
/// Numerical fields and audits must reproduce exactly; wall time is metadata.
pub fn same_solution(a: &Report, b: &Report) -> Result<bool> {
    let mut a = a.clone();
    let mut b = b.clone();
    if let Some(r) = &mut a.linear {
        r.elapsed_seconds = 0.;
    }
    if let Some(r) = &mut b.linear {
        r.elapsed_seconds = 0.;
    }
    for r in a
        .refinement_solves
        .iter_mut()
        .chain(&mut b.refinement_solves)
    {
        r.elapsed_seconds = 0.;
    }
    for balance in a.balances.iter_mut().chain(&mut b.balances) {
        for r in &mut balance.linear_solves {
            r.elapsed_seconds = 0.;
        }
    }
    Ok(serde_json::to_vec(&a)? == serde_json::to_vec(&b)?)
}

pub(crate) fn linear_solve(
    a: &CsrMatrix,
    b: &[f64],
    s: &Settings,
) -> Result<(Vec<f64>, LinearReport)> {
    match &s.linear_solver {
        LinearSolver::SparseGmresIlu0 { options } => Ok(linalg::gmres_ilu0(
            a,
            b,
            s.relative_linear_tolerance,
            options,
        )?),
        LinearSolver::DenseLu => {
            let dimension = b.len();
            let start = Instant::now();
            let mut dense = DMatrix::zeros(dimension, dimension);
            for i in 0..dimension {
                for j in a.row_ptr()[i]..a.row_ptr()[i + 1] {
                    dense[(i, a.columns()[j])] = a.values()[j];
                }
            }
            let x = dense
                .lu()
                .solve(&DVector::from_vec(b.to_vec()))
                .ok_or_else(|| Error::Solve("singular hydraulic matrix".into()))?;
            let value = a.matvec(x.as_slice())?;
            let residual = norm(&value.iter().zip(b).map(|(a, b)| a - b).collect::<Vec<_>>())
                / norm(b).max(1e-30);
            if !residual.is_finite() || residual > s.relative_linear_tolerance {
                return Err(Error::Solve(
                    "hydraulic dense solve true residual exceeded".into(),
                ));
            }
            Ok((
                x.as_slice().to_vec(),
                LinearReport {
                    method: "dense_lu".into(),
                    dimension,
                    stored_nonzeros: dimension * dimension,
                    iterations: 1,
                    restarts: 0,
                    true_relative_residual: residual,
                    matrix_storage_bytes: dimension * dimension * 8,
                    workspace_storage_bytes: (dimension * dimension + 2 * dimension) * 8,
                    elapsed_seconds: start.elapsed().as_secs_f64(),
                },
            ))
        }
    }
}
fn component_residual(
    k: &CsrMatrix,
    delta: &[f64],
    rhs: &[f64],
    nodes: &[usize],
    free: &[usize],
) -> Result<(f64, Vec<f64>)> {
    let action = k.matvec(delta)?;
    let correction: Vec<_> = free.iter().map(|&n| rhs[n] - action[n]).collect();
    let scale = norm(&nodes.iter().map(|&n| action[n]).collect::<Vec<_>>())
        .max(norm(&nodes.iter().map(|&n| rhs[n]).collect::<Vec<_>>()))
        .max(1e-30);
    Ok((norm(&correction) / scale, correction))
}

pub fn solve(
    mesh: &Mesh,
    case: &Case,
    active: &[bool],
    material_ids: Option<&[u32]>,
) -> Result<Report> {
    mesh.validate()?;
    case.settings.linear_solver.validate()?;
    let s = &case.settings;
    if case.schema != "soil-fem-darcy/1"
        || case.units != "m-s"
        || active.len() != mesh.tetrahedra.len()
        || !active.iter().any(|v| *v)
        || !s.relative_linear_tolerance.is_finite()
        || !(1e-13..=1e-6).contains(&s.relative_linear_tolerance)
        || !s.pressure_head_tolerance_m.is_finite()
        || !(0. ..=1e-5).contains(&s.pressure_head_tolerance_m)
        || s.maximum_free_heads == 0
        || s.maximum_free_heads > s.linear_solver.dimension_limit()
    {
        return Err(Error::Input(
            "invalid steady Darcy schema, units, active set or settings".into(),
        ));
    }
    let ids = material_ids
        .filter(|v| !v.is_empty())
        .unwrap_or(&mesh.material_ids);
    if ids.len() != active.len() {
        return Err(Error::Input(
            "hydraulic material map length mismatch".into(),
        ));
    }
    let mut materials = BTreeMap::new();
    for m in &case.materials {
        if materials
            .insert(m.id, (tensor(m)?, m.volume_source_per_s))
            .is_some()
        {
            return Err(Error::Input("duplicate hydraulic material id".into()));
        }
    }
    let owners = boundary::surface_owners(mesh);
    let mut parent: Vec<_> = (0..active.len()).collect();
    for values in owners.values() {
        let live: Vec<_> = values.iter().copied().filter(|&i| active[i]).collect();
        if live.len() == 2 {
            let a = root(&mut parent, live[0]);
            let b = root(&mut parent, live[1]);
            parent[b] = a;
        }
    }
    let mut components = BTreeMap::new();
    let mut element_component = vec![usize::MAX; active.len()];
    let mut node_component = vec![None; mesh.nodes.len()];
    for i in 0..active.len() {
        if active[i] {
            let r = root(&mut parent, i);
            let next = components.len();
            let c = *components.entry(r).or_insert(next);
            element_component[i] = c;
            for n in mesh.cell_nodes(i) {
                if node_component[n].is_some_and(|v| v != c) {
                    return Err(Error::Input("hydraulic regions touch only at nodes or edges; split hydraulic nodes or provide a common face".into()));
                }
                node_component[n] = Some(c);
            }
            if !materials.contains_key(&ids[i]) {
                return Err(Error::Input(
                    "active region has no hydraulic material".into(),
                ));
            }
        }
    }
    let mut head_values = vec![None; mesh.nodes.len()];
    let mut reference_heads = vec![None; components.len()];
    let mut head_faces = BTreeSet::new();
    for h in &case.heads {
        if !h.head_m.is_finite()
            || !h
                .gradient
                .iter()
                .chain(&h.reference_point_m)
                .all(|x| x.is_finite())
        {
            return Err(Error::Input("non-finite prescribed head".into()));
        }
        if mesh.face_sets.get(&h.face_set).is_some_and(Vec::is_empty) {
            return Err(Error::Input("empty hydraulic head face set".into()));
        }
        for key in mesh
            .face_sets
            .get(&h.face_set)
            .ok_or_else(|| Error::Input("unknown hydraulic head face set".into()))?
        {
            let mut key = *key;
            key.sort_unstable();
            let live: Vec<_> = owners
                .get(&key)
                .ok_or_else(|| Error::Input("unknown hydraulic face".into()))?
                .iter()
                .copied()
                .filter(|&i| active[i])
                .collect();
            if live.len() != 1 {
                return Err(Error::Input(
                    "head condition is not on an exposed active face".into(),
                ));
            }
            head_faces.insert(key);
            let component = element_component[live[0]];
            let datum = *reference_heads[component].get_or_insert(h.head_m);
            for n in mesh.face_nodes(live[0], &key)? {
                let value = (h.head_m - datum)
                    + (0..3)
                        .map(|a| h.gradient[a] * (mesh.nodes[n][a] - h.reference_point_m[a]))
                        .sum::<f64>();
                if !value.is_finite()
                    || head_values[n].is_some_and(|v: f64| {
                        (v - value).abs() > 1e-11 * (1. + v.abs().max(value.abs()))
                    })
                {
                    return Err(Error::Input(
                        "conflicting or non-finite nodal hydraulic heads".into(),
                    ));
                }
                head_values[n] = Some(value);
            }
        }
    }
    let anchored: BTreeSet<_> = head_values
        .iter()
        .enumerate()
        .filter_map(|(i, h)| h.and(node_component[i]))
        .collect();
    if anchored.len() != components.len() {
        return Err(Error::Input("every face-connected hydraulic region needs a prescribed head; no arbitrary gauge or pure-Neumann fallback".into()));
    }
    let reference_heads = reference_heads
        .into_iter()
        .map(|v| v.ok_or_else(|| Error::Input("missing hydraulic reference head".into())))
        .collect::<Result<Vec<_>>>()?;
    let mut source = vec![0.; mesh.nodes.len()];
    let mut prescribed = vec![0.; mesh.nodes.len()];
    let mut source_components = vec![0.; components.len()];
    let mut absolute_source_components = vec![0.; components.len()];
    let mut k = CsrMatrix::zeros(CsrPattern::from_connectivity(
        mesh.nodes.len(),
        (0..active.len())
            .filter(|&i| active[i])
            .map(|i| mesh.cell_nodes(i).to_vec()),
    )?);
    let mut point_count = 0;
    for i in 0..active.len() {
        if active[i] {
            let (conductivity, rate) = materials[&ids[i]];
            let nodes = mesh.cell_nodes(i);
            for q in mesh.integration(i)? {
                point_count += 1;
                if point_count > 10_000_000 {
                    return Err(Error::Input("hydraulic quadrature budget exceeded".into()));
                }
                let g = gradients(&q);
                source_components[element_component[i]] += rate * q.weight_m3;
                absolute_source_components[element_component[i]] += rate.abs() * q.weight_m3;
                for (a, &n) in nodes.iter().enumerate() {
                    source[n] += rate * q.shape[a] * q.weight_m3;
                    for (b, &m) in nodes.iter().enumerate() {
                        k.add(n, m, q.weight_m3 * g[a].dot(&(conductivity * g[b])))?;
                    }
                }
            }
        }
    }
    let mut flux_faces = BTreeMap::new();
    let mut neumann_components = vec![0.; components.len()];
    let mut absolute_neumann_components = vec![0.; components.len()];
    for f in &case.fluxes {
        if mesh.face_sets.get(f.face_set()).is_some_and(Vec::is_empty) {
            return Err(Error::Input("empty hydraulic flux face set".into()));
        }
        for key in mesh
            .face_sets
            .get(f.face_set())
            .ok_or_else(|| Error::Input("unknown hydraulic flux face set".into()))?
        {
            let mut key = *key;
            key.sort_unstable();
            let live: Vec<_> = owners
                .get(&key)
                .ok_or_else(|| Error::Input("unknown hydraulic flux face".into()))?
                .iter()
                .copied()
                .filter(|&i| active[i])
                .collect();
            if live.len() != 1 || head_faces.contains(&key) || flux_faces.contains_key(&key) {
                return Err(Error::Input(
                    "hydraulic face has conflicting or inactive boundary conditions".into(),
                ));
            }
            let nodes = mesh.face_nodes(live[0], &key)?;
            let positions: Vec<_> = nodes.iter().map(|&n| mesh.nodes[n]).collect();
            let mut total = 0.;
            for q in elements::face_points(&positions)? {
                let normal = Vector3::from(q.area_vector);
                let out = match f {
                    Flux::Outward { flux_m_s, .. } => *flux_m_s * normal.norm(),
                    Flux::Vector {
                        flux_m_s,
                        gradient_per_s,
                        reference_point_m,
                        ..
                    } => Vector3::from(std::array::from_fn(|a| {
                        flux_m_s[a]
                            + (0..3)
                                .map(|b| {
                                    gradient_per_s[a][b] * (q.position[b] - reference_point_m[b])
                                })
                                .sum::<f64>()
                    }))
                    .dot(&normal),
                };
                if !out.is_finite() {
                    return Err(Error::Input("non-finite hydraulic boundary flux".into()));
                }
                total += out;
                absolute_neumann_components[element_component[live[0]]] += out.abs();
                for (a, &n) in nodes.iter().enumerate() {
                    prescribed[n] += q.shape[a] * out;
                }
            }
            flux_faces.insert(key, total);
            neumann_components[element_component[live[0]]] += total;
        }
    }
    if !source
        .iter()
        .chain(&source_components)
        .chain(&prescribed)
        .chain(&neumann_components)
        .chain(&absolute_source_components)
        .chain(&absolute_neumann_components)
        .all(|v| v.is_finite())
    {
        return Err(Error::Solve(
            "hydraulic source or boundary integration overflow".into(),
        ));
    }
    let mut delta: Vec<_> = head_values.iter().map(|h| h.unwrap_or(0.)).collect();
    let fixed_action = k.matvec(&delta)?;
    let rhs: Vec<_> = source.iter().zip(&prescribed).map(|(a, b)| a - b).collect();
    let free: Vec<_> = node_component
        .iter()
        .enumerate()
        .filter(|(i, c)| c.is_some() && head_values[*i].is_none())
        .map(|(i, _)| i)
        .collect();
    if free.len() > s.maximum_free_heads {
        return Err(Error::Input("hydraulic free head limit exceeded".into()));
    }
    let mut component_nodes = vec![Vec::new(); components.len()];
    let mut component_free = vec![Vec::new(); components.len()];
    for (n, c) in node_component.iter().enumerate() {
        if let Some(c) = c {
            component_nodes[*c].push(n);
            if head_values[n].is_none() {
                component_free[*c].push(n);
            }
        }
    }
    let mut mapping = vec![None; mesh.nodes.len()];
    for nodes in &component_free {
        for (i, &n) in nodes.iter().enumerate() {
            mapping[n] = Some(i);
        }
    }
    let mut component_solves = vec![Vec::new(); components.len()];
    let mut component_relative = vec![0.; components.len()];
    let mut linear = None;
    let mut refinement_solves = Vec::new();
    for c in 0..components.len() {
        let nodes = &component_nodes[c];
        let free = &component_free[c];
        let reduced = if free.is_empty() {
            None
        } else {
            let mut a = CsrMatrix::zeros(CsrPattern::from_connectivity(
                free.len(),
                (0..active.len())
                    .filter(|&i| active[i] && element_component[i] == c)
                    .map(|i| {
                        mesh.cell_nodes(i)
                            .iter()
                            .filter_map(|&n| mapping[n])
                            .collect()
                    }),
            )?);
            for (i, &n) in free.iter().enumerate() {
                for j in k.row_ptr()[n]..k.row_ptr()[n + 1] {
                    if let Some(col) = mapping[k.columns()[j]] {
                        a.add(i, col, k.values()[j])?;
                    }
                }
            }
            let b: Vec<_> = free.iter().map(|&n| rhs[n] - fixed_action[n]).collect();
            let (x, r) = linear_solve(&a, &b, s)?;
            for (&n, &value) in free.iter().zip(&x) {
                delta[n] = value;
            }
            if components.len() == 1 {
                linear = Some(r.clone());
            }
            component_solves[c].push(r);
            Some(a)
        };
        for _ in 0..3 {
            let (relative, correction_rhs) = component_residual(&k, &delta, &rhs, nodes, free)?;
            if relative <= s.relative_linear_tolerance {
                break;
            }
            let a = reduced
                .as_ref()
                .ok_or_else(|| Error::Solve("hydraulic residual without free heads".into()))?;
            let (correction, report) = linear_solve(a, &correction_rhs, s)?;
            for (&n, &value) in free.iter().zip(&correction) {
                delta[n] += value;
            }
            component_solves[c].push(report.clone());
            refinement_solves.push(report);
        }
        let (relative, _) = component_residual(&k, &delta, &rhs, nodes, free)?;
        if !relative.is_finite() || relative > s.relative_linear_tolerance * 1.01 {
            return Err(Error::Solve(format!(
                "hydraulic component {c} original free-row residual exceeded: {relative:e}"
            )));
        }
        component_relative[c] = relative;
    }
    let action = k.matvec(&delta)?;
    let residual: Vec<_> = action.iter().zip(&rhs).map(|(a, b)| a - b).collect();
    let scale = norm(&action).max(norm(&rhs)).max(1e-30);
    let relative = norm(&free.iter().map(|&n| residual[n]).collect::<Vec<_>>()) / scale;
    if !relative.is_finite() || relative > s.relative_linear_tolerance * 1.01 {
        return Err(Error::Solve(
            "original hydraulic free-row residual exceeded".into(),
        ));
    }
    let field = HeadField {
        relative_head_m: delta.clone(),
        datum_head_m: node_component
            .iter()
            .map(|c| c.map_or(0., |c| reference_heads[c]))
            .collect(),
        head_m: delta
            .iter()
            .enumerate()
            .map(|(i, v)| node_component[i].map_or(0., |c| reference_heads[c] + v))
            .collect(),
        defined_nodes: node_component.iter().map(Option::is_some).collect(),
    };
    field.validate(mesh, active)?;
    let dirichlet: Vec<_> = head_values
        .iter()
        .enumerate()
        .map(|(i, v)| v.map(|_| -residual[i]))
        .collect();
    let mut points = Vec::new();
    let mut minimum = f64::INFINITY;
    let mut bound = f64::INFINITY;
    for i in 0..active.len() {
        if active[i] {
            let nodes = mesh.cell_nodes(i);
            let (conductivity, _) = materials[&ids[i]];
            let anchor = field.relative_head_m[nodes[0]];
            let pressure: Vec<_> = nodes
                .iter()
                .map(|&n| field.datum_head_m[n] - mesh.nodes[n][2] + field.relative_head_m[n])
                .collect();
            minimum = minimum.min(pressure.iter().copied().fold(f64::INFINITY, f64::min));
            let mut coefficients = pressure[..4].to_vec();
            if nodes.len() == 10 {
                for (e, [a, b]) in elements::EDGES.into_iter().enumerate() {
                    let value = pressure[a] + 2. * (pressure[e + 4] - pressure[a])
                        - (pressure[b] - pressure[a]) / 2.;
                    if !value.is_finite() {
                        return Err(Error::Solve(
                            "quadratic pressure certificate overflow".into(),
                        ));
                    }
                    coefficients.push(value);
                }
            }
            bound = bound.min(coefficients.iter().copied().fold(f64::INFINITY, f64::min));
            for (j, q) in mesh.integration(i)?.into_iter().enumerate() {
                let g = gradients(&q);
                let gradient = nodes.iter().zip(&g).fold(Vector3::zeros(), |a, (&n, g)| {
                    a + g * (field.relative_head_m[n] - anchor)
                });
                let head = field.datum_head_m[nodes[0]]
                    + anchor
                    + nodes
                        .iter()
                        .zip(&q.shape)
                        .map(|(&n, w)| (field.relative_head_m[n] - anchor) * w)
                        .sum::<f64>();
                let p = field.pressure_head(&nodes, &q.shape, q.position[2])?;
                minimum = minimum.min(p);
                points.push(Point {
                    element_index: i,
                    local_index: j,
                    position_m: q.position,
                    weight_m3: q.weight_m3,
                    head_m: head,
                    pressure_head_m: p,
                    head_gradient: gradient.into(),
                    darcy_flux_m_s: (-conductivity * gradient).into(),
                });
            }
        }
    }
    let mut boundary_fluxes = Vec::new();
    let mut gradient_outflow = vec![0.; components.len()];
    for (face, all_owners) in &owners {
        let live: Vec<_> = all_owners.iter().copied().filter(|&i| active[i]).collect();
        if live.len() != 1 {
            continue;
        }
        let owner = live[0];
        let nodes = mesh.face_nodes(owner, face)?;
        let positions: Vec<_> = nodes.iter().map(|&n| mesh.nodes[n]).collect();
        let cell = mesh.cell_nodes(owner);
        let cell_positions: Vec<_> = cell.iter().map(|&n| mesh.nodes[n]).collect();
        let (conductivity, _) = materials[&ids[owner]];
        let anchor = field.relative_head_m[cell[0]];
        let mut flux = 0.;
        for q in elements::face_points(&positions)? {
            let mut bary = [0.; 4];
            for a in 0..3 {
                let mut w = q.shape[a];
                if nodes.len() == 6 {
                    for (j, [b, c]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
                        if a == b || a == c {
                            w += q.shape[j + 3] / 2.;
                        }
                    }
                }
                let v = cell[..4]
                    .iter()
                    .position(|&n| n == nodes[a])
                    .ok_or_else(|| Error::Input("hydraulic face/owner mismatch".into()))?;
                bary[v] = w;
            }
            let g = elements::scalar_gradient(&cell_positions, [bary[1], bary[2], bary[3]])?;
            let grad = cell.iter().zip(g).fold(Vector3::zeros(), |v, (&n, g)| {
                v + Vector3::from(g) * (field.relative_head_m[n] - anchor)
            });
            flux += (-conductivity * grad).dot(&Vector3::from(q.area_vector));
        }
        gradient_outflow[element_component[owner]] += flux;
        boundary_fluxes.push(FaceFlux {
            corners: *face,
            owner,
            boundary_kind: if head_faces.contains(face) {
                "head"
            } else if flux_faces.contains_key(face) {
                "prescribed_flux"
            } else {
                "no_flow"
            }
            .into(),
            gradient_outflow_m3_s: flux,
            prescribed_outflow_m3_s: if head_faces.contains(face) {
                None
            } else {
                Some(*flux_faces.get(face).unwrap_or(&0.))
            },
        });
    }
    let mut balances = Vec::new();
    for c in 0..components.len() {
        let out = dirichlet
            .iter()
            .enumerate()
            .filter(|(n, _)| node_component[*n] == Some(c))
            .map(|(_, v)| v.unwrap_or(0.))
            .sum::<f64>();
        let den = absolute_source_components[c]
            + absolute_neumann_components[c]
            + dirichlet
                .iter()
                .enumerate()
                .filter(|(n, _)| node_component[*n] == Some(c))
                .map(|(_, v)| v.unwrap_or(0.).abs())
                .sum::<f64>();
        let error = (out + neumann_components[c] - source_components[c]).abs() / den.max(1e-30);
        if !den.is_finite() || !error.is_finite() || error > s.relative_linear_tolerance * 20. {
            return Err(Error::Solve(
                "hydraulic component water balance exceeded".into(),
            ));
        }
        balances.push(Balance {
            reference_head_m: reference_heads[c],
            relative_free_head_residual: component_relative[c],
            linear_solves: component_solves[c].clone(),
            absolute_source_m3_s: absolute_source_components[c],
            absolute_neumann_flow_m3_s: absolute_neumann_components[c],
            component: c,
            elements: element_component.iter().filter(|&&v| v == c).count(),
            source_m3_s: source_components[c],
            neumann_outflow_m3_s: neumann_components[c],
            dirichlet_reaction_outflow_m3_s: out,
            relative_nodal_balance_error: error,
            gradient_boundary_outflow_m3_s: gradient_outflow[c],
        });
    }
    if !minimum.is_finite()
        || !bound.is_finite()
        || !points.iter().all(|p| {
            p.head_m.is_finite()
                && p.pressure_head_m.is_finite()
                && p.head_gradient
                    .iter()
                    .chain(&p.darcy_flux_m_s)
                    .all(|v| v.is_finite())
        })
        || !boundary_fluxes
            .iter()
            .all(|p| p.gradient_outflow_m3_s.is_finite())
        || !balances.iter().all(|b| {
            [
                b.source_m3_s,
                b.absolute_source_m3_s,
                b.absolute_neumann_flow_m3_s,
                b.neumann_outflow_m3_s,
                b.dirichlet_reaction_outflow_m3_s,
                b.gradient_boundary_outflow_m3_s,
            ]
            .iter()
            .all(|v| v.is_finite())
        })
    {
        return Err(Error::Solve(
            "non-finite hydraulic output; no complete report can be serialized".into(),
        ));
    }
    let complete = bound >= -s.pressure_head_tolerance_m;
    Ok(Report {
        refinement_solves,
        schema: "soil-fem-darcy-result/1".into(),
        complete,
        failure_reason: if complete {
            None
        } else {
            Some(if minimum< -s.pressure_head_tolerance_m {"negative pressure head in saturated-only request; free-surface/Richards solve required"}else{"nonnegative pressure head not certified by quadratic Bernstein bounds"}.into())
        },
        mesh_sha256: crate::checkpoint::hash_json(mesh)?,
        case_sha256: crate::checkpoint::hash_json(case)?,
        active_elements: active.to_vec(),
        material_ids: ids.to_vec(),
        field,
        points,
        boundary_fluxes,
        dirichlet_nodal_outflow_m3_s: dirichlet,
        balances,
        relative_free_head_residual: relative,
        minimum_sampled_pressure_head_m: minimum,
        bernstein_pressure_head_lower_bound_m: bound,
        linear,
    })
}
