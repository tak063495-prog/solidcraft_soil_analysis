//! Rigid-component kinematics, separate from force equilibrium and stability.
//! Face-connected solids each carry six modes. Shared edge/vertex DOFs are
//! imposed as equality rows, retaining genuine hinge mechanisms.
use crate::{Error, Result, boundary, mesh::Mesh};
use nalgebra::{DMatrix, DVector, Matrix3, SMatrix, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

type Motion = SMatrix<f64, 3, 6>;
struct Component {
    center: Vector3<f64>,
    radius: f64,
}
pub struct Geometry {
    origin: Vector3<f64>,
    positions: Vec<[f64; 3]>,
    components: Vec<Component>,
    element_component: Vec<Option<usize>>,
    primary: Vec<Option<usize>>,
    ties: Vec<(usize, usize, usize)>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub continuum_components: usize,
    pub rigid_modes: usize,
    pub restrained_modes: usize,
    pub unrestrained_modes: usize,
}
pub struct Space<'a> {
    geometry: &'a Geometry,
    basis: Vec<DVector<f64>>,
    weak_rows: Vec<DVector<f64>>,
    weak_basis: Vec<DVector<f64>>,
}
fn root(parents: &mut [usize], mut i: usize) -> usize {
    while parents[i] != i {
        parents[i] = parents[parents[i]];
        i = parents[i];
    }
    i
}
impl Geometry {
    pub fn new(mesh: &Mesh, active: &[bool]) -> Result<Self> {
        mesh.validate()?;
        if active.len() != mesh.tetrahedra.len()
            || !mesh.nodes.iter().flatten().all(|x| x.is_finite())
            || (0..mesh.tetrahedra.len())
                .any(|e| mesh.cell_nodes(e).iter().any(|&n| n >= mesh.nodes.len()))
        {
            return Err(Error::Input(
                "invalid rigid-component geometry/layout".into(),
            ));
        }
        let mut parents: Vec<_> = (0..active.len()).collect();
        for owners in boundary::surface_owners(mesh).values() {
            if owners.len() > 2 {
                return Err(Error::Input("nonmanifold rigid-component face".into()));
            }
            if owners.len() == 2 && active[owners[0]] && active[owners[1]] {
                let a = root(&mut parents, owners[0]);
                let b = root(&mut parents, owners[1]);
                parents[b] = a;
            }
        }
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for (e, &enabled) in active.iter().enumerate() {
            if enabled {
                groups.entry(root(&mut parents, e)).or_default().push(e);
            }
        }
        if groups.is_empty() || groups.len() > 500 {
            return Err(Error::Input(
                "rigid-mode audit requires 1..=500 face-connected solids".into(),
            ));
        }
        let mut components = Vec::new();
        let common_origin = Vector3::from(mesh.nodes[0]);
        let positions: Vec<[f64; 3]> = mesh
            .nodes
            .iter()
            .map(|&p| (Vector3::from(p) - common_origin).into())
            .collect();
        let mut element_component = vec![None; active.len()];
        let mut primary = vec![None; mesh.nodes.len()];
        let mut ties = Vec::new();
        for elements in groups.values() {
            let index = components.len();
            let mut nodes = BTreeSet::new();
            for &e in elements {
                element_component[e] = Some(index);
                nodes.extend(mesh.cell_nodes(e));
            }
            let origin = Vector3::from(
                positions[*nodes
                    .first()
                    .ok_or_else(|| Error::Input("empty rigid component".into()))?],
            );
            let center = origin
                + nodes
                    .iter()
                    .map(|&n| Vector3::from(positions[n]) - origin)
                    .sum::<Vector3<f64>>()
                    / nodes.len() as f64;
            let radius = nodes
                .iter()
                .map(|&n| (Vector3::from(positions[n]) - center).norm())
                .fold(0., f64::max);
            if !radius.is_finite() || radius <= 0. {
                return Err(Error::Input("degenerate rigid-component radius".into()));
            }
            for &n in &nodes {
                if let Some(first) = primary[n] {
                    ties.push((n, first, index));
                } else {
                    primary[n] = Some(index);
                }
            }
            components.push(Component { center, radius });
        }
        Ok(Self {
            origin: common_origin,
            positions,
            components,
            element_component,
            primary,
            ties,
        })
    }
    pub fn dimension(&self) -> usize {
        6 * self.components.len()
    }
    fn motion(&self, component: usize, position: [f64; 3]) -> Motion {
        let c = &self.components[component];
        let r = (Vector3::from(position) - c.center) / c.radius;
        Motion::from_fn(|a, b| {
            if b < 3 {
                f64::from(a == b)
            } else {
                let mut e = Vector3::zeros();
                e[b - 3] = 1.;
                e.cross(&r)[a]
            }
        })
    }
    fn node_row(&self, node: usize, axis: usize) -> DVector<f64> {
        let mut row = DVector::zeros(self.dimension());
        if let Some(c) = self.primary[node] {
            let m = self.motion(c, self.positions[node]);
            for k in 0..6 {
                row[6 * c + k] = m[(axis, k)];
            }
        }
        row
    }
    fn relative_row(
        &self,
        owner_a: usize,
        owner_b: usize,
        position_a: [f64; 3],
        position_b: [f64; 3],
        direction: [f64; 3],
    ) -> Result<DVector<f64>> {
        let a = self
            .element_component
            .get(owner_a)
            .copied()
            .flatten()
            .ok_or_else(|| Error::Input("inactive/unknown rigid owner a".into()))?;
        let b = self
            .element_component
            .get(owner_b)
            .copied()
            .flatten()
            .ok_or_else(|| Error::Input("inactive/unknown rigid owner b".into()))?;
        let direction = Vector3::from(direction);
        let ma = self.motion(a, (Vector3::from(position_a) - self.origin).into());
        let mb = self.motion(b, (Vector3::from(position_b) - self.origin).into());
        let mut row = DVector::zeros(self.dimension());
        for k in 0..6 {
            row[6 * a + k] -= direction.dot(&ma.column(k));
            row[6 * b + k] += direction.dot(&mb.column(k));
        }
        Ok(row)
    }
    pub fn constrained(&self, rows: &[Vec<(usize, f64)>]) -> Result<Space<'_>> {
        if rows.len() != 3 * self.positions.len() {
            return Err(Error::Input("rigid constraint transform mismatch".into()));
        }
        // This projector is valid for node-local orthonormal nullspaces.
        // A nonlocal/tied-coordinate transform needs a different reduction.
        let mut owners = BTreeMap::new();
        for node in 0..self.positions.len() {
            let mut columns = BTreeMap::<usize, Vector3<f64>>::new();
            for axis in 0..3 {
                for &(index, value) in &rows[3 * node + axis] {
                    if !value.is_finite()
                        || index >= rows.len()
                        || owners
                            .insert(index, node)
                            .is_some_and(|other| other != node)
                    {
                        return Err(Error::Input(
                            "rigid audit requires finite node-local constraints".into(),
                        ));
                    }
                    columns.entry(index).or_default()[axis] += value;
                }
            }
            for (i, a) in &columns {
                for (j, b) in &columns {
                    if (a.dot(b) - f64::from(i == j)).abs() > 1e-12 {
                        return Err(Error::Input(
                            "rigid audit requires orthonormal constraint columns".into(),
                        ));
                    }
                }
            }
        }
        let mut space = Space {
            geometry: self,
            basis: Vec::new(),
            weak_rows: Vec::new(),
            weak_basis: Vec::new(),
        };
        for (n, component) in self.primary.iter().enumerate() {
            if component.is_none() {
                continue;
            }
            let projection = Matrix3::from_fn(|a, b| {
                f64::from(a == b)
                    - rows[3 * n + a]
                        .iter()
                        .map(|&(i, c)| {
                            c * rows[3 * n + b]
                                .iter()
                                .filter(|&&(j, _)| j == i)
                                .map(|&(_, d)| d)
                                .sum::<f64>()
                        })
                        .sum::<f64>()
            });
            for a in 0..3 {
                let mut r = DVector::zeros(self.dimension());
                for b in 0..3 {
                    r += self.node_row(n, b) * projection[(a, b)];
                }
                space.add(r)?;
            }
        }
        for &(node, a, b) in &self.ties {
            let ma = self.motion(a, self.positions[node]);
            let mb = self.motion(b, self.positions[node]);
            for axis in 0..3 {
                let mut r = DVector::zeros(self.dimension());
                for k in 0..6 {
                    r[6 * a + k] -= ma[(axis, k)];
                    r[6 * b + k] += mb[(axis, k)];
                }
                space.add(r)?;
            }
        }
        Ok(space)
    }
    /// Map component null motions to orthonormal free-coordinate columns.
    pub(crate) fn reduced_modes(
        &self,
        modes: &[DVector<f64>],
        rows: &[Vec<(usize, f64)>],
        count: usize,
    ) -> Result<DMatrix<f64>> {
        if rows.len() != 3 * self.positions.len()
            || modes
                .iter()
                .any(|m| m.len() != self.dimension() || !m.iter().all(|v| v.is_finite()))
            || rows
                .iter()
                .flatten()
                .any(|&(i, c)| i >= count || !c.is_finite())
        {
            return Err(Error::Input(
                "rigid predictor transform/mode mismatch".into(),
            ));
        }
        if modes.len() > 128 || count.saturating_mul(modes.len()) > 20_000_000 {
            return Err(Error::Solve(
                "rigid predictor workspace budget exceeded".into(),
            ));
        }
        let mut basis: Vec<DVector<f64>> = Vec::new();
        for mode in modes {
            let mut v = DVector::zeros(count);
            for (node, component) in self.primary.iter().enumerate() {
                let Some(c) = component else { continue };
                let physical = self.motion(*c, self.positions[node]) * mode.rows(6 * c, 6);
                for a in 0..3 {
                    for &(i, w) in &rows[3 * node + a] {
                        v[i] += w * physical[a];
                    }
                }
            }
            for _ in 0..2 {
                for q in &basis {
                    let dot = q.dot(&v);
                    v -= q * dot;
                }
            }
            let norm = v.norm();
            if norm > 1e-10 {
                basis.push(v / norm);
            }
        }
        Ok(if basis.is_empty() {
            DMatrix::zeros(count, 0)
        } else {
            DMatrix::from_columns(&basis)
        })
    }
}

/// Pick independent free coordinates spanning the current Newton null modes.
/// These only fix zero increments in a temporary predictor equation.
pub(crate) struct NewtonGauges {
    pub coordinates: Vec<usize>,
    pub modes: DMatrix<f64>,
}
pub(crate) fn gauge_coordinates(
    rigid: &DMatrix<f64>,
    action: impl Fn(&DVector<f64>) -> Result<DVector<f64>>,
    scale: f64,
) -> Result<NewtonGauges> {
    if rigid.ncols() == 0 {
        return Ok(NewtonGauges {
            coordinates: vec![],
            modes: DMatrix::zeros(rigid.nrows(), 0),
        });
    }
    if !scale.is_finite() || scale <= 0. {
        return Err(Error::Solve(
            "invalid tangent scale for contact predictor".into(),
        ));
    }
    let mut columns = Vec::new();
    for column in rigid.column_iter() {
        columns.push(action(&column.into_owned())?);
    }
    let az = DMatrix::from_columns(&columns);
    let svd = az.svd(false, true);
    let vt = svd
        .v_t
        .ok_or_else(|| Error::Solve("missing contact predictor singular vectors".into()))?;
    let null: Vec<_> = svd
        .singular_values
        .iter()
        .enumerate()
        .filter(|(_, s)| **s <= scale * 1e-10)
        .map(|(i, _)| vt.row(i).transpose())
        .collect();
    if null.is_empty() {
        return Ok(NewtonGauges {
            coordinates: vec![],
            modes: DMatrix::zeros(rigid.nrows(), 0),
        });
    }
    let modes = rigid * DMatrix::from_columns(&null);
    let mut selected = Vec::new();
    let mut basis: Vec<DVector<f64>> = Vec::new();
    for _ in 0..modes.ncols() {
        let mut best = None;
        let mut maximum = 0.;
        for i in 0..modes.nrows() {
            if selected.contains(&i) {
                continue;
            }
            let mut row = modes.row(i).transpose();
            for _ in 0..2 {
                for q in &basis {
                    let dot = q.dot(&row);
                    row -= q * dot;
                }
            }
            let norm = row.norm();
            if norm > maximum {
                maximum = norm;
                best = Some((i, row));
            }
        }
        let (i, row) = best
            .filter(|_| maximum > 1e-10)
            .ok_or_else(|| Error::Solve("could not select independent Newton gauges".into()))?;
        selected.push(i);
        basis.push(row / maximum);
    }
    Ok(NewtonGauges {
        coordinates: selected,
        modes,
    })
}
impl Space<'_> {
    pub fn full(&self) -> bool {
        self.basis.len() + self.weak_basis.len() == self.geometry.dimension()
    }
    pub fn add(&mut self, mut row: DVector<f64>) -> Result<()> {
        if row.len() != self.geometry.dimension() || !row.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("invalid rigid support row".into()));
        }
        if self.full() {
            return Ok(());
        }
        let original = row.norm();
        if original <= 1e-10 {
            return Ok(());
        }
        row /= original;
        for _ in 0..2 {
            for q in &self.basis {
                let dot = q.dot(&row);
                row -= q * dot;
            }
        }
        let norm = row.norm();
        if norm >= f64::EPSILON.sqrt() {
            self.basis.push(row / norm);
            // A nearly collinear early quadrature row must not be scaled to
            // unit length before later, well resolved rows arrive. Otherwise
            // roundoff in that weak pivot can invent global restraints.
            let pending = std::mem::take(&mut self.weak_rows);
            self.weak_basis.clear();
            for mut weak in pending {
                for _ in 0..2 {
                    for q in &self.basis {
                        let dot = q.dot(&weak);
                        weak -= q * dot;
                    }
                }
                self.add_weak(weak);
            }
        } else if norm > 1e-10 {
            self.add_weak(row);
        }
        Ok(())
    }
    fn add_weak(&mut self, row: DVector<f64>) {
        let mut independent = row.clone();
        for _ in 0..2 {
            for q in &self.weak_basis {
                let dot = q.dot(&independent);
                independent -= q * dot;
            }
        }
        let norm = independent.norm();
        if norm > 1e-10 {
            self.weak_rows.push(row);
            self.weak_basis.push(independent / norm);
        }
    }
    pub fn relative(
        &mut self,
        owner_a: usize,
        owner_b: usize,
        position_a: [f64; 3],
        position_b: [f64; 3],
        direction: [f64; 3],
    ) -> Result<()> {
        if !self.full() {
            self.add(
                self.geometry
                    .relative_row(owner_a, owner_b, position_a, position_b, direction)?,
            )?;
        }
        Ok(())
    }
    pub fn physical_row(&mut self, coefficients: &[(usize, f64)]) -> Result<()> {
        if self.full() {
            return Ok(());
        }
        let mut row = DVector::zeros(self.geometry.dimension());
        for &(g, c) in coefficients {
            if g >= 3 * self.geometry.positions.len() || !c.is_finite() {
                return Err(Error::Input("invalid structural rigid row".into()));
            }
            row += self.geometry.node_row(g / 3, g % 3) * c;
        }
        self.add(row)
    }
    pub fn report(&self) -> Report {
        Report {
            continuum_components: self.geometry.components.len(),
            rigid_modes: self.geometry.dimension(),
            restrained_modes: self.basis.len() + self.weak_basis.len(),
            unrestrained_modes: self.geometry.dimension()
                - self.basis.len()
                - self.weak_basis.len(),
        }
    }
    pub fn null_modes(&self) -> Vec<DVector<f64>> {
        let mut basis = self.basis.clone();
        basis.extend(self.weak_basis.iter().cloned());
        let mut free = Vec::new();
        for i in 0..self.geometry.dimension() {
            let mut v = DVector::zeros(self.geometry.dimension());
            v[i] = 1.;
            for _ in 0..2 {
                for q in &basis {
                    let dot = q.dot(&v);
                    v -= q * dot;
                }
            }
            let norm = v.norm();
            if norm > 1e-10 {
                let q = v / norm;
                basis.push(q.clone());
                free.push(q);
            }
        }
        free
    }
    pub fn require_supported(&self) -> Result<Report> {
        let report = self.report();
        if !self.full() {
            return Err(Error::Solve(format!(
                "JOINT rigid mechanism: {} unrestrained modes in {} continuum components",
                report.unrestrained_modes, report.continuum_components
            )));
        }
        Ok(report)
    }
}
