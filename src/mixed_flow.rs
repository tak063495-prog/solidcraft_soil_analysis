//! Locally conservative hybrid RT0 geometry on isoparametric tetrahedra.
//! Four integrated outward face fluxes. Contravariant Piola mapping preserves
//! face integrals also on curved Tet10. Pressure/head is cellwise constant.
use crate::{Error, Result, elements, mesh::Mesh};
use nalgebra::{Matrix3, Matrix4, Vector3};

#[derive(Debug, Clone)]
pub struct Point {
    pub natural: [f64; 3],
    pub position_m: [f64; 3],
    pub weight_m3: f64,
    pub basis_per_m2: [Vector3<f64>; 4],
    pub divergence_per_m3: f64,
}
#[derive(Debug, Clone)]
pub struct Cell {
    pub element_index: usize,
    pub volume_m3: f64,
    pub points: Vec<Point>,
    /// Integral of z times a unit-flux RT0 normal trace (reference face mean).
    pub face_elevation_m: [f64; 4],
    pub face_area_m2: [f64; 4],
}
impl Cell {
    pub fn from_mesh(mesh: &Mesh, index: usize) -> Result<Self> {
        mesh.validate()?;
        Self::new(mesh, index)
    }
    pub(crate) fn new(mesh: &Mesh, index: usize) -> Result<Self> {
        if index >= mesh.tetrahedra.len() {
            return Err(Error::Input("unknown mixed-flow element".into()));
        }
        let coordinates: Vec<_> = mesh
            .cell_nodes(index)
            .iter()
            .map(|&n| mesh.nodes[n])
            .collect();
        let ref_vertices = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let mut points = Vec::new();
        for &(r, w) in elements::volume_rule(mesh.quadrature_order.unwrap_or(3))? {
            let (x, j) = elements::tet_mapping(&coordinates, r)?;
            let det = j.determinant();
            points.push(Point {
                natural: r,
                position_m: x,
                weight_m3: w * det,
                basis_per_m2: std::array::from_fn(|i| {
                    j * (Vector3::from(r) - Vector3::from(ref_vertices[i])) * (2. / det)
                }),
                divergence_per_m3: 6. / det,
            });
        }
        let volume_m3: f64 = points.iter().map(|p| p.weight_m3).sum();
        let mut face_elevation_m = [0.; 4];
        let mut face_area_m2 = [0.; 4];
        for i in 0..4 {
            let key = elements::FACES[i].map(|j| mesh.tetrahedra[index][j]);
            let p: Vec<_> = mesh
                .face_nodes(index, &key)?
                .iter()
                .map(|&n| mesh.nodes[n])
                .collect();
            face_area_m2[i] = elements::face_points(&p)?
                .iter()
                .map(|q| Vector3::from(q.area_vector).norm())
                .sum();
            // Quadratic elevation integrates exactly at the three face points.
            for l in [
                [2. / 3., 1. / 6., 1. / 6.],
                [1. / 6., 2. / 3., 1. / 6.],
                [1. / 6., 1. / 6., 2. / 3.],
            ] {
                let r = std::array::from_fn(|a| {
                    (0..3)
                        .map(|k| l[k] * ref_vertices[elements::FACES[i][k]][a])
                        .sum()
                });
                face_elevation_m[i] += elements::tet_mapping(&coordinates, r)?.0[2] / 3.;
            }
        }
        if !volume_m3.is_finite()
            || volume_m3 <= 0.
            || face_area_m2.iter().any(|a| !a.is_finite() || *a <= 0.)
        {
            return Err(Error::Input("invalid mixed-flow geometry".into()));
        }
        Ok(Self {
            element_index: index,
            volume_m3,
            points,
            face_elevation_m,
            face_area_m2,
        })
    }
    /// Mass matrix for inverse hydraulic conductivity at every quadrature point.
    pub fn mass(&self, inverse_conductivity: &[Matrix3<f64>]) -> Result<Matrix4<f64>> {
        let m = self.weighted_mass(inverse_conductivity)?;
        if m.cholesky().is_none() {
            return Err(Error::Solve(
                "non-positive RT0 hydraulic mass matrix".into(),
            ));
        }
        Ok(m)
    }
    /// The same bilinear integral, also for signed conductivity derivatives.
    pub fn weighted_mass(&self, inverse_conductivity: &[Matrix3<f64>]) -> Result<Matrix4<f64>> {
        if inverse_conductivity.len() != self.points.len()
            || inverse_conductivity
                .iter()
                .any(|m| !m.iter().all(|v| v.is_finite()))
        {
            return Err(Error::Input("invalid mixed conductivity layout".into()));
        }
        let mut m = Matrix4::<f64>::zeros();
        for (q, k) in self.points.iter().zip(inverse_conductivity) {
            for i in 0..4 {
                for j in 0..4 {
                    m[(i, j)] += q.weight_m3 * q.basis_per_m2[i].dot(&(k * q.basis_per_m2[j]));
                }
            }
        }
        if !m.iter().all(|v| v.is_finite()) {
            return Err(Error::Solve(
                "non-positive/overflowing RT0 hydraulic mass matrix".into(),
            ));
        }
        Ok(m)
    }
}
