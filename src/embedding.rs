//! Variationally consistent point embedding into Tet4/Tet10 host translations.
use crate::{Error, Result, elements, mesh::Mesh};
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Point {
    pub position_m: [f64; 3],
    pub body: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Binding {
    pub position_m: [f64; 3],
    pub host_elements: Vec<usize>,
    pub weights: Vec<(usize, f64)>,
    pub maximum_coordinate_residual_m: f64,
}

fn mapping(
    mesh: &Mesh,
    nodes: &[usize],
    r: Vector3<f64>,
) -> (Vector3<f64>, Matrix3<f64>, Vec<f64>) {
    let (shape, gradient) = if nodes.len() == 10 {
        let (n, g) = elements::tet10_shape(r.into());
        (n.to_vec(), g.to_vec())
    } else {
        (
            vec![1. - r.sum(), r[0], r[1], r[2]],
            vec![[-1.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        )
    };
    let mut x = Vector3::zeros();
    let mut jacobian = Matrix3::zeros();
    for (i, &node) in nodes.iter().enumerate() {
        let p = Vector3::from(mesh.nodes[node]);
        x += p * shape[i];
        jacobian += p * Vector3::from(gradient[i]).transpose();
    }
    (x, jacobian, shape)
}

fn inverse(mesh: &Mesh, nodes: &[usize], point: Vector3<f64>, tolerance: f64) -> Option<Vec<f64>> {
    let p: Vec<_> = nodes
        .iter()
        .map(|&n| Vector3::from(mesh.nodes[n]))
        .collect();
    // Bernstein control points bound the curved isoparametric volume. Bounds
    // of only the ten Lagrange nodes are not sufficient for quadratic maps.
    let mut controls = p[..4].to_vec();
    if nodes.len() == 10 {
        for (i, [a, b]) in elements::EDGES.into_iter().enumerate() {
            controls.push(2. * p[i + 4] - (p[a] + p[b]) / 2.);
        }
    }
    for a in 0..3 {
        let lo = controls.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min);
        let hi = controls
            .iter()
            .map(|p| p[a])
            .fold(f64::NEG_INFINITY, f64::max);
        if point[a] < lo - tolerance || point[a] > hi + tolerance {
            return None;
        }
    }
    let chord = Matrix3::from_columns(&[p[1] - p[0], p[2] - p[0], p[3] - p[0]]);
    let affine = chord.lu().solve(&(point - p[0]))?;
    let scale = (p[1] - p[0])
        .norm()
        .max((p[2] - p[0]).norm())
        .max((p[3] - p[0]).norm());
    let natural_tolerance = (tolerance / scale).max(1e-12);
    let starts = [
        affine,
        Vector3::repeat(0.25),
        Vector3::repeat(0.05),
        Vector3::new(0.85, 0.05, 0.05),
        Vector3::new(0.05, 0.85, 0.05),
        Vector3::new(0.05, 0.05, 0.85),
    ];
    for mut r in starts {
        for _ in 0..40 {
            let (x, j, shape) = mapping(mesh, nodes, r);
            let residual = x - point;
            if residual.norm() <= tolerance {
                if r.iter().all(|&a| a >= -natural_tolerance) && r.sum() <= 1. + natural_tolerance {
                    return Some(shape);
                }
                break;
            }
            let Some(delta) = j.lu().solve(&residual) else {
                break;
            };
            let mut alpha = 1.;
            let mut accepted = false;
            for _ in 0..16 {
                let candidate = r - alpha * delta;
                if mapping(mesh, nodes, candidate).0.metric_distance(&point) < residual.norm() {
                    r = candidate;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                break;
            }
        }
    }
    None
}

impl Binding {
    /// Mesh quality must be validated by the caller once for the whole model.
    pub fn locate(mesh: &Mesh, point: &Point, tolerance: f64) -> Result<Self> {
        if !point.position_m.iter().all(|x| x.is_finite())
            || !tolerance.is_finite()
            || tolerance <= 0.
            || tolerance > 1e-5
        {
            return Err(Error::Input("invalid embedded point or tolerance".into()));
        }
        let bodies: BTreeMap<_, _> = mesh.regions.iter().map(|r| (r.id, &r.body)).collect();
        if let Some(body) = &point.body
            && !bodies.values().any(|b| *b == body)
        {
            return Err(Error::Input(
                "embedded point host body does not exist".into(),
            ));
        }
        let target = Vector3::from(point.position_m);
        let mut candidates = Vec::new();
        for ei in 0..mesh.tetrahedra.len() {
            if let Some(body) = &point.body
                && bodies.get(&mesh.region_ids[ei]).copied() != Some(body)
            {
                continue;
            }
            let nodes = mesh.cell_nodes(ei);
            if let Some(shape) = inverse(mesh, &nodes, target, tolerance) {
                let weights: BTreeMap<_, _> = nodes
                    .into_iter()
                    .zip(shape)
                    .filter(|(_, v)| v.abs() > 1e-12)
                    .collect();
                candidates.push((ei, weights));
            }
        }
        let Some((_, first)) = candidates.first() else {
            return Err(Error::Input(format!(
                "embedded point {:?} has no host volume",
                point.position_m
            )));
        };
        for (_, candidate) in &candidates {
            if first
                .iter()
                .any(|(n, w)| (w - candidate.get(n).unwrap_or(&0.)).abs() > 1e-8)
                || candidate
                    .iter()
                    .any(|(n, w)| (w - first.get(n).unwrap_or(&0.)).abs() > 1e-8)
            {
                return Err(Error::Input("ambiguous embedded point across nonconforming/overlapping hosts; specify a host body".into()));
            }
        }
        let weights: Vec<_> = first.iter().map(|(&n, &w)| (n, w)).collect();
        let position = weights.iter().fold(Vector3::zeros(), |x, &(n, w)| {
            x + Vector3::from(mesh.nodes[n]) * w
        });
        let residual = (position - target).norm();
        if (weights.iter().map(|(_, w)| w).sum::<f64>() - 1.).abs() > 1e-10
            || residual > tolerance * 1.01
        {
            return Err(Error::Input(
                "embedded partition of unity/coordinate consistency failed".into(),
            ));
        }
        Ok(Self {
            position_m: point.position_m,
            host_elements: candidates.iter().map(|(ei, _)| *ei).collect(),
            weights,
            maximum_coordinate_residual_m: residual,
        })
    }
    pub fn displacement(&self, u: &[[f64; 3]]) -> Vector3<f64> {
        self.weights
            .iter()
            .fold(Vector3::zeros(), |x, &(n, w)| x + Vector3::from(u[n]) * w)
    }
    pub fn is_active(&self, active_elements: &[bool]) -> bool {
        self.host_elements.iter().any(|&ei| active_elements[ei])
    }
    pub fn distribute(&self, value: Vector3<f64>, nodal: &mut [[f64; 3]]) {
        for &(node, weight) in &self.weights {
            for a in 0..3 {
                nodal[node][a] += weight * value[a];
            }
        }
    }
}
