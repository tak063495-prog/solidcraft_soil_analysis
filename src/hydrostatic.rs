//! Hydrostatic face traction on the actual Tri3/Tri6 isoparametric surface.
//! Quadratic height and oriented Jacobian Bernstein bounds control the error
//! caused by replacing the wet domain with a clipped planar reference patch.
//! Polynomial integrands are integrated on that domain without flattening CAD.
use crate::{Error, Result, elements};
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, collections::BinaryHeap, sync::OnceLock};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Integration {
    pub relative_tolerance: f64,
    pub absolute_tolerance_kn: f64,
    pub maximum_cells: usize,
}
impl Default for Integration {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-11,
            absolute_tolerance_kn: 1e-11,
            maximum_cells: 50_000,
        }
    }
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FaceCase {
    pub schema: String,
    pub points_m: Vec<[f64; 3]>,
    pub water_level_m: f64,
    pub water_unit_weight_kn_m3: f64,
    #[serde(default)]
    pub integration: Integration,
}
pub fn run_face(case: FaceCase) -> Result<FaceLoad> {
    if case.schema != "soil-fem-hydrostatic-face/1" {
        return Err(Error::Input("invalid hydrostatic face schema".into()));
    }
    face_load(
        &case.points_m,
        case.water_level_m,
        case.water_unit_weight_kn_m3,
        case.integration,
    )
}
#[derive(Debug, Serialize)]
pub struct FaceLoad {
    pub nodal_force_kn: Vec<[f64; 3]>,
    /// L2 bound for all nodal force components from wet-domain approximation.
    /// Floating point roundoff is additional; this is not an FE spatial bound.
    pub domain_error_bound_kn: f64,
    pub permitted_domain_error_kn: f64,
    pub visited_cells: usize,
    pub maximum_depth: usize,
}
#[derive(Clone)]
struct Cell {
    vertices: [[f64; 2]; 3],
    depth: usize,
    depths: [f64; 3],
    error: f64,
    wet: bool,
    dry: bool,
}
impl PartialEq for Cell {
    fn eq(&self, b: &Self) -> bool {
        self.error == b.error
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for Cell {
    fn cmp(&self, b: &Self) -> Ordering {
        self.error.total_cmp(&b.error)
    }
}

fn mapping(points: &[[f64; 3]], r: [f64; 2], level: f64) -> (Vec<f64>, f64, Vector3<f64>) {
    let (n, g) = elements::tri_shape(r, points.len() == 6);
    let origin = Vector3::from(points[0]);
    let mut depth = level - points[0][2];
    let mut du = Vector3::zeros();
    let mut dv = Vector3::zeros();
    for i in 0..points.len() {
        let offset = Vector3::from(points[i]) - origin;
        depth -= offset[2] * n[i];
        du += offset * g[i][0];
        dv += offset * g[i][1];
    }
    (n, depth, du.cross(&dv))
}
fn midpoint(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [(a[0] + b[0]) / 2., (a[1] + b[1]) / 2.]
}
fn determinant(v: [[f64; 2]; 3]) -> f64 {
    (v[1][0] - v[0][0]) * (v[2][1] - v[0][1]) - (v[1][1] - v[0][1]) * (v[2][0] - v[0][0])
}
fn info(points: &[[f64; 3]], level: f64, gamma: f64, v: [[f64; 2]; 3], depth: usize) -> Cell {
    let corners = v.map(|r| mapping(points, r, level));
    let depths = std::array::from_fn(|i| corners[i].1);
    let mut coefficients = depths.to_vec();
    let mut jacobian_bound = corners.iter().map(|q| q.2.norm()).fold(0., f64::max);
    let mut deviation: f64 = 0.;
    for [i, j] in [[0, 1], [1, 2], [2, 0]] {
        let m = mapping(points, midpoint(v[i], v[j]), level);
        let b = 2. * m.1 - (depths[i] + depths[j]) / 2.;
        deviation = deviation.max((b - (depths[i] + depths[j]) / 2.).abs());
        coefficients.push(b);
        jacobian_bound = jacobian_bound.max((2. * m.2 - (corners[i].2 + corners[j].2) / 2.).norm());
    }
    let wet = coefficients.iter().all(|d| *d >= 0.);
    let dry = coefficients.iter().all(|d| *d <= 0.);
    let error = if wet || dry {
        0.
    } else {
        gamma * deviation * determinant(v).abs() / 2.
            * jacobian_bound
            * (points.len() as f64).sqrt()
    };
    Cell {
        vertices: v,
        depth,
        depths,
        error,
        wet,
        dry,
    }
}

/// Sufficient whole-surface depth bounds, including curved edge extrema.
pub fn depth_bounds(points: &[[f64; 3]], level: f64) -> Result<[f64; 2]> {
    if ![3, 6].contains(&points.len())
        || !level.is_finite()
        || !points.iter().flatten().all(|x| x.is_finite())
    {
        return Err(Error::Input(
            "hydrostatic face requires finite Tri3/Tri6 and level".into(),
        ));
    }
    let d = [
        mapping(points, [0., 0.], level).1,
        mapping(points, [1., 0.], level).1,
        mapping(points, [0., 1.], level).1,
    ];
    let mut c = d.to_vec();
    for (r, [i, j]) in [
        ([0.5, 0.], [0, 1]),
        ([0.5, 0.5], [1, 2]),
        ([0., 0.5], [2, 0]),
    ] {
        c.push(2. * mapping(points, r, level).1 - (d[i] + d[j]) / 2.);
    }
    Ok([
        c.iter().copied().fold(f64::INFINITY, f64::min),
        c.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    ])
}

fn integrate_triangle(
    points: &[[f64; 3]],
    level: f64,
    gamma: f64,
    v: [[f64; 2]; 3],
    force: &mut [[f64; 3]],
) -> Result<()> {
    static RULE: OnceLock<Vec<([f64; 2], f64)>> = OnceLock::new();
    let rule = RULE.get_or_init(|| {
        // Four-point Gauss-Jacobi x Gauss-Legendre integrates degree <=7.
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0.5, 0., 0.],
            [0.5, 0.5, 0.],
            [0., 0.5, 0.],
        ];
        elements::face_points(&p)
            .map(|q| {
                q.into_iter()
                    .map(|q| ([q.position[0], q.position[1]], q.area_vector[2]))
                    .collect()
            })
            .unwrap_or_default()
    });
    if rule.is_empty() {
        return Err(Error::Input(
            "hydrostatic integration rule unavailable".into(),
        ));
    }
    let jac = determinant(v);
    if jac < 0. || !jac.is_finite() {
        return Err(Error::Input("invalid hydrostatic reference patch".into()));
    }
    if jac == 0. {
        return Ok(());
    }
    for &(r, w) in rule {
        let a = std::array::from_fn(|i| {
            v[0][i] + r[0] * (v[1][i] - v[0][i]) + r[1] * (v[2][i] - v[0][i])
        });
        let (n, d, normal) = mapping(points, a, level);
        for i in 0..points.len() {
            for axis in 0..3 {
                force[i][axis] -= gamma * d * n[i] * normal[axis] * w * jac;
            }
        }
    }
    Ok(())
}
fn integrate_cell(
    points: &[[f64; 3]],
    level: f64,
    gamma: f64,
    c: &Cell,
    force: &mut [[f64; 3]],
) -> Result<()> {
    if c.dry {
        return Ok(());
    }
    let mut polygon = Vec::new();
    if c.wet {
        polygon.extend(c.vertices);
    } else {
        for i in 0..3 {
            let j = (i + 1) % 3;
            let a = c.depths[i];
            let b = c.depths[j];
            if a >= 0. {
                polygon.push(c.vertices[i]);
            }
            if (a > 0. && b < 0.) || (a < 0. && b > 0.) {
                let t = a / (a - b);
                polygon.push(std::array::from_fn(|axis| {
                    c.vertices[i][axis] + t * (c.vertices[j][axis] - c.vertices[i][axis])
                }));
            }
        }
    }
    for i in 1..polygon.len().saturating_sub(1) {
        integrate_triangle(
            points,
            level,
            gamma,
            [polygon[0], polygon[i], polygon[i + 1]],
            force,
        )?;
    }
    Ok(())
}

pub fn face_load(
    points: &[[f64; 3]],
    level: f64,
    gamma: f64,
    settings: Integration,
) -> Result<FaceLoad> {
    let bounds = depth_bounds(points, level)?;
    if !gamma.is_finite()
        || gamma <= 0.
        || !settings.relative_tolerance.is_finite()
        || !(1e-13..=1e-5).contains(&settings.relative_tolerance)
        || !settings.absolute_tolerance_kn.is_finite()
        || settings.absolute_tolerance_kn <= 0.
        || !(1..=200_000).contains(&settings.maximum_cells)
    {
        return Err(Error::Input(
            "invalid hydrostatic integration controls".into(),
        ));
    }
    // Validate the original surface Jacobian before any dry-domain shortcut.
    elements::face_points(points)?;
    let root = info(points, level, gamma, [[0., 0.], [1., 0.], [0., 1.]], 0);
    let mut jac: f64 = 0.;
    for q in elements::face_points(points)? {
        jac += Vector3::from(q.area_vector).norm();
    }
    let permitted = settings.absolute_tolerance_kn
        + settings.relative_tolerance
            * gamma
            * bounds[1].max(0.)
            * jac
            * (points.len() as f64).sqrt();
    if !jac.is_finite() || !permitted.is_finite() || !root.error.is_finite() {
        return Err(Error::Input(
            "hydrostatic area or error scale overflow".into(),
        ));
    }
    let mut force = vec![[0.; 3]; points.len()];
    let mut heap = BinaryHeap::new();
    let mut error = root.error;
    let mut visited = 1;
    let mut maximum_depth = 0;
    if root.wet || root.dry {
        integrate_cell(points, level, gamma, &root, &mut force)?;
    } else {
        heap.push(root);
    }
    while error > permitted {
        let c = heap
            .pop()
            .ok_or_else(|| Error::Solve("hydrostatic error accounting failed".into()))?;
        if visited + 4 > settings.maximum_cells || c.depth >= 22 {
            return Err(Error::Solve(
                "hydrostatic wet-domain integration budget exceeded".into(),
            ));
        }
        error = (error - c.error).max(0.);
        let [a, b, d] = c.vertices;
        let ab = midpoint(a, b);
        let bd = midpoint(b, d);
        let da = midpoint(d, a);
        for vertices in [[a, ab, da], [ab, b, bd], [da, bd, d], [ab, bd, da]] {
            let child = info(points, level, gamma, vertices, c.depth + 1);
            visited += 1;
            maximum_depth = maximum_depth.max(child.depth);
            if child.wet || child.dry {
                integrate_cell(points, level, gamma, &child, &mut force)?;
            } else {
                error += child.error;
                heap.push(child);
            }
        }
    }
    for c in heap {
        integrate_cell(points, level, gamma, &c, &mut force)?;
    }
    if !force.iter().flatten().all(|f| f.is_finite()) || !error.is_finite() {
        return Err(Error::Solve("non-finite hydrostatic face force".into()));
    }
    Ok(FaceLoad {
        nodal_force_kn: force,
        domain_error_bound_kn: error,
        permitted_domain_error_kn: permitted,
        visited_cells: visited,
        maximum_depth,
    })
}
