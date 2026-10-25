//! Error-bounded integral of max(n_z,0) on an actual Tri3/Tri6 surface.
//! The oriented Jacobian is quadratic. Bernstein signs identify exact patches;
//! linear clipping on ambiguous patches has a quadratic-deviation area bound.
use crate::{Error, Result, elements};
use nalgebra::Vector3;
use std::{cmp::Ordering, collections::BinaryHeap};

#[derive(Debug)]
pub struct Projection {
    pub area_m2: f64,
    pub domain_error_bound_m2: f64,
    pub visited_cells: usize,
}
#[derive(Clone)]
struct Cell {
    v: [[f64; 2]; 3],
    estimate: f64,
    error: f64,
    depth: usize,
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
fn mid(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [(a[0] + b[0]) / 2., (a[1] + b[1]) / 2.]
}
fn determinant(v: [[f64; 2]; 3]) -> f64 {
    ((v[1][0] - v[0][0]) * (v[2][1] - v[0][1]) - (v[2][0] - v[0][0]) * (v[1][1] - v[0][1])).abs()
}
fn value(points: &[[f64; 3]], r: [f64; 2]) -> f64 {
    let (_, gradient) = elements::tri_shape(r, points.len() == 6);
    let mut du = Vector3::zeros();
    let mut dv = Vector3::zeros();
    for (p, g) in points.iter().zip(gradient) {
        let x = Vector3::from(*p) - Vector3::from(points[0]);
        du += x * g[0];
        dv += x * g[1];
    }
    du.cross(&dv)[2]
}
fn integral(points: &[[f64; 3]], v: [[f64; 2]; 3]) -> f64 {
    determinant(v) / 6.
        * (value(points, mid(v[0], v[1]))
            + value(points, mid(v[1], v[2]))
            + value(points, mid(v[2], v[0])))
}
fn cell(points: &[[f64; 3]], v: [[f64; 2]; 3], depth: usize) -> Cell {
    let q = v.map(|r| value(points, r));
    let mut coefficients = q.to_vec();
    let mut deviation: f64 = 0.;
    for [i, j] in [[0, 1], [1, 2], [2, 0]] {
        let b = 2. * value(points, mid(v[i], v[j])) - (q[i] + q[j]) / 2.;
        deviation = deviation.max((b - (q[i] + q[j]) / 2.).abs());
        coefficients.push(b);
    }
    if coefficients.iter().all(|q| *q >= 0.) {
        return Cell {
            v,
            estimate: integral(points, v),
            error: 0.,
            depth,
        };
    }
    if coefficients.iter().all(|q| *q <= 0.) {
        return Cell {
            v,
            estimate: 0.,
            error: 0.,
            depth,
        };
    }
    let mut polygon = vec![];
    for i in 0..3 {
        let j = (i + 1) % 3;
        if q[i] >= 0. {
            polygon.push(v[i]);
        }
        if (q[i] > 0. && q[j] < 0.) || (q[i] < 0. && q[j] > 0.) {
            let t = q[i] / (q[i] - q[j]);
            polygon.push(std::array::from_fn(|a| v[i][a] + t * (v[j][a] - v[i][a])));
        }
    }
    let estimate = (1..polygon.len().saturating_sub(1))
        .map(|i| integral(points, [polygon[0], polygon[i], polygon[i + 1]]))
        .sum();
    Cell {
        v,
        estimate,
        error: deviation * determinant(v) / 2.,
        depth,
    }
}
pub fn upward_area(points: &[[f64; 3]]) -> Result<Projection> {
    let scale = elements::face_points(points)?
        .iter()
        .map(|q| Vector3::from(q.area_vector).norm())
        .sum::<f64>();
    let first = cell(points, [[0., 0.], [1., 0.], [0., 1.]], 0);
    let mut area = first.estimate;
    let mut error = first.error;
    let mut heap = BinaryHeap::from([first]);
    let mut count = 1;
    let tolerance = 1e-11 * scale;
    while error > tolerance {
        let c = heap
            .pop()
            .ok_or_else(|| Error::Solve("surface projection lost error heap".into()))?;
        if count + 4 > 100_000 || c.depth >= 30 {
            return Err(Error::Solve(
                "surface projection integration budget exhausted".into(),
            ));
        }
        let [a, b, d] = c.v;
        let ab = mid(a, b);
        let bd = mid(b, d);
        let da = mid(d, a);
        area -= c.estimate;
        error = (error - c.error).max(0.);
        for v in [[a, ab, da], [ab, b, bd], [da, bd, d], [ab, bd, da]] {
            let child = cell(points, v, c.depth + 1);
            area += child.estimate;
            error += child.error;
            heap.push(child);
            count += 1;
        }
    }
    if !area.is_finite() || area < -tolerance {
        return Err(Error::Solve("invalid upward projected area".into()));
    }
    Ok(Projection {
        area_m2: area.max(0.),
        domain_error_bound_m2: error,
        visited_cells: count,
    })
}

/// Smooth surface area: require agreement between successively higher rules.
/// This is a quadrature convergence estimate, not a certified interval bound.
pub fn surface_area(points: &[[f64; 3]]) -> Result<(f64, f64)> {
    let mut previous = elements::face_points(points)?
        .iter()
        .map(|q| Vector3::from(q.area_vector).norm())
        .sum::<f64>();
    for order in [8, 12, 16, 24] {
        let area = elements::face_points_order(points, order)?
            .iter()
            .map(|q| Vector3::from(q.area_vector).norm())
            .sum::<f64>();
        let error = (area - previous).abs();
        if error <= 1e-12 * area {
            return Ok((area, error));
        }
        previous = area;
    }
    Err(Error::Solve(
        "curved surface area quadrature did not converge".into(),
    ))
}
