//! Reference geometry and the linearisation of a unit interpolated thickness line.
use crate::plates::mindlin;
use crate::{Error, Result};
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub type V = Vector3<f64>;
pub const NODES: [[f64; 2]; 6] = [
    [0., 0.],
    [1., 0.],
    [0., 1.],
    [0.5, 0.],
    [0.5, 0.5],
    [0., 0.5],
];
const EDGES: [(usize, usize); 3] = [(0, 1), (1, 2), (2, 0)];

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Face {
    pub id: String,
    pub group: String,
    pub nodes: [usize; 6],
    pub offset_m: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Director {
    pub group: String,
    pub node: usize,
    pub direction: [f64; 3],
}

pub fn shape(p: [f64; 2]) -> ([f64; 6], [[f64; 6]; 2]) {
    let [r, s] = p;
    let l = 1. - r - s;
    (
        mindlin::shape([l, r, s]),
        [
            [1. - 4. * l, 4. * r - 1., 0., 4. * (l - r), 4. * s, -4. * s],
            [1. - 4. * l, 0., 4. * s - 1., -4. * r, 4. * r, 4. * (l - s)],
        ],
    )
}

pub fn interp(weights: &[f64; 6], nodes: &[V; 6]) -> V {
    weights.iter().zip(nodes).map(|(n, v)| v * *n).sum()
}

#[derive(Clone)]
pub struct Geometry {
    pub x: [V; 6],
    pub directors: [V; 6],
    pub e: f64,
}

pub struct Surface {
    pub n: [f64; 6],
    pub dn: [[f64; 6]; 2],
    pub tangent: [V; 2],
    pub d: V,
    pub dd: [V; 2],
    pub length: f64,
    pub length_gradient: [f64; 2],
}

impl Geometry {
    pub fn surface(&self, p: [f64; 2]) -> Result<Surface> {
        let (n, dn) = shape(p);
        let tangent = dn.map(|w| interp(&w, &self.x));
        let v = interp(&n, &self.directors);
        let length = v.norm();
        if !(length.is_finite() && length > 64. * f64::EPSILON) {
            return Err(Error::Input("curved Tri6 director cancellation".into()));
        }
        let d = v / length;
        let dv = dn.map(|w| interp(&w, &self.directors));
        let length_gradient = dv.map(|v| d.dot(&v));
        let dd = std::array::from_fn(|k| (dv[k] - d * length_gradient[k]) / length);
        Ok(Surface {
            n,
            dn,
            tangent,
            d,
            dd,
            length,
            length_gradient,
        })
    }

    pub fn basis(&self, s: &Surface, z: f64) -> Matrix3<f64> {
        Matrix3::from_columns(&[
            s.tangent[0] + (z - self.e) * s.dd[0],
            s.tangent[1] + (z - self.e) * s.dd[1],
            s.d,
        ])
    }
}

impl Surface {
    /// Variation of the normalized director and its two parameter derivatives.
    pub fn variation(&self, i: usize, rotation: V) -> (V, [V; 2]) {
        let projected = rotation - self.d * self.d.dot(&rotation);
        let delta = projected * (self.n[i] / self.length);
        let grad = std::array::from_fn(|k| {
            (projected * self.dn[k][i]
                - (self.dd[k] * self.d.dot(&rotation) + self.d * self.dd[k].dot(&rotation))
                    * self.n[i]
                - delta * self.length_gradient[k])
                / self.length
        });
        (delta, grad)
    }
}

fn unit(v: V) -> Result<V> {
    v.try_normalize(f64::MIN_POSITIVE)
        .filter(|u| u.iter().all(|a| a.is_finite()))
        .ok_or_else(|| Error::Input("degenerate curved Tri6 surface".into()))
}

/// Faces declared smooth share a director at every common node. Accumulate
/// one-sided unit normals in face ID order, weighted by physical face area.
pub fn geometries(
    nodes: &[[f64; 3]],
    faces: &[Face],
    explicit: &[Director],
) -> Result<Vec<Geometry>> {
    let mut order: Vec<_> = (0..faces.len()).collect();
    order.sort_by_key(|&i| &faces[i].id);
    let mut ids = BTreeSet::new();
    let corners: BTreeSet<_> = faces
        .iter()
        .flat_map(|f| f.nodes[..3].iter().copied())
        .collect();
    let mut triangles = BTreeSet::new();
    let mut midpoints = BTreeMap::new();
    let mut incident = BTreeMap::<(&str, usize), BTreeSet<usize>>::new();
    let mut adjacent = BTreeMap::<(usize, usize), Vec<usize>>::new();
    let mut normals = BTreeMap::<(&str, usize), Vec<(V, f64)>>::new();
    let mut edges = BTreeMap::<(usize, usize), Vec<(usize, usize, usize)>>::new();
    for i in order {
        let f = &faces[i];
        if f.id.is_empty()
            || f.group.is_empty()
            || !ids.insert(&f.id)
            || f.nodes.iter().any(|&j| j >= nodes.len())
            || f.nodes.into_iter().collect::<BTreeSet<_>>().len() != 6
            || !f.offset_m.is_finite()
        {
            return Err(Error::Input("invalid curved Tri6 face or nodes".into()));
        }
        let mut key = [f.nodes[0], f.nodes[1], f.nodes[2]];
        key.sort_unstable();
        if !triangles.insert(key) {
            return Err(Error::Input("repeated curved Tri6 face".into()));
        }
        let x = f.nodes.map(|j| V::from(nodes[j]));
        let normal = |p| {
            let (_, dn) = shape(p);
            interp(&dn[0], &x).cross(&interp(&dn[1], &x))
        };
        let area = 0.5
            * mindlin::quadrature()
                .map(|(l, w)| w * normal([l[1], l[2]]).norm())
                .sum::<f64>();
        for (j, p) in NODES.into_iter().enumerate() {
            incident
                .entry((&f.group, f.nodes[j]))
                .or_default()
                .insert(i);
            normals
                .entry((&f.group, f.nodes[j]))
                .or_default()
                .push((unit(normal(p))?, area));
        }
        for (k, (a, b)) in EDGES.into_iter().enumerate() {
            let (a, b) = (f.nodes[a], f.nodes[b]);
            let key = (a.min(b), a.max(b));
            let middle = f.nodes[k + 3];
            if corners.contains(&middle)
                || midpoints
                    .insert(middle, key)
                    .is_some_and(|prior| prior != key)
            {
                return Err(Error::Input("aliased curved Tri6 midpoint".into()));
            }
            let owners = edges.entry((a.min(b), a.max(b))).or_default();
            for &(other, first, mid) in owners.iter() {
                if faces[other].group == f.group {
                    for node in [a, b, middle] {
                        adjacent.entry((node, i)).or_default().push(other);
                        adjacent.entry((node, other)).or_default().push(i);
                    }
                }
                if mid != f.nodes[k + 3] || (faces[other].group == f.group && first == a) {
                    return Err(Error::Input(
                        "inconsistent curved Tri6 shared edge or orientation".into(),
                    ));
                }
            }
            owners.push((i, a, f.nodes[k + 3]));
            if owners.len() > 2 {
                return Err(Error::Input("nonmanifold curved Tri6 edge".into()));
            }
        }
    }
    for ((_, node), fan) in incident {
        let mut seen = BTreeSet::new();
        let mut pending = vec![*fan.first().unwrap()];
        while let Some(e) = pending.pop() {
            if seen.insert(e)
                && let Some(neighbors) = adjacent.get(&(node, e))
            {
                pending.extend(neighbors);
            }
        }
        if seen != fan {
            return Err(Error::Input(
                "disconnected smooth shell fan at shared node".into(),
            ));
        }
    }
    let mut overrides = BTreeMap::new();
    for d in explicit {
        let key = (d.group.as_str(), d.node);
        let v = V::from(d.direction);
        if !normals.contains_key(&key)
            || (v.norm() - 1.).abs() > 64. * f64::EPSILON
            || overrides.insert(key, unit(v)?).is_some()
        {
            return Err(Error::Input(
                "invalid or repeated explicit curved Tri6 director".into(),
            ));
        }
    }
    let mut directors = BTreeMap::new();
    let mut arms = BTreeMap::<usize, V>::new();
    for (key, sides) in normals {
        let d = match overrides.get(&key) {
            Some(d) => *d,
            None => unit(sides.iter().map(|(n, a)| n * *a).sum())?,
        };
        if sides.iter().any(|(n, _)| d.dot(n) <= 0.) {
            return Err(Error::Input("opposed curved Tri6 directors".into()));
        }
        directors.insert(key, d);
    }
    faces
        .iter()
        .map(|f| {
            let ds = f.nodes.map(|j| directors[&(f.group.as_str(), j)]);
            for (j, d) in f.nodes.iter().zip(ds) {
                let arm = f.offset_m * d;
                if let Some(other) = arms.insert(*j, arm)
                    && (other - arm).norm() > 64. * f64::EPSILON * other.norm().max(arm.norm())
                {
                    return Err(Error::Input(
                        "inconsistent signed curved Tri6 attachment arms".into(),
                    ));
                }
            }
            Ok(Geometry {
                x: f.nodes.map(|j| V::from(nodes[j])),
                directors: ds,
                e: f.offset_m,
            })
        })
        .collect()
}

impl Geometry {
    /// A sufficient injectivity certificate on the entire convex reference
    /// prism. In affine coordinates Y=A(r,s,z), require sym(dX/dY)>0.
    /// Integrating this inequality along the segment between any two points
    /// proves distinct images. Subdivision sharpens interval bounds on d.
    pub fn certify(&self, t: f64) -> Result<()> {
        let mid = self.surface([1. / 3.; 2])?;
        let a = self.basis(&mid, 0.);
        let inv = a
            .try_inverse()
            .filter(|_| a.determinant() > 0.)
            .ok_or_else(|| Error::Input("singular reference shell prism".into()))?;
        self.certify_cell([[0., 0.], [1., 0.], [0., 1.]], t, &inv, 0)
    }
    fn certify_cell(
        &self,
        tri: [[f64; 2]; 3],
        t: f64,
        inv: &Matrix3<f64>,
        depth: u32,
    ) -> Result<()> {
        let center = std::array::from_fn(|k| tri.iter().map(|v| v[k]).sum::<f64>() / 3.);
        let s = self.surface(center)?;
        let v0 = interp(&s.n, &self.directors);
        let dv0 = s.dn.map(|w| interp(&w, &self.directors));
        let corners = tri.map(|p| shape(p).1);
        let mut radius = 0_f64;
        let mut da = [0_f64; 2];
        let mut dv = [0_f64; 2];
        for (p, grad) in tri.iter().zip(corners) {
            for k in 0..2 {
                da[k] = da[k].max((interp(&grad[k], &self.x) - s.tangent[k]).norm());
                dv[k] = dv[k].max((interp(&grad[k], &self.directors) - dv0[k]).norm());
            }
            let diff = [p[0] - center[0], p[1] - center[1]];
            // Directional derivative of a quadratic is affine; its norm
            // over the subtriangle is bounded by its three vertex values.
            for g in corners {
                radius = radius.max(
                    (interp(&g[0], &self.directors) * diff[0]
                        + interp(&g[1], &self.directors) * diff[1])
                        .norm(),
                );
            }
        }
        let lower = v0.norm() - radius;
        if lower > 0. {
            let dd = radius / lower;
            let deriv = std::array::from_fn::<_, 2, _>(|k| {
                dv[k] / lower + dv0[k].norm() * (2. * dd / lower + radius / (lower * v0.norm()))
            });
            let certified = [-t / 2., t / 2.].into_iter().all(|z| {
                let jac = self.basis(&s, z) * inv;
                let sym = (jac + jac.transpose()) * 0.5;
                let error = (da[0] + (z - self.e).abs() * deriv[0]) * inv.row(0).norm()
                    + (da[1] + (z - self.e).abs() * deriv[1]) * inv.row(1).norm()
                    + dd * inv.row(2).norm();
                sym.symmetric_eigenvalues().min() > error + 128. * f64::EPSILON * (1. + jac.norm())
            });
            if certified {
                return Ok(());
            }
        }
        if depth == 6 {
            return Err(Error::Input(
                "curved shell thickness prism cannot be certified injective".into(),
            ));
        }
        let [a, b, c] = tri;
        let middle = |a: [f64; 2], b: [f64; 2]| [(a[0] + b[0]) / 2., (a[1] + b[1]) / 2.];
        let (ab, bc, ca) = (middle(a, b), middle(b, c), middle(c, a));
        for child in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
            self.certify_cell(child, t, inv, depth + 1)?;
        }
        Ok(())
    }
}
