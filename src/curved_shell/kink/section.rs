//! One-sided section at a point on an owned edge, integrated through
//! thickness in the directed line frame. No cross-crease extrapolation.
use super::{Joint, Line, invalid};
use crate::Result;
use crate::curved_shell::{
    Geometry, Layer, Q36, Resultant, Shear, Stretch, Tying, V, geometry::NODES, path::Case,
};
use nalgebra::{DVector, Matrix3, Vector2};
use serde::{Deserialize, Serialize};

pub(super) struct Probe {
    group: String,
    face: String,
    nodes: [usize; 6],
    e: f64,
    arc: f64,
    outward: i8,
    axes: Matrix3<f64>,
    rs: [f64; 2],
    layers: Vec<(f64, f64, Stretch, Shear)>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Cut {
    pub group: String,
    pub face: String,
    pub arc_m: f64,
    /// Cut outward co-normal = sign * (director cross tangent).
    pub outward_sign: i8,
    pub section: Resultant,
}
impl Line {
    pub(super) fn probes(&self, case: &Case, joint: &Joint, g: &[Geometry]) -> Result<Vec<Probe>> {
        let mut probes = Vec::new();
        for w in &joint.wings {
            let sec = joint.section(case, &w.group);
            let nu = sec.poisson;
            let d = Matrix3::new(1., nu, 0., nu, 1., 0., 0., 0., (1. - nu) / 2.)
                * (sec.young_kpa / (1. - nu * nu));
            for i in (0..self.arc.len() - 2).step_by(2) {
                let node = w.edge_nodes[i + 1];
                let (f, geo) = case
                    .faces
                    .iter()
                    .zip(g)
                    .find(|(f, _)| f.group == w.group && f.nodes[3..].contains(&node))
                    .unwrap();
                let rs = NODES[f.nodes.iter().position(|j| *j == node).unwrap()];
                let s = geo.surface(rs)?;
                if self.axis.dot(&s.d).abs() > 1e-10 {
                    return Err(invalid("joint axis outside section tangent plane"));
                }
                let axes = Matrix3::from_columns(&[self.axis, s.d.cross(&self.axis), s.d]);
                let mid = geo.basis(&s, 0.);
                let area = mid.column(0).cross(&mid.column(1)).norm();
                let side = axes.column(1).dot(
                    &(geo.x[0] + geo.x[1] + geo.x[2]
                        - 3. * geo.x[f.nodes.iter().position(|j| *j == node).unwrap()]),
                );
                let mut layers = Vec::new();
                for (z, weight) in super::gauss(-sec.thickness_m / 2., sec.thickness_m / 2.) {
                    let (b, h) = Layer::new(geo, z)?.operators(geo, rs, &axes, Tying::default())?;
                    layers.push((
                        z,
                        weight * geo.basis(&s, z).determinant() / area,
                        d * b,
                        h * (5. * sec.young_kpa / (12. * (1. + nu))),
                    ));
                }
                probes.push(Probe {
                    group: w.group.clone(),
                    face: f.id.clone(),
                    nodes: f.nodes,
                    e: geo.e,
                    arc: self.arc[i + 1],
                    outward: if side > 0. { -1 } else { 1 },
                    axes,
                    rs,
                    layers,
                });
            }
        }
        Ok(probes)
    }
}
impl Probe {
    pub(super) fn cut(&self, delta: &DVector<f64>) -> Cut {
        let local = Q36::from_fn(|i, _| delta[6 * self.nodes[i / 6] + i % 6]);
        let mut n = V::zeros();
        let mut m = V::zeros();
        let mut q = Vector2::zeros();
        for (z, w, s, h) in &self.layers {
            let dn = (s * local) * *w;
            n += dn;
            m += dn * *z;
            q += (h * local) * *w;
        }
        Cut {
            group: self.group.clone(),
            face: self.face.clone(),
            arc_m: self.arc,
            outward_sign: self.outward,
            section: Resultant {
                rs: self.rs,
                mid_axes: self.axes.into(),
                n_kn_m: n.into(),
                m_kn: m.into(),
                m_attachment_kn: (m - self.e * n).into(),
                q_kn_m: q.into(),
            },
        }
    }
}
