//! Small-rotation, zero-width joint on explicitly paired straight Tri6 edges.
mod reduction;
mod section;
use super::{
    geometry::{Geometry, V},
    path::Case,
};
use crate::{Error, Result, plates::Section};
use nalgebra::{DMatrix, DVector};
pub(crate) use reduction::Reduction;
use reduction::Row;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub shell: Case,
    pub joint: Joint,
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Joint {
    pub id: String,
    pub rotation: Rotation,
    /// End, midpoint, end, midpoint, ... in increasing arc length.
    pub line_m: Vec<[f64; 3]>,
    pub wings: Vec<Wing>,
    #[serde(default)]
    pub loads: BTreeMap<usize, Vec<LineLoad>>,
    /// One separating plane for each pair of trimmed material domains.
    pub separators: Vec<Separator>,
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Separator {
    pub negative: String,
    pub positive: String,
    pub point_m: [f64; 3],
    pub normal: [f64; 3],
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Rotation {
    Rigid,
    Hinge,
    Spring { stiffness_kn_m_per_rad_m: f64 },
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Wing {
    pub group: String,
    pub edge_nodes: Vec<usize>,
    pub section: Option<Section>,
}
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LineLoad {
    pub id: String,
    pub wing: String,
    pub interval_m: [f64; 2],
    pub force_kn_m: [f64; 3],
    pub couple_kn_m_per_m: [f64; 3],
}
pub(crate) struct Line {
    probes: Vec<section::Probe>,
    points: Vec<V>,
    axis: V,
    arc: Vec<f64>,
    // Original material (attachment) arm, directed from the installed line.
    arms: Vec<Vec<V>>,
}
fn invalid(what: &str) -> Error {
    Error::Input(format!("offset kink joint: {what}"))
}
fn vector(q: &[f64], j: usize) -> V {
    V::from_column_slice(&q[j..j + 3])
}
fn close(a: V, b: V) -> bool {
    (a - b).norm() <= 256. * f64::EPSILON * (1. + a.norm() + b.norm())
}
fn quadratic(x: f64) -> [f64; 3] {
    [
        (1. - x) * (1. - 2. * x),
        4. * x * (1. - x),
        x * (2. * x - 1.),
    ]
}
fn gauss(a: f64, b: f64) -> [(f64, f64); 3] {
    [
        (-(0.6_f64).sqrt(), 5. / 9.),
        (0., 8. / 9.),
        ((0.6_f64).sqrt(), 5. / 9.),
    ]
    .map(|(x, w)| (a + (1. + x) * (b - a) / 2., w * (b - a) / 2.))
}
impl Joint {
    pub(crate) fn section<'a>(&'a self, case: &'a Case, group: &str) -> &'a Section {
        self.wings
            .iter()
            .find(|w| w.group == group)
            .and_then(|w| w.section.as_ref())
            .unwrap_or(&case.section)
    }
    pub(crate) fn shared_tie_error(&self, case: &Case, soil: usize, a: usize, b: usize) -> Error {
        // If both points lie on this line, measure the rank of equality of
        // their motions *after* imposing the joint kinematics.
        let station = self
            .wings
            .iter()
            .find_map(|w| w.edge_nodes.iter().position(|n| *n == a));
        if let Some(i) = station
            && self.wings.iter().any(|w| w.edge_nodes.get(i) == Some(&b))
            && let (Some(point), Some(first), Some(last)) =
                (self.line_m.get(i), self.line_m.first(), self.line_m.last())
            && let (Some(sa), Some(sb)) = (case.nodes.get(a), case.nodes.get(b))
        {
            let t = (V::from(*last) - V::from(*first)).normalize();
            let r = V::from(*sb) - V::from(*point);
            let mut matrix = DMatrix::zeros(3, 4);
            for k in 0..3 {
                matrix.set_column(k, &V::ith(k, 1.).cross(&(V::from(*sb) - V::from(*sa))));
            }
            if !matches!(self.rotation, Rotation::Rigid) {
                matrix.set_column(3, &t.cross(&r));
            }
            let rank = matrix.svd(false, false).rank(128. * f64::EPSILON);
            return invalid(&format!(
                "{} soil {soil}, wing nodes {a}/{b}: shared tie adds rank {rank}/3 after joint constraints",
                self.id
            ));
        }
        invalid(&format!(
            "{} repeated soil {soil} binding for nodes {a}/{b} outside a paired station",
            self.id
        ))
    }

    fn spring(&self) -> f64 {
        match self.rotation {
            Rotation::Spring {
                stiffness_kn_m_per_rad_m,
            } => stiffness_kn_m_per_rad_m,
            _ => 0.,
        }
    }
    pub(crate) fn install(&self, case: &Case, geo: &[Geometry], q: &[f64]) -> Result<Line> {
        let count = self.line_m.len();
        if self.id.trim().is_empty()
            || count < 3
            || count.is_multiple_of(2)
            || self.wings.len() < 2
            || (!matches!(self.rotation, Rotation::Rigid) && self.wings.len() != 2)
            || matches!(self.rotation, Rotation::Spring { .. })
                && !(self.spring().is_finite() && self.spring() > 0.)
        {
            return Err(invalid("invalid line layout or rotation law"));
        }
        let mut owned = BTreeSet::new();
        let mut groups = BTreeSet::new();
        let mut points = Vec::new();
        let mut arms = Vec::new();
        for w in &self.wings {
            if !groups.insert(&w.group) || w.edge_nodes.len() != count {
                return Err(invalid("repeated wing or missing stations"));
            }
            let mut wing_arms = Vec::new();
            for (i, &n) in w.edge_nodes.iter().enumerate() {
                if n >= case.nodes.len()
                    || !owned.insert(n)
                    || case
                        .faces
                        .iter()
                        .any(|f| f.group != w.group && f.nodes.contains(&n))
                {
                    return Err(invalid(&format!(
                        "{} node {n}: dependent or shared wing DOFs",
                        w.group
                    )));
                }
                let s = V::from(case.nodes[n]);
                let arm = s - V::from(self.line_m[i]);
                let theta = vector(q, n * 6 + 3);
                let current = s + vector(q, n * 6);
                let j = current - arm - theta.cross(&arm);
                if points.len() == i {
                    points.push(j);
                } else if !close(points[i], j) {
                    return Err(invalid(
                        "wing installation positions disagree at the joint line",
                    ));
                }
                wing_arms.push(current - j);
            }
            for triple in w.edge_nodes.windows(3).step_by(2) {
                let matches = case
                    .faces
                    .iter()
                    .filter(|f| {
                        f.group == w.group
                            && (0..3).any(|i| {
                                let [a, b] = [f.nodes[i], f.nodes[(i + 1) % 3]];
                                f.nodes[i + 3] == triple[1]
                                    && ((a == triple[0] && b == triple[2])
                                        || (a == triple[2] && b == triple[0]))
                            })
                    })
                    .count();
                if matches != 1 {
                    return Err(invalid("each ordered edge must have one boundary owner"));
                }
            }
            arms.push(wing_arms);
        }
        let axis = (points[count - 1] - points[0])
            .try_normalize(f64::MIN_POSITIVE)
            .ok_or_else(|| invalid("zero length line"))?;
        let arc: Vec<_> = points.iter().map(|p| axis.dot(&(p - points[0]))).collect();
        for i in 0..count {
            if !close(points[i], points[0] + axis * arc[i]) || (i > 0 && arc[i] <= arc[i - 1]) {
                return Err(invalid("nonstraight or unordered joint line"));
            }
            if i % 2 == 1
                && (!close(points[i], (points[i - 1] + points[i + 1]) / 2.)
                    || arms
                        .iter()
                        .any(|w| !close(w[i], w[i - 1]) || !close(w[i], w[i + 1])))
            {
                return Err(invalid("nonconforming straight edge or varying arm"));
            }
        }
        if case.faces.iter().any(|f| !groups.contains(&f.group)) {
            return Err(invalid(
                "shell group without an explicit trimmed wing owner",
            ));
        }
        self.check_domains(case, geo)?;
        for (event, loads) in &self.loads {
            if !matches!(
                event.checked_sub(1).and_then(|i| case.steps.get(i)),
                Some(crate::offset_composite::path::Step::Load { .. })
            ) {
                return Err(invalid("line load outside a load event"));
            }
            let mut ids = BTreeSet::new();
            for l in loads {
                let [a, b] = l.interval_m;
                if l.id.is_empty()
                    || !ids.insert(&l.id)
                    || !groups.contains(&l.wing)
                    || !(0. <= a
                        && a < b
                        && b <= arc[count - 1] + 128. * f64::EPSILON * arc[count - 1])
                    || l.force_kn_m
                        .iter()
                        .chain(&l.couple_kn_m_per_m)
                        .any(|x| !x.is_finite())
                {
                    return Err(invalid(
                        "invalid load owner, interval or duplicate line load ID",
                    ));
                }
            }
        }
        let mut line = Line {
            points,
            axis,
            arc,
            arms,
            probes: Vec::new(),
        };
        line.probes = line.probes(case, self, geo)?;
        Ok(line)
    }
    fn check_domains(&self, case: &Case, geo: &[Geometry]) -> Result<()> {
        let mut pairs = BTreeSet::new();
        for plane in &self.separators {
            let n = V::from(plane.normal);
            if (n.norm() - 1.).abs() > 1e-12
                || !n.norm().is_finite()
                || plane.point_m.iter().any(|v| !v.is_finite())
                || plane.negative == plane.positive
                || [&plane.negative, &plane.positive]
                    .iter()
                    .any(|id| !self.wings.iter().any(|w| &w.group == *id))
                || !pairs.insert((
                    plane.negative.as_str().min(&plane.positive),
                    plane.negative.as_str().max(&plane.positive),
                ))
            {
                return Err(invalid("invalid material-domain separator"));
            }
            for (f, g) in case.faces.iter().zip(geo) {
                let sign = if f.group == plane.negative {
                    -1.
                } else if f.group == plane.positive {
                    1.
                } else {
                    continue;
                };
                let half = self.section(case, &f.group).thickness_m / 2.;
                // The degree-two Bernstein control hull encloses the whole
                // attachment triangle. Bound a varying unit director by one.
                let mut hull = g.x;
                for i in 0..3 {
                    hull[3 + i] = 2. * g.x[3 + i] - (g.x[i] + g.x[(i + 1) % 3]) / 2.;
                }
                let flat = g.directors.iter().all(|d| *d == g.directors[0]);
                let dn = sign * n.dot(&g.directors[0]);
                let lower = if flat {
                    -g.e * dn - half * dn.abs()
                } else {
                    -(half + g.e.abs())
                };
                if hull.iter().any(|p| {
                    sign * n.dot(&(p - V::from(plane.point_m))) + lower
                        < -128. * f64::EPSILON * (1. + p.norm())
                }) {
                    return Err(invalid("overlap or uncertified trimmed wing volume"));
                }
            }
        }
        if pairs.len() != self.wings.len() * (self.wings.len() - 1) / 2 {
            return Err(invalid("missing trimmed volume owner separator"));
        }
        Ok(())
    }
}

impl Line {
    pub(crate) fn reduction(&self, joint: &Joint, free: &[usize]) -> Result<Reduction> {
        let t = self.axis;
        let unit = V::ith(t.iamin(), 1.);
        let a = t.cross(&unit).normalize();
        let b = t.cross(&a);
        let axes = if matches!(joint.rotation, Rotation::Rigid) {
            vec![V::x(), V::y(), V::z()]
        } else {
            vec![a, b]
        };
        let mut rows = Vec::new();
        for i in 0..self.points.len() {
            for side in 1..joint.wings.len() {
                let owner = format!(
                    "{} station {i} nodes {} and {}",
                    joint.id, joint.wings[0].edge_nodes[i], joint.wings[side].edge_nodes[i]
                );
                for k in 0..3 {
                    let mut row = Row::new();
                    for (s, sign) in [(0, -1.), (side, 1.)] {
                        let n = joint.wings[s].edge_nodes[i] * 6;
                        row.insert(n + k, sign);
                        // u_J = u_s - theta x (s-J).
                        for h in 0..3 {
                            row.insert(n + 3 + h, -sign * V::ith(h, 1.).cross(&self.arms[s][i])[k]);
                        }
                    }
                    rows.push((owner.clone(), row));
                }
                for axis in &axes {
                    let mut row = Row::new();
                    for (s, sign) in [(0, -1.), (side, 1.)] {
                        let n = joint.wings[s].edge_nodes[i] * 6;
                        for k in 0..3 {
                            row.insert(n + 3 + k, sign * axis[k]);
                        }
                    }
                    rows.push((owner.clone(), row));
                }
            }
        }
        Reduction::new(rows, free)
    }
    fn angle_row(&self, joint: &Joint, i: usize, h: [f64; 3]) -> Vec<(usize, f64)> {
        let mut row = Vec::new();
        for (s, sign) in [(0, -1.), (1, 1.)] {
            for (j, w) in h.iter().enumerate() {
                let n = joint.wings[s].edge_nodes[i + j];
                for k in 0..3 {
                    row.push((6 * n + 3 + k, sign * w * self.axis[k]));
                }
            }
        }
        row
    }
    pub(crate) fn add_spring(&self, joint: &Joint, k: &mut DMatrix<f64>) {
        if joint.spring() == 0. {
            return;
        }
        for i in (0..self.arc.len() - 2).step_by(2) {
            let len = self.arc[i + 2] - self.arc[i];
            for (x, w) in gauss(0., 1.) {
                let row = self.angle_row(joint, i, quadratic(x));
                for &(a, v) in &row {
                    for &(b, u) in &row {
                        k[(a, b)] += joint.spring() * len * w * v * u;
                    }
                }
            }
        }
    }
    pub(crate) fn spring_force(&self, joint: &Joint, delta: &DVector<f64>) -> (DVector<f64>, f64) {
        let mut f = DVector::zeros(delta.len());
        let mut energy = 0.;
        if joint.spring() != 0. {
            for i in (0..self.arc.len() - 2).step_by(2) {
                let len = self.arc[i + 2] - self.arc[i];
                let phi: [f64; 3] = std::array::from_fn(|j| {
                    let spin = |s: usize| {
                        vector(delta.as_slice(), 6 * joint.wings[s].edge_nodes[i + j] + 3)
                    };
                    self.axis.dot(&(spin(1) - spin(0)))
                });
                for (s, w) in gauss(0., 1.) {
                    let h = quadratic(s);
                    let angle = V::from(h).dot(&V::from(phi));
                    let moment = joint.spring() * angle * w * len;
                    energy += 0.5 * moment * angle;
                    for (k, sign) in [(0, -1.), (1, 1.)] {
                        for (j, v) in h.iter().enumerate() {
                            let n = joint.wings[k].edge_nodes[i + j];
                            for d in 0..3 {
                                f[6 * n + 3 + d] += sign * v * moment * self.axis[d];
                            }
                        }
                    }
                }
            }
        }
        (f, energy)
    }
    pub(crate) fn loads(&self, joint: &Joint, event: usize, f: &mut DVector<f64>) {
        for l in joint.loads.get(&event).into_iter().flatten() {
            let side = joint.wings.iter().position(|w| w.group == l.wing).unwrap();
            for i in (0..self.arc.len() - 2).step_by(2) {
                let start = self.arc[i];
                let end = self.arc[i + 2];
                let a = start.max(l.interval_m[0]);
                let b = end.min(l.interval_m[1]);
                if a >= b {
                    continue;
                }
                for (s, w) in gauss(a, b) {
                    for (j, h) in quadratic((s - start) / (end - start))
                        .into_iter()
                        .enumerate()
                    {
                        let n = joint.wings[side].edge_nodes[i + j];
                        let force = V::from(l.force_kn_m) * w * h;
                        let m = V::from(l.couple_kn_m_per_m) * w * h
                            - self.arms[side][i + j].cross(&force);
                        for k in 0..3 {
                            f[n * 6 + k] += force[k];
                            f[n * 6 + 3 + k] += m[k];
                        }
                    }
                }
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Station {
    pub group: String,
    pub node: usize,
    pub arc_m: f64,
    pub line_axis: [f64; 3],
    pub joint_m: [f64; 3],
    pub attachment_m: [f64; 3],
    pub middle_m: [f64; 3],
    pub director: [f64; 3],
    pub u_joint_m: [f64; 3],
    pub u_attachment_m: [f64; 3],
    pub u_middle_m: [f64; 3],
    pub rotation_rad: [f64; 3],
    pub tangent_rotation_rad: [f64; 3],
    pub normal_rotation_rad: f64,
    pub force_kn: Option<[f64; 3]>,
    pub force_density_kn_m: Option<[f64; 3]>,
    pub moment_density_kn_m_per_m: Option<[f64; 3]>,
    pub moment_at_joint_kn_m: Option<[f64; 3]>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub installation_relative_angle_rad: Vec<f64>,
    pub relative_angle_increment_rad: Vec<f64>,
    pub line_load_ids: Vec<String>,
    pub id: String,
    pub stations: Vec<Station>,
    pub cuts: Vec<section::Cut>,
    pub spring_energy_kn_m: f64,
}

impl Line {
    pub(crate) fn state(
        &self,
        joint: &Joint,
        g: &super::path::Generation,
        q: &DVector<f64>,
        transfer: Option<&DVector<f64>>,
        energy: f64,
    ) -> State {
        let mut stations = Vec::new();
        for (s, wing) in joint.wings.iter().enumerate() {
            for (i, &n) in wing.edge_nodes.iter().enumerate() {
                let e = g
                    .elements
                    .iter()
                    .find(|e| e.face.group == wing.group && e.face.nodes.contains(&n))
                    .unwrap();
                let j = e.face.nodes.iter().position(|v| *v == n).unwrap();
                let d = V::from(e.directors[j]);
                let mid_arm = -e.face.offset_m * d;
                let arm = self.arms[s][i];
                let theta = vector(q.as_slice(), 6 * n + 3);
                let rot = theta - vector(&g.installation, 6 * n + 3);
                let u = vector(q.as_slice(), 6 * n);
                let s0 = self.points[i] + arm;
                let initial_arm = V::from(e.initial_directors[j]) * (-e.face.offset_m);
                let original_j = V::from(joint.line_m[i]);
                let force = transfer.map(|f| vector(f.as_slice(), 6 * n));
                let moment = transfer.map(|f| {
                    vector(f.as_slice(), 6 * n + 3) + arm.cross(&vector(f.as_slice(), 6 * n))
                });
                stations.push(Station {
                    group: wing.group.clone(),
                    node: n,
                    arc_m: self.arc[i],
                    line_axis: self.axis.into(),
                    joint_m: self.points[i].into(),
                    attachment_m: s0.into(),
                    middle_m: (s0 + mid_arm).into(),
                    director: d.into(),
                    u_joint_m: (self.points[i] - original_j + u
                        - vector(&g.installation, 6 * n)
                        - rot.cross(&arm))
                    .into(),
                    u_attachment_m: u.into(),
                    u_middle_m: (u + mid_arm - initial_arm + rot.cross(&mid_arm)).into(),
                    rotation_rad: theta.into(),
                    tangent_rotation_rad: (theta - d * theta.dot(&d)).into(),
                    normal_rotation_rad: theta.dot(&d),
                    force_kn: force.map(Into::into),
                    force_density_kn_m: None,
                    moment_density_kn_m_per_m: None,
                    moment_at_joint_kn_m: moment.map(Into::into),
                });
            }
        }
        if transfer.is_some() {
            let n = self.arc.len();
            let mut mass = DMatrix::zeros(n, n);
            let coefficients = [[4., 2., -1.], [2., 16., 2.], [-1., 2., 4.]];
            for i in (0..n - 2).step_by(2) {
                for (a, row) in coefficients.iter().enumerate() {
                    for (b, w) in row.iter().enumerate() {
                        mass[(i + a, i + b)] += w * (self.arc[i + 2] - self.arc[i]) / 30.;
                    }
                }
            }
            // Consistent line-mass projection of nodal resultants (a density
            // representation of the same force, never an additional load).
            let factor = mass.cholesky().expect("positive straight segment lengths");
            for wing in stations.chunks_mut(n) {
                let rhs = DMatrix::from_fn(n, 6, |i, j| {
                    if j < 3 {
                        wing[i].force_kn.unwrap()[j]
                    } else {
                        wing[i].moment_at_joint_kn_m.unwrap()[j - 3]
                    }
                });
                let density = factor.solve(&rhs);
                for (i, w) in wing.iter_mut().enumerate() {
                    w.force_density_kn_m = Some(std::array::from_fn(|j| density[(i, j)]));
                    w.moment_density_kn_m_per_m =
                        Some(std::array::from_fn(|j| density[(i, j + 3)]));
                }
            }
        }
        let angles = |u: &[f64]| {
            if joint.wings.len() != 2 {
                return Vec::new();
            }
            (0..self.points.len())
                .map(|i| {
                    let spin = |w: usize| vector(u, 6 * joint.wings[w].edge_nodes[i] + 3);
                    self.axis.dot(&(spin(1) - spin(0)))
                })
                .collect::<Vec<_>>()
        };
        let reference = angles(&g.installation);
        let increment = angles(q.as_slice())
            .into_iter()
            .zip(&reference)
            .map(|(a, b)| a - b)
            .collect();
        State {
            installation_relative_angle_rad: reference,
            relative_angle_increment_rad: increment,
            line_load_ids: Vec::new(),
            id: joint.id.clone(),
            cuts: self
                .probes
                .iter()
                .map(|p| p.cut(&(q - DVector::from_column_slice(&g.installation))))
                .collect(),
            stations,
            spring_energy_kn_m: energy,
        }
    }
}
