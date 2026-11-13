//! Matching Tri3/Tri6 interfaces; only declared Tet10 crack-tip edges may share nodes.
#[path = "joints_global_step.rs"]
mod global_step;
pub use global_step::GlobalStepControl;

use crate::{
    Error, Result, boundary, elements,
    joint_law::{Material, PointState},
    mesh::Mesh,
};
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
#[path = "joints_birth.rs"]
mod birth;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub id: String,
    pub side_a: String,
    pub side_b: String,
    pub material: Material,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub birth_policy: Option<BirthPolicy>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub water: Option<crate::joint_water::Policy>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub crack_tip: Option<CrackTip>,
}

/// Explicit, immutable topology certificate for a degenerate Tet10 crack front.
/// Edges are [endpoint, endpoint, intermediate node], in displacement IDs.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CrackTip {
    pub edges: Vec<[usize; 3]>,
    pub faces: Vec<OwnedFace>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OwnedFace {
    pub face_a: [usize; 3],
    pub face_b: [usize; 3],
    pub owner_a: usize,
    pub owner_b: usize,
}
fn face_key(mut face: [usize; 3]) -> [usize; 3] {
    face.sort_unstable();
    face
}
impl CrackTip {
    fn validate(
        &self,
        mesh: &Mesh,
        fa: &[[usize; 3]],
        fb: &[[usize; 3]],
    ) -> Result<BTreeSet<usize>> {
        let invalid =
            || Error::Input("JOINT invalid crack-tip edge/surface/owner declaration".into());
        if mesh.tetrahedra10.is_empty() || self.edges.is_empty() || self.faces.len() != fa.len() {
            return Err(invalid());
        }
        let owners = boundary::surface_owners(mesh);
        let mut a_faces = BTreeSet::new();
        let mut b_faces = BTreeSet::new();
        let mut surface_edges = BTreeMap::<[usize; 2], (usize, usize)>::new();
        for f in &self.faces {
            let a = face_key(f.face_a);
            let b = face_key(f.face_b);
            if !a_faces.insert(a)
                || !b_faces.insert(b)
                || f.owner_a == f.owner_b
                || owners.get(&a).map(Vec::as_slice) != Some(&[f.owner_a])
                || owners.get(&b).map(Vec::as_slice) != Some(&[f.owner_b])
            {
                return Err(invalid());
            }
            let nodes = mesh.face_nodes(f.owner_a, &a)?;
            for (k, [i, j]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
                let mut edge = [nodes[i], nodes[j]];
                edge.sort_unstable();
                let entry = surface_edges.entry(edge).or_insert((0, nodes[k + 3]));
                entry.0 += 1;
                if entry.0 > 2 || entry.1 != nodes[k + 3] {
                    return Err(invalid());
                }
            }
        }
        if a_faces != fa.iter().copied().map(face_key).collect()
            || b_faces != fb.iter().copied().map(face_key).collect()
        {
            return Err(invalid());
        }
        let mut unique = BTreeSet::new();
        let mut degree = BTreeMap::<usize, usize>::new();
        let mut allowed = BTreeSet::new();
        for &[v, w, mid] in &self.edges {
            let mut edge = [v, w];
            edge.sort_unstable();
            if v == w
                || mid == v
                || mid == w
                || !unique.insert(edge)
                || surface_edges.get(&edge) != Some(&(1, mid))
            {
                return Err(invalid());
            }
            for n in [v, w] {
                *degree.entry(n).or_default() += 1;
            }
            allowed.extend([v, w, mid]);
            // Reconstruct the pre-split fan using physical glued faces and
            // declared contact pairs. A disconnected/branched fan is invalid.
            let star: BTreeSet<_> = mesh
                .tetrahedra
                .iter()
                .enumerate()
                .filter(|(_, t)| t.contains(&v) && t.contains(&w))
                .map(|(e, _)| e)
                .collect();
            let mut links: BTreeMap<usize, BTreeSet<usize>> =
                star.iter().map(|&e| (e, BTreeSet::new())).collect();
            let mut connect = |a, b| {
                links.get_mut(&a).ok_or_else(invalid)?.insert(b);
                links.get_mut(&b).ok_or_else(invalid)?.insert(a);
                Ok::<_, Error>(())
            };
            for (face, os) in &owners {
                if face.contains(&v) && face.contains(&w) && os.len() == 2 {
                    connect(os[0], os[1])?;
                }
            }
            for f in &self.faces {
                if f.face_a.contains(&v) && f.face_a.contains(&w) {
                    if !f.face_b.contains(&v) || !f.face_b.contains(&w) {
                        return Err(invalid());
                    }
                    connect(f.owner_a, f.owner_b)?;
                }
            }
            let mut stack = vec![*star.first().ok_or_else(invalid)?];
            let mut seen = BTreeSet::new();
            while let Some(e) = stack.pop() {
                if seen.insert(e) {
                    stack.extend(&links[&e]);
                }
            }
            if seen != star || links.values().any(|l| l.len() > 2) {
                return Err(invalid());
            }
        }
        if degree.values().any(|&d| d > 2) {
            return Err(invalid());
        }
        Ok(allowed)
    }
    fn matches(&self, ka: [usize; 3], kb: [usize; 3], oa: usize, ob: usize) -> bool {
        self.faces.iter().any(|f| {
            face_key(f.face_a) == ka
                && face_key(f.face_b) == kb
                && f.owner_a == oa
                && f.owner_b == ob
        })
    }
}

#[cfg(test)]
mod predictor_tests {
    use super::*;
    #[test]
    fn cold_elastic_prediction_matches_a_closed_sticking_physical_increment() -> Result<()> {
        let mesh: Mesh =
            serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
        let case: crate::stages::StageCase =
            serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
        let model = Model::new(&mesh, case.joint_system.unwrap())?;
        let upper: BTreeSet<_> = mesh
            .region_ids
            .iter()
            .enumerate()
            .filter(|(_, r)| **r == 2)
            .flat_map(|(e, _)| mesh.cell_nodes(e))
            .collect();
        let active = vec![true; mesh.tetrahedra.len()];
        let old = model.initial_state()?;
        let matrix = model.cold_activation_tangent(&old, &active)?.unwrap();
        let u: Vec<_> = (0..mesh.nodes.len())
            .map(|n| {
                if upper.contains(&n) {
                    [1e-6, -0.5e-6, -1e-5]
                } else {
                    [0.; 3]
                }
            })
            .collect();
        let actual = model.evaluate(&old, &u, &active, None, false)?;
        let mut predicted = vec![0.; 3 * u.len()];
        for (g, h, k) in matrix {
            predicted[g] += k * u[h / 3][h % 3];
        }
        let error = predicted
            .iter()
            .zip(&actual.internal_kn)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(error < 1e-12);
        assert!(
            model
                .cold_activation_tangent(&actual.state, &active)?
                .is_none()
        );
        assert!(
            old.points
                .iter()
                .all(|p| p.jump_m == [0.; 3] && p.traction_kpa == [0.; 3])
        );
        Ok(())
    }
    #[test]
    fn neutral_predictor_returns_an_overshot_sliding_trial_to_the_elastic_interior() -> Result<()> {
        let mesh: Mesh =
            serde_json::from_str(include_str!("../examples/joint-column-tet10-mesh.json"))?;
        let case: crate::stages::StageCase =
            serde_json::from_str(include_str!("../examples/stages-joint-column.json"))?;
        let model = Model::new(&mesh, case.joint_system.unwrap())?;
        let upper: BTreeSet<_> = mesh
            .region_ids
            .iter()
            .enumerate()
            .filter(|(_, r)| **r == 2)
            .flat_map(|(e, _)| mesh.cell_nodes(e))
            .collect();
        let mut rows = vec![vec![]; 3 * mesh.nodes.len()];
        let mut count = 0;
        for &n in &upper {
            for a in 0..3 {
                rows[3 * n + a].push((count, 1.));
                count += 1;
            }
        }
        let u: Vec<_> = (0..mesh.nodes.len())
            .map(|n| {
                if upper.contains(&n) {
                    [0.001, 0., -0.0001]
                } else {
                    [0.; 3]
                }
            })
            .collect();
        let active = vec![true; mesh.tetrahedra.len()];
        let old = model.initial_state()?;
        let current = model.evaluate(&old, &u, &active, None, false)?.state;
        let geometry = crate::rigidity::Geometry::new(&mesh, &active)?;
        let mut space = geometry.constrained(&rows)?;
        model.add_rigid_support(&current, &active, 1e-8, &mut space)?;
        assert_eq!(space.report().unrestrained_modes, 1);
        let modes = geometry.reduced_modes(&space.null_modes(), &rows, count)?;
        let mut residual = DVector::zeros(count);
        for p in &model.points {
            for &(node, c) in &p.coefficients {
                for &(r, w) in &rows[3 * node] {
                    residual[r] -= 10. * p.info.weight_m2 * c * w;
                }
            }
        }
        let eta = model.neutral_elastic_predictor(
            (&old, &current),
            &active,
            &rows,
            &modes,
            &DVector::zeros(count),
            modes.transpose() * residual,
        )?;
        let correction = &modes * eta;
        let mut candidate = u.clone();
        for (g, terms) in rows.iter().enumerate() {
            candidate[g / 3][g % 3] += terms.iter().map(|&(r, w)| w * correction[r]).sum::<f64>();
        }
        let corrected = model
            .evaluate(&old, &candidate, &active, None, false)?
            .state;
        for p in corrected.points {
            assert!((p.traction_kpa[0] - 20.).abs() < 1e-8);
            assert!(p.accumulated_plastic_slip_m < 1e-12);
        }
        assert!(
            current
                .points
                .iter()
                .all(|p| p.accumulated_plastic_slip_m > 0.)
        );
        assert!(old.points.iter().all(|p| p.jump_m == [0.; 3]));
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    /// Global equilibrium cutback at endpoint contact-status changes; absent is legacy.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub global_step_control: Option<GlobalStepControl>,
    pub schema: String,
    pub geometric_tolerance_m: f64,
    pub maximum_increment_m: f64,
    pub tangent_step_m: f64,
    pub joints: Vec<Definition>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub structure_bodies: Vec<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BirthPolicy {
    TransferForce,
    StressFree,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Generation {
    pub generation: u64,
    pub points: Vec<(usize, PointState)>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Lifecycle {
    pub joint_id: String,
    pub policy: Option<BirthPolicy>,
    pub active: bool,
    pub tied: bool,
    pub generation: u64,
    pub archived: Vec<Generation>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lifecycle: Vec<Lifecycle>,
    pub points: Vec<PointState>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointInfo {
    pub joint_id: String,
    pub face_index: usize,
    pub position_m: [f64; 3],
    pub normal: [f64; 3],
    pub weight_m2: f64,
    pub nodes_a: Vec<usize>,
    pub nodes_b: Vec<usize>,
    pub shape: Vec<f64>,
    pub owner_a: usize,
    pub owner_b: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointReport {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hydraulic: Option<crate::joint_water::PointReport>,
    pub index: usize,
    pub active: bool,
    pub gap_m: f64,
    pub opening_m: f64,
    pub normal_pressure_kpa: f64,
    pub shear_traction_kpa: [f64; 3],
    pub elastic_energy_kpa_m: f64,
    pub contact_reference_offset_m: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub water: Option<crate::joint_water::Field>,
    pub points: Vec<PointReport>,
    pub resultant_internal_force_kn: [f64; 3],
    pub resultant_internal_moment_kn_m: [f64; 3],
    pub elastic_energy_kn_m: f64,
    pub plastic_dissipation_kn_m: f64,
    pub numerical_dissipation_kn_m: f64,
    pub endpoint_work_kn_m: f64,
}
pub struct Evaluation {
    pub state: State,
    pub internal_kn: Vec<f64>,
    pub tangent: Vec<(usize, usize, f64)>,
    pub report: Report,
}
struct Point {
    info: PointInfo,
    material: usize,
    coefficients: Vec<(usize, f64)>,
}
pub struct Model {
    pub configuration: Configuration,
    points: Vec<Point>,
    positions: Vec<[f64; 3]>,
    element_count: usize,
    counterparts: BTreeMap<[usize; 3], usize>,
    pub face_keys: BTreeSet<[usize; 3]>,
    pub owner_pairs: Vec<[usize; 2]>,
}

impl Model {
    pub fn new(mesh: &Mesh, configuration: Configuration) -> Result<Self> {
        mesh.validate()?;
        let c = &configuration;
        if c.schema != "soil-fem-joints/1"
            || c.joints.is_empty()
            || !c.geometric_tolerance_m.is_finite()
            || !(1e-12..=1e-5).contains(&c.geometric_tolerance_m)
            || !c.maximum_increment_m.is_finite()
            || c.maximum_increment_m <= 0.
            || !c.tangent_step_m.is_finite()
            || c.tangent_step_m <= 0.
            || c.tangent_step_m > c.maximum_increment_m / 10.
        {
            return Err(Error::Input("invalid JOINT configuration".into()));
        }
        // The front is a property of the complete contact surface, even when
        // callers assign different JOINT IDs/materials to adjacent patches.
        let mut tip_edges = BTreeSet::new();
        let mut tip_degree = BTreeMap::<usize, usize>::new();
        for tip in c.joints.iter().filter_map(|d| d.crack_tip.as_ref()) {
            for &[a, b, _] in &tip.edges {
                let mut edge = [a, b];
                edge.sort_unstable();
                if !tip_edges.insert(edge) {
                    return Err(Error::Input(
                        "duplicate crack-tip edge across JOINTs".into(),
                    ));
                }
                for n in edge {
                    *tip_degree.entry(n).or_default() += 1;
                }
            }
        }
        if tip_degree.values().any(|&degree| degree > 2) {
            return Err(Error::Input("branching crack tip across JOINTs".into()));
        }
        if c.joints.iter().any(|d| d.birth_policy.is_some())
            && c.joints.iter().any(|d| d.crack_tip.is_some())
            && (c.joints.len() != 1 || c.joints[0].water.is_some())
        {
            return Err(Error::Input(
                "declared-tip JOINT BIRTH is limited to one dry interface".into(),
            ));
        }
        let owners = boundary::surface_owners(mesh);
        let mut ids = BTreeSet::new();
        let mut face_keys = BTreeSet::new();
        let mut points = Vec::new();
        let mut pairs = Vec::new();
        let mut counterparts = BTreeMap::new();
        for (mi, d) in c.joints.iter().enumerate() {
            d.material.validate()?;
            if d.id.is_empty() || !ids.insert(&d.id) || d.side_a == d.side_b {
                return Err(Error::Input(
                    "JOINT requires unique ID and different side sets".into(),
                ));
            }
            let fa = mesh
                .face_sets
                .get(&d.side_a)
                .filter(|f| !f.is_empty())
                .ok_or_else(|| Error::Input("JOINT side-a face set missing/empty".into()))?;
            let fb = mesh
                .face_sets
                .get(&d.side_b)
                .filter(|f| !f.is_empty())
                .ok_or_else(|| Error::Input("JOINT side-b face set missing/empty".into()))?;
            if fa.len() != fb.len() {
                return Err(Error::Input("JOINT matching face counts differ".into()));
            }
            let allowed_tip_nodes = match &d.crack_tip {
                Some(tip) => tip.validate(mesh, fa, fb)?,
                None => BTreeSet::new(),
            };
            let mut used = BTreeSet::new();
            for (fi, key_a) in fa.iter().enumerate() {
                let mut ka = *key_a;
                ka.sort_unstable();
                let oa = owners.get(&ka).filter(|o| o.len() == 1).ok_or_else(|| {
                    Error::Input("JOINT sides require separated displacement nodes".into())
                })?[0];
                let na = mesh.face_nodes(oa, &ka)?;
                let mut found = Vec::new();
                for (bi, key_b) in fb.iter().enumerate() {
                    // A declared tip must reject every geometric ambiguity,
                    // including candidates paired with an earlier face.
                    if used.contains(&bi) && d.crack_tip.is_none() {
                        continue;
                    }
                    let mut kb = *key_b;
                    kb.sort_unstable();
                    let ob = owners.get(&kb).filter(|o| o.len() == 1).ok_or_else(|| {
                        Error::Input("JOINT side is not an exterior side of its volume".into())
                    })?[0];
                    if oa == ob {
                        continue;
                    }
                    let nb = mesh.face_nodes(ob, &kb)?;
                    if nb.len() != na.len() {
                        continue;
                    }
                    let mut permutation = Vec::new();
                    for &n in &na {
                        let matches: Vec<_> = nb
                            .iter()
                            .enumerate()
                            .filter(|(_, b)| {
                                (Vector3::from(mesh.nodes[n]) - Vector3::from(mesh.nodes[**b]))
                                    .norm()
                                    <= c.geometric_tolerance_m
                            })
                            .map(|(i, _)| i)
                            .collect();
                        if matches.len() != 1 {
                            break;
                        }
                        permutation.push(matches[0]);
                    }
                    if permutation.len() == na.len()
                        && permutation.iter().collect::<BTreeSet<_>>().len() == na.len()
                    {
                        found.push((bi, kb, ob, nb, permutation));
                    }
                }
                if found.len() != 1 {
                    return Err(Error::Input(
                        "JOINT opposing face geometry missing or ambiguous".into(),
                    ));
                }
                let (bi, kb, ob, nb, permutation) = found.pop().unwrap();
                if permutation[..3].iter().any(|&i| i >= 3)
                    || na
                        .iter()
                        .enumerate()
                        .any(|(i, &n)| (n == nb[permutation[i]]) != allowed_tip_nodes.contains(&n))
                    || na.iter().enumerate().all(|(i, &n)| n == nb[permutation[i]])
                    || d.crack_tip
                        .as_ref()
                        .is_some_and(|tip| !tip.matches(ka, kb, oa, ob))
                {
                    return Err(Error::Input(
                        "JOINT interpolation, declared-tip sharing or owner pair mismatch".into(),
                    ));
                }
                let inversions = (0..3)
                    .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
                    .filter(|&(i, j)| permutation[i] > permutation[j])
                    .count();
                if inversions % 2 != 1 {
                    return Err(Error::Input("JOINT side normals must oppose".into()));
                }
                if na.len() == 6 {
                    for (i, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
                        let pair = BTreeSet::from([permutation[a], permutation[b]]);
                        let edge = [[0, 1], [1, 2], [2, 0]]
                            .iter()
                            .position(|e| BTreeSet::from(*e) == pair)
                            .unwrap();
                        if permutation[i + 3] != edge + 3 {
                            return Err(Error::Input(
                                "JOINT quadratic edge ordering mismatch".into(),
                            ));
                        }
                    }
                }
                if !face_keys.insert(ka) || !face_keys.insert(kb) {
                    return Err(Error::Input("JOINT face is selected more than once".into()));
                }
                used.insert(bi);
                pairs.push([oa, ob]);
                counterparts.insert(ka, ob);
                counterparts.insert(kb, oa);
                let nb: Vec<_> = permutation.iter().map(|&i| nb[i]).collect();
                let pa: Vec<_> = na.iter().map(|&n| mesh.nodes[n]).collect();
                for q in elements::face_points(&pa)? {
                    let area = Vector3::from(q.area_vector);
                    let weight = area.norm();
                    let normal = area / weight;
                    let mut coefficients = BTreeMap::<usize, f64>::new();
                    for ((&a, &b), &w) in na.iter().zip(&nb).zip(&q.shape) {
                        *coefficients.entry(a).or_default() -= w;
                        *coefficients.entry(b).or_default() += w;
                    }
                    // Exact +/- Ni cancellation removes shared tip DOFs from
                    // force, tangent and connectivity, including inactive history.
                    for &n in &allowed_tip_nodes {
                        if coefficients.get(&n).is_some_and(|&v| v != 0.) {
                            return Err(Error::Input(
                                "JOINT tip jump coefficient did not cancel".into(),
                            ));
                        }
                    }
                    coefficients.retain(|_, w| *w != 0.);
                    points.push(Point {
                        info: PointInfo {
                            joint_id: d.id.clone(),
                            face_index: fi,
                            position_m: q.position,
                            normal: normal.into(),
                            weight_m2: weight,
                            nodes_a: na.clone(),
                            nodes_b: nb.clone(),
                            shape: q.shape,
                            owner_a: oa,
                            owner_b: ob,
                        },
                        material: mi,
                        coefficients: coefficients.into_iter().collect(),
                    });
                }
            }
        }
        Ok(Self {
            configuration,
            points,
            positions: mesh.nodes.clone(),
            element_count: mesh.tetrahedra.len(),
            counterparts,
            face_keys,
            owner_pairs: pairs,
        })
    }
    pub fn point_count(&self) -> usize {
        self.points.len()
    }
    pub fn counterpart(&self, face: &[usize; 3]) -> Option<usize> {
        self.counterparts.get(face).copied()
    }
    pub fn shared_tip_nodes(&self) -> BTreeSet<usize> {
        self.configuration
            .joints
            .iter()
            .filter_map(|d| d.crack_tip.as_ref())
            .flat_map(|tip| tip.edges.iter().flatten().copied())
            .collect()
    }
    fn validate_layout(&self, state: &State, u: &[[f64; 3]], active: &[bool]) -> Result<()> {
        if state.points.len() != self.points.len()
            || u.len() != self.positions.len()
            || active.len() != self.element_count
            || !u.iter().flatten().all(|v| v.is_finite())
        {
            return Err(Error::Input(
                "JOINT history/displacement/active layout mismatch".into(),
            ));
        }
        self.validate_lifecycle(state)?;
        Ok(())
    }
    pub fn layout(&self) -> Vec<PointInfo> {
        self.points.iter().map(|p| p.info.clone()).collect()
    }
    pub fn initial_state(&self) -> Result<State> {
        Ok(State {
            lifecycle: if self
                .configuration
                .joints
                .iter()
                .any(|d| d.birth_policy.is_some())
            {
                self.configuration
                    .joints
                    .iter()
                    .map(|d| Lifecycle {
                        joint_id: d.id.clone(),
                        policy: d.birth_policy,
                        active: d.birth_policy.is_none(),
                        tied: d.birth_policy.is_some(),
                        generation: 0,
                        archived: vec![],
                    })
                    .collect()
            } else {
                vec![]
            },
            points: self
                .points
                .iter()
                .map(|p| PointState::initial(p.info.normal))
                .collect::<Result<_>>()?,
        })
    }
    /// Exact assembled jump coefficients, for topology/work audits.
    pub fn jump_coefficients(&self) -> Vec<Vec<(usize, f64)>> {
        self.points.iter().map(|p| p.coefficients.clone()).collect()
    }
    pub fn connectivity(&self) -> Vec<Vec<usize>> {
        self.points
            .iter()
            .map(|p| {
                p.coefficients
                    .iter()
                    .flat_map(|&(n, _)| (0..3).map(move |a| 3 * n + a))
                    .collect()
            })
            .collect()
    }
    fn contact_active(&self, state: &State, p: &Point, active: &[bool]) -> bool {
        Self::active(p, active) && state.lifecycle.get(p.material).is_none_or(|l| l.active)
    }
    fn active(p: &Point, active: &[bool]) -> bool {
        active[p.info.owner_a] && active[p.info.owner_b]
    }
    fn jump(p: &Point, u: &[[f64; 3]]) -> Vector3<f64> {
        p.coefficients
            .iter()
            .fold(Vector3::zeros(), |v, &(n, w)| v + w * Vector3::from(u[n]))
    }
    pub fn partitions(&self, state: &State, u: &[[f64; 3]], active: &[bool]) -> Result<Vec<usize>> {
        self.validate_layout(state, u, active)?;
        if state.points.len() != self.points.len() || u.len() != self.positions.len() {
            return Err(Error::Input(
                "JOINT state/displacement length mismatch".into(),
            ));
        }
        self.points
            .iter()
            .zip(&state.points)
            .map(|(p, s)| {
                if self.contact_active(state, p, active) {
                    self.configuration.joints[p.material].material.partition(
                        s,
                        Self::jump(p, u).into(),
                        self.configuration.maximum_increment_m,
                    )
                } else {
                    Ok(1)
                }
            })
            .collect()
    }
    pub fn validate_state(&self, state: &State, u: &[[f64; 3]], active: &[bool]) -> Result<()> {
        self.validate_layout(state, u, active)?;
        self.validate_lifecycle(state)?;
        for (p, s) in self.points.iter().zip(&state.points) {
            self.configuration.joints[p.material]
                .material
                .validate_state(s)?;
            if (Vector3::from(s.normal) - Vector3::from(p.info.normal)).norm() > 1e-10
                || (self.contact_active(state, p, active)
                    && (Self::jump(p, u) - Vector3::from(s.jump_m)).norm() > 1e-10)
            {
                return Err(Error::Input("JOINT geometry/jump cache mismatch".into()));
            }
        }
        Ok(())
    }
    /// Sustained, one-sided loading stiffness. The central-difference tangent
    /// at zero pressure or at the friction limit is not a support certificate.
    pub fn add_rigid_support(
        &self,
        state: &State,
        active: &[bool],
        pressure_tolerance_kpa: f64,
        space: &mut crate::rigidity::Space<'_>,
    ) -> Result<()> {
        if state.points.len() != self.points.len()
            || active.len() != self.element_count
            || !pressure_tolerance_kpa.is_finite()
            || pressure_tolerance_kpa < 1e-10
        {
            return Err(Error::Input("JOINT rigid support layout mismatch".into()));
        }
        for (p, s) in self.points.iter().zip(&state.points) {
            if space.full() {
                break;
            }
            if !self.contact_active(state, p, active) {
                continue;
            }
            let n = Vector3::from(p.info.normal);
            let q = Vector3::from(s.traction_kpa);
            let pressure = -n.dot(&q);
            if pressure <= pressure_tolerance_kpa {
                continue;
            }
            let mut add = |direction: Vector3<f64>| {
                // Use the same nodal jump operator as force assembly. Inferring
                // rigid jumps from rounded world quadrature positions can
                // invent constraints on global rotation at large coordinates.
                let row: Vec<_> = p
                    .coefficients
                    .iter()
                    .flat_map(|&(node, w)| (0..3).map(move |a| (3 * node + a, w * direction[a])))
                    .collect();
                space.physical_row(&row)
            };
            add(n)?;
            let cap = self.configuration.joints[p.material]
                .material
                .friction_coefficient
                * pressure;
            let shear_tolerance = pressure_tolerance_kpa
                * (1.
                    + self.configuration.joints[p.material]
                        .material
                        .friction_coefficient);
            if cap <= shear_tolerance {
                continue;
            }
            let tau = q + n * pressure;
            if tau.norm() < cap - shear_tolerance.max(cap * 1e-8) {
                let axis = if n[0].abs() < 0.8 {
                    Vector3::x()
                } else {
                    Vector3::y()
                };
                let tangent = n.cross(&axis).normalize();
                add(tangent)?;
                add(n.cross(&tangent))?;
            } else if tau.norm() > 0. {
                // At incipient/sliding Coulomb contact, the radial loading
                // direction has no restoring stiffness (psi = 0).
                add(n.cross(&tau).normalize())?;
            }
        }
        Ok(())
    }
    /// An elastic predictor restricted to current Newton-neutral rigid modes.
    /// It can return an overshot sliding trial to the friction cone interior.
    /// This contributes no force, stiffness, energy or history to equilibrium.
    pub(crate) fn neutral_elastic_predictor(
        &self,
        history: (&State, &State),
        active: &[bool],
        rows: &[Vec<(usize, f64)>],
        modes: &DMatrix<f64>,
        range: &DVector<f64>,
        mut rhs: DVector<f64>,
    ) -> Result<DVector<f64>> {
        let (old, current) = history;
        if old.points.len() != self.points.len()
            || current.points.len() != self.points.len()
            || active.len() != self.element_count
            || rows.len() != 3 * self.positions.len()
            || modes.nrows() != range.len()
            || modes.ncols() != rhs.len()
            || rows
                .iter()
                .flatten()
                .any(|&(i, c)| i >= range.len() || !c.is_finite())
        {
            return Err(Error::Input(
                "neutral contact predictor layout mismatch".into(),
            ));
        }
        let mut matrix = DMatrix::zeros(rhs.len(), rhs.len());
        for (i, p) in self.points.iter().enumerate() {
            if !self.contact_active(current, p, active) {
                continue;
            }
            let n = Vector3::from(p.info.normal);
            let s = &current.points[i];
            let o = &old.points[i];
            let jump = Vector3::from(s.jump_m);
            let old_jump = Vector3::from(o.jump_m);
            let gap = n.dot(&s.elastic_jump());
            let old_gap = n.dot(&o.elastic_jump());
            if gap > 0. {
                continue;
            }
            let fraction = if old_gap > 0. {
                (gap / (gap - old_gap)).clamp(0., 1.)
            } else {
                1.
            };
            let material = &self.configuration.joints[p.material].material;
            let projection = Matrix3::identity() - n * n.transpose();
            let k = projection * (material.tangential_stiffness_kpa_m * fraction);
            let unlimited = if old_gap > 0. {
                Vector3::zeros()
            } else {
                projection * Vector3::from(o.traction_kpa)
            } + k * (jump - old_jump);
            let extra = unlimited - projection * Vector3::from(s.traction_kpa);
            let mut c = DMatrix::zeros(3, rhs.len());
            let mut du = Vector3::zeros();
            for &(node, shape) in &p.coefficients {
                for a in 0..3 {
                    for &(r, w) in &rows[3 * node + a] {
                        du[a] += shape * w * range[r];
                        for j in 0..rhs.len() {
                            c[(a, j)] += shape * w * modes[(r, j)];
                        }
                    }
                }
            }
            rhs -= c.transpose() * (extra + k * du) * p.info.weight_m2;
            matrix += c.transpose() * k * c * p.info.weight_m2;
        }
        if !matrix.iter().chain(rhs.iter()).all(|v| v.is_finite()) {
            return Err(Error::Solve("nonfinite neutral elastic predictor".into()));
        }
        let eig = matrix.symmetric_eigen();
        let scale = eig.eigenvalues.iter().map(|v| v.abs()).fold(0., f64::max);
        let mut eta = DVector::zeros(rhs.len());
        if scale > 0. {
            for (i, &value) in eig.eigenvalues.iter().enumerate() {
                if value > scale * 1e-10 {
                    let column = eig.eigenvectors.column(i);
                    eta += column * (column.dot(&rhs) / value);
                }
            }
        }
        Ok(eta)
    }
    /// At an unstressed, exactly touching interface, the Coulomb law has no
    /// single linear derivative. A bilateral elastic trial predicts closure;
    /// its force is never used to certify equilibrium or commit history.
    pub(crate) fn cold_activation_tangent(
        &self,
        state: &State,
        active: &[bool],
    ) -> Result<Option<Vec<(usize, usize, f64)>>> {
        if state.points.len() != self.points.len() || active.len() != self.element_count {
            return Err(Error::Input(
                "cold contact predictor layout mismatch".into(),
            ));
        }
        let mut found = false;
        for (p, s) in self.points.iter().zip(&state.points) {
            if self.contact_active(state, p, active) {
                found = true;
                if s.elastic_jump() != Vector3::zeros() || s.traction_kpa != [0.; 3] {
                    return Ok(None);
                }
            }
        }
        if !found {
            return Ok(None);
        }
        let mut matrix = Vec::new();
        for p in &self.points {
            if !self.contact_active(state, p, active) {
                continue;
            }
            let material = &self.configuration.joints[p.material].material;
            let n = Vector3::from(p.info.normal);
            let normal = n * n.transpose();
            let d = normal * material.normal_stiffness_kpa_m
                + (Matrix3::identity() - normal) * material.tangential_stiffness_kpa_m;
            for &(node, c) in &p.coefficients {
                for &(other, v) in &p.coefficients {
                    for a in 0..3 {
                        for b in 0..3 {
                            let k = p.info.weight_m2 * c * v * d[(a, b)];
                            if k != 0. {
                                matrix.push((3 * node + a, 3 * other + b, k));
                            }
                        }
                    }
                }
            }
        }
        Ok(Some(matrix))
    }
    pub fn evaluate(
        &self,
        state: &State,
        u: &[[f64; 3]],
        active: &[bool],
        partitions: Option<&[usize]>,
        tangent: bool,
    ) -> Result<Evaluation> {
        self.validate_layout(state, u, active)?;
        if state.points.len() != self.points.len()
            || u.len() != self.positions.len()
            || partitions.is_some_and(|p| p.len() != self.points.len())
        {
            return Err(Error::Input("JOINT evaluation layout mismatch".into()));
        }
        let mut next = state.clone();
        let mut internal = vec![0.; 3 * self.positions.len()];
        let mut matrix = Vec::new();
        let mut reports = Vec::new();
        let mut energy = 0.;
        let mut plastic = 0.;
        let mut numerical = 0.;
        let mut work = 0.;
        for (i, p) in self.points.iter().enumerate() {
            let is_active = self.contact_active(state, p, active);
            let material = &self.configuration.joints[p.material].material;
            let old = &state.points[i];
            // Removed-side histories remain archived and must still be valid.
            material.validate_state(old)?;
            if (Vector3::from(old.normal) - Vector3::from(p.info.normal)).norm() > 1e-10 {
                return Err(Error::Input("JOINT normal/history mismatch".into()));
            }
            if is_active {
                let jump = Self::jump(p, u);
                let count = if let Some(counts) = partitions {
                    counts[i]
                } else {
                    material.partition(old, jump.into(), self.configuration.maximum_increment_m)?
                };
                next.points[i] = material.integrate_partition(old, jump.into(), count)?;
                let q = Vector3::from(next.points[i].traction_kpa);
                for &(n, c) in &p.coefficients {
                    for a in 0..3 {
                        internal[3 * n + a] += p.info.weight_m2 * c * q[a];
                    }
                }
                if tangent {
                    let d = material.tangent(
                        old,
                        jump.into(),
                        count,
                        self.configuration.tangent_step_m,
                    )?;
                    for &(n, c) in &p.coefficients {
                        for &(m, v) in &p.coefficients {
                            for a in 0..3 {
                                for b in 0..3 {
                                    let k = p.info.weight_m2 * c * v * d[(a, b)];
                                    if k != 0. {
                                        matrix.push((3 * n + a, 3 * m + b, k));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            let s = &next.points[i];
            let n = Vector3::from(p.info.normal);
            let q = Vector3::from(s.traction_kpa);
            let gn = s.elastic_jump().dot(&n);
            let normal = q.dot(&n);
            let tau = q - normal * n;
            let elastic = material.energy(s);
            if is_active {
                energy += elastic * p.info.weight_m2;
                plastic += s.plastic_dissipation_kpa_m * p.info.weight_m2;
                numerical += s.numerical_dissipation_kpa_m * p.info.weight_m2;
                work += s.endpoint_work_kpa_m * p.info.weight_m2;
            }
            let tangential = Vector3::from(s.jump_m) - gn * n;
            reports.push(PointReport {
                hydraulic: None,
                index: i,
                active: is_active,
                gap_m: gn,
                opening_m: gn.max(0.),
                normal_pressure_kpa: -normal,
                shear_traction_kpa: tau.into(),
                elastic_energy_kpa_m: elastic,
                contact_reference_offset_m: (tangential
                    - tau / material.tangential_stiffness_kpa_m
                    - Vector3::from(s.plastic_slip_m))
                .into(),
            });
        }
        let mut f = Vector3::zeros();
        let mut moment = Vector3::zeros();
        for (n, x) in self.positions.iter().enumerate() {
            let q = Vector3::new(internal[3 * n], internal[3 * n + 1], internal[3 * n + 2]);
            f += q;
            moment += Vector3::from(*x).cross(&q);
        }
        Ok(Evaluation {
            state: next,
            internal_kn: internal,
            tangent: matrix,
            report: Report {
                water: None,
                points: reports,
                resultant_internal_force_kn: f.into(),
                resultant_internal_moment_kn_m: moment.into(),
                elastic_energy_kn_m: energy,
                plastic_dissipation_kn_m: plastic,
                numerical_dissipation_kn_m: numerical,
                endpoint_work_kn_m: work,
            },
        })
    }
}
