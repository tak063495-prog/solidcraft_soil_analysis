//! Small-strain tension-only tendons: unbonded free span and perfectly bonded
//! segmented anchorage. Force-controlled jacking is released by exact lock-off.
use crate::{
    Error, Result,
    embedding::{Binding, Point},
    mesh::Mesh,
};
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub id: String,
    pub young_kpa: f64,
    pub area_m2: f64,
    pub head: Point,
    /// First point is the start of grout. Every successive segment is bonded.
    pub bonded_path: Vec<Point>,
    pub initially_active: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    pub schema: String,
    pub embedding_tolerance_m: f64,
    pub anchors: Vec<Definition>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnchorState {
    pub active: bool,
    pub reference_extensions_m: Vec<f64>,
    pub force_control_kn: Option<f64>,
    pub last_lock_off_force_kn: f64,
    pub segment_forces_kn: Vec<f64>,
    pub segment_taut: Vec<bool>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub anchors: Vec<AnchorState>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnchorResult {
    pub id: String,
    pub point_positions_m: Vec<[f64; 3]>,
    pub point_displacement_m: Vec<[f64; 3]>,
    /// Forces exerted by the tendon on the soil, opposite to FEM internal force.
    pub point_soil_force_kn: Vec<[f64; 3]>,
    pub free_length_m: f64,
    pub bonded_length_m: f64,
    pub segment_lengths_m: Vec<f64>,
    pub segment_extension_m: Vec<f64>,
    pub segment_forces_kn: Vec<f64>,
    pub segment_taut: Vec<bool>,
    pub force_control_active: bool,
    pub elastic_energy_kn_m: f64,
    pub bindings: Vec<Binding>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub anchors: Vec<AnchorResult>,
    pub resultant_internal_force_kn: [f64; 3],
    pub resultant_internal_moment_kn_m: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct Control {
    pub forces_kn: BTreeMap<String, f64>,
    pub lock_off: bool,
}
pub struct Evaluation {
    pub state: State,
    pub internal_kn: Vec<f64>,
    pub tangent: Vec<(usize, usize, f64)>,
    pub report: Report,
}
struct Anchor {
    bindings: Vec<Binding>,
    axis: Vector3<f64>,
    lengths: Vec<f64>,
    stiffness: Vec<f64>,
    coefficients: Vec<Vec<(usize, f64)>>,
}
pub struct Model {
    pub configuration: Configuration,
    anchors: Vec<Anchor>,
    positions: Vec<[f64; 3]>,
    element_count: usize,
}

impl Model {
    pub fn new(mesh: &Mesh, configuration: Configuration) -> Result<Self> {
        if configuration.schema != "soil-fem-anchors/1"
            || configuration.anchors.is_empty()
            || !configuration.embedding_tolerance_m.is_finite()
            || configuration.embedding_tolerance_m <= 0.
            || configuration.embedding_tolerance_m > 1e-5
        {
            return Err(Error::Input("invalid anchor configuration".into()));
        }
        mesh.validate()?;
        let mut ids = BTreeSet::new();
        let mut anchors = Vec::new();
        for d in &configuration.anchors {
            if d.id.is_empty()
                || !ids.insert(&d.id)
                || !d.young_kpa.is_finite()
                || d.young_kpa <= 0.
                || !d.area_m2.is_finite()
                || d.area_m2 <= 0.
                || d.bonded_path.len() < 2
            {
                return Err(Error::Input(
                    "anchor requires unique ID, positive EA and at least two bonded points".into(),
                ));
            }
            let points: Vec<_> = std::iter::once(&d.head).chain(&d.bonded_path).collect();
            let bindings = points
                .iter()
                .map(|p| Binding::locate(mesh, p, configuration.embedding_tolerance_m))
                .collect::<Result<Vec<_>>>()?;
            let first = Vector3::from(d.head.position_m);
            let last = Vector3::from(points.last().unwrap().position_m);
            let chord = last - first;
            let total_length = chord.norm();
            if total_length <= 1e-9 {
                return Err(Error::Input("zero anchor length".into()));
            }
            let axis = chord / total_length;
            let mut lengths = Vec::new();
            let mut coefficients = Vec::new();
            for pair in bindings.windows(2) {
                let delta = Vector3::from(pair[1].position_m) - Vector3::from(pair[0].position_m);
                let length = delta.dot(&axis);
                if length <= 1e-9
                    || (delta - length * axis).norm() > configuration.embedding_tolerance_m
                {
                    return Err(Error::Input(
                        "anchor points must be ordered along one straight CAD axis".into(),
                    ));
                }
                lengths.push(length);
                let mut c = BTreeMap::<usize, f64>::new();
                for (point, sign) in [(0, -1.), (1, 1.)] {
                    for &(node, weight) in &pair[point].weights {
                        for a in 0..3 {
                            *c.entry(3 * node + a).or_default() += sign * weight * axis[a];
                        }
                    }
                }
                coefficients.push(c.into_iter().filter(|(_, v)| v.abs() > 1e-15).collect());
            }
            let stiffness = lengths
                .iter()
                .map(|l| d.young_kpa * d.area_m2 / l)
                .collect();
            anchors.push(Anchor {
                bindings,
                axis,
                lengths,
                stiffness,
                coefficients,
            });
        }
        Ok(Self {
            configuration,
            anchors,
            positions: mesh.nodes.clone(),
            element_count: mesh.tetrahedra.len(),
        })
    }
    pub fn initial_state(&self) -> State {
        State {
            anchors: self
                .anchors
                .iter()
                .zip(&self.configuration.anchors)
                .map(|(a, d)| AnchorState {
                    active: d.initially_active,
                    reference_extensions_m: vec![0.; a.lengths.len()],
                    force_control_kn: None,
                    last_lock_off_force_kn: 0.,
                    segment_forces_kn: vec![0.; a.lengths.len()],
                    segment_taut: vec![true; a.lengths.len()],
                })
                .collect(),
        }
    }
    pub fn validate_state(&self, state: &State) -> Result<()> {
        if state.anchors.len() != self.anchors.len() {
            return Err(Error::Input("anchor history count mismatch".into()));
        }
        for (s, a) in state.anchors.iter().zip(&self.anchors) {
            if s.reference_extensions_m.len() != a.lengths.len()
                || s.segment_forces_kn.len() != a.lengths.len()
                || s.segment_taut.len() != a.lengths.len()
                || !s.reference_extensions_m.iter().all(|v| v.is_finite())
                || !s
                    .segment_forces_kn
                    .iter()
                    .chain([&s.last_lock_off_force_kn])
                    .all(|v| v.is_finite() && *v >= 0.)
                || s.force_control_kn.is_some_and(|v| !v.is_finite() || v < 0.)
            {
                return Err(Error::Input(
                    "invalid anchor reference/control/force history".into(),
                ));
            }
        }
        Ok(())
    }
    pub fn connectivity(&self) -> Vec<Vec<usize>> {
        self.anchors
            .iter()
            .flat_map(|a| {
                a.coefficients
                    .iter()
                    .map(|c| c.iter().map(|(g, _)| *g).collect())
            })
            .collect()
    }
    fn selection(&self, names: &[String]) -> Result<Vec<usize>> {
        if names.is_empty() || names.iter().collect::<BTreeSet<_>>().len() != names.len() {
            return Err(Error::Input(
                "anchor selection must be nonempty and unique".into(),
            ));
        }
        names
            .iter()
            .map(|name| {
                self.configuration
                    .anchors
                    .iter()
                    .position(|a| &a.id == name)
                    .ok_or_else(|| Error::Input("unknown anchor ID".into()))
            })
            .collect()
    }
    pub fn change(
        &self,
        state: &State,
        u: &[[f64; 3]],
        names: &[String],
        install: bool,
    ) -> Result<State> {
        self.validate_state(state)?;
        if u.len() != self.positions.len() {
            return Err(Error::Input("anchor displacement size mismatch".into()));
        }
        let mut next = state.clone();
        for i in self.selection(names)? {
            let s = &mut next.anchors[i];
            let a = &self.anchors[i];
            if s.active == install || s.force_control_kn.is_some() {
                return Err(Error::Input(
                    "invalid anchor activation or active jack".into(),
                ));
            }
            s.active = install;
            if install {
                s.reference_extensions_m = a
                    .coefficients
                    .iter()
                    .map(|c| c.iter().map(|&(g, v)| v * u[g / 3][g % 3]).sum())
                    .collect();
                s.segment_forces_kn.fill(0.);
                s.segment_taut.fill(true);
                s.last_lock_off_force_kn = 0.;
            }
        }
        Ok(next)
    }
    pub fn control(
        &self,
        origin: &State,
        names: &[String],
        target: f64,
        progress: f64,
    ) -> Result<Control> {
        self.validate_state(origin)?;
        if !target.is_finite()
            || target < 0.
            || !progress.is_finite()
            || !(0. ..=1.).contains(&progress)
        {
            return Err(Error::Input("invalid anchor jack target/progress".into()));
        }
        let mut forces = BTreeMap::new();
        for i in self.selection(names)? {
            let s = &origin.anchors[i];
            if !s.active || s.force_control_kn.is_some() {
                return Err(Error::Input(
                    "jacking requires an installed, locked anchor".into(),
                ));
            }
            forces.insert(
                self.configuration.anchors[i].id.clone(),
                s.segment_forces_kn[0] + progress * (target - s.segment_forces_kn[0]),
            );
        }
        Ok(Control {
            forces_kn: forces,
            lock_off: progress == 1.,
        })
    }
    pub fn add_rigid_support(
        &self,
        state: &State,
        space: &mut crate::rigidity::Space<'_>,
    ) -> Result<()> {
        self.validate_state(state)?;
        for (a, s) in self.anchors.iter().zip(&state.anchors) {
            if !s.active {
                continue;
            }
            for (j, coefficients) in a.coefficients.iter().enumerate() {
                // A force-controlled segment supplies load, not stiffness.
                // A zero-force tendon can shorten freely into slackness.
                if s.segment_forces_kn[j] > 1e-10 && !(j == 0 && s.force_control_kn.is_some()) {
                    space.physical_row(coefficients)?;
                }
            }
        }
        Ok(())
    }
    pub fn evaluate(
        &self,
        state: &State,
        u: &[[f64; 3]],
        active: &[bool],
        control: Option<&Control>,
    ) -> Result<Evaluation> {
        self.validate_state(state)?;
        if u.len() != self.positions.len()
            || active.len() != self.element_count
            || !u.iter().flatten().all(|v| v.is_finite())
        {
            return Err(Error::Input(
                "anchor displacement/host size or nonfinite values".into(),
            ));
        }
        if let Some(c) = control {
            for (name, force) in &c.forces_kn {
                let i = self.selection(std::slice::from_ref(name))?[0];
                if !state.anchors[i].active || !force.is_finite() || *force < 0. {
                    return Err(Error::Input("invalid controlled anchor".into()));
                }
            }
        }
        let mut next = state.clone();
        let mut internal = vec![0.; self.positions.len() * 3];
        let mut tangent = Vec::new();
        let mut results = Vec::new();
        for (i, a) in self.anchors.iter().enumerate() {
            let old = &state.anchors[i];
            if !old.active {
                continue;
            }
            if a.bindings.iter().any(|p| !p.is_active(active)) {
                return Err(Error::Input(
                    "active anchor lost its embedded host; remove anchor before excavation".into(),
                ));
            }
            let id = &self.configuration.anchors[i].id;
            let prescribed = control
                .and_then(|c| c.forces_kn.get(id).copied())
                .or(old.force_control_kn);
            let s = &mut next.anchors[i];
            let mut point_forces = vec![Vector3::zeros(); a.bindings.len()];
            let mut extensions = Vec::new();
            let mut energy = 0.;
            for (j, c) in a.coefficients.iter().enumerate() {
                let extension = c.iter().map(|&(g, v)| v * u[g / 3][g % 3]).sum::<f64>();
                extensions.push(extension);
                let elastic = a.stiffness[j] * (extension - old.reference_extensions_m[j]);
                let (force, k, taut) = if let (0, Some(force)) = (j, prescribed) {
                    (force, 0., force > 0.)
                } else {
                    (
                        elastic.max(0.),
                        if elastic >= 0. { a.stiffness[j] } else { 0. },
                        elastic >= 0.,
                    )
                };
                s.segment_forces_kn[j] = force;
                s.segment_taut[j] = taut;
                energy += force * force / (2. * a.stiffness[j]);
                for &(g, v) in c {
                    internal[g] += force * v;
                    if k != 0. {
                        for &(h, w) in c {
                            tangent.push((g, h, k * v * w));
                        }
                    }
                }
                point_forces[j] += force * a.axis;
                point_forces[j + 1] -= force * a.axis;
            }
            s.force_control_kn = prescribed;
            if control.is_some_and(|c| c.lock_off && c.forces_kn.contains_key(id)) {
                let force = s.segment_forces_kn[0];
                // At zero tension the stress-free length is not unique. Do
                // not shorten an already slack tendon without positive work.
                s.reference_extensions_m[0] = if force == 0. {
                    old.reference_extensions_m[0].max(extensions[0])
                } else {
                    extensions[0] - force / a.stiffness[0]
                };
                s.force_control_kn = None;
                s.last_lock_off_force_kn = force;
                if force == 0. {
                    s.segment_taut[0] = extensions[0] >= s.reference_extensions_m[0];
                }
            }
            results.push(AnchorResult {
                id: id.clone(),
                point_positions_m: a.bindings.iter().map(|b| b.position_m).collect(),
                point_displacement_m: a
                    .bindings
                    .iter()
                    .map(|b| b.displacement(u).into())
                    .collect(),
                point_soil_force_kn: point_forces.iter().map(|v| (*v).into()).collect(),
                free_length_m: a.lengths[0],
                bonded_length_m: a.lengths[1..].iter().sum(),
                segment_lengths_m: a.lengths.clone(),
                segment_extension_m: extensions,
                segment_forces_kn: s.segment_forces_kn.clone(),
                segment_taut: s.segment_taut.clone(),
                force_control_active: s.force_control_kn.is_some(),
                elastic_energy_kn_m: energy,
                bindings: a.bindings.clone(),
            });
        }
        let mut force = Vector3::zeros();
        let mut moment = Vector3::zeros();
        for (n, x) in self.positions.iter().enumerate() {
            let f = Vector3::new(internal[3 * n], internal[3 * n + 1], internal[3 * n + 2]);
            force += f;
            moment += Vector3::from(*x).cross(&f);
        }
        Ok(Evaluation {
            state: next,
            internal_kn: internal,
            tangent,
            report: Report {
                anchors: results,
                resultant_internal_force_kn: force.into(),
                resultant_internal_moment_kn_m: moment.into(),
            },
        })
    }
}
