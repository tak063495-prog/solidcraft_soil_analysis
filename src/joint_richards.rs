//! RT0 traces on matching displacement faces. Each trace belongs to an
//! element and its opposite local vertex, independent of the displacement IDs.
use super::*;
use crate::{elements, richards};

/// Ordered [physical seconds, total head in metres] knots of a reservoir.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct HeadHistory(pub Vec<[f64; 2]>);
impl HeadHistory {
    pub(crate) fn validate(&self, end: f64) -> Result<()> {
        if self.0.len() < 2
            || !self.0.iter().flatten().all(|v| v.is_finite())
            || self.0.windows(2).any(|w| w[0][0] >= w[1][0])
            || self.0[0][0] > 0.
            || self.0.last().unwrap()[0] < end
        {
            return Err(Error::Input(
                "reservoir head table must cover all physical times with ordered finite knots"
                    .into(),
            ));
        }
        Ok(())
    }
    pub(crate) fn at(&self, time: f64) -> Result<f64> {
        let i = self.0.partition_point(|v| v[0] < time);
        if let Some(&[t, h]) = self.0.get(i)
            && t == time
        {
            return Ok(h);
        }
        if i == 0 || i >= self.0.len() {
            return Err(Error::Input("reservoir time outside table".into()));
        }
        let [t0, h0] = self.0[i - 1];
        let [t1, h1] = self.0[i];
        Ok(h0 + (time - t0) / (t1 - t0) * (h1 - h0))
    }
}
#[derive(Clone)]
pub(crate) struct Wall {
    pub corners: [usize; 3],
    pub element: usize,
    pub local: usize,
}

#[derive(Clone)]
pub(crate) struct Pair {
    pub a: Wall,
    pub b: Wall,
    pub policy: Policy,
    pub history: Option<HeadHistory>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Trace {
    pub face: usize,
    pub owner: usize,
    pub local_face: usize,
    pub total_head_m: f64,
    pub pressure_head_m: f64,
    pub signed_pressure_kpa: f64,
    pub outward_flow_m3_s: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PointTraces {
    pub joint_id: String,
    pub face_index: usize,
    pub time_s: f64,
    pub elevation_m: f64,
    pub side_a: Option<Trace>,
    pub side_b: Option<Trace>,
    pub applied_pressure_kpa: f64,
}

impl Model {
    pub(crate) fn empty_field(&self) -> Result<Field> {
        self.field(
            vec![0.; self.layout.len()],
            self.layout
                .iter()
                .map(|p| {
                    let di = self
                        .configuration
                        .joints
                        .iter()
                        .position(|d| d.id == p.joint_id)
                        .unwrap();
                    self.source(di)
                })
                .collect(),
            self.faces
                .iter()
                .map(|f| vec![[0.; 3]; f.nodes_a.len()])
                .collect(),
        )
    }

    pub(crate) fn richards_pairs(&self, active: &[bool]) -> Result<Vec<Pair>> {
        self.soil_domain(active)?;
        let wall = |owner: usize, nodes: &[usize]| -> Result<Wall> {
            let mut corners = [nodes[0], nodes[1], nodes[2]];
            corners.sort_unstable();
            let local = elements::FACES
                .iter()
                .position(|v| {
                    let mut key = v.map(|i| self.mesh.tetrahedra[owner][i]);
                    key.sort_unstable();
                    key == corners
                })
                .ok_or_else(|| Error::Input("wet JOINT face/owner mismatch".into()))?;
            Ok(Wall {
                corners,
                element: owner,
                local,
            })
        };
        self.faces
            .iter()
            .filter(|f| active[f.owner_a] && active[f.owner_b])
            .map(|f| {
                Ok(Pair {
                    a: wall(f.owner_a, &f.nodes_a)?,
                    b: wall(f.owner_b, &f.nodes_b)?,
                    history: self
                        .histories
                        .get(&self.configuration.joints[f.definition].id)
                        .cloned(),
                    policy: self.configuration.joints[f.definition]
                        .water
                        .clone()
                        .ok_or_else(|| Error::Input("missing Richards JOINT policy".into()))?,
                })
            })
            .collect()
    }

    /// Water-filled void volume at the common RT0 trace. The normal jump is
    /// positive in opening. Negative aperture is outside this constitutive law.
    pub(crate) fn saturated_volumes(
        &self,
        b0: f64,
        joint: &joints::State,
        traces: &[PointTraces],
    ) -> Result<BTreeMap<usize, f64>> {
        let mut indices = BTreeMap::new();
        for t in traces {
            let side = t
                .side_a
                .as_ref()
                .ok_or_else(|| Error::Input("missing wet joint wall".into()))?;
            indices.insert((t.joint_id.as_str(), t.face_index), side);
        }
        if self.layout.len() != joint.points.len() {
            return Err(Error::Input("invalid joint storage point layout".into()));
        }
        let mut volumes = BTreeMap::<usize, f64>::new();
        for (p, state) in self.layout.iter().zip(&joint.points) {
            let wall = indices
                .get(&(p.joint_id.as_str(), p.face_index))
                .ok_or_else(|| Error::Input("missing internal storage face".into()))?;
            if p.nodes_a
                .iter()
                .any(|&node| wall.total_head_m < self.mesh.nodes[node][2])
            {
                return Err(Error::Solve(
                    "joint aperture storage requires flooded faces".into(),
                ));
            }
            let aperture = b0 + Vector3::from(p.normal).dot(&Vector3::from(state.jump_m));
            if !aperture.is_finite() || aperture <= 0. {
                return Err(Error::Solve("nonpositive hydraulic aperture".into()));
            }
            *volumes.entry(wall.face).or_default() += p.weight_m2 * aperture;
        }
        Ok(volumes)
    }

    pub(crate) fn richards_field(
        &self,
        active: &[bool],
        state: &richards::State,
        faces: &[richards::Face],
        datum: f64,
        time: f64,
    ) -> Result<(Field, Vec<PointTraces>)> {
        if time != state.time_s
            || !time.is_finite()
            || faces.len() != state.relative_face_head_m.len()
            || state.relative_cell_head_m.len() != self.soil.len()
        {
            return Err(Error::Input(
                "Richards JOINT endpoint/trace layout mismatch".into(),
            ));
        }
        let cells: Vec<_> = state
            .relative_cell_head_m
            .iter()
            .enumerate()
            .filter_map(|(i, h)| h.map(|_| i))
            .collect();
        let mut indices = BTreeMap::new();
        for (fi, f) in faces.iter().enumerate() {
            for &[ci, j] in &f.owners {
                let owner = *cells
                    .get(ci)
                    .ok_or_else(|| Error::Input("invalid RT0 owner".into()))?;
                if indices.insert((owner, j), fi).is_some() {
                    return Err(Error::Input("duplicate RT0 wall".into()));
                }
            }
        }
        let trace = |w: &Wall, z: f64| -> Result<Option<Trace>> {
            if !self.is_soil(w.element) {
                return Ok(None);
            }
            let face = *indices
                .get(&(w.element, w.local))
                .ok_or_else(|| Error::Input("missing soil RT0 trace".into()))?;
            let relative = state.relative_face_head_m[face];
            let psi = (datum - z) + relative;
            let outward = state
                .cell_outflow_m3_s
                .get(w.element)
                .copied()
                .flatten()
                .ok_or_else(|| Error::Input("missing RT0 owner outflow".into()))?[w.local];
            Ok(Some(Trace {
                face,
                owner: w.element,
                local_face: w.local,
                total_head_m: datum + relative,
                pressure_head_m: psi,
                signed_pressure_kpa: self.gamma * psi,
                outward_flow_m3_s: outward,
            }))
        };
        let mut field = self.empty_field()?;
        let mut points = Vec::new();
        let pairs = self.richards_pairs(active)?;
        for ((fi, f), pair) in self
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| active[f.owner_a] && active[f.owner_b])
            .zip(pairs)
        {
            let head = match &pair.policy {
                Policy::Continuous => {
                    let a = trace(&pair.a, 0.)?
                        .ok_or_else(|| Error::Input("missing continuous trace".into()))?;
                    let b = trace(&pair.b, 0.)?
                        .ok_or_else(|| Error::Input("missing continuous trace".into()))?;
                    if a.face != b.face {
                        return Err(Error::Input(
                            "CONTINUOUS requires one shared RT0 trace".into(),
                        ));
                    }
                    a.total_head_m
                }
                Policy::Soil { soil_side } => {
                    trace(
                        if matches!(soil_side, Side::A) {
                            &pair.a
                        } else {
                            &pair.b
                        },
                        0.,
                    )?
                    .ok_or_else(|| Error::Input("missing declared soil trace".into()))?
                    .total_head_m
                }
                Policy::Prescribed { level_m } => {
                    if let Some(h) = &pair.history {
                        h.at(time)?
                    } else {
                        level_m.ok_or_else(|| {
                            Error::Input("explicit Richards reservoir head required".into())
                        })?
                    }
                }
            };
            let xyz: Vec<_> = f.nodes_a.iter().map(|&n| self.mesh.nodes[n]).collect();
            field.face_nodal_force_a_kn[fi] =
                hydrostatic::face_load(&xyz, head, self.gamma, Default::default())?.nodal_force_kn;
            for (i, p) in self
                .layout
                .iter()
                .enumerate()
                .filter(|(_, p)| p.nodes_a == f.nodes_a && p.owner_a == f.owner_a)
            {
                let z = p.position_m[2];
                let pressure = self.gamma * (head - z).max(0.);
                field.physical_pressure_kpa[i] = pressure;
                points.push(PointTraces {
                    joint_id: p.joint_id.clone(),
                    face_index: p.face_index,
                    time_s: time,
                    elevation_m: z,
                    side_a: trace(&pair.a, z)?,
                    side_b: trace(&pair.b, z)?,
                    applied_pressure_kpa: pressure,
                });
            }
        }
        Ok((
            self.field(
                field.physical_pressure_kpa,
                field.sources,
                field.face_nodal_force_a_kn,
            )?,
            points,
        ))
    }
}
