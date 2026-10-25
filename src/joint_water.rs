//! Physical interface pressure is a pair of external forces, separate from
//! effective contact springs. Displacement nodes never merge here.
use crate::{Error, Result, darcy, ground::Ground, hydrostatic, joints, mesh::Mesh};
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Side {
    A,
    B,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Policy {
    Soil {
        soil_side: Side,
    },
    Continuous,
    /// None follows the current hydrostatic/exterior-water level.
    Prescribed {
        level_m: Option<f64>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointReport {
    pub physical_pore_pressure_kpa: f64,
    pub effective_contact_pressure_kpa: f64,
    pub total_normal_pressure_kpa: f64,
    pub friction_limit_kpa: f64,
    pub source: String,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Field {
    pub physical_pressure_kpa: Vec<f64>,
    pub sources: Vec<String>,
    pub nodal_load_kn: Vec<f64>,
    /// Per matching face, forces on side a in its node order. The opposing
    /// forces are reconstructed, so a partial DEATH cannot retain stale loads.
    pub face_nodal_force_a_kn: Vec<Vec<[f64; 3]>>,
    pub resultant_force_kn: [f64; 3],
    pub resultant_moment_kn_m: [f64; 3],
}
struct Face {
    definition: usize,
    nodes_a: Vec<usize>,
    nodes_b: Vec<usize>,
    owner_a: usize,
    owner_b: usize,
}
pub struct Model {
    mesh: Mesh,
    configuration: joints::Configuration,
    layout: Vec<joints::PointInfo>,
    faces: Vec<Face>,
    soil: Vec<bool>,
    gamma: f64,
}
impl Model {
    pub fn new(mesh: &Mesh, contact: &joints::Model, ground: &Ground) -> Result<Self> {
        let configuration = contact.configuration.clone();
        if !ground.initialization.water_unit_weight_kn_m3.is_finite()
            || ground.initialization.water_unit_weight_kn_m3 <= 0.
        {
            return Err(Error::Input(
                "wet JOINT requires finite positive water unit weight".into(),
            ));
        }
        if configuration.joints.iter().any(|d| d.water.is_none()) {
            return Err(Error::Input(
                "wet JOINT requires an explicit pressure policy on every interface".into(),
            ));
        }
        let structures: BTreeSet<_> = configuration.structure_bodies.iter().collect();
        if structures.len() != configuration.structure_bodies.len()
            || structures
                .iter()
                .any(|b| !mesh.regions.iter().any(|r| &r.body == *b))
        {
            return Err(Error::Input(
                "unknown or duplicate impermeable structure body".into(),
            ));
        }
        let mut soil = Vec::with_capacity(mesh.tetrahedra.len());
        for &region in &mesh.region_ids {
            let r = mesh
                .regions
                .iter()
                .find(|r| r.id == region)
                .ok_or_else(|| Error::Input("missing JOINT body".into()))?;
            let is_soil = !structures.contains(&r.body);
            if !is_soil {
                let m = ground
                    .materials
                    .iter()
                    .find(|m| m.id == r.material_id)
                    .ok_or_else(|| Error::Input("missing structure material".into()))?;
                if !matches!(m.model, crate::soil::Model::MohrCoulomb(_)) {
                    return Err(Error::Input("impermeable total-stress continuum structure requires an MC material, not a soil compression-hardening model".into()));
                }
            }
            soil.push(is_soil);
        }
        let layout = contact.layout();
        let mut faces = Vec::new();
        let mut visited = BTreeSet::new();
        for p in &layout {
            let di = configuration
                .joints
                .iter()
                .position(|d| d.id == p.joint_id)
                .ok_or_else(|| Error::Input("missing joint definition".into()))?;
            let d = &configuration.joints[di];
            match d
                .water
                .as_ref()
                .ok_or_else(|| Error::Input("missing pressure policy".into()))?
            {
                Policy::Continuous if !(soil[p.owner_a] && soil[p.owner_b]) => {
                    return Err(Error::Input(
                        "CONTINUOUS requires soil on both sides".into(),
                    ));
                }
                Policy::Soil { soil_side } => {
                    let expected = matches!(soil_side, Side::A);
                    if soil[p.owner_a] != expected || soil[p.owner_b] == expected {
                        return Err(Error::Input("SOIL requires the declared soil side and one impermeable structure side".into()));
                    }
                }
                Policy::Prescribed { level_m: Some(h) } if !h.is_finite() => {
                    return Err(Error::Input(
                        "non-finite prescribed JOINT water level".into(),
                    ));
                }
                _ => {}
            }
            if visited.insert((di, p.face_index)) {
                faces.push(Face {
                    definition: di,
                    nodes_a: p.nodes_a.clone(),
                    nodes_b: p.nodes_b.clone(),
                    owner_a: p.owner_a,
                    owner_b: p.owner_b,
                });
            }
        }
        Ok(Self {
            mesh: mesh.clone(),
            configuration,
            layout,
            faces,
            soil,
            gamma: ground.initialization.water_unit_weight_kn_m3,
        })
    }
    pub fn is_soil(&self, element: usize) -> bool {
        self.soil[element]
    }
    pub fn soil_domain(&self, active: &[bool]) -> Result<Vec<bool>> {
        if active.len() != self.soil.len() {
            return Err(Error::Input("wet JOINT domain length mismatch".into()));
        }
        Ok(active
            .iter()
            .zip(&self.soil)
            .map(|(a, s)| *a && *s)
            .collect())
    }
    fn level(&self, definition: usize, current: Option<f64>) -> Option<f64> {
        match self.configuration.joints[definition].water.as_ref() {
            Some(Policy::Prescribed { level_m: Some(h) }) => Some(*h),
            _ => current,
        }
    }
    fn source(&self, definition: usize) -> String {
        match self.configuration.joints[definition].water.as_ref() {
            Some(Policy::Soil { .. }) => "soil",
            Some(Policy::Continuous) => "continuous",
            _ => "prescribed",
        }
        .into()
    }
    /// Hydrostatics use adaptive partial-submergence integration on the actual
    /// curved face, rather than integrating clipped GP samples at fixed order.
    pub fn hydrostatic(&self, active: &[bool], current: Option<f64>) -> Result<Field> {
        self.soil_domain(active)?;
        if current.is_some_and(|h| !h.is_finite()) {
            return Err(Error::Input("non-finite JOINT water level".into()));
        }
        let mut forces: Vec<_> = self
            .faces
            .iter()
            .map(|f| vec![[0.; 3]; f.nodes_a.len()])
            .collect();
        for (fi, f) in self.faces.iter().enumerate() {
            if !(active[f.owner_a] && active[f.owner_b]) {
                continue;
            }
            let h = self.level(f.definition, current);
            let points: Vec<_> = f.nodes_a.iter().map(|&n| self.mesh.nodes[n]).collect();
            if matches!(
                self.configuration.joints[f.definition].water,
                Some(Policy::Prescribed { .. })
            ) && (self.soil[f.owner_a] || self.soil[f.owner_b])
            {
                let same = match (h, current) {
                    (Some(a), Some(b)) => {
                        (a - b).abs()
                            <= 1e-11 * (1. + (a - points[0][2]).abs().max((b - points[0][2]).abs()))
                                + 16. * f64::EPSILON * a.abs().max(b.abs())
                    }
                    (None, None) => true,
                    _ => false,
                };
                let wet_at = |level: Option<f64>| -> Result<bool> {
                    Ok(match level {
                        Some(v) => hydrostatic::depth_bounds(&points, v)?[1] > 0.,
                        None => false,
                    })
                };
                if !same && (wet_at(h)? || wet_at(current)?) {
                    return Err(Error::Input("prescribed JOINT pressure conflicts with the adjacent hydrostatic soil; specify a consistent level or solve a hydraulic head field".into()));
                }
            }
            if let Some(h) = h {
                let integrated = hydrostatic::face_load(
                    &points,
                    h,
                    self.gamma,
                    hydrostatic::Integration::default(),
                )?;
                forces[fi] = integrated.nodal_force_kn;
            }
        }
        let mut pressures = Vec::with_capacity(self.layout.len());
        let mut sources = Vec::with_capacity(self.layout.len());
        for p in &self.layout {
            let di = self
                .configuration
                .joints
                .iter()
                .position(|d| d.id == p.joint_id)
                .ok_or_else(|| Error::Input("missing interface definition".into()))?;
            pressures.push(if active[p.owner_a] && active[p.owner_b] {
                self.level(di, current)
                    .map_or(0., |h| self.gamma * (h - p.position_m[2]).max(0.))
            } else {
                0.
            });
            sources.push(self.source(di));
        }
        self.field(pressures, sources, forces)
    }
    fn field(
        &self,
        physical_pressure_kpa: Vec<f64>,
        sources: Vec<String>,
        face_nodal_force_a_kn: Vec<Vec<[f64; 3]>>,
    ) -> Result<Field> {
        if physical_pressure_kpa.len() != self.layout.len()
            || sources.len() != self.layout.len()
            || face_nodal_force_a_kn.len() != self.faces.len()
            || face_nodal_force_a_kn
                .iter()
                .zip(&self.faces)
                .any(|(q, f)| q.len() != f.nodes_a.len())
        {
            return Err(Error::Input("JOINT water force layout mismatch".into()));
        }
        for (p, source) in self.layout.iter().zip(&sources) {
            let di = self
                .configuration
                .joints
                .iter()
                .position(|d| d.id == p.joint_id)
                .ok_or_else(|| Error::Input("missing JOINT pressure source".into()))?;
            if *source != self.source(di) {
                return Err(Error::Input(
                    "stored JOINT pressure source differs from its definition".into(),
                ));
            }
        }
        let mut nodal_load_kn = vec![0.; 3 * self.mesh.nodes.len()];
        for (f, forces) in self.faces.iter().zip(&face_nodal_force_a_kn) {
            for (j, force) in forces.iter().enumerate() {
                for a in 0..3 {
                    nodal_load_kn[3 * f.nodes_a[j] + a] += force[a];
                    nodal_load_kn[3 * f.nodes_b[j] + a] -= force[a];
                }
            }
        }
        if !physical_pressure_kpa
            .iter()
            .chain(&nodal_load_kn)
            .all(|v| v.is_finite())
            || physical_pressure_kpa.iter().any(|p| *p < 0.)
        {
            return Err(Error::Solve("JOINT water force overflow".into()));
        }
        let mut force = Vector3::zeros();
        let mut moment = Vector3::zeros();
        for (i, x) in self.mesh.nodes.iter().enumerate() {
            let f = Vector3::new(
                nodal_load_kn[3 * i],
                nodal_load_kn[3 * i + 1],
                nodal_load_kn[3 * i + 2],
            );
            force += f;
            moment += Vector3::from(*x).cross(&f);
        }
        Ok(Field {
            physical_pressure_kpa,
            sources,
            nodal_load_kn,
            face_nodal_force_a_kn,
            resultant_force_kn: force.into(),
            resultant_moment_kn_m: moment.into(),
        })
    }
    pub fn annotate(&self, report: &mut joints::Report, field: &Field) -> Result<()> {
        if report.points.len() != self.layout.len()
            || field.physical_pressure_kpa.len() != self.layout.len()
            || field.sources.len() != self.layout.len()
        {
            return Err(Error::Input(
                "JOINT hydraulic report layout mismatch".into(),
            ));
        }
        for (i, p) in report.points.iter_mut().enumerate() {
            let di = self
                .configuration
                .joints
                .iter()
                .position(|d| d.id == self.layout[i].joint_id)
                .ok_or_else(|| Error::Input("missing joint definition".into()))?;
            let effective = if p.active { p.normal_pressure_kpa } else { 0. };
            let water = field.physical_pressure_kpa[i];
            p.hydraulic = Some(PointReport {
                physical_pore_pressure_kpa: water,
                effective_contact_pressure_kpa: effective,
                total_normal_pressure_kpa: effective + water,
                friction_limit_kpa: self.configuration.joints[di].material.friction_coefficient
                    * effective.max(0.),
                source: field.sources[i].clone(),
            });
        }
        report.water = Some(field.clone());
        Ok(())
    }

    pub fn retain(&self, field: &Field, active: &[bool]) -> Result<Field> {
        self.soil_domain(active)?;
        let checked = self.field(
            field.physical_pressure_kpa.clone(),
            field.sources.clone(),
            field.face_nodal_force_a_kn.clone(),
        )?;
        if &checked != field || field.physical_pressure_kpa.iter().any(|v| *v < 0.) {
            return Err(Error::Input(
                "inconsistent stored JOINT pressure/force field".into(),
            ));
        }
        let mut pressures = field.physical_pressure_kpa.clone();
        let mut forces = field.face_nodal_force_a_kn.clone();
        for (i, p) in self.layout.iter().enumerate() {
            if !(active[p.owner_a] && active[p.owner_b]) {
                pressures[i] = 0.;
            }
        }
        for (f, q) in self.faces.iter().zip(&mut forces) {
            if !(active[f.owner_a] && active[f.owner_b]) {
                q.fill([0.; 3]);
            }
        }
        self.field(pressures, field.sources.clone(), forces)
    }
    pub fn interpolate(
        &self,
        from: &Field,
        to: &Field,
        active: &[bool],
        progress: f64,
    ) -> Result<Field> {
        if !progress.is_finite() || !(0. ..=1.).contains(&progress) {
            return Err(Error::Input(
                "invalid JOINT water interpolation progress".into(),
            ));
        }
        let a = self.retain(from, active)?;
        let b = self.retain(to, active)?;
        if progress == 0. {
            return Ok(a);
        }
        if progress == 1. {
            return Ok(b);
        }
        if a.sources != b.sources {
            return Err(Error::Input(
                "JOINT pressure source cannot change implicitly".into(),
            ));
        }
        let pressures = a
            .physical_pressure_kpa
            .iter()
            .zip(&b.physical_pressure_kpa)
            .map(|(a, b)| a + (b - a) * progress)
            .collect();
        let forces = a
            .face_nodal_force_a_kn
            .iter()
            .zip(&b.face_nodal_force_a_kn)
            .map(|(a, b)| {
                a.iter()
                    .zip(b)
                    .map(|(a, b)| std::array::from_fn(|i| a[i] + (b[i] - a[i]) * progress))
                    .collect()
            })
            .collect();
        self.field(pressures, b.sources, forces)
    }

    /// Tie only matching active soil-interface head nodes. A separate compact
    /// mesh makes the mapping and the actual hydraulic geometry auditable.
    pub fn solve_steady(
        &self,
        case: &darcy::Case,
        active: &[bool],
        material_ids: Option<&[u32]>,
        current_level: Option<f64>,
    ) -> Result<Flow> {
        let soil_active = self.soil_domain(active)?;
        let mut parent: Vec<_> = (0..self.mesh.nodes.len()).collect();
        fn root(p: &mut [usize], mut n: usize) -> usize {
            while p[n] != n {
                p[n] = p[p[n]];
                n = p[n];
            }
            n
        }
        let mut flow_case = case.clone();
        let mut prescribed_faces = BTreeMap::<String, (f64, Vec<[usize; 3]>)>::new();
        for f in &self.faces {
            if !(active[f.owner_a] && active[f.owner_b]) {
                continue;
            }
            let definition = &self.configuration.joints[f.definition];
            match definition
                .water
                .as_ref()
                .ok_or_else(|| Error::Input("missing hydraulic policy".into()))?
            {
                Policy::Continuous => {
                    for (&a, &b) in f.nodes_a.iter().zip(&f.nodes_b) {
                        let ra = root(&mut parent, a);
                        let rb = root(&mut parent, b);
                        parent[ra.max(rb)] = ra.min(rb);
                    }
                }
                Policy::Prescribed { .. } => {
                    let h=self.level(f.definition,current_level).ok_or_else(||Error::Input("PRESCRIBED steady JOINT needs an explicit current or fixed water level".into()))?;
                    for (owner, side, nodes) in
                        [(f.owner_a, "a", &f.nodes_a), (f.owner_b, "b", &f.nodes_b)]
                    {
                        if self.soil[owner] {
                            let name = format!("__joint_water_{}_{}", f.definition, side);
                            if self.mesh.face_sets.contains_key(&name) {
                                return Err(Error::Input(
                                    "reserved hydraulic JOINT face-set name collision".into(),
                                ));
                            }
                            prescribed_faces
                                .entry(name)
                                .or_insert((h, vec![]))
                                .1
                                .push([nodes[0], nodes[1], nodes[2]]);
                        }
                    }
                }
                Policy::Soil { .. } => {}
            }
        }
        let mut mesh = self.mesh.clone();
        let mut compact = BTreeMap::new();
        let mut coordinates = Vec::new();
        let mut mapping = Vec::new();
        for n in 0..parent.len() {
            let r = root(&mut parent, n);
            let next = compact.len();
            let id = *compact.entry(r).or_insert_with(|| {
                coordinates.push(self.mesh.nodes[r]);
                next
            });
            mapping.push(id);
        }
        mesh.nodes = coordinates;
        for cell in &mut mesh.tetrahedra {
            for n in cell {
                *n = mapping[*n];
            }
        }
        for cell in &mut mesh.tetrahedra10 {
            for n in cell {
                *n = mapping[*n];
            }
        }
        for set in mesh.node_sets.values_mut() {
            *set = set
                .iter()
                .map(|&n| mapping[n])
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
        }
        for set in mesh.face_sets.values_mut() {
            let mut unique = BTreeSet::new();
            *set = set
                .iter()
                .filter_map(|f| {
                    let face = f.map(|n| mapping[n]);
                    let mut key = face;
                    key.sort_unstable();
                    unique.insert(key).then_some(face)
                })
                .collect();
        }
        for (name, (h, faces)) in prescribed_faces {
            mesh.face_sets.insert(
                name.clone(),
                faces.into_iter().map(|f| f.map(|n| mapping[n])).collect(),
            );
            flow_case.heads.push(darcy::Head {
                face_set: name,
                head_m: h,
                gradient: [0.; 3],
                reference_point_m: [0.; 3],
            });
        }
        mesh.validate()?;
        let solution = darcy::solve(&mesh, &flow_case, &soil_active, material_ids)?;
        let field = &solution.field;
        let original_field = darcy::HeadField {
            relative_head_m: mapping.iter().map(|&n| field.relative_head_m[n]).collect(),
            datum_head_m: mapping.iter().map(|&n| field.datum_head_m[n]).collect(),
            head_m: mapping.iter().map(|&n| field.head_m[n]).collect(),
            defined_nodes: mapping.iter().map(|&n| field.defined_nodes[n]).collect(),
        };
        original_field.validate(&self.mesh, &soil_active)?;
        Ok(Flow {
            original_mesh_sha256: crate::checkpoint::hash_json(&self.mesh)?,
            original_case_sha256: crate::checkpoint::hash_json(case)?,
            hydraulic_mesh: mesh,
            hydraulic_case: flow_case,
            original_node_to_hydraulic_node: mapping,
            original_field,
            solution,
        })
    }

    pub fn from_heads(
        &self,
        active: &[bool],
        field: &darcy::HeadField,
        current_level: Option<f64>,
    ) -> Result<Field> {
        field.validate(&self.mesh, &self.soil_domain(active)?)?;
        let pressure = |definition: usize,
                        owner_a: usize,
                        owner_b: usize,
                        nodes_a: &[usize],
                        nodes_b: &[usize],
                        shape: &[f64],
                        z: f64|
         -> Result<f64> {
            let pa = if self.soil[owner_a] {
                Some(self.gamma * field.pressure_head(nodes_a, shape, z)?)
            } else {
                None
            };
            let pb = if self.soil[owner_b] {
                Some(self.gamma * field.pressure_head(nodes_b, shape, z)?)
            } else {
                None
            };
            let tolerance = self.gamma * (2. * self.configuration.geometric_tolerance_m + 1e-9);
            if pa.into_iter().chain(pb).any(|p| p < -tolerance) {
                return Err(Error::Input("wet saturated JOINT cannot replace a free-surface solve by clipping negative head".into()));
            }
            let p = match self.configuration.joints[definition]
                .water
                .as_ref()
                .ok_or_else(|| Error::Input("missing pressure policy".into()))?
            {
                Policy::Soil { soil_side: Side::A } => pa,
                Policy::Soil { soil_side: Side::B } => pb,
                Policy::Continuous => {
                    if (pa.ok_or_else(|| Error::Input("missing side-a soil head".into()))?
                        - pb.ok_or_else(|| Error::Input("missing side-b soil head".into()))?)
                    .abs()
                        > tolerance
                    {
                        return Err(Error::Input(
                            "CONTINUOUS JOINT has inconsistent paired heads".into(),
                        ));
                    }
                    pa
                }
                Policy::Prescribed { .. } => {
                    let p = self
                        .level(definition, current_level)
                        .map_or(0., |h| self.gamma * (h - z).max(0.));
                    if pa
                        .into_iter()
                        .chain(pb)
                        .any(|other| (other - p).abs() > tolerance)
                    {
                        return Err(Error::Input(
                            "PRESCRIBED JOINT pressure conflicts with an adjacent soil head".into(),
                        ));
                    }
                    Some(p)
                }
            };
            Ok(
                p.ok_or_else(|| Error::Input("missing JOINT pressure source".into()))?
                    .max(0.),
            )
        };
        let mut forces: Vec<_> = self
            .faces
            .iter()
            .map(|f| vec![[0.; 3]; f.nodes_a.len()])
            .collect();
        for (fi, f) in self.faces.iter().enumerate() {
            if !(active[f.owner_a] && active[f.owner_b]) {
                continue;
            }
            let coordinates: Vec<_> = f.nodes_a.iter().map(|&n| self.mesh.nodes[n]).collect();
            if !self.soil[f.owner_a] && !self.soil[f.owner_b] {
                if let Some(h) = self.level(f.definition, current_level) {
                    let value = hydrostatic::face_load(
                        &coordinates,
                        h,
                        self.gamma,
                        hydrostatic::Integration::default(),
                    )?;
                    forces[fi] = value.nodal_force_kn;
                }
            } else {
                for q in crate::elements::face_points(&coordinates)? {
                    let p = pressure(
                        f.definition,
                        f.owner_a,
                        f.owner_b,
                        &f.nodes_a,
                        &f.nodes_b,
                        &q.shape,
                        q.position[2],
                    )?;
                    for (i, &n) in q.shape.iter().enumerate() {
                        for (a, value) in forces[fi][i].iter_mut().enumerate() {
                            let v = -n * p * q.area_vector[a];
                            *value += v;
                        }
                    }
                }
            }
        }
        let mut pressures = Vec::new();
        let mut sources = Vec::new();
        for p in &self.layout {
            let di = self
                .configuration
                .joints
                .iter()
                .position(|d| d.id == p.joint_id)
                .ok_or_else(|| Error::Input("missing pressure definition".into()))?;
            pressures.push(if active[p.owner_a] && active[p.owner_b] {
                pressure(
                    di,
                    p.owner_a,
                    p.owner_b,
                    &p.nodes_a,
                    &p.nodes_b,
                    &p.shape,
                    p.position_m[2],
                )?
            } else {
                0.
            });
            sources.push(self.source(di));
        }
        self.field(pressures, sources, forces)
    }
}

/// A hydraulic mesh and explicit map keep the mechanical split mesh unchanged.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Flow {
    pub original_mesh_sha256: String,
    pub original_case_sha256: String,
    pub hydraulic_mesh: Mesh,
    pub hydraulic_case: darcy::Case,
    pub original_node_to_hydraulic_node: Vec<usize>,
    pub original_field: darcy::HeadField,
    pub solution: darcy::Report,
}
pub fn same_flow(a: &Option<Flow>, b: &Option<Flow>) -> Result<bool> {
    Ok(match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.original_mesh_sha256 == b.original_mesh_sha256
                && a.original_case_sha256 == b.original_case_sha256
                && crate::checkpoint::hash_json(&a.hydraulic_mesh)?
                    == crate::checkpoint::hash_json(&b.hydraulic_mesh)?
                && crate::checkpoint::hash_json(&a.hydraulic_case)?
                    == crate::checkpoint::hash_json(&b.hydraulic_case)?
                && a.original_node_to_hydraulic_node == b.original_node_to_hydraulic_node
                && a.original_field == b.original_field
                && darcy::same_solution(&a.solution, &b.solution)?
        }
        _ => false,
    })
}
