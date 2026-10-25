use crate::{Error, Result};
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Region {
    pub id: u32,
    pub body: String,
    pub material_id: u32,
    pub bounds_m: [[f64; 3]; 2],
    pub volume_m3: f64,
    pub source_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mesh {
    pub schema: String,
    pub units: String,
    pub axis: String,
    pub crs: String,
    pub origin_m: [f64; 3],
    pub nodes: Vec<[f64; 3]>,
    pub tetrahedra: Vec<[usize; 4]>,
    /// For mesh/3, one VTK-ordered quadratic cell per corner tetrahedron.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tetrahedra10: Vec<[usize; 10]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub quadrature_order: Option<usize>,
    pub material_ids: Vec<u32>,
    pub node_sets: BTreeMap<String, Vec<usize>>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub face_sets: BTreeMap<String, Vec<[usize; 3]>>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub regions: Vec<Region>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub region_ids: Vec<u32>,
    /// Boolean-fragment membership of CAD excavation tools; not material bodies.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub excavation_sets: BTreeMap<String, Vec<usize>>,
}

pub fn geometry(points: &[[f64; 3]; 4]) -> Result<(f64, [Vector3<f64>; 4])> {
    let p = points.map(|p| Vector3::from_row_slice(&p));
    let j = Matrix3::from_columns(&[p[1] - p[0], p[2] - p[0], p[3] - p[0]]);
    let det = j.determinant();
    let edge = (p[1] - p[0])
        .norm()
        .max((p[2] - p[0]).norm())
        .max((p[3] - p[0]).norm());
    if !det.is_finite() || det <= 1e-12 * edge.powi(3) {
        return Err(Error::Input("inverted or degenerate tetrahedron".into()));
    }
    let inverse = j
        .try_inverse()
        .ok_or_else(|| Error::Input("singular tetrahedron".into()))?;
    let g1 = inverse.row(0).transpose();
    let g2 = inverse.row(1).transpose();
    let g3 = inverse.row(2).transpose();
    Ok((det / 6.0, [-g1 - g2 - g3, g1, g2, g3]))
}

impl Mesh {
    /// Exact elevation of straight-sided cells, for verification/imported
    /// linear meshes. CAD-curved edge nodes must instead come from the mesher.
    pub fn elevate_straight(&self) -> Result<Self> {
        self.validate()?;
        if !self.tetrahedra10.is_empty() || self.schema != "soil-fem-mesh/2" {
            return Err(Error::Input(
                "elevation requires linear body-resolved mesh/2".into(),
            ));
        }
        let mut result = self.clone();
        result.schema = "soil-fem-mesh/3".into();
        let mut edges = BTreeMap::new();
        for tet in &self.tetrahedra {
            let mut cell = [0; 10];
            cell[..4].copy_from_slice(tet);
            for (k, [i, j]) in crate::elements::EDGES.into_iter().enumerate() {
                let mut edge = [tet[i], tet[j]];
                edge.sort_unstable();
                let id = if let Some(&id) = edges.get(&edge) {
                    id
                } else {
                    let id = result.nodes.len();
                    result.nodes.push(std::array::from_fn(|a| {
                        (self.nodes[edge[0]][a] + self.nodes[edge[1]][a]) / 2.
                    }));
                    edges.insert(edge, id);
                    id
                };
                cell[k + 4] = id;
            }
            result.tetrahedra10.push(cell);
        }
        for (name, set) in &self.node_sets {
            let original: BTreeSet<_> = set.iter().copied().collect();
            let elevated = result
                .node_sets
                .get_mut(name)
                .ok_or_else(|| Error::Input("lost node set during elevation".into()))?;
            elevated.extend(
                edges
                    .iter()
                    .filter(|(edge, _)| original.contains(&edge[0]) && original.contains(&edge[1]))
                    .map(|(_, id)| *id),
            );
            elevated.sort_unstable();
            elevated.dedup();
        }
        result.validate()?;
        Ok(result)
    }
    pub fn cell_nodes(&self, index: usize) -> Vec<usize> {
        if self.tetrahedra10.is_empty() {
            self.tetrahedra[index].to_vec()
        } else {
            self.tetrahedra10[index].to_vec()
        }
    }
    pub fn integration(&self, index: usize) -> Result<Vec<crate::elements::VolumePoint>> {
        let points: Vec<_> = self
            .cell_nodes(index)
            .iter()
            .map(|&i| self.nodes[i])
            .collect();
        crate::elements::volume_points(&points, self.quadrature_order.unwrap_or(3))
    }
    pub fn integration_layout(&self) -> Result<Vec<crate::elements::PointInfo>> {
        let mut result = Vec::new();
        for i in 0..self.tetrahedra.len() {
            for (k, q) in self.integration(i)?.into_iter().enumerate() {
                result.push(crate::elements::PointInfo {
                    element_index: i,
                    local_index: k,
                    natural_coordinates: q.natural,
                    position_m: q.position,
                    weight_m3: q.weight_m3,
                    shape_functions: q.shape,
                });
            }
        }
        Ok(result)
    }
    /// Outward local reference ordering; face keys themselves may be sorted.
    pub fn face_nodes(&self, owner: usize, key: &[usize; 3]) -> Result<Vec<usize>> {
        let tet = self.tetrahedra[owner];
        let mut sorted_key = *key;
        sorted_key.sort_unstable();
        for local in crate::elements::FACES {
            let vertices = local.map(|a| tet[a]);
            let mut sorted = vertices;
            sorted.sort_unstable();
            if sorted != sorted_key {
                continue;
            }
            let mut nodes = vertices.to_vec();
            if !self.tetrahedra10.is_empty() {
                for [a, b] in [[0, 1], [1, 2], [2, 0]] {
                    let va = local[a];
                    let vb = local[b];
                    let edge = crate::elements::EDGES
                        .iter()
                        .position(|&[i, j]| (i == va && j == vb) || (j == va && i == vb))
                        .ok_or_else(|| Error::Input("unknown Tet10 edge".into()))?;
                    nodes.push(self.tetrahedra10[owner][edge + 4]);
                }
            }
            return Ok(nodes);
        }
        Err(Error::Input(
            "face is not part of the owning tetrahedron".into(),
        ))
    }
    pub fn points(&self, tet: &[usize; 4]) -> [[f64; 3]; 4] {
        tet.map(|i| self.nodes[i])
    }

    pub fn validate(&self) -> Result<f64> {
        if !["soil-fem-mesh/1", "soil-fem-mesh/2", "soil-fem-mesh/3"]
            .contains(&self.schema.as_str())
            || self.units != "m-kN-kPa"
            || self.axis != "z_up"
        {
            return Err(Error::Input("schema/units/axis mismatch".into()));
        }
        if self.crs.is_empty() || self.nodes.len() < 4 || self.tetrahedra.is_empty() {
            return Err(Error::Input(
                "CRS, nodes and volume elements required".into(),
            ));
        }
        if !self
            .origin_m
            .iter()
            .chain(self.nodes.iter().flatten())
            .all(|x| x.is_finite())
        {
            return Err(Error::Input("non-finite coordinates".into()));
        }
        if self.material_ids.len() != self.tetrahedra.len() || self.material_ids.contains(&0) {
            return Err(Error::Input(
                "one positive material ID per tetrahedron required".into(),
            ));
        }
        if self.schema == "soil-fem-mesh/3" {
            if self.tetrahedra10.len() != self.tetrahedra.len()
                || ![3, 4, 5].contains(&self.quadrature_order.unwrap_or(3))
            {
                return Err(Error::Input(
                    "mesh/3 requires matching Tet10 connectivity and valid quadrature".into(),
                ));
            }
        } else if !self.tetrahedra10.is_empty() || self.quadrature_order.is_some() {
            return Err(Error::Input(
                "quadratic cells/quadrature require mesh/3".into(),
            ));
        }
        if self.schema != "soil-fem-mesh/1" {
            if self.regions.is_empty() || self.region_ids.len() != self.tetrahedra.len() {
                return Err(Error::Input("region metadata required for mesh/2".into()));
            }
            let mut ids = BTreeMap::new();
            let mut names = BTreeSet::new();
            for region in &self.regions {
                if region.id == 0
                    || region.material_id == 0
                    || region.body.is_empty()
                    || ids.insert(region.id, region.material_id).is_some()
                    || !names.insert(&region.body)
                    || !region.volume_m3.is_finite()
                    || region.volume_m3 <= 0.
                    || !region.bounds_m.iter().flatten().all(|x| x.is_finite())
                    || (0..3).any(|i| region.bounds_m[0][i] >= region.bounds_m[1][i])
                    || region.source_sha256.len() != 64
                    || !region.source_sha256.chars().all(|c| c.is_ascii_hexdigit())
                {
                    return Err(Error::Input("invalid/duplicate CAD region metadata".into()));
                }
            }
            let used_regions: BTreeSet<_> = self.region_ids.iter().copied().collect();
            if used_regions.len() != ids.len()
                || self
                    .region_ids
                    .iter()
                    .zip(&self.material_ids)
                    .any(|(r, m)| ids.get(r) != Some(m))
            {
                return Err(Error::Input(
                    "element region/material mismatch or empty region".into(),
                ));
            }
        } else if !self.regions.is_empty() || !self.region_ids.is_empty() {
            return Err(Error::Input("region metadata requires mesh/2".into()));
        }
        let mut used = BTreeSet::<usize>::new();
        for (name, indices) in &self.excavation_sets {
            let unique: BTreeSet<_> = indices.iter().collect();
            if name.is_empty()
                || self.regions.iter().any(|r| r.body == *name)
                || indices.is_empty()
                || unique.len() != indices.len()
                || indices.iter().any(|&i| i >= self.tetrahedra.len())
            {
                return Err(Error::Input("invalid CAD excavation element set".into()));
            }
        }
        let mut unique = BTreeSet::new();
        let mut faces = BTreeMap::<[usize; 3], usize>::new();
        let mut volume = 0.0;
        let mut edges = BTreeMap::<[usize; 2], usize>::new();
        let mut midpoint_edges = BTreeMap::<usize, [usize; 2]>::new();
        let corners: BTreeSet<_> = self.tetrahedra.iter().flatten().copied().collect();
        for (ei, tet) in self.tetrahedra.iter().enumerate() {
            if tet.iter().any(|&i| i >= self.nodes.len()) {
                return Err(Error::Input("element references unknown node".into()));
            }
            let mut key = *tet;
            key.sort_unstable();
            if key.windows(2).any(|w| w[0] == w[1]) || !unique.insert(key) {
                return Err(Error::Input("repeated node or duplicate element".into()));
            }
            for mut face in [
                [tet[0], tet[1], tet[2]],
                [tet[0], tet[1], tet[3]],
                [tet[0], tet[2], tet[3]],
                [tet[1], tet[2], tet[3]],
            ] {
                face.sort_unstable();
                let count = faces.entry(face).or_default();
                *count += 1;
                if *count > 2 {
                    return Err(Error::Input("non-manifold volume mesh".into()));
                }
            }
            used.extend(tet);
            if self.tetrahedra10.is_empty() {
                volume += geometry(&self.points(tet))?.0;
            } else {
                let cell = self.tetrahedra10[ei];
                if cell[..4] != tet[..]
                    || cell.iter().any(|&i| i >= self.nodes.len())
                    || cell.iter().copied().collect::<BTreeSet<_>>().len() != 10
                {
                    return Err(Error::Input(
                        "invalid Tet10 node/corner connectivity".into(),
                    ));
                }
                for (k, [i, j]) in crate::elements::EDGES.into_iter().enumerate() {
                    let mut edge = [tet[i], tet[j]];
                    edge.sort_unstable();
                    let middle = cell[k + 4];
                    if corners.contains(&middle)
                        || edges.insert(edge, middle).is_some_and(|old| old != middle)
                        || midpoint_edges
                            .insert(middle, edge)
                            .is_some_and(|old| old != edge)
                    {
                        return Err(Error::Input(
                            "nonconforming/aliased Tet10 edge nodes".into(),
                        ));
                    }
                }
                used.extend(cell);
                volume += self
                    .integration(ei)?
                    .iter()
                    .map(|q| q.weight_m3)
                    .sum::<f64>();
            }
        }
        if used.len() != self.nodes.len() {
            return Err(Error::Input("orphan nodes".into()));
        }
        if self
            .node_sets
            .values()
            .flatten()
            .any(|&i| i >= self.nodes.len())
        {
            return Err(Error::Input("node set references unknown node".into()));
        }
        for (name, nodes) in &self.node_sets {
            if name.is_empty() || nodes.iter().collect::<BTreeSet<_>>().len() != nodes.len() {
                return Err(Error::Input("invalid/duplicate node set membership".into()));
            }
        }
        for (name, set) in &self.face_sets {
            let unique: BTreeSet<_> = set
                .iter()
                .map(|f| {
                    let mut f = *f;
                    f.sort_unstable();
                    f
                })
                .collect();
            if name.is_empty()
                || unique.len() != set.len()
                || unique.iter().any(|f| !faces.contains_key(f))
            {
                return Err(Error::Input("invalid/duplicate face set membership".into()));
            }
        }
        if !volume.is_finite() {
            return Err(Error::Input("non-finite mesh volume".into()));
        }
        Ok(volume)
    }
}
