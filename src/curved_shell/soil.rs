//! Eliminate conforming attachment translations by assigning the same row.
//! Remaining soil nodes follow the shell unknowns in ascending mesh-node order.
use super::path::Case;
use crate::{
    Error, Result,
    mesh::Mesh,
    offset_composite::{
        Face,
        path::{self, Cell, Soil},
    },
};
use nalgebra::{DMatrix, DVector, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Attachment {
    pub mesh: Mesh,
    pub materials: BTreeMap<u32, Soil>,
    /// Every attached shell face owns one complete exterior Tet10 face.
    pub owners: BTreeMap<String, Face>,
    pub exposed_side: BTreeMap<String, i8>,
    pub tie_at_physical_surface: bool,
}
pub(super) struct Domain {
    pub stiffness: Option<DMatrix<f64>>,
    pub cells: Vec<Cell>,
    pub size: usize,
    pub rows: Vec<[usize; 3]>,
}
impl Domain {
    pub fn new(case: &Case) -> Result<Self> {
        let mut size = 6 * case.nodes.len();
        let Some(input) = &case.soil else {
            return Ok(Self {
                stiffness: None,
                cells: vec![],
                size,
                rows: vec![],
            });
        };
        let mesh = &input.mesh;
        mesh.validate()?;
        if mesh.tetrahedra10.len() != mesh.tetrahedra.len()
            || mesh.tetrahedra.len() > 24
            || input.owners.len() != case.faces.len()
            || input.exposed_side.len() != case.faces.len()
        {
            return Err(Error::Input("curved shell attachment needs complete exterior faces and a dry Tet10 patch of at most 24 cells".into()));
        }
        let mut bindings = BTreeMap::new();
        let mut reverse = BTreeMap::new();
        let mut used = BTreeSet::new();
        for f in &case.faces {
            let owner = input
                .owners
                .get(&f.id)
                .ok_or_else(|| Error::Input("curved shell soil face owner missing".into()))?;
            owner.validate(mesh)?;
            let sigma = *input
                .exposed_side
                .get(&f.id)
                .ok_or_else(|| Error::Input("missing curved shell exposed side".into()))?;
            if ![-1, 1].contains(&sigma)
                || (input.tie_at_physical_surface
                    && f.offset_m != -f64::from(sigma) * case.section.thickness_m / 2.)
            {
                return Err(Error::Input(
                    "curved shell physical attachment side and offset disagree".into(),
                ));
            }
            let exterior = mesh.face_nodes(
                owner.cell,
                &[owner.nodes[0], owner.nodes[1], owner.nodes[2]],
            )?;
            let corner = exterior
                .iter()
                .take(3)
                .position(|&n| n == owner.nodes[0])
                .unwrap();
            if (exterior[(corner + 1) % 3] == owner.nodes[1]) != (sigma == 1) {
                return Err(Error::Input(
                    "curved shell exposed side points into soil owner".into(),
                ));
            }
            let mut key = owner.nodes;
            key.sort_unstable();
            if !used.insert(key) {
                return Err(Error::Input("duplicate curved shell soil face".into()));
            }
            for j in 0..3 {
                let ij = [owner.nodes[j], owner.nodes[(j + 1) % 3]];
                let cell = mesh.tetrahedra10[owner.cell];
                let edge = crate::elements::EDGES
                    .iter()
                    .position(|&[a, b]| ij.contains(&cell[a]) && ij.contains(&cell[b]));
                if edge.is_none_or(|e| cell[e + 4] != owner.nodes[j + 3]) {
                    return Err(Error::Input(
                        "curved shell soil edge midpoint mismatch".into(),
                    ));
                }
            }
            for (&shell, &soil) in f.nodes.iter().zip(&owner.nodes) {
                if shell >= case.nodes.len()
                    || mesh.nodes[soil] != case.nodes[shell]
                    || bindings.insert(soil, shell).is_some_and(|n| n != shell)
                    || reverse.insert(shell, soil).is_some_and(|n| n != soil)
                {
                    return Err(Error::Input(
                        "curved shell must coincide with one full quadratic soil face".into(),
                    ));
                }
            }
        }
        let rows: Vec<_> = (0..mesh.nodes.len())
            .map(|i| {
                let base = if let Some(s) = bindings.get(&i) {
                    6 * s
                } else {
                    let b = size;
                    size += 3;
                    b
                };
                [base, base + 1, base + 2]
            })
            .collect();
        let (k, cells) = path::assemble_soil(mesh, &input.materials, size, &rows)?;
        Ok(Self {
            stiffness: Some(k),
            cells,
            size,
            rows,
        })
    }
    pub fn force(&self, q: &DVector<f64>) -> DVector<f64> {
        self.stiffness
            .as_ref()
            .map_or_else(|| DVector::zeros(self.size), |k| k * q)
    }
    pub fn stress(&self, q: &DVector<f64>) -> Vec<Vec<[f64; 6]>> {
        self.cells.iter().map(|c| c.stresses(q)).collect()
    }
    pub fn translations(&self, q: &DVector<f64>) -> Vec<[f64; 3]> {
        self.rows
            .iter()
            .map(|r| Vector3::new(q[r[0]], q[r[1]], q[r[2]]).into())
            .collect()
    }
}
