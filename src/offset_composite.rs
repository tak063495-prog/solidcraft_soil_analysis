//! Dry, conforming offset plate on a Tet10 face. Translations on side A
//! are eliminated exactly. Only the plate has rotational unknowns.
use crate::plates::{Section, offset};
use crate::{Error, Result, elements, joint_law, mesh::Mesh};
use nalgebra::{DMatrix, DVector, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub mod path;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Face {
    pub cell: usize,
    /// Vertices then midsides 01, 12, 20 in the plate's node order.
    pub nodes: [usize; 6],
}

impl Face {
    fn validate(&self, mesh: &Mesh) -> Result<()> {
        let key = [self.nodes[0], self.nodes[1], self.nodes[2]];
        if self.cell >= mesh.tetrahedra.len() || mesh.tetrahedra10.is_empty() {
            return Err(Error::Input(
                "offset composite requires a Tet10 owner".into(),
            ));
        }
        let actual = mesh.face_nodes(self.cell, &key)?;
        let members = |v: &[usize]| v.iter().copied().collect::<BTreeSet<_>>();
        if members(&actual) != members(&self.nodes) {
            return Err(Error::Input(
                "offset composite face nodes disagree with owner".into(),
            ));
        }
        // Only exterior volume faces may supply a physical side.
        if mesh
            .tetrahedra
            .iter()
            .filter(|tet| key.iter().all(|n| tet.contains(n)))
            .count()
            != 1
        {
            return Err(Error::Input(
                "offset composite requires an exterior owner face".into(),
            ));
        }
        Ok(())
    }

    fn other_corner(&self, mesh: &Mesh) -> Vector3<f64> {
        let index = mesh.tetrahedra[self.cell]
            .iter()
            .find(|n| !self.nodes.contains(n))
            .unwrap();
        mesh.nodes[*index].into()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Contact {
    pub id: String,
    pub soil_b: Face,
    pub material: joint_law::Material,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Plate {
    pub id: String,
    pub section: Section,
    pub soil_a: Face,
    pub offset_m: f64,
    /// +1 or -1 relative to the ordered plate normal.
    pub exposed_side: i8,
    pub tie_at_physical_surface: bool,
    pub contacts: Vec<Contact>,
}

pub struct Point {
    pub weight: f64,
    pub position: [f64; 3],
    pub shape: Vec<f64>,
    /// Maps [six A translations & rotations, six B translations] to u_B-u_c.
    pub jump: DMatrix<f64>,
}
pub struct Patch {
    pub plate: offset::Kernel,
    pub normal: Vector3<f64>,
    pub arm: Vector3<f64>,
    pub points: Vec<Point>,
}
impl Patch {
    pub fn new(mesh: &Mesh, plate: &Plate) -> Result<Self> {
        mesh.validate()?;
        plate.soil_a.validate(mesh)?;
        if plate.id.is_empty() || ![-1, 1].contains(&plate.exposed_side) {
            return Err(Error::Input(
                "offset composite needs a plate ID and signed exposed side".into(),
            ));
        }
        let coordinates = plate.soil_a.nodes.map(|i| mesh.nodes[i]);
        let kernel = offset::Kernel::new(coordinates, &plate.section, plate.offset_m)?;
        let axis = kernel.neutral.axes.row(2).transpose().into_owned();
        let normal = axis * f64::from(plate.exposed_side);
        let half = plate.section.thickness_m / 2.;
        if plate.tie_at_physical_surface && plate.offset_m != -f64::from(plate.exposed_side) * half
        {
            return Err(Error::Input(
                "declared physical tie must lie on the opposite plate surface".into(),
            ));
        }
        let arm = normal * half - axis * plate.offset_m;
        let contact_x = coordinates.map(|x| <[f64; 3]>::from(Vector3::from(x) + arm));
        let x0 = Vector3::from(coordinates[0]);
        if normal.dot(&(plate.soil_a.other_corner(mesh) - x0)) >= 0. {
            return Err(Error::Input(
                "offset composite exposed side points into tied owner".into(),
            ));
        }
        let mut owners = BTreeSet::new();
        let mut ids = BTreeSet::new();
        for c in &plate.contacts {
            c.soil_b.validate(mesh)?;
            c.material.validate()?;
            let mut region = c.soil_b.nodes;
            region.sort_unstable();
            if c.soil_b.cell == plate.soil_a.cell
                || !owners.insert((c.soil_b.cell, region))
                || c.id.is_empty()
                || !ids.insert(&c.id)
            {
                return Err(Error::Input(format!(
                    "plate {} contact {} duplicates/self-pairs face {:?}, normal and tangential components",
                    plate.id, c.id, c.soil_b
                )));
            }
            let scale = kernel.neutral.area_m2.sqrt();
            if normal.dot(&(c.soil_b.other_corner(mesh) - Vector3::from(contact_x[0]))) <= 0.
                || c.soil_b.nodes.iter().zip(contact_x).any(|(&b, x)| {
                    (Vector3::from(mesh.nodes[b]) - Vector3::from(x)).norm()
                        > 64. * f64::EPSILON * scale
                })
            {
                return Err(Error::Input(
                    "offset composite requires a conforming oriented physical face pair".into(),
                ));
            }
        }
        let points = elements::face_points(&contact_x)?
            .into_iter()
            .map(|gp| {
                let mut jump = DMatrix::zeros(3, 54);
                for (i, &s) in gp.shape.iter().enumerate() {
                    for a in 0..3 {
                        jump[(a, 6 * i + a)] = -s;
                        jump[(a, 36 + 3 * i + a)] = s;
                        let direction = Vector3::ith(a, 1.);
                        jump.set_column(6 * i + 3 + a, &(-s * direction.cross(&arm)));
                    }
                }
                Point {
                    weight: Vector3::from(gp.area_vector).norm(),
                    position: gp.position,
                    shape: gp.shape,
                    jump,
                }
            })
            .collect();
        Ok(Self {
            plate: kernel,
            normal,
            arm,
            points,
        })
    }

    pub fn assemble_contact(
        &self,
        material: &joint_law::Material,
        old: &[joint_law::PointState],
        d: &DVector<f64>,
    ) -> Result<ContactTrial> {
        if old.len() != self.points.len() || d.len() != 54 {
            return Err(Error::Input("offset contact state layout mismatch".into()));
        }
        let mut trial = ContactTrial {
            points: Vec::new(),
            force: DVector::zeros(54),
            tangent: DMatrix::zeros(54, 54),
        };
        for (p, s) in self.points.iter().zip(old) {
            let u = &p.jump * d;
            let jump = [u[0], u[1], u[2]];
            let state = material.integrate_partition(s, jump, 1)?;
            let tangent = material.tangent(s, jump, 1, 1e-8)?;
            trial.force += p.jump.transpose() * Vector3::from(state.traction_kpa) * p.weight;
            trial.tangent += p.jump.transpose() * tangent * &p.jump * p.weight;
            trial.points.push(state);
        }
        Ok(trial)
    }
}

pub struct ContactTrial {
    pub points: Vec<joint_law::PointState>,
    pub force: DVector<f64>,
    pub tangent: DMatrix<f64>,
}
