//! Small-rotation shell on the actual attachment surface, with strain tying
//! in the centroid basis (Kim & Bathe 2009, Fig. 1 and equations 10–13).
pub mod geometry;
pub mod path;
pub mod soil;
use crate::plates::{Section, mindlin};
use crate::{Error, Result};
use geometry::{Geometry, Surface, V};
use nalgebra::{Matrix3, SMatrix, SVector};
use serde::{Deserialize, Serialize};

pub type K36 = SMatrix<f64, 36, 36>;
pub type Q36 = SVector<f64, 36>;
type Row = SMatrix<f64, 1, 36>;
type Strain = SMatrix<f64, 6, 36>;
type Stretch = SMatrix<f64, 3, 36>;
type Shear = SMatrix<f64, 2, 36>;
const CENTER: [f64; 2] = [1. / 3.; 2];
const EDGES: [([f64; 2], [f64; 2]); 3] = [
    ([0., 0.], [1., 0.]),
    ([0., 0.], [0., 1.]),
    ([1., 0.], [-1., 1.]),
];

#[derive(Clone, Copy)]
pub struct Tying {
    pub membrane: bool,
    pub shear: bool,
}
impl Default for Tying {
    fn default() -> Self {
        Self {
            membrane: true,
            shear: true,
        }
    }
}

fn inverse(m: Matrix3<f64>) -> Result<Matrix3<f64>> {
    if !m.determinant().is_finite() || m.determinant() <= 0. {
        return Err(Error::Input(
            "inverted curved Tri6 thickness mapping".into(),
        ));
    }
    m.try_inverse()
        .ok_or_else(|| Error::Input("singular curved Tri6 thickness mapping".into()))
}

// Tensor components rr, ss, rs, rz, sz, zz, without engineering factors.
fn components(m: &Matrix3<f64>) -> SVector<f64, 6> {
    SVector::from_row_slice(&[
        m[(0, 0)],
        m[(1, 1)],
        m[(0, 1)],
        m[(0, 2)],
        m[(1, 2)],
        m[(2, 2)],
    ])
}
fn tensor(v: impl std::ops::Index<usize, Output = f64>) -> Matrix3<f64> {
    Matrix3::new(v[0], v[2], v[3], v[2], v[1], v[4], v[3], v[4], v[5])
}
fn transform(b: &Strain, map: &Matrix3<f64>) -> Strain {
    let mut out = Strain::zeros();
    for (j, v) in b.column_iter().enumerate() {
        out.set_column(j, &components(&(map.transpose() * tensor(v) * map)));
    }
    out
}
fn displacement(s: &Surface, geo: &Geometry, z: f64) -> Stretch {
    let mut u = Stretch::zeros();
    for (i, d) in geo.directors.iter().enumerate() {
        for a in 0..3 {
            let axis = V::ith(a, 1.);
            u.set_column(6 * i + a, &(axis * s.n[i]));
            u.set_column(
                6 * i + 3 + a,
                &(s.variation(i, axis.cross(d)).0 * (z - geo.e)),
            );
        }
    }
    u
}
fn raw(geo: &Geometry, p: [f64; 2], z: f64, centroid: &Matrix3<f64>) -> Result<Strain> {
    let s = geo.surface(p)?;
    let basis = geo.basis(&s, z);
    let map = inverse(basis)? * centroid;
    let mut out = Strain::zeros();
    for (i, d) in geo.directors.iter().enumerate() {
        for a in 0..3 {
            let axis = V::ith(a, 1.);
            let (tilt, dt) = s.variation(i, axis.cross(d));
            let du = Matrix3::from_columns(&[axis * s.dn[0][i], axis * s.dn[1][i], V::zeros()]);
            let rot = Matrix3::from_columns(&[dt[0] * (z - geo.e), dt[1] * (z - geo.e), tilt]);
            for (col, grad) in [(6 * i + a, du), (6 * i + a + 3, rot)] {
                let x = centroid.transpose() * grad * map;
                out.set_column(col, &components(&((x + x.transpose()) * 0.5)));
            }
        }
    }
    Ok(out)
}
fn shear_space([r, s]: [f64; 2]) -> SMatrix<f64, 2, 8> {
    SMatrix::from_row_slice(&[
        1.,
        r,
        s,
        0.,
        0.,
        0.,
        s * r,
        s * s,
        0.,
        0.,
        0.,
        1.,
        r,
        s,
        -r * r,
        -r * s,
    ])
}
struct Layer {
    z: f64,
    center: Matrix3<f64>,
    inv_center: Matrix3<f64>,
    normal: [Stretch; 3],
    transverse: SMatrix<f64, 8, 36>,
}
impl Layer {
    fn new(geo: &Geometry, z: f64) -> Result<Self> {
        let center = geo.basis(&geo.surface(CENTER)?, z);
        let mut normal = [Stretch::zeros(); 3];
        let mut shear_matrix = SMatrix::<f64, 8, 8>::zeros();
        let mut shear_samples = SMatrix::<f64, 8, 36>::zeros();
        let a = (1. - 1. / 3_f64.sqrt()) / 2.;
        let interior = [[a, 1. - 2. * a], [1. - 2. * a, a], [a, a]];
        for (edge, (origin, t)) in EDGES.into_iter().enumerate() {
            let mut points = Matrix3::zeros();
            let mut samples = Stretch::zeros();
            let along = |v| [origin[0] + v * t[0], origin[1] + v * t[1]];
            for (j, p) in [along(a), along(1. - a), interior[edge]]
                .into_iter()
                .enumerate()
            {
                let b = raw(geo, p, z, &center)?;
                points
                    .row_mut(j)
                    .copy_from(&V::new(1., p[0], p[1]).transpose());
                samples.row_mut(j).copy_from(
                    &(b.row(0) * t[0].powi(2)
                        + b.row(1) * t[1].powi(2)
                        + b.row(2) * (2. * t[0] * t[1])),
                );
                if j < 2 {
                    let space = shear_space(p);
                    shear_matrix
                        .row_mut(2 * edge + j)
                        .copy_from(&(space.row(0) * t[0] + space.row(1) * t[1]));
                    shear_samples
                        .row_mut(2 * edge + j)
                        .copy_from(&(b.row(3) * t[0] + b.row(4) * t[1]));
                }
            }
            normal[edge] = points
                .lu()
                .solve(&samples)
                .ok_or_else(|| Error::Input("singular membrane tying system".into()))?;
        }
        shear_matrix
            .fixed_rows_mut::<2>(6)
            .copy_from(&shear_space(CENTER));
        shear_samples
            .fixed_rows_mut::<2>(6)
            .copy_from(&raw(geo, CENTER, z, &center)?.fixed_rows::<2>(3));
        let transverse = shear_matrix
            .lu()
            .solve(&shear_samples)
            .ok_or_else(|| Error::Input("singular transverse tying system".into()))?;
        Ok(Self {
            z,
            center,
            inv_center: inverse(center)?,
            normal,
            transverse,
        })
    }
    fn operators(
        &self,
        geo: &Geometry,
        p: [f64; 2],
        axes: &Matrix3<f64>,
        tying: Tying,
    ) -> Result<(Stretch, Shear)> {
        let mut b = raw(geo, p, self.z, &self.center)?;
        if tying.membrane {
            let h = V::new(1., p[0], p[1]).transpose();
            let [rr, ss, diag] = self.normal.map(|c| h * c);
            b.row_mut(0).copy_from(&rr);
            b.row_mut(1).copy_from(&ss);
            b.row_mut(2).copy_from(&((rr + ss - diag) * 0.5));
        }
        if tying.shear {
            b.fixed_rows_mut::<2>(3)
                .copy_from(&(shear_space(p) * self.transverse));
        }
        let b = transform(&b, &(self.inv_center * axes));
        let mut inplane = b.fixed_rows::<3>(0).into_owned();
        inplane.row_mut(2).scale_mut(2.);
        Ok((inplane, b.fixed_rows::<2>(3) * 2.))
    }
}

fn spin(s: &Surface) -> Result<Row> {
    let [a, b] = s.tangent;
    let n = a.cross(&b).normalize();
    let inv = inverse(Matrix3::from_columns(&[a, b, n]))?;
    let mut row = Row::zeros();
    for i in 0..6 {
        for k in 0..3 {
            let axis = V::ith(k, 1.);
            let grad =
                Matrix3::from_columns(&[axis * s.dn[0][i], axis * s.dn[1][i], V::zeros()]) * inv;
            let dn = -grad.transpose() * n;
            let curl = V::new(
                grad[(2, 1)] - grad[(1, 2)],
                grad[(0, 2)] - grad[(2, 0)],
                grad[(1, 0)] - grad[(0, 1)],
            );
            row[6 * i + k] = -s.d.dot(&(n.cross(&dn) + n * (0.5 * n.dot(&curl))));
            row[6 * i + 3 + k] = s.n[i] * s.d[k];
        }
    }
    // The drilling rotation is the interpolated *full* rotation about d.
    // Projection onto each nodal axis before interpolation would lose rigid spin.
    Ok(row)
}

pub struct Kernel {
    pub stiffness: K36,
    pub drill: K36,
    pub traction: Stretch,
    pub gravity: Q36,
    pub samples: Vec<Point>,
    elastic: Matrix3<f64>,
    shear_modulus: f64,
    drill_points: Vec<(Row, f64)>,
}
pub struct Point {
    pub p: [f64; 2],
    pub z: f64,
    pub weight: f64,
    pub membrane: Stretch,
    pub transverse: Shear,
    pub spin: Row,
    pub axes: Matrix3<f64>,
    pub jacobian: f64,
    pub mid_area_weight: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Field {
    pub rs: [f64; 2],
    pub z_m: f64,
    pub jacobian: f64,
    pub weight: f64,
    pub strain: [f64; 3],
    pub shear_strain: [f64; 2],
    pub stress_kpa: [f64; 3],
    pub shear_kpa: [f64; 2],
    pub drill_rad: f64,
    pub material_axes: [[f64; 3]; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Resultant {
    pub rs: [f64; 2],
    pub mid_axes: [[f64; 3]; 3],
    pub n_kn_m: [f64; 3],
    pub m_kn: [f64; 3],
    pub m_attachment_kn: [f64; 3],
    pub q_kn_m: [f64; 2],
}
impl Kernel {
    pub fn new(geo: &Geometry, section: &Section, tying: Tying) -> Result<Self> {
        section.validate()?;
        let t = section.thickness_m;
        if !geo.e.is_finite() || geo.e.abs() > t / 2. {
            return Err(Error::Input("invalid signed curved shell offset".into()));
        }
        geo.certify(t)?;
        let nu = section.poisson;
        let elastic = Matrix3::new(1., nu, 0., nu, 1., 0., 0., 0., (1. - nu) / 2.)
            * (section.young_kpa / (1. - nu * nu));
        let g = section.young_kpa / (2. * (1. + nu));
        let mut k = Self {
            stiffness: K36::zeros(),
            drill: K36::zeros(),
            traction: Stretch::zeros(),
            gravity: Q36::zeros(),
            samples: vec![],
            elastic,
            shear_modulus: g * 5. / 6.,
            drill_points: vec![],
        };
        for (l, w) in mindlin::quadrature() {
            let p = [l[1], l[2]];
            let s = geo.surface(p)?;
            let mid = geo.basis(&s, 0.);
            let area = mid.column(0).cross(&mid.column(1)).norm() * w / 2.;
            let drill = spin(&s)?;
            k.drill_points
                .push((drill, section.drilling_factor * g * t * area));
            k.drill += drill.transpose() * drill * (section.drilling_factor * g * t * area);
            k.traction += displacement(&s, geo, 0.) * area;
        }
        for (xi, wz) in [
            (-(0.6_f64).sqrt(), 5. / 9.),
            (0., 8. / 9.),
            ((0.6_f64).sqrt(), 5. / 9.),
        ] {
            let z = xi * t / 2.;
            let layer = Layer::new(geo, z)?;
            for (l, w) in mindlin::quadrature() {
                let p = [l[1], l[2]];
                let s = geo.surface(p)?;
                let basis = geo.basis(&s, z);
                let jacobian = basis.determinant();
                let weight = w * wz * t * jacobian / 4.;
                let er = basis.column(1).cross(&s.d).normalize();
                let axes = Matrix3::from_columns(&[er, s.d.cross(&er), s.d]);
                let (membrane, transverse) = layer.operators(geo, p, &axes, tying)?;
                k.stiffness += (membrane.transpose() * elastic * membrane
                    + transverse.transpose() * transverse * (g * 5. / 6.))
                    * weight;
                k.gravity -= displacement(&s, geo, z).row(2).transpose()
                    * (section.unit_weight_kn_m3 * weight);
                k.samples.push(Point {
                    p,
                    z,
                    weight,
                    membrane,
                    transverse,
                    spin: spin(&s)?,
                    axes,
                    jacobian,
                    mid_area_weight: {
                        let m = geo.basis(&s, 0.);
                        m.column(0).cross(&m.column(1)).norm() * w / 2.
                    },
                });
            }
        }
        k.stiffness += k.drill;
        Ok(k)
    }

    /// Recover nodal forces from stresses and drilling work at the quadrature points.
    pub fn internal(&self, q: &Q36) -> Q36 {
        let mut f = Q36::zeros();
        for p in &self.samples {
            let stress = self.elastic * (p.membrane * q);
            let shear = p.transverse * q * self.shear_modulus;
            f += (p.membrane.transpose() * stress + p.transverse.transpose() * shear) * p.weight;
        }
        for (spin, scale) in &self.drill_points {
            f += spin.transpose() * ((spin * q)[0] * scale);
        }
        f
    }

    /// Integrate Cauchy stresses on a common material frame, per unit
    /// midsurface area. Use that same area and axes for both moments.
    pub fn resultants(&self, q: &Q36, e: f64) -> Vec<Resultant> {
        (0..mindlin::POINTS.len())
            .map(|i| {
                let mid = &self.samples[mindlin::POINTS.len() + i];
                let mut force = V::zeros();
                let mut moment = V::zeros();
                let mut shear = nalgebra::Vector2::zeros();
                for p in self.samples[i..].iter().step_by(mindlin::POINTS.len()) {
                    let s = self.elastic * (p.membrane * q);
                    let h = (p.transverse * q) * self.shear_modulus;
                    let stress = Matrix3::new(s.x, s.z, h.x, s.z, s.y, h.y, h.x, h.y, 0.);
                    let rotate = mid.axes.transpose() * p.axes;
                    let stress =
                        rotate * stress * rotate.transpose() * (p.weight / mid.mid_area_weight);
                    let inplane = V::new(stress[(0, 0)], stress[(1, 1)], stress[(0, 1)]);
                    force += inplane;
                    moment += inplane * p.z;
                    shear += nalgebra::Vector2::new(stress[(0, 2)], stress[(1, 2)]);
                }
                Resultant {
                    rs: mid.p,
                    mid_axes: mid.axes.into(),
                    n_kn_m: force.into(),
                    m_kn: moment.into(),
                    m_attachment_kn: (moment - e * force).into(),
                    q_kn_m: shear.into(),
                }
            })
            .collect()
    }
    pub fn fields(&self, q: &Q36) -> Vec<Field> {
        self.samples
            .iter()
            .map(|p| Field {
                rs: p.p,
                stress_kpa: (self.elastic * (p.membrane * q)).into(),
                shear_kpa: (p.transverse * q * self.shear_modulus).into(),
                z_m: p.z,
                jacobian: p.jacobian,
                weight: p.weight,
                strain: (p.membrane * q).into(),
                shear_strain: (p.transverse * q).into(),
                drill_rad: (p.spin * q)[0],
                material_axes: p.axes.into(),
            })
            .collect()
    }
}

/// Dead load on a specified layer. Pressure is positive into the selected
/// exterior side. The pointwise area vector and the kinematic variation are
/// evaluated on the same layer, including rotation of the offset arm.
pub fn surface_action(geo: &Geometry, t: f64, z: f64, traction: V, pressure: f64) -> Result<Q36> {
    if !(t.is_finite()
        && t > 0.
        && z.is_finite()
        && z.abs() <= t / 2.
        && traction.iter().all(|v| v.is_finite())
        && pressure.is_finite())
        || (pressure != 0. && z.abs() != t / 2.)
    {
        return Err(Error::Input(
            "invalid curved Tri6 load layer or pressure side".into(),
        ));
    }
    let mut action = Q36::zeros();
    for (l, w) in mindlin::quadrature() {
        let s = geo.surface([l[1], l[2]])?;
        let m = geo.basis(&s, z);
        let area = m.column(0).cross(&m.column(1));
        let force = traction * area.norm() - pressure * z.signum() * area;
        action += displacement(&s, geo, z).transpose() * force * (w / 2.);
    }
    Ok(action)
}
