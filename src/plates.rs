//! Small-strain, flat CST membrane / discrete Kirchhoff triangle bending.
//!
//! DKT rotations are quadratic. Edge tangential rotations follow the midpoint
//! derivative of the cubic Hermite deflection; normal rotations are linear.
//! A declared drilling penalty ties the interpolated normal rotation to the
//! membrane spin. Its energy is reported separately: it is not plate bending.
use crate::{Error, Result};
use nalgebra::{Matrix3, SMatrix, SVector, Vector2, Vector3};
use serde::{Deserialize, Serialize};

pub type M18 = SMatrix<f64, 18, 18>;
pub type V18 = SVector<f64, 18>;
type B3 = SMatrix<f64, 3, 18>;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Section {
    pub id: String,
    pub young_kpa: f64,
    pub poisson: f64,
    pub thickness_m: f64,
    pub unit_weight_kn_m3: f64,
    /// Dimensionless multiplier of G t in the membrane-spin penalty.
    pub drilling_factor: f64,
}

impl Section {
    pub fn validate(&self) -> Result<()> {
        crate::elastic::constitutive(self.young_kpa, self.poisson)?;
        if self.id.is_empty()
            || !self.thickness_m.is_finite()
            || self.thickness_m <= 0.
            || !self.unit_weight_kn_m3.is_finite()
            || self.unit_weight_kn_m3 < 0.
            || !self.drilling_factor.is_finite()
            || !(1e-8..=0.01).contains(&self.drilling_factor)
        {
            return Err(Error::Input(
                "invalid plate section/thickness/drilling factor".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Member {
    pub id: String,
    pub nodes: [String; 3],
    pub section_id: String,
    pub initially_active: bool,
    #[serde(default)]
    pub constant_traction_kpa: [f64; 3],
    #[serde(default)]
    pub boundary_traction_kpa: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sample {
    pub barycentric: [f64; 3],
    /// [epsilon_xx, epsilon_yy, gamma_xy], tension-positive.
    pub membrane_strain: [f64; 3],
    /// [beta_x,x, beta_y,y, beta_x,y + beta_y,x] [1/m].
    pub curvature_per_m: [f64; 3],
    /// [Nxx, Nyy, Nxy] [kN/m], tension-positive.
    pub membrane_force_kn_m: [f64; 3],
    /// [Mxx, Myy, Mxy] [kN m/m], positive with the curvature above.
    pub bending_moment_kn: [f64; 3],
    pub drilling_mismatch_rad: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemberResult {
    pub id: String,
    pub area_m2: f64,
    pub local_axes: [[f64; 3]; 3],
    pub samples: Vec<Sample>,
    pub membrane_energy_kn_m: f64,
    pub bending_energy_kn_m: f64,
    pub drilling_penalty_energy_kn_m: f64,
}

pub struct Kernel {
    pub area_m2: f64,
    pub axes: Matrix3<f64>,
    pub transform: M18,
    pub stiffness: M18,
    membrane: Matrix3<f64>,
    bending: Matrix3<f64>,
    drilling: f64,
    gradient: [Vector2<f64>; 3],
    /// Quadratic rotation interpolation at corners and edge midpoints.
    beta: [SMatrix<f64, 2, 18>; 6],
    bm: B3,
}

const POINTS: [[f64; 3]; 3] = [
    [2. / 3., 1. / 6., 1. / 6.],
    [1. / 6., 2. / 3., 1. / 6.],
    [1. / 6., 1. / 6., 2. / 3.],
];

impl Kernel {
    pub fn new(points: [[f64; 3]; 3], section: &Section) -> Result<Self> {
        section.validate()?;
        if !points.iter().flatten().all(|v| v.is_finite()) {
            return Err(Error::Input("non-finite plate geometry".into()));
        }
        let x = points.map(Vector3::from);
        let edge = x[1] - x[0];
        let normal = edge.cross(&(x[2] - x[0]));
        let scale = edge
            .norm()
            .max((x[2] - x[0]).norm())
            .max((x[2] - x[1]).norm());
        if scale <= 1e-12 || normal.norm() <= 1e-10 * scale * scale {
            return Err(Error::Input("degenerate/ill-shaped plate triangle".into()));
        }
        let ex = edge.normalize();
        let ez = normal.normalize();
        let ey = ez.cross(&ex);
        let axes = Matrix3::from_rows(&[ex.transpose(), ey.transpose(), ez.transpose()]);
        let xy = x.map(|p| Vector2::new(ex.dot(&(p - x[0])), ey.dot(&(p - x[0]))));
        let twice_area = normal.norm();
        let gradient = std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            Vector2::new(xy[j][1] - xy[k][1], xy[k][0] - xy[j][0]) / twice_area
        });
        let mut transform = M18::zeros();
        for n in 0..3 {
            for part in 0..2 {
                for a in 0..3 {
                    for b in 0..3 {
                        transform[(6 * n + 3 * part + a, 6 * n + 3 * part + b)] = axes[(a, b)];
                    }
                }
            }
        }
        let mut beta = [SMatrix::<f64, 2, 18>::zeros(); 6];
        let mut bm = B3::zeros();
        for i in 0..3 {
            beta[i][(0, 6 * i + 4)] = 1.;
            beta[i][(1, 6 * i + 3)] = -1.;
            bm[(0, 6 * i)] = gradient[i][0];
            bm[(1, 6 * i + 1)] = gradient[i][1];
            bm[(2, 6 * i)] = gradient[i][1];
            bm[(2, 6 * i + 1)] = gradient[i][0];
        }
        for (e, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            let delta = xy[j] - xy[i];
            let length = delta.norm();
            let tangent = delta / length;
            let sum = beta[i] + beta[j];
            beta[3 + e] = sum * 0.5 - tangent * (tangent.transpose() * sum) * 0.75;
            for a in 0..2 {
                beta[3 + e][(a, 6 * i + 2)] += 1.5 * tangent[a] / length;
                beta[3 + e][(a, 6 * j + 2)] -= 1.5 * tangent[a] / length;
            }
        }
        let nu = section.poisson;
        let c = Matrix3::new(1., nu, 0., nu, 1., 0., 0., 0., (1. - nu) / 2.)
            * (section.young_kpa / (1. - nu * nu));
        let mut kernel = Self {
            area_m2: twice_area / 2.,
            axes,
            transform,
            stiffness: M18::zeros(),
            membrane: c * section.thickness_m,
            bending: c * (section.thickness_m.powi(3) / 12.),
            drilling: section.drilling_factor * section.young_kpa / (2. * (1. + nu))
                * section.thickness_m,
            gradient,
            beta,
            bm,
        };
        let mut local = kernel.bm.transpose() * kernel.membrane * kernel.bm * kernel.area_m2;
        for l in POINTS {
            let (bb, bd) = kernel.operators(l);
            local += (bb.transpose() * kernel.bending * bb + bd.transpose() * bd * kernel.drilling)
                * (kernel.area_m2 / 3.);
        }
        kernel.stiffness = kernel.transform.transpose() * local * kernel.transform;
        if !kernel.stiffness.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("plate stiffness overflow".into()));
        }
        Ok(kernel)
    }

    fn operators(&self, l: [f64; 3]) -> (B3, SMatrix<f64, 1, 18>) {
        let mut derivative = [SMatrix::<f64, 2, 18>::zeros(); 2];
        for (i, value) in l.iter().enumerate() {
            for (a, d) in derivative.iter_mut().enumerate() {
                *d += self.beta[i] * ((4. * value - 1.) * self.gradient[i][a]);
            }
        }
        for (e, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            for (a, d) in derivative.iter_mut().enumerate() {
                *d += self.beta[3 + e]
                    * (4. * (l[i] * self.gradient[j][a] + l[j] * self.gradient[i][a]));
            }
        }
        let mut bb = B3::zeros();
        bb.row_mut(0).copy_from(&derivative[0].row(0));
        bb.row_mut(1).copy_from(&derivative[1].row(1));
        bb.row_mut(2)
            .copy_from(&(derivative[1].row(0) + derivative[0].row(1)));
        let mut bd = SMatrix::<f64, 1, 18>::zeros();
        for (i, shape) in l.into_iter().enumerate() {
            bd[(0, 6 * i)] = 0.5 * self.gradient[i][1];
            bd[(0, 6 * i + 1)] = -0.5 * self.gradient[i][0];
            bd[(0, 6 * i + 5)] = shape;
        }
        (bb, bd)
    }

    /// Uniform dead traction in global axes, using linear transverse load work.
    pub fn uniform_load(&self, traction: [f64; 3]) -> V18 {
        V18::from_fn(|i, _| {
            if i % 6 < 3 {
                traction[i % 6] * self.area_m2 / 3.
            } else {
                0.
            }
        })
    }

    pub fn recover(&self, id: String, global_relative: &V18) -> MemberResult {
        let d = self.transform * global_relative;
        let strain = self.bm * d;
        let force = self.membrane * strain;
        let mut bending_energy = 0.;
        let mut drilling_energy = 0.;
        let samples = POINTS
            .into_iter()
            .map(|l| {
                let (bb, bd) = self.operators(l);
                let curvature = bb * d;
                let moment = self.bending * curvature;
                let drilling = (bd * d)[0];
                bending_energy += curvature.dot(&moment) * self.area_m2 / 6.;
                drilling_energy += drilling * drilling * self.drilling * self.area_m2 / 6.;
                Sample {
                    barycentric: l,
                    membrane_strain: strain.into(),
                    curvature_per_m: curvature.into(),
                    membrane_force_kn_m: force.into(),
                    bending_moment_kn: moment.into(),
                    drilling_mismatch_rad: drilling,
                }
            })
            .collect();
        MemberResult {
            id,
            area_m2: self.area_m2,
            local_axes: std::array::from_fn(|i| std::array::from_fn(|j| self.axes[(i, j)])),
            samples,
            membrane_energy_kn_m: strain.dot(&force) * self.area_m2 / 2.,
            bending_energy_kn_m: bending_energy,
            drilling_penalty_energy_kn_m: drilling_energy,
        }
    }
}
