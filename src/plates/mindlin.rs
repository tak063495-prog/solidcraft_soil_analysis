//! Straight, planar quadratic Reissner–Mindlin triangle with MITC-type
//! covariant shear interpolation into the eight-dimensional Nedelec space.
//! See docs/PLATES.ja.md for the fixed tying functionals and energy contract.
use super::{MemberResult, Sample, Section};
use crate::{Error, Result};
use nalgebra::{DMatrix, DVector, Matrix2, Matrix3, SMatrix, Vector2, Vector3};
type B = SMatrix<f64, 3, 36>;
type H = SMatrix<f64, 2, 36>;
type D = SMatrix<f64, 1, 36>;

pub const POINTS: [[f64; 3]; 6] = [
    [0.816847572980459, 0.091576213509771, 0.091576213509771],
    [0.091576213509771, 0.816847572980459, 0.091576213509771],
    [0.091576213509771, 0.091576213509771, 0.816847572980459],
    [0.108103018168070, 0.445948490915965, 0.445948490915965],
    [0.445948490915965, 0.108103018168070, 0.445948490915965],
    [0.445948490915965, 0.445948490915965, 0.108103018168070],
];
use std::ops::AddAssign;

/// Same ordering as the Tet10 face: vertices, then edges 01,12,20.
pub fn shape(l: [f64; 3]) -> [f64; 6] {
    [
        l[0] * (2. * l[0] - 1.),
        l[1] * (2. * l[1] - 1.),
        l[2] * (2. * l[2] - 1.),
        4. * l[0] * l[1],
        4. * l[1] * l[2],
        4. * l[2] * l[0],
    ]
}
fn basis(r: f64, s: f64) -> SMatrix<f64, 2, 8> {
    SMatrix::from_row_slice(&[
        1.,
        r,
        s,
        0.,
        0.,
        0.,
        -r * s,
        -s * s,
        0.,
        0.,
        0.,
        1.,
        r,
        s,
        r * r,
        r * s,
    ])
}
pub struct Kernel {
    pub area_m2: f64,
    pub axes: Matrix3<f64>,
    pub stiffness: DMatrix<f64>,
    transform: SMatrix<f64, 36, 36>,
    gradient: [Vector2<f64>; 3],
    jacobian: Matrix2<f64>,
    shear_coefficients: SMatrix<f64, 8, 36>,
    membrane: Matrix3<f64>,
    bending: Matrix3<f64>,
    shear: f64,
    drilling: f64,
}
impl Kernel {
    pub fn new(points: [[f64; 3]; 6], section: &Section) -> Result<Self> {
        section.validate()?;
        let base = super::Kernel::new([points[0], points[1], points[2]], section)?;
        let x = points.map(Vector3::from);
        let scale = (x[1] - x[0])
            .norm()
            .max((x[2] - x[0]).norm())
            .max((x[2] - x[1]).norm());
        if !points.iter().flatten().all(|v| v.is_finite()) {
            return Err(Error::Input("non-finite Mindlin Tri6 geometry".into()));
        }
        for (e, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            if (x[3 + e] - (x[i] + x[j]) * 0.5).norm() > 64. * f64::EPSILON * scale {
                return Err(Error::Input("Mindlin Tri6 requires straight edges and a planar reference surface; curved/nonplanar/offset attachment unsupported".into()));
            }
        }
        let xy = x.map(|p| {
            let q = base.axes * (p - x[0]);
            Vector2::new(q[0], q[1])
        });
        let jacobian = Matrix2::from_columns(&[xy[1], xy[2]]);
        let mut transform = SMatrix::<f64, 36, 36>::zeros();
        for block in 0..12 {
            transform
                .fixed_view_mut::<3, 3>(3 * block, 3 * block)
                .copy_from(&base.axes);
        }
        let nu = section.poisson;
        let c = Matrix3::new(1., nu, 0., nu, 1., 0., 0., 0., (1. - nu) / 2.)
            * (section.young_kpa / (1. - nu * nu));
        let g = section.young_kpa / (2. * (1. + nu));
        let mut kernel = Self {
            area_m2: base.area_m2,
            axes: base.axes,
            stiffness: DMatrix::zeros(36, 36),
            transform,
            gradient: base.gradient,
            jacobian,
            shear_coefficients: SMatrix::zeros(),
            membrane: c * section.thickness_m,
            bending: c * (section.thickness_m.powi(3) / 12.),
            shear: (5. / 6.) * g * section.thickness_m,
            drilling: section.drilling_factor * g * section.thickness_m,
        };
        let mut tying = SMatrix::<f64, 8, 8>::zeros();
        let mut raw = SMatrix::<f64, 8, 36>::zeros();
        let vertices = [
            Vector2::new(0., 0.),
            Vector2::new(1., 0.),
            Vector2::new(0., 1.),
        ];
        for (e, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            let tangent = vertices[j] - vertices[i];
            for (q, z) in [(1. - 1. / 3_f64.sqrt()) / 2., (1. + 1. / 3_f64.sqrt()) / 2.]
                .into_iter()
                .enumerate()
            {
                let p = vertices[i] + tangent * z;
                let h = kernel.raw_shear([1. - p[0] - p[1], p[0], p[1]]);
                tying
                    .row_mut(2 * e + q)
                    .copy_from(&(tangent.transpose() * basis(p[0], p[1])));
                raw.row_mut(2 * e + q)
                    .copy_from(&(tangent.transpose() * jacobian.transpose() * h));
            }
        }
        // Interior moments (area-normalised); degree-two exact three-point rule.
        for l in super::POINTS {
            let h = jacobian.transpose() * kernel.raw_shear(l) / 3.;
            let b = basis(l[1], l[2]) / 3.;
            for a in 0..2 {
                tying.row_mut(6 + a).add_assign(b.row(a));
                raw.row_mut(6 + a).add_assign(h.row(a));
            }
        }
        kernel.shear_coefficients = tying
            .lu()
            .solve(&raw)
            .ok_or_else(|| Error::Input("singular MITC tying functionals".into()))?;
        let mut local = SMatrix::<f64, 36, 36>::zeros();
        for (l, w) in quadrature() {
            let (bm, bb, bs, bd) = kernel.operators(l);
            local += (bm.transpose() * kernel.membrane * bm
                + bb.transpose() * kernel.bending * bb
                + bs.transpose() * bs * kernel.shear
                + bd.transpose() * bd * kernel.drilling)
                * (w * kernel.area_m2);
        }
        let k = transform.transpose() * local * transform;
        kernel.stiffness = DMatrix::from_column_slice(36, 36, k.as_slice());
        if !kernel.stiffness.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("Mindlin stiffness overflow".into()));
        }
        Ok(kernel)
    }
    fn derivatives(&self, l: [f64; 3]) -> [Vector2<f64>; 6] {
        let mut d = [Vector2::zeros(); 6];
        for i in 0..3 {
            d[i] = self.gradient[i] * (4. * l[i] - 1.);
        }
        for (e, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            d[3 + e] = (self.gradient[j] * l[i] + self.gradient[i] * l[j]) * 4.;
        }
        d
    }
    fn raw_shear(&self, l: [f64; 3]) -> H {
        let mut h = H::zeros();
        let n = shape(l);
        for (i, d) in self.derivatives(l).iter().enumerate() {
            h[(0, 6 * i + 2)] = d[0];
            h[(1, 6 * i + 2)] = d[1];
            h[(0, 6 * i + 4)] = n[i];
            h[(1, 6 * i + 3)] = -n[i];
        }
        h
    }
    fn operators(&self, l: [f64; 3]) -> (B, B, H, D) {
        let mut bm = B::zeros();
        let mut bb = B::zeros();
        let mut bd = D::zeros();
        let n = shape(l);
        for (i, d) in self.derivatives(l).iter().enumerate() {
            bm[(0, 6 * i)] = d[0];
            bm[(1, 6 * i + 1)] = d[1];
            bm[(2, 6 * i)] = d[1];
            bm[(2, 6 * i + 1)] = d[0];
            bb[(0, 6 * i + 4)] = d[0];
            bb[(1, 6 * i + 3)] = -d[1];
            bb[(2, 6 * i + 4)] = d[1];
            bb[(2, 6 * i + 3)] = -d[0];
            bd[(0, 6 * i)] = d[1] / 2.;
            bd[(0, 6 * i + 1)] = -d[0] / 2.;
            bd[(0, 6 * i + 5)] = n[i];
        }
        let bs = self.jacobian.transpose().try_inverse().unwrap()
            * basis(l[1], l[2])
            * self.shear_coefficients;
        (bm, bb, bs, bd)
    }
    /// Global translations on the same isoparametric reference face.
    pub fn translation(&self, l: [f64; 3], global: &DVector<f64>) -> [f64; 3] {
        let n = shape(l);
        std::array::from_fn(|a| (0..6).map(|i| n[i] * global[6 * i + a]).sum())
    }
    pub fn uniform_load(&self, traction: [f64; 3]) -> DVector<f64> {
        DVector::from_fn(36, |i, _| {
            if i / 6 >= 3 && i % 6 < 3 {
                traction[i % 6] * self.area_m2 / 3.
            } else {
                0.
            }
        })
    }
    pub fn recover(&self, id: String, global: &DVector<f64>) -> MemberResult {
        let d = self.transform * nalgebra::SVector::<f64, 36>::from_column_slice(global.as_slice());
        let mut energy = [0.; 4];
        let mut samples = Vec::new();
        for (l, w) in quadrature() {
            let (bm, bb, bs, bd) = self.operators(l);
            let strain = bm * d;
            let curvature = bb * d;
            let gamma = bs * d;
            let drill = (bd * d)[0];
            let force = self.membrane * strain;
            let moment = self.bending * curvature;
            let shear = gamma * self.shear;
            for (a, v) in [
                strain.dot(&force),
                curvature.dot(&moment),
                gamma.dot(&shear),
                drill * drill * self.drilling,
            ]
            .into_iter()
            .enumerate()
            {
                energy[a] += v * w * self.area_m2 / 2.;
            }
            samples.push(Sample {
                barycentric: l,
                membrane_strain: strain.into(),
                curvature_per_m: curvature.into(),
                membrane_force_kn_m: force.into(),
                bending_moment_kn: moment.into(),
                drilling_mismatch_rad: drill,
                transverse_shear_strain: gamma.into(),
                transverse_shear_force_kn_m: shear.into(),
            });
        }
        MemberResult {
            id,
            area_m2: self.area_m2,
            local_axes: std::array::from_fn(|a| std::array::from_fn(|b| self.axes[(a, b)])),
            samples,
            membrane_energy_kn_m: energy[0],
            bending_energy_kn_m: energy[1],
            transverse_shear_energy_kn_m: energy[2],
            drilling_penalty_energy_kn_m: energy[3],
        }
    }
}
fn quadrature() -> impl Iterator<Item = ([f64; 3], f64)> {
    POINTS.into_iter().enumerate().map(|(i, l)| {
        (
            l,
            if i < 3 {
                0.109951743655322
            } else {
                0.223381589678011
            },
        )
    })
}
