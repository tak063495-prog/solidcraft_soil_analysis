//! Explicit flat Tri6 offset scope. Unknown translations and structural node
//! positions belong to the attachment surface. All three rotations are global.
//! x_p = x_s - e n; u_p = u_s - theta cross (e n) = u_s + e [n]x theta.
//! A maps attachment to neutral kinematics; K_s=A^T K_p A, f_s=A^T f_p.
use super::{MemberResult, Section, mindlin};
use crate::{Error, Result};
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

pub struct Kernel {
    pub neutral: mindlin::Kernel,
    pub stiffness: DMatrix<f64>,
    pub offset_m: f64,
    transform: DMatrix<f64>,
}
impl Kernel {
    pub fn new(points: [[f64; 3]; 6], section: &Section, e: f64) -> Result<Self> {
        if !e.is_finite() || e.abs() > section.thickness_m / 2. {
            return Err(Error::Input(
                "offset scope requires finite |e| <= thickness/2".into(),
            ));
        }
        // Translation of a flat plane does not change its operators. Building
        // at x_s preserves the exact e=0 arithmetic of the existing kernel.
        let neutral = mindlin::Kernel::new(points, section)?;
        let n = neutral.axes.row(2).transpose();
        let cross = Matrix3::new(0., -n[2], n[1], n[2], 0., -n[0], -n[1], n[0], 0.) * e;
        let mut transform = DMatrix::identity(36, 36);
        for i in 0..6 {
            transform
                .view_mut((6 * i, 6 * i + 3), (3, 3))
                .copy_from(&cross);
        }
        let stiffness = if e == 0. {
            neutral.stiffness.clone()
        } else {
            transform.transpose() * &neutral.stiffness * &transform
        };
        Ok(Self {
            neutral,
            stiffness,
            offset_m: e,
            transform,
        })
    }
    pub fn neutral_kinematics(&self, surface: &DVector<f64>) -> DVector<f64> {
        if self.offset_m == 0. {
            surface.clone()
        } else {
            &self.transform * surface
        }
    }
    pub fn translation(&self, l: [f64; 3], surface: &DVector<f64>) -> [f64; 3] {
        self.neutral
            .translation(l, &self.neutral_kinematics(surface))
    }
    pub fn attachment_translation(&self, l: [f64; 3], surface: &DVector<f64>) -> [f64; 3] {
        let mut u = Vector3::from(self.translation(l, surface));
        let shape = mindlin::shape(l);
        let theta =
            Vector3::from_fn(|a, _| (0..6).map(|i| shape[i] * surface[6 * i + 3 + a]).sum());
        u += theta.cross(&(self.neutral.axes.row(2).transpose() * self.offset_m));
        u.into()
    }
    /// Dead traction and self-weight act on the neutral plane; their moment
    /// about x_s is included by the same virtual-work transformation.
    pub fn uniform_load(&self, traction: [f64; 3]) -> DVector<f64> {
        let f = self.neutral.uniform_load(traction);
        if self.offset_m == 0. {
            f
        } else {
            self.transform.transpose() * f
        }
    }
    pub fn recover(&self, id: String, surface: &DVector<f64>) -> MemberResult {
        self.neutral.recover(id, &self.neutral_kinematics(surface))
    }
}
