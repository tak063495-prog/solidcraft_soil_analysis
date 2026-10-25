//! Dry unilateral Coulomb interface, psi=0, with two tangential slip directions.
//! Global vectors preserve history independently of the chosen tangent basis.
use crate::{Error, Result};
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub normal_stiffness_kpa_m: f64,
    pub tangential_stiffness_kpa_m: f64,
    pub friction_coefficient: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    Open,
    Stick,
    Slip,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointState {
    pub normal: [f64; 3],
    /// u_b - u_a; positive normal component opens the interface.
    pub jump_m: [f64; 3],
    /// Conjugate to jump: tension-positive normal, restoring shear on side b.
    pub traction_kpa: [f64; 3],
    pub plastic_slip_m: [f64; 3],
    pub accumulated_plastic_slip_m: f64,
    pub plastic_dissipation_kpa_m: f64,
    pub numerical_dissipation_kpa_m: f64,
    /// Backward-Euler endpoint work, separate from physical frictional loss.
    pub endpoint_work_kpa_m: f64,
    pub status: Status,
}
impl PointState {
    pub fn initial(normal: [f64; 3]) -> Result<Self> {
        let n = Vector3::from(normal);
        if !n.iter().all(|v| v.is_finite()) || n.norm() < 1e-12 {
            return Err(Error::Input("invalid JOINT normal".into()));
        }
        Ok(Self {
            normal: (n / n.norm()).into(),
            jump_m: [0.; 3],
            traction_kpa: [0.; 3],
            plastic_slip_m: [0.; 3],
            accumulated_plastic_slip_m: 0.,
            plastic_dissipation_kpa_m: 0.,
            numerical_dissipation_kpa_m: 0.,
            endpoint_work_kpa_m: 0.,
            status: Status::Stick,
        })
    }
    pub fn rotated(&self, r: &Matrix3<f64>) -> Result<Self> {
        if (r.transpose() * r - Matrix3::identity()).norm() > 1e-10
            || (r.determinant() - 1.).abs() > 1e-10
        {
            return Err(Error::Input(
                "JOINT history requires a proper rotation".into(),
            ));
        }
        let mut s = self.clone();
        s.normal = (r * Vector3::from(self.normal)).into();
        s.jump_m = (r * Vector3::from(self.jump_m)).into();
        s.traction_kpa = (r * Vector3::from(self.traction_kpa)).into();
        s.plastic_slip_m = (r * Vector3::from(self.plastic_slip_m)).into();
        Ok(s)
    }
}
impl Material {
    pub fn validate(&self) -> Result<()> {
        if !self.normal_stiffness_kpa_m.is_finite()
            || self.normal_stiffness_kpa_m <= 0.
            || !self.tangential_stiffness_kpa_m.is_finite()
            || self.tangential_stiffness_kpa_m <= 0.
            || !self.friction_coefficient.is_finite()
            || self.friction_coefficient < 0.
        {
            return Err(Error::Input(
                "JOINT requires finite positive stiffnesses and nonnegative friction".into(),
            ));
        }
        Ok(())
    }
    pub fn energy(&self, s: &PointState) -> f64 {
        let n = Vector3::from(s.normal);
        let q = Vector3::from(s.traction_kpa);
        let normal = q.dot(&n);
        let shear = q - normal * n;
        normal * normal / (2. * self.normal_stiffness_kpa_m)
            + shear.norm_squared() / (2. * self.tangential_stiffness_kpa_m)
    }
    pub fn validate_state(&self, s: &PointState) -> Result<()> {
        self.validate()?;
        if !s
            .normal
            .iter()
            .chain(&s.jump_m)
            .chain(&s.traction_kpa)
            .chain(&s.plastic_slip_m)
            .chain([
                &s.accumulated_plastic_slip_m,
                &s.plastic_dissipation_kpa_m,
                &s.numerical_dissipation_kpa_m,
                &s.endpoint_work_kpa_m,
            ])
            .all(|v| v.is_finite())
            || s.accumulated_plastic_slip_m < 0.
            || s.plastic_dissipation_kpa_m < 0.
            || s.numerical_dissipation_kpa_m < 0.
        {
            return Err(Error::Input("invalid JOINT history".into()));
        }
        let n = Vector3::from(s.normal);
        let jump = Vector3::from(s.jump_m);
        let q = Vector3::from(s.traction_kpa);
        let gn = n.dot(&jump);
        let qn = n.dot(&q);
        let tau = q - qn * n;
        let normal = self.normal_stiffness_kpa_m * gn.min(0.);
        let scale = 1. + normal.abs() + q.norm();
        let slip = Vector3::from(s.plastic_slip_m);
        if (n.norm() - 1.).abs() > 1e-10
            || (qn - normal).abs() > 1e-9 * scale
            || tau.norm() > self.friction_coefficient * (-normal) + 1e-9 * scale
            || slip.dot(&n).abs() > 1e-10 * (1. + slip.norm())
            || slip.norm()
                > s.accumulated_plastic_slip_m + 1e-10 * (1. + s.accumulated_plastic_slip_m)
            || (gn > 0. && (q.norm() > 1e-9 * scale || s.status != Status::Open))
            || (gn < 0. && s.status == Status::Open)
        {
            return Err(Error::Input(
                "JOINT stress/contact/slip history is inconsistent".into(),
            ));
        }
        let work = self.energy(s) + s.plastic_dissipation_kpa_m + s.numerical_dissipation_kpa_m;
        if (work - s.endpoint_work_kpa_m).abs()
            > 1e-8 * (1. + work.abs() + s.endpoint_work_kpa_m.abs())
        {
            return Err(Error::Input(
                "JOINT work/energy history is inconsistent".into(),
            ));
        }
        Ok(())
    }
    fn closed(&self, old: &PointState, target: Vector3<f64>) -> Result<PointState> {
        let n = Vector3::from(old.normal);
        let before = Vector3::from(old.jump_m);
        let delta = target - before;
        let qold = Vector3::from(old.traction_kpa);
        let tauold = qold - qold.dot(&n) * n;
        let dt = delta - delta.dot(&n) * n;
        let trial = tauold + self.tangential_stiffness_kpa_m * dt;
        let normal = self.normal_stiffness_kpa_m * target.dot(&n).min(0.);
        let limit = self.friction_coefficient * (-normal);
        let norm = trial.norm();
        let mut s = old.clone();
        let (tau, plastic) = if norm > limit && norm > 0. {
            (
                trial * (limit / norm),
                trial * ((norm - limit) / (self.tangential_stiffness_kpa_m * norm)),
            )
        } else {
            (trial, Vector3::zeros())
        };
        let q = normal * n + tau;
        s.jump_m = target.into();
        s.traction_kpa = q.into();
        s.status = if plastic.norm() > 0. {
            Status::Slip
        } else {
            Status::Stick
        };
        s.plastic_slip_m = (Vector3::from(old.plastic_slip_m) + plastic).into();
        s.accumulated_plastic_slip_m += plastic.norm();
        let dissipation = tau.dot(&plastic).max(0.);
        let work = q.dot(&delta);
        let numerical = work - (self.energy(&s) - self.energy(old)) - dissipation;
        let tol = 1e-10 * (1. + work.abs() + self.energy(old));
        if !numerical.is_finite() || numerical < -tol {
            return Err(Error::Material("JOINT discrete work balance failed".into()));
        }
        s.plastic_dissipation_kpa_m += dissipation;
        s.numerical_dissipation_kpa_m += numerical.max(0.);
        s.endpoint_work_kpa_m += work;
        Ok(s)
    }
    fn segment(&self, old: &PointState, target: Vector3<f64>) -> Result<PointState> {
        let n = Vector3::from(old.normal);
        let before = Vector3::from(old.jump_m);
        let g0 = n.dot(&before);
        let g1 = n.dot(&target);
        if g0 <= 0. && g1 <= 0. {
            return self.closed(old, target);
        }
        if g0 > 0. && g1 > 0. {
            let mut s = old.clone();
            s.jump_m = target.into();
            s.traction_kpa = [0.; 3];
            s.status = Status::Open;
            return Ok(s);
        }
        let crossing = before + (target - before) * (g0 / (g0 - g1));
        // Remove normal roundoff so a crossing is exactly on the contact plane.
        let crossing = crossing - n * crossing.dot(&n);
        if g0 <= 0. {
            let mut s = self.closed(old, crossing)?;
            s.jump_m = target.into();
            s.traction_kpa = [0.; 3];
            s.status = Status::Open;
            Ok(s)
        } else {
            let mut s = old.clone();
            s.jump_m = crossing.into();
            s.traction_kpa = [0.; 3];
            s.status = Status::Stick;
            self.closed(&s, target)
        }
    }
    pub fn partition(
        &self,
        old: &PointState,
        target: [f64; 3],
        maximum_increment_m: f64,
    ) -> Result<usize> {
        self.validate_state(old)?;
        if !target.iter().all(|v| v.is_finite())
            || !maximum_increment_m.is_finite()
            || maximum_increment_m <= 0.
        {
            return Err(Error::Input("invalid JOINT jump/increment".into()));
        }
        let count = ((Vector3::from(target) - Vector3::from(old.jump_m)).norm()
            / maximum_increment_m)
            .ceil()
            .max(1.);
        if count > 50000. {
            return Err(Error::Material(
                "JOINT increment partition budget exceeded".into(),
            ));
        }
        Ok(count as usize)
    }
    pub fn integrate_partition(
        &self,
        old: &PointState,
        target: [f64; 3],
        count: usize,
    ) -> Result<PointState> {
        self.validate_state(old)?;
        if count == 0 || count > 50000 || !target.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("invalid JOINT partition".into()));
        }
        let start = Vector3::from(old.jump_m);
        if target == old.jump_m {
            return Ok(old.clone());
        }
        let delta = Vector3::from(target) - start;
        let mut s = old.clone();
        for k in 1..=count {
            s = self.segment(
                &s,
                if k == count {
                    Vector3::from(target)
                } else {
                    start + delta * (k as f64 / count as f64)
                },
            )?;
        }
        self.validate_state(&s)?;
        Ok(s)
    }
    pub fn integrate(
        &self,
        old: &PointState,
        target: [f64; 3],
        maximum_increment_m: f64,
    ) -> Result<PointState> {
        self.integrate_partition(
            old,
            target,
            self.partition(old, target, maximum_increment_m)?,
        )
    }
    pub fn tangent(
        &self,
        old: &PointState,
        target: [f64; 3],
        count: usize,
        h: f64,
    ) -> Result<Matrix3<f64>> {
        if !h.is_finite() || h <= 0. {
            return Err(Error::Input("invalid JOINT tangent step".into()));
        }
        let mut d = Matrix3::zeros();
        for a in 0..3 {
            let mut plus = target;
            let mut minus = target;
            plus[a] += h;
            minus[a] -= h;
            let p = self.integrate_partition(old, plus, count)?;
            let m = self.integrate_partition(old, minus, count)?;
            d.set_column(
                a,
                &((Vector3::from(p.traction_kpa) - Vector3::from(m.traction_kpa)) / (2. * h)),
            );
        }
        Ok(d)
    }
}
