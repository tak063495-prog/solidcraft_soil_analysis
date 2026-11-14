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
pub struct BirthReference {
    pub jump_m: [f64; 3],
    pub initial_traction_kpa: [f64; 3],
    pub elastic_offset_m: [f64; 3],
    pub initial_energy_kpa_m: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointState {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub birth: Option<BirthReference>,
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
            birth: None,
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
    pub fn elastic_jump(&self) -> Vector3<f64> {
        Vector3::from(self.jump_m) - self.elastic_offset()
    }
    fn elastic_offset(&self) -> Vector3<f64> {
        self.birth
            .as_ref()
            .map_or(Vector3::zeros(), |b| Vector3::from(b.elastic_offset_m))
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
        if let Some(b) = &mut s.birth {
            b.jump_m = (r * Vector3::from(b.jump_m)).into();
            b.initial_traction_kpa = (r * Vector3::from(b.initial_traction_kpa)).into();
            b.elastic_offset_m = (r * Vector3::from(b.elastic_offset_m)).into();
        }
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
    /// Fresh generation. Initial effective traction is admitted, never clipped.
    pub fn birth(
        &self,
        normal: [f64; 3],
        jump: [f64; 3],
        traction: [f64; 3],
    ) -> Result<PointState> {
        self.validate()?;
        let mut s = PointState::initial(normal)?;
        let n = Vector3::from(s.normal);
        let q = Vector3::from(traction);
        let qn = n.dot(&q);
        let shear = q - qn * n;
        if !jump.iter().chain(&traction).all(|v| v.is_finite())
            || qn > 0.
            || shear.norm() > self.friction_coefficient * (-qn)
        {
            return Err(Error::Input(
                "JOINT BIRTH initial effective force violates no-tension/Coulomb cone".into(),
            ));
        }
        let elastic =
            qn * n / self.normal_stiffness_kpa_m + shear / self.tangential_stiffness_kpa_m;
        s.jump_m = jump;
        s.traction_kpa = traction;
        s.birth = Some(BirthReference {
            jump_m: jump,
            initial_traction_kpa: traction,
            elastic_offset_m: (Vector3::from(jump) - elastic).into(),
            initial_energy_kpa_m: self.energy(&s),
        });
        self.validate_state(&s)?;
        Ok(s)
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
        let jump = s.elastic_jump();
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
        let initial_energy = if let Some(b) = &s.birth {
            let q0 = Vector3::from(b.initial_traction_kpa);
            let qn0 = n.dot(&q0);
            let qt0 = q0 - qn0 * n;
            let e0 = qn0 * qn0 / (2. * self.normal_stiffness_kpa_m)
                + qt0.norm_squared() / (2. * self.tangential_stiffness_kpa_m);
            let elastic0 =
                qn0 * n / self.normal_stiffness_kpa_m + qt0 / self.tangential_stiffness_kpa_m;
            if !b
                .jump_m
                .iter()
                .chain(&b.initial_traction_kpa)
                .chain(&b.elastic_offset_m)
                .chain([&b.initial_energy_kpa_m])
                .all(|v| v.is_finite())
                || qn0 > 0.
                || qt0.norm() > self.friction_coefficient * (-qn0)
                || (Vector3::from(b.jump_m) - elastic0 - Vector3::from(b.elastic_offset_m)).norm()
                    > 1e-10
                || (b.initial_energy_kpa_m - e0).abs() > 1e-10 * (1. + e0)
            {
                return Err(Error::Input(
                    "JOINT BIRTH initial reference/energy is inconsistent".into(),
                ));
            }
            e0
        } else {
            0.
        };
        let work = self.energy(s) + s.plastic_dissipation_kpa_m + s.numerical_dissipation_kpa_m;
        if (work - initial_energy - s.endpoint_work_kpa_m).abs()
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
        let normal = self.normal_stiffness_kpa_m * (target - old.elastic_offset()).dot(&n).min(0.);
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
        let offset = old.elastic_offset();
        let g0 = n.dot(&(before - offset));
        let g1 = n.dot(&(target - offset));
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
        let crossing = crossing - n * (crossing - offset).dot(&n);
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
