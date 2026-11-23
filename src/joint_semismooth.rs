//! Active-branch derivative of the unilateral Coulomb projection along a fixed partition.
use super::*;
impl Material {
    pub fn semismooth_tangent(
        &self,
        old: &PointState,
        target: [f64; 3],
        count: usize,
    ) -> Result<Matrix3<f64>> {
        if count == 0 || count > 50000 {
            return Err(Error::Input("invalid JOINT active-set partition".into()));
        }
        let n = Vector3::from(old.normal);
        let plane = Matrix3::identity() - n * n.transpose();
        let start = Vector3::from(old.jump_m);
        let full = Vector3::from(target) - start;
        let dt = plane * full / count as f64;
        let dt_jac = plane / count as f64;
        let mut tau = plane * Vector3::from(old.traction_kpa);
        let mut jac = Matrix3::zeros();
        let mut g0 = n.dot(&(start - old.elastic_offset()));
        let mut grad0 = Vector3::zeros();
        for j in 1..=count {
            let t = j as f64 / count as f64;
            let end = if j == count {
                Vector3::from(target)
            } else {
                start + t * full
            };
            let g1 = n.dot(&(end - old.elastic_offset()));
            let grad1 = t * n;
            if g1 > 0. {
                tau.fill(0.);
                jac.fill(0.);
            } else {
                let mut fraction = 1.;
                let mut df = Vector3::zeros();
                if g0 > 0. {
                    // Only the closed part of a crossing contributes tangential slip.
                    tau.fill(0.);
                    jac.fill(0.);
                    let dg = g1 - g0;
                    fraction = g1 / dg;
                    df = (g1 * grad0 - g0 * grad1) / (dg * dg);
                }
                let trial = tau + self.tangential_stiffness_kpa_m * fraction * dt;
                let trial_jac = jac
                    + self.tangential_stiffness_kpa_m * (fraction * dt_jac + dt * df.transpose());
                let cap = -self.friction_coefficient * self.normal_stiffness_kpa_m * g1;
                let norm = trial.norm();
                if norm > cap && norm > 0. {
                    let dir = trial / norm;
                    let grad_cap = -self.friction_coefficient * self.normal_stiffness_kpa_m * grad1;
                    tau = cap * dir;
                    jac = (cap / norm) * (Matrix3::identity() - dir * dir.transpose()) * trial_jac
                        + dir * grad_cap.transpose();
                } else {
                    tau = trial;
                    jac = trial_jac;
                }
            }
            g0 = g1;
            grad0 = grad1;
        }
        if g0 <= 0. {
            jac += self.normal_stiffness_kpa_m * n * n.transpose();
        }
        Ok(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_sticking_sliding_and_reclosing_branches_have_their_own_jacobians() {
        let mat = Material {
            normal_stiffness_kpa_m: 8.,
            tangential_stiffness_kpa_m: 4.,
            friction_coefficient: 0.5,
        };
        let old = PointState::initial([0., 0., 1.]).unwrap();
        assert_eq!(
            mat.semismooth_tangent(&old, [0., 0., 0.25], 4).unwrap(),
            Matrix3::zeros()
        );
        assert_eq!(
            mat.semismooth_tangent(&old, [0.125, 0., -1.], 4).unwrap(),
            Matrix3::from_diagonal(&Vector3::new(4., 4., 8.))
        );
        assert_eq!(
            mat.semismooth_tangent(&old, [1., 0., -0.25], 1).unwrap(),
            Matrix3::new(0., 0., -4., 0., 1., 0., 0., 0., 8.)
        );
        let open = mat.integrate_partition(&old, [0., 0., 0.25], 1).unwrap();
        assert_eq!(
            mat.semismooth_tangent(&open, [0.25, 0., -0.25], 1).unwrap(),
            Matrix3::new(2., 0., -1., 0., 2., 0., 0., 0., 8.)
        );
    }
}
