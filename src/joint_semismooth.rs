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
        if self.friction_coefficient == 0. {
            // The frictionless law is identically zero in both tangential directions.
            let gap = n.dot(&(Vector3::from(target) - old.elastic_offset()));
            return Ok(if gap <= 0. {
                self.normal_stiffness_kpa_m * n * n.transpose()
            } else {
                Matrix3::zeros()
            });
        }
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
            } else if g0 > 0. && g1 == 0. {
                // Both the trial shear and the Coulomb cap vanish at recontact.
                // Compare their slopes per unit penetration on the closed side.
                tau.fill(0.);
                let trial_increment = self.tangential_stiffness_kpa_m * dt;
                let trial_norm = trial_increment.norm();
                let cap_slope = self.friction_coefficient * self.normal_stiffness_kpa_m;
                let closed_slope = if trial_norm > cap_slope * g0 {
                    cap_slope * (trial_increment / trial_norm)
                } else {
                    trial_increment / g0
                };
                jac = -closed_slope * grad1.transpose();
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

    fn traction(mat: &Material, old: &PointState, target: [f64; 3], count: usize) -> Vector3<f64> {
        mat.integrate_partition(old, target, count)
            .unwrap()
            .traction_kpa
            .into()
    }

    fn assert_exact_recontact(shear: f64, expected_x: f64, status: Status) {
        let mat = Material {
            normal_stiffness_kpa_m: 8.,
            tangential_stiffness_kpa_m: 4.,
            friction_coefficient: 0.5,
        };
        let initial = PointState::initial([0., 0., 1.]).unwrap();
        let open = mat
            .integrate_partition(&initial, [0., 0., 0.00025], 1)
            .unwrap();
        assert_eq!(open.status, Status::Open);
        let h = 1e-11;
        let target = [shear, 0., 0.];
        let expected = Vector3::new(expected_x, 0., 8.);
        for count in [1, 4] {
            let at = traction(&mat, &open, target, count);
            let inside = mat
                .integrate_partition(&open, [shear, 0., -h], count)
                .unwrap();
            let outside = mat
                .integrate_partition(&open, [shear, 0., h], count)
                .unwrap();
            assert_eq!(inside.status, status);
            assert_eq!(outside.status, Status::Open);
            assert_eq!(at, Vector3::zeros());
            let closed = (at - Vector3::from(inside.traction_kpa)) / h;
            let opened = (Vector3::from(outside.traction_kpa) - at) / h;
            assert!((closed - expected).norm() < 1e-6, "{closed:?}");
            assert_eq!(opened, Vector3::zeros());

            let actual = mat.semismooth_tangent(&open, target, count).unwrap();
            assert_eq!(actual.column(0), Vector3::zeros());
            assert_eq!(actual.column(1), Vector3::zeros());
            assert!((actual.column(2) - expected).norm() < 1e-12, "{actual:?}");
            assert!((actual.column(2) - closed).norm() < 1e-6);
        }
    }

    #[test]
    fn exact_recontact_uses_sliding_closed_limit() {
        assert_exact_recontact(0.0005, -4., Status::Slip);
    }

    #[test]
    fn exact_recontact_uses_sticking_closed_limit() {
        assert_exact_recontact(0.000125, -2., Status::Stick);
    }

    #[test]
    fn frictionless_opt_in_has_only_the_normal_derivative() {
        let initial = PointState::initial([0., 0., 1.]).unwrap();
        let mut mat = Material {
            normal_stiffness_kpa_m: 8.,
            tangential_stiffness_kpa_m: 4.,
            friction_coefficient: 0.,
        };
        let open = mat
            .integrate_partition(&initial, [0., 0., 0.25], 1)
            .unwrap();
        let closed = Matrix3::from_diagonal(&Vector3::new(0., 0., 8.));
        for count in [1, 4] {
            for shear in [0., 1e-6] {
                let target = [shear, 0., -1.];
                assert!(
                    (mat.tangent(&initial, target, count, 1e-6).unwrap() - closed).norm() < 1e-8
                );
                assert_eq!(
                    mat.semismooth_tangent(&initial, target, count).unwrap(),
                    closed
                );
                assert_eq!(
                    mat.semismooth_tangent(&open, [shear, 0., 0.], count)
                        .unwrap(),
                    closed
                );
                assert_eq!(
                    mat.semismooth_tangent(&open, [shear, 0., 0.125], count)
                        .unwrap(),
                    Matrix3::zeros()
                );
            }
        }
        mat.friction_coefficient = 0.5;
        for count in [1, 4] {
            assert_eq!(
                mat.semismooth_tangent(&initial, [0., 0., -1.], count)
                    .unwrap(),
                Matrix3::from_diagonal(&Vector3::new(4., 4., 8.))
            );
        }
    }

    #[test]
    fn smooth_recontact_matches_the_existing_centered_tangent() {
        let mat = Material {
            normal_stiffness_kpa_m: 8.,
            tangential_stiffness_kpa_m: 4.,
            friction_coefficient: 0.5,
        };
        let initial = PointState::initial([0., 0., 1.]).unwrap();
        let open = mat
            .integrate_partition(&initial, [0., 0., 0.25], 1)
            .unwrap();
        for shear in [0.125, 0.75] {
            let target = [shear, 0., -0.3125];
            for count in [1, 4] {
                let reference = mat.tangent(&open, target, count, 1e-7).unwrap();
                let actual = mat.semismooth_tangent(&open, target, count).unwrap();
                assert!(
                    (actual - reference).norm() < 1e-7,
                    "{actual:?} {reference:?}"
                );
            }
        }
    }
}
