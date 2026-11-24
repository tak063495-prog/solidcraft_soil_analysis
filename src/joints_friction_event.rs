//! Locations on an equilibrated half-step chord. They propose load coordinates only.
use crate::joint_law::{Material, PointState, Status};
use nalgebra::Vector3;

pub(super) fn fraction(old: &PointState, new: &PointState, mat: &Material) -> Option<f64> {
    let n = Vector3::from(old.normal);
    let a = old.elastic_jump().dot(&n);
    let b = new.elastic_jump().dot(&n);
    if (a > 0.) != (b > 0.) {
        return Some(a / (a - b));
    }
    if a > 0. || old.status == new.status {
        return None;
    }
    if new.status == Status::Stick && old.status == Status::Slip {
        // Endpoint tractions lie on or inside the cone. Bracket the change of activity.
        return Some(0.5);
    }
    if new.status != Status::Slip || old.status != Status::Stick {
        return None;
    }
    let traction = Vector3::from(old.traction_kpa);
    let shear = traction - n * traction.dot(&n);
    let dj = Vector3::from(new.jump_m) - Vector3::from(old.jump_m);
    let dt = mat.tangential_stiffness_kpa_m * (dj - n * dj.dot(&n));
    let cap = mat.friction_coefficient * mat.normal_stiffness_kpa_m;
    let margin = |t: f64| (shear + dt * t).norm() + cap * (a + t * (b - a));
    if margin(0.) >= 0. || margin(1.) <= 0. {
        return None;
    }
    let (mut lo, mut hi) = (0., 1.);
    for _ in 0..48 {
        let mid = (lo + hi) * 0.5;
        if margin(mid) > 0. {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some((lo + hi) * 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn friction_crossing_uses_the_trial_traction_before_return() {
        let mat = Material {
            normal_stiffness_kpa_m: 8.,
            tangential_stiffness_kpa_m: 4.,
            friction_coefficient: 0.5,
        };
        let mut old = PointState::initial([0., 0., 1.]).unwrap();
        old.jump_m[2] = -0.001;
        old.traction_kpa[2] = -0.008;
        let mut end = old.clone();
        end.jump_m[0] = 0.004;
        end.traction_kpa[0] = 0.004;
        end.status = Status::Slip;
        assert!((fraction(&old, &end, &mat).unwrap() - 0.25).abs() < 1e-14);
        assert!(fraction(&old, &old, &mat).is_none());
        end.jump_m[2] = 0.003;
        end.status = Status::Open;
        assert_eq!(fraction(&old, &end, &mat), Some(0.25));
    }
}
