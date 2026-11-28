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

/// With a closed affine gap and straight tangential jump, the elastic trial
/// norm minus Coulomb capacity is convex. Starting on the cone and ending
/// inside it means unloading begins at the beginning of this chord.
pub(super) fn unloads_from_start(old: &PointState, new: &PointState, mat: &Material) -> bool {
    if old.status != Status::Slip || new.status != Status::Stick {
        return false;
    }
    let n = Vector3::from(old.normal);
    let end_gap = new.elastic_jump().dot(&n);
    if old.elastic_jump().dot(&n) > 0. || end_gap > 0. {
        return false;
    }
    let tau = Vector3::from(old.traction_kpa);
    let change = Vector3::from(new.jump_m) - Vector3::from(old.jump_m);
    let trial =
        tau - n * tau.dot(&n) + mat.tangential_stiffness_kpa_m * (change - n * change.dot(&n));
    trial.norm() <= -mat.friction_coefficient * mat.normal_stiffness_kpa_m * end_gap
}

/// At opening the tangent can change discontinuously. Use the pair of equilibria
/// on the same side if its extrapolated zero lies in the sign-changing half.
pub(super) fn opening_fraction(gaps: [f64; 3], half: usize) -> Option<f64> {
    let [a, b, c] = gaps;
    if gaps[half] > 0. || gaps[half + 1] <= 0. {
        return None;
    }
    let root = if half == 0 {
        1. - b / (c - b)
    } else {
        -b / (b - a)
    };
    (0. < root && root < 1.).then_some(root)
}

/// Root on one sign-changing half of the quadratic through the three global
/// equilibria. The proposed load coordinate is subsequently solved again.
pub(super) fn quadratic_gap_fraction(gaps: [f64; 3], half: usize) -> Option<f64> {
    let [a, b, c] = gaps;
    if (gaps[half] > 0.) == (gaps[half + 1] > 0.) {
        return None;
    }
    let mut lo = 0.;
    let mut hi = 1.;
    for _ in 0..48 {
        let s = (lo + hi) * 0.5;
        let t = half as f64 + s;
        let value = a + t * (b - a) + 0.5 * t * (t - 1.) * (a - 2. * b + c);
        if (value > 0.) == (gaps[half] > 0.) {
            lo = s;
        } else {
            hi = s;
        }
    }
    Some((lo + hi) * 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn branch_extrapolation_resolves_the_opening_kink() {
        // Closed slope 1, open slope 3, a zero at t = 0.75 or t = 1.25.
        assert_eq!(opening_fraction([-0.75, 0.75, 3.75], 0), Some(0.75));
        assert_eq!(opening_fraction([-1.25, -0.25, 2.25], 1), Some(0.25));
        assert_eq!(opening_fraction([-1., 1., 1.], 0), None);
    }

    #[test]
    fn gap_crossing_selects_the_root_in_the_requested_half() {
        // (t - 0.25)(t - 1.5) sampled at t = 0, 1, 2.
        let gaps = [0.375, -0.375, 0.875];
        assert!((quadratic_gap_fraction(gaps, 0).unwrap() - 0.25).abs() < 1e-14);
        assert!((quadratic_gap_fraction(gaps, 1).unwrap() - 0.5).abs() < 1e-14);
        assert_eq!(quadratic_gap_fraction([1., 2., 3.], 0), None);
    }

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
