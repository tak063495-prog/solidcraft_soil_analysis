//! Constant-G MCC with an exactly integrated elastic path up to first yield.
use super::*;

impl Parameters {
    pub fn integrate_with_elastic_split(
        &self,
        committed: &State,
        strain: [f64; 6],
        max_substep: f64,
        shear_kpa: f64,
    ) -> Result<Update> {
        self.integrate_steps(
            committed,
            strain,
            max_substep,
            Some(shear_kpa),
            ElasticStep::Split,
        )
    }

    /// Preserve the outer strain grid and the exact elastic leg. Bound the
    /// strain one-norm only on the remainder integrated by the plastic return.
    pub fn integrate_with_plastic_refinement(
        &self,
        committed: &State,
        strain: [f64; 6],
        max_substep: f64,
        shear_kpa: f64,
        plastic_width: f64,
    ) -> Result<Update> {
        if !(1e-7..=1e-3).contains(&plastic_width) {
            return Err(Error::Input(
                "MCC plastic width must be in [1e-7,1e-3]".into(),
            ));
        }
        self.integrate_steps(
            committed,
            strain,
            max_substep,
            Some(shear_kpa),
            ElasticStep::RefinePlastic(plastic_width),
        )
    }

    pub(super) fn elastic_split_substep(
        &self,
        old: &State,
        de: [f64; 6],
        g: f64,
        plastic_width: Option<f64>,
    ) -> Result<(State, usize, usize)> {
        let ev = de[..3].iter().sum::<f64>();
        let p = invariants(&old.effective_stress_compression_positive_kpa).0;
        let v = 1. + old.void_ratio;
        let predictor = |fraction: f64| {
            // v1 = v0 exp(-dev), and kappa ln(p1/p0) = v0 - v1.
            let dv = v * (-fraction * ev).exp_m1();
            let dp = p * (-dv / self.kappa).exp_m1();
            let mut state = old.clone();
            state.void_ratio += dv;
            for (i, s) in state
                .effective_stress_compression_positive_kpa
                .iter_mut()
                .enumerate()
            {
                *s += if i < 3 {
                    dp + fraction * (2. * g) * (de[i] - ev / 3.)
                } else {
                    fraction * g * de[i]
                };
            }
            state
        };
        let trial = predictor(1.);
        let yield_value = self.normalized_yield(&trial);
        if !yield_value.is_finite() || trial.void_ratio <= 0. {
            return Err(Error::Material(
                "invalid analytic MCC elastic predictor".into(),
            ));
        }
        if yield_value <= 1e-11 {
            return Ok((trial, 0, 1));
        }
        let plastic_leg = |start: &State, strain: [f64; 6]| -> Result<(State, usize, usize)> {
            let width = strain.iter().map(|x| x.abs()).sum::<f64>();
            let steps = plastic_width.map_or(1, |h| (width / h).ceil().max(1.) as usize);
            let increment = strain.map(|x| x / steps as f64);
            let mut state = start.clone();
            let mut plastic = 0;
            for _ in 0..steps {
                let (next, yielded) = self.substep_with_shear(&state, increment, Some(g))?;
                state = next;
                plastic += usize::from(yielded);
            }
            Ok((state, plastic, steps))
        };
        // A parent already on the loading surface has no initial elastic leg.
        // On reversal a step can first enter the interior and then yield again.
        let mut lo = 0.;
        let mut hi = 1.;
        if self.normalized_yield(old) >= -1e-11 {
            if self.outward_elastic_direction(old, de, g) {
                return plastic_leg(old, de);
            }
            let mut probe = 0.5;
            while probe > f64::EPSILON && self.normalized_yield(&predictor(probe)) >= 0. {
                probe *= 0.5;
            }
            if self.normalized_yield(&predictor(probe)) >= 0. {
                return plastic_leg(old, de);
            }
            lo = probe;
        }
        for _ in 0..52 {
            let mid = (lo + hi) / 2.;
            if self.normalized_yield(&predictor(mid)) <= 0. {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let boundary = predictor(lo);
        plastic_leg(&boundary, de.map(|x| (1. - lo) * x))
            .map(|(state, plastic, count)| (state, plastic, count + 1))
    }

    fn outward_elastic_direction(&self, old: &State, de: [f64; 6], g: f64) -> bool {
        let s = &old.effective_stress_compression_positive_kpa;
        let p = invariants(s).0;
        let ev: f64 = de[..3].iter().sum();
        let dev_dot_de: f64 = (0..3)
            .map(|i| (s[i] - p) * de[i] + s[i + 3] * de[i + 3])
            .sum();
        let dp = (1. + old.void_ratio) * p * ev / self.kappa;
        6. * g * dev_dot_de
            + self.critical_state_m.powi(2) * (2. * p - old.preconsolidation_pressure_kpa) * dp
            >= 0.
    }
}

#[cfg(test)]
mod regression {
    use super::super::tests::{isotropic, parameters};
    use super::*;

    #[test]
    fn constant_shear_elastic_pressure_and_crossing() -> Result<()> {
        let law = parameters();
        let start = isotropic(5., 40.)?;
        let strain = [0.0002, 0.0002, 0.0002, 0., 0., 0.];
        let u = law.integrate_with_elastic_split(&start, strain, 0.001, 250.)?;
        let (p, _) = invariants(&u.state.effective_stress_compression_positive_kpa);
        assert_eq!((u.substeps, u.plastic_substeps), (1, 0));
        assert!((law.kappa * (p / 5.).ln() + u.state.void_ratio - start.void_ratio).abs() < 1e-14);

        let near_yield = isotropic(5., 5.00005)?;
        let shear = [0., 0., 0., 0.0002, 0., 0.];
        let u = law.integrate_with_elastic_split(&near_yield, shear, 0.001, 250.)?;
        assert_eq!((u.substeps, u.plastic_substeps), (2, 1));
        assert!(law.normalized_yield(&u.state).abs() < 1e-10);
        // For pure shear p is constant. The closed-form shear stress at yield
        // gives an independent parent for the identical backward Euler return.
        let tau_y = law.critical_state_m
            * (5. * (near_yield.preconsolidation_pressure_kpa - 5.) / 3.).sqrt();
        let mut boundary = near_yield.clone();
        boundary.effective_stress_compression_positive_kpa[3] = tau_y;
        let remainder = [0., 0., 0., shear[3] - tau_y / 250., 0., 0.];
        let manual = law.integrate_with_shear(&boundary, remainder, 0.001, Some(250.))?;
        for (split, exact) in u
            .state
            .effective_stress_compression_positive_kpa
            .iter()
            .zip(manual.state.effective_stress_compression_positive_kpa)
        {
            assert!((split - exact).abs() < 1e-13);
        }
        Ok(())
    }
}
