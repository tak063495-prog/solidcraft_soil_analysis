//! Constant-G MCC with an exactly integrated elastic path up to first yield.
use super::*;

const YIELD_TOLERANCE: f64 = 1e-11;

// The normalized yield function along the elastic path is a quadratic in the
// strain fraction (the deviatoric term) plus P(P - 1), where P = p / pc.
struct ElasticYield {
    pressure_ratio: f64,
    specific_volume: f64,
    kappa: f64,
    volumetric_strain: f64,
    deviatoric_coefficients: [f64; 3],
}

#[derive(Clone, Copy)]
struct YieldPoint {
    fraction: f64,
    pressure_ratio: f64,
    log_pressure_slope: f64,
    value: f64,
    slope: f64,
}

impl ElasticYield {
    fn new(law: &Parameters, old: &State, de: [f64; 6], g: f64, ev: f64, p: f64) -> Self {
        let pc = old.preconsolidation_pressure_kpa;
        let mut path = Self {
            pressure_ratio: p / pc,
            specific_volume: 1. + old.void_ratio,
            kappa: law.kappa,
            volumetric_strain: ev,
            deviatoric_coefficients: [0.; 3],
        };
        for (i, stress) in old
            .effective_stress_compression_positive_kpa
            .into_iter()
            .enumerate()
        {
            let (start, change, weight) = if i < 3 {
                (stress - p, 2. * g * (de[i] - ev / 3.), 1.5)
            } else {
                (stress, g * de[i], 3.)
            };
            let a = start / pc / law.critical_state_m;
            let b = change / pc / law.critical_state_m;
            path.deviatoric_coefficients[0] += weight * a * a;
            path.deviatoric_coefficients[1] += 2. * weight * a * b;
            path.deviatoric_coefficients[2] += weight * b * b;
        }
        path
    }

    fn at(&self, fraction: f64) -> YieldPoint {
        let volume_exponent = -fraction * self.volumetric_strain;
        let pressure_ratio = self.pressure_ratio
            * (-self.specific_volume * volume_exponent.exp_m1() / self.kappa).exp();
        let log_pressure_slope =
            self.volumetric_strain * self.specific_volume * volume_exponent.exp() / self.kappa;
        let [constant, linear, quadratic] = self.deviatoric_coefficients;
        YieldPoint {
            fraction,
            pressure_ratio,
            log_pressure_slope,
            value: constant
                + fraction * (linear + fraction * quadratic)
                + pressure_ratio * (pressure_ratio - 1.),
            slope: linear
                + 2. * fraction * quadratic
                + (2. * pressure_ratio - 1.) * pressure_ratio * log_pressure_slope,
        }
    }

    fn curvature_bounds(&self, left: YieldPoint, right: YieldPoint) -> (f64, f64) {
        // P and |d ln(P)/dt| are monotone. In F'', the pressure contribution is
        // P * W * ((4P - 1) * W - (2P - 1) * |ev|), where W = |d ln(P)/dt|.
        // The bilinear factor takes its extrema at corners of the enclosing
        // rectangle. P * W is nonnegative, so bound the product by its sign.
        let pressures = [left.pressure_ratio, right.pressure_ratio];
        let rates = [
            left.log_pressure_slope.abs(),
            right.log_pressure_slope.abs(),
        ];
        let ev = self.volumetric_strain.abs();
        let mut lower = f64::INFINITY;
        let mut upper = f64::NEG_INFINITY;
        for p in pressures {
            for w in rates {
                let factor = (4. * p - 1.) * w - (2. * p - 1.) * ev;
                lower = lower.min(factor);
                upper = upper.max(factor);
            }
        }
        let min_product = pressures[0].min(pressures[1]) * rates[0].min(rates[1]);
        let max_product = pressures[0].max(pressures[1]) * rates[0].max(rates[1]);
        let shear = 2. * self.deviatoric_coefficients[2];
        (
            shear + lower.min(0.) * max_product + lower.max(0.) * min_product,
            shear + upper.min(0.) * min_product + upper.max(0.) * max_product,
        )
    }

    // Return an inside-to-outside bracket on the first excursion above level.
    // Sign -1 searches for entry into the interior from a positive residual.
    fn first_above(
        &self,
        left: YieldPoint,
        right: YieldPoint,
        level: f64,
        sign: f64,
    ) -> Option<(YieldPoint, YieldPoint)> {
        let left_value = sign * left.value;
        let right_value = sign * right.value;
        if left_value > level {
            return Some((left, left));
        }
        let (mut lower, mut upper) = self.curvature_bounds(left, right);
        if sign < 0. {
            (lower, upper) = (-upper, -lower);
        }
        let width = right.fraction - left.fraction;
        // A lower bound on F'' bounds how far F can rise above its chord.
        if left_value.max(right_value) - lower.min(0.) * width * width / 8. <= level {
            return None;
        }
        let left_slope = sign * left.slope;
        if lower >= 0. || left_slope + lower.min(0.) * width >= 0. || width <= f64::EPSILON {
            return (right_value > level).then_some((left, right));
        }
        if upper <= 0. {
            // A concave interval has one maximum. Use its rising leg for
            // the bracket, even if the descending leg ends inside.
            let mut peak = right;
            if sign * right.slope < 0. {
                if left_slope <= 0. {
                    return None;
                }
                let mut rising = left;
                for _ in 0..52 {
                    let mid = self.at((rising.fraction + peak.fraction) / 2.);
                    if sign * mid.slope > 0. {
                        rising = mid;
                    } else {
                        peak = mid;
                    }
                }
            }
            return (sign * peak.value > level).then_some((left, peak));
        }
        if left_slope + upper.max(0.) * width <= 0. {
            return None;
        }
        let middle = self.at((left.fraction + right.fraction) / 2.);
        self.first_above(left, middle, level, sign)
            .or_else(|| self.first_above(middle, right, level, sign))
    }

    fn first_yield(&self, start: YieldPoint, end: YieldPoint) -> Option<(f64, f64)> {
        // Preserve the loading-surface convention for an outward endpoint.
        if end.value > YIELD_TOLERANCE && start.value >= -YIELD_TOLERANCE && start.slope >= 0. {
            return Some((0., end.fraction));
        }
        let mut inside = start;
        if start.value > 0. && start.slope < 0. {
            // The stored parent can be slightly exterior within the state
            // tolerance. Find the inward crossing before looking for reload.
            if let Some((_, inward)) = self.first_above(start, end, 0., -1.) {
                inside = inward;
            } else {
                return (end.value > YIELD_TOLERANCE).then_some((0., end.fraction));
            }
        }
        let (_, witness) = self.first_above(inside, end, YIELD_TOLERANCE, 1.)?;
        if inside.fraction == 0. && start.value >= 0. && start.slope >= 0. {
            return Some((0., witness.fraction));
        }
        let (mut lo, mut hi) = self.first_above(inside, witness, 0., 1.)?;
        for _ in 0..52 {
            let mid = self.at((lo.fraction + hi.fraction) / 2.);
            if mid.value <= 0. {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some((lo.fraction, witness.fraction))
    }
}

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
        let stress = old.effective_stress_compression_positive_kpa;
        let p = invariants(&stress).0;
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
        let path = ElasticYield::new(self, old, de, g, ev, p);
        let start = path.at(0.);
        let end = path.at(1.);
        let Some((first, witness)) = path.first_yield(start, end) else {
            return Ok((trial, 0, 1));
        };
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
        let boundary = predictor(first);
        let elastic_count = usize::from(first > 0.);
        let (state, plastic, count) = plastic_leg(&boundary, de.map(|x| (1. - first) * x))?;
        if plastic > 0 || end.value > YIELD_TOLERANCE || witness <= first {
            return Ok((state, plastic, count + elastic_count));
        }

        // The full remainder can also hide an excursion from the frozen
        // elastic trial used by the return. Advance to an exterior point
        // first, then continue from its updated state through the tail.
        let (peak, plastic, count) = plastic_leg(&boundary, de.map(|x| (witness - first) * x))?;
        if plastic == 0 {
            return Err(Error::Material(
                "MCC elastic crossing requires a smaller plastic increment".into(),
            ));
        }
        let (state, tail_plastic, tail_count) = plastic_leg(&peak, de.map(|x| (1. - witness) * x))?;
        Ok((
            state,
            plastic + tail_plastic,
            elastic_count + count + tail_count,
        ))
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
    fn interior_excursion(
        shear_offset: f64,
        duration: f64,
    ) -> Result<(Parameters, State, [f64; 6])> {
        let law = Parameters {
            kappa: 0.05,
            ..parameters()
        };
        let mut start = isotropic(1., 40.)?;
        start.void_ratio = 1.;
        start.effective_stress_compression_positive_kpa[3] = 4.326655757054095 + shear_offset;
        let strain = [
            0.00055 / 3.,
            0.00055 / 3.,
            0.00055 / 3.,
            0.00018638173592571168,
            0.,
            0.,
        ]
        .map(|x| x * duration);
        Ok((law, start, strain))
    }

    #[test]
    fn first_yield_between_elastic_endpoints() -> Result<()> {
        let (law, start, strain) = interior_excursion(0., 1.)?;
        let path = ElasticYield::new(&law, &start, strain, 250., 0.00055, 1.);
        let (begin, end) = (path.at(0.), path.at(1.));
        assert!(begin.value < 0. && end.value < 0. && path.at(0.5).value > 0.);
        let (first, _) = path.first_yield(begin, end).expect("interior yield");
        // First root from an independent 60-digit pressure/shear calculation.
        assert!((first - 0.0266576652546463).abs() < 1e-9);

        let update = |width| law.integrate_with_elastic_split(&start, strain, width, 250.);
        let coarse = update(0.001)?;
        assert_eq!((coarse.substeps, coarse.plastic_substeps), (2, 1));
        let refined = update(0.00001)?;
        let fine = update(0.0000001)?;
        let stress_error = |u: &Update| {
            u.state
                .effective_stress_compression_positive_kpa
                .iter()
                .zip(fine.state.effective_stress_compression_positive_kpa)
                .map(|(actual, reference)| (actual - reference).abs())
                .sum::<f64>()
        };
        let pc_error = |u: &Update| {
            (u.state.preconsolidation_pressure_kpa - fine.state.preconsolidation_pressure_kpa).abs()
        };
        let strain_error = |u: &Update| {
            u.state
                .plastic_strain_compression_positive
                .iter()
                .zip(fine.state.plastic_strain_compression_positive)
                .map(|(actual, reference)| (actual - reference).abs())
                .sum::<f64>()
        };
        assert!(stress_error(&refined) < stress_error(&coarse));
        assert!(pc_error(&refined) < pc_error(&coarse));
        assert!(strain_error(&refined) < strain_error(&coarse));
        Ok(())
    }

    #[test]
    fn narrow_excursion_can_miss_the_midpoint() -> Result<()> {
        let (law, start, strain) = interior_excursion(-0.000049, 1.3)?;
        let path = ElasticYield::new(&law, &start, strain, 250., 0.00055 * 1.3, 1.);
        let (begin, end) = (path.at(0.), path.at(1.));
        assert!(begin.value < 0. && end.value < 0. && path.at(0.5).value < 0.);
        let (first, _) = path
            .first_yield(begin, end)
            .expect("off-center interior yield");
        assert!((first - 0.3355505054253187).abs() < 1e-8);
        let u = law.integrate_with_elastic_split(&start, strain, 0.001, 250.)?;
        assert!(u.plastic_substeps > 0);
        Ok(())
    }

    #[test]
    fn unloading_crosses_the_opposite_side_before_return() -> Result<()> {
        let law = parameters();
        let mut start = isotropic(5., 5.00005)?;
        let tau =
            law.critical_state_m * (5. * (start.preconsolidation_pressure_kpa - 5.) / 3.).sqrt();
        start.effective_stress_compression_positive_kpa[3] = tau;
        let mut strain = [0.; 6];
        strain[3] = -3. * tau / 250.;
        let result = law.integrate_with_elastic_split(&start, strain, 0.001, 250.)?;

        start.effective_stress_compression_positive_kpa[3] = -tau;
        strain[3] = -tau / 250.;
        let expected = law.integrate_with_shear(&start, strain, 0.001, Some(250.))?;
        assert_eq!((result.substeps, result.plastic_substeps), (2, 1));
        for (actual, reference) in result
            .state
            .effective_stress_compression_positive_kpa
            .iter()
            .zip(expected.state.effective_stress_compression_positive_kpa)
        {
            assert!((actual - reference).abs() < 1e-13);
        }
        Ok(())
    }
}
