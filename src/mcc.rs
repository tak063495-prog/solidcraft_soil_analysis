//! Modified Cam-Clay material point. Compression-positive effective stress and
//! compression-positive engineering strain: xx,yy,zz,xy,yz,zx. No pressure clipping.
use crate::{Error, Result};
use nalgebra::{Matrix4, Vector4};
use serde::{Deserialize, Serialize};

#[path = "mcc_sensitivity.rs"]
mod sensitivity;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {
    pub critical_state_m: f64,
    pub lambda: f64,
    pub kappa: f64,
    pub poisson: f64,
    pub reference_pressure_kpa: f64,
    pub specific_volume_ncl_at_reference: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub effective_stress_compression_positive_kpa: [f64; 6],
    pub void_ratio: f64,
    pub preconsolidation_pressure_kpa: f64,
    pub plastic_strain_compression_positive: [f64; 6],
}

#[derive(Debug, Serialize)]
pub struct Update {
    pub state: State,
    pub substeps: usize,
    pub plastic_substeps: usize,
}

pub fn invariants(s: &[f64; 6]) -> (f64, f64) {
    let p = (s[0] + s[1] + s[2]) / 3.;
    let ss = (s[0] - p).powi(2)
        + (s[1] - p).powi(2)
        + (s[2] - p).powi(2)
        + 2. * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2));
    (p, (1.5 * ss).sqrt())
}

impl Parameters {
    /// Differentiate the converged local return and propagate sensitivities
    /// through the complete fixed partition. The primal update is unchanged.
    /// At a switching point this selects the branch taken by the primal update.
    pub fn implicit_stress_tangent_partition(
        &self,
        committed: &State,
        strain: [f64; 6],
        steps: usize,
    ) -> Result<(Update, crate::elastic::Tensor)> {
        sensitivity::integrate(self, committed, strain, steps)
    }
    /// Fixed partition for differentiating the complete increment map. In
    /// particular, finite-difference perturbations must not change ceil(nsteps).
    pub fn integrate_partition(
        &self,
        committed: &State,
        strain: [f64; 6],
        steps: usize,
    ) -> Result<Update> {
        self.validate_state(committed)?;
        if steps == 0 || steps > 100_000 || !strain.iter().all(|x| x.is_finite()) {
            return Err(Error::Input("invalid fixed MCC partition".into()));
        }
        let mut state = committed.clone();
        let mut plastic = 0;
        for _ in 0..steps {
            let (next, yielded) = self.substep(&state, strain.map(|x| x / steps as f64))?;
            state = next;
            plastic += usize::from(yielded);
        }
        self.validate_state(&state)?;
        Ok(Update {
            state,
            substeps: steps,
            plastic_substeps: plastic,
        })
    }

    /// Numerical algorithmic tangent of the *whole* material increment, without
    /// symmetrization. At a switching surface this is a central generalized
    /// derivative; global line search still accepts only the material residual.
    pub fn stress_tangent(
        &self,
        committed: &State,
        strain: [f64; 6],
        maximum_substep: f64,
        h: f64,
    ) -> Result<(Update, crate::elastic::Tensor)> {
        if !maximum_substep.is_finite()
            || !(1e-7..=1e-3).contains(&maximum_substep)
            || !h.is_finite()
            || !(1e-9..=1e-5).contains(&h)
        {
            return Err(Error::Input("invalid MCC tangent step controls".into()));
        }
        let update = self.integrate(committed, strain, maximum_substep)?;
        self.stress_tangent_partition(committed, strain, update.substeps, h)
    }

    pub fn stress_tangent_partition(
        &self,
        committed: &State,
        strain: [f64; 6],
        steps: usize,
        h: f64,
    ) -> Result<(Update, crate::elastic::Tensor)> {
        if !h.is_finite() || !(1e-9..=1e-5).contains(&h) {
            return Err(Error::Input("invalid MCC tangent perturbation".into()));
        }
        let update = self.integrate_partition(committed, strain, steps)?;
        let mut tangent = crate::elastic::Tensor::zeros();
        for j in 0..6 {
            let mut plus = strain;
            plus[j] += h;
            let mut minus = strain;
            minus[j] -= h;
            let a = self
                .integrate_partition(committed, plus, update.substeps)?
                .state;
            let b = self
                .integrate_partition(committed, minus, update.substeps)?
                .state;
            for i in 0..6 {
                tangent[(i, j)] = (a.effective_stress_compression_positive_kpa[i]
                    - b.effective_stress_compression_positive_kpa[i])
                    / (2. * h);
            }
        }
        if !tangent.iter().all(|x| x.is_finite()) {
            return Err(Error::Material("non-finite MCC tangent".into()));
        }
        Ok((update, tangent))
    }

    pub fn validate(&self) -> Result<()> {
        if ![
            self.critical_state_m,
            self.lambda,
            self.kappa,
            self.poisson,
            self.reference_pressure_kpa,
            self.specific_volume_ncl_at_reference,
        ]
        .iter()
        .all(|x| x.is_finite())
            || self.critical_state_m <= 0.
            || self.kappa <= 0.
            || self.lambda <= self.kappa
            || self.poisson <= -1.
            || self.poisson >= 0.5
            || self.reference_pressure_kpa <= 0.
            || self.specific_volume_ncl_at_reference <= 1.
        {
            return Err(Error::Input(
                "MCC requires finite M>0, lambda>kappa>0, -1<nu<0.5, reference p>0 and v>1".into(),
            ));
        }
        Ok(())
    }

    pub fn reference_void_ratio(&self, p: f64, pc: f64) -> Result<f64> {
        self.validate()?;
        if !p.is_finite() || !pc.is_finite() || p <= 0. || pc <= 0. {
            return Err(Error::Input(
                "MCC initial mean pressure and pc must be positive".into(),
            ));
        }
        let e = self.specific_volume_ncl_at_reference
            - self.lambda * (pc / self.reference_pressure_kpa).ln()
            + self.kappa * (pc / p).ln()
            - 1.;
        if !e.is_finite() || e <= 0. {
            return Err(Error::Input(
                "MCC reference state implies non-positive void ratio".into(),
            ));
        }
        Ok(e)
    }

    pub fn moduli(&self, state: &State) -> Result<(f64, f64)> {
        let p = invariants(&state.effective_stress_compression_positive_kpa).0;
        let k = (1. + state.void_ratio) * p / self.kappa;
        let g = 3. * k * (1. - 2. * self.poisson) / (2. * (1. + self.poisson));
        if !k.is_finite() || !g.is_finite() || k <= 0. || g <= 0. {
            return Err(Error::Material("invalid MCC elastic moduli".into()));
        }
        Ok((k, g))
    }

    pub fn normalized_yield(&self, state: &State) -> f64 {
        let (p, q) = invariants(&state.effective_stress_compression_positive_kpa);
        let pc = state.preconsolidation_pressure_kpa;
        ((q / pc).powi(2) + self.critical_state_m.powi(2) * (p / pc) * (p / pc - 1.))
            / self.critical_state_m.powi(2)
    }

    pub fn validate_state(&self, state: &State) -> Result<()> {
        self.validate()?;
        let (p, q) = invariants(&state.effective_stress_compression_positive_kpa);
        let yield_value = self.normalized_yield(state);
        if !state
            .effective_stress_compression_positive_kpa
            .iter()
            .chain(&state.plastic_strain_compression_positive)
            .all(|x| x.is_finite())
            || !state.void_ratio.is_finite()
            || state.void_ratio <= 0.
            || !state.preconsolidation_pressure_kpa.is_finite()
            || state.preconsolidation_pressure_kpa <= 0.
            || !p.is_finite()
            || !q.is_finite()
            || p <= 0.
            || !yield_value.is_finite()
            || yield_value > 1e-8
        {
            return Err(Error::Input(
                "inadmissible MCC state: require p>0, pc>0, e>0 and stress inside yield surface"
                    .into(),
            ));
        }
        self.moduli(state)?;
        Ok(())
    }

    /// Immutable trial update. A failed call cannot alter the committed state.
    /// Substeps freeze K,G and hardening coefficient at each substep's start;
    /// backward Euler return is solved in p,q,ln(pc),plastic multiplier.
    pub fn integrate(
        &self,
        committed: &State,
        strain: [f64; 6],
        max_substep: f64,
    ) -> Result<Update> {
        self.validate_state(committed)?;
        if !strain.iter().all(|x| x.is_finite())
            || !max_substep.is_finite()
            || !(1e-7..=1e-3).contains(&max_substep)
        {
            return Err(Error::Input(
                "finite strain and MCC max_substep in [1e-7,1e-3] required".into(),
            ));
        }
        let magnitude = strain.iter().map(|x| x.abs()).sum::<f64>();
        let required = (magnitude / max_substep).ceil().max(1.);
        if required > 100_000. {
            return Err(Error::Material("MCC substep budget exceeded".into()));
        }
        let mut steps = required as usize;
        let mut last = String::new();
        while steps <= 100_000 {
            let mut state = committed.clone();
            let mut plastic = 0;
            let increment = strain.map(|x| x / steps as f64);
            let mut failed = false;
            for _ in 0..steps {
                match self.substep(&state, increment) {
                    Ok((next, yielded)) => {
                        state = next;
                        plastic += usize::from(yielded);
                    }
                    Err(e) => {
                        failed = true;
                        last = e.to_string();
                        break;
                    }
                }
            }
            if !failed {
                self.validate_state(&state)?;
                return Ok(Update {
                    state,
                    substeps: steps,
                    plastic_substeps: plastic,
                });
            }
            steps *= 2;
        }
        Err(Error::Material(format!(
            "MCC return failed without committing state: {last}"
        )))
    }

    fn substep(&self, old: &State, de: [f64; 6]) -> Result<(State, bool)> {
        let (k, g) = self.moduli(old)?;
        let ev = de[0] + de[1] + de[2];
        let mut trial = old.effective_stress_compression_positive_kpa;
        for i in 0..3 {
            trial[i] += k * ev + 2. * g * (de[i] - ev / 3.);
        }
        for i in 3..6 {
            trial[i] += g * de[i];
        }
        let (pt, qt) = invariants(&trial);
        let mut next = old.clone();
        next.void_ratio = (1. + old.void_ratio) * (-ev).exp() - 1.;
        if pt <= 0.
            || !pt.is_finite()
            || !qt.is_finite()
            || !next.void_ratio.is_finite()
            || next.void_ratio <= 0.
        {
            return Err(Error::Material(
                "MCC tensile mean stress or invalid specific volume; no clipping".into(),
            ));
        }
        next.effective_stress_compression_positive_kpa = trial;
        if self.normalized_yield(&next) <= 1e-11 {
            return Ok((next, false));
        }
        let scale = old.preconsolidation_pressure_kpa.max(pt).max(qt);
        let p_trial = pt / scale;
        let q_trial = qt / scale;
        let c_old = old.preconsolidation_pressure_kpa / scale;
        let kn = k / scale;
        let gn = g / scale;
        let a = (1. + old.void_ratio) / (self.lambda - self.kappa);
        let m2 = self.critical_state_m.powi(2);
        let evaluate = |x: &Vector4<f64>| {
            let (p, q, c, l) = (x[0], x[1], x[2].exp(), x[3]);
            let h = m2 * (2. * p - c);
            Vector4::new(
                p - p_trial + kn * l * h,
                q - q_trial + 6. * gn * l * q,
                x[2] - c_old.ln() - a * l * h,
                q * q + m2 * p * (p - c),
            )
        };
        let mut x = Vector4::new(p_trial, q_trial, c_old.ln(), 0.);
        let mut success = false;
        for _ in 0..60 {
            let r = evaluate(&x);
            if r.norm() < 1e-11 {
                success = true;
                break;
            }
            let (p, q, c, l) = (x[0], x[1], x[2].exp(), x[3]);
            let h = m2 * (2. * p - c);
            let jac = Matrix4::from_row_slice(&[
                1. + 2. * kn * l * m2,
                0.,
                -kn * l * m2 * c,
                kn * h,
                0.,
                1. + 6. * gn * l,
                0.,
                6. * gn * q,
                -2. * a * l * m2,
                0.,
                1. + a * l * m2 * c,
                -a * h,
                m2 * (2. * p - c),
                2. * q,
                -m2 * p * c,
                0.,
            ]);
            let correction = jac
                .lu()
                .solve(&(-r))
                .ok_or_else(|| Error::Material("singular MCC local Jacobian".into()))?;
            let mut accepted = false;
            let mut alpha = 1.;
            for _ in 0..30 {
                let candidate = x + alpha * correction;
                if candidate.iter().all(|v| v.is_finite())
                    && candidate[0] > 0.
                    && candidate[1] >= 0.
                    && candidate[3] >= 0.
                    && candidate[2].abs() < 100.
                    && evaluate(&candidate).norm() < r.norm()
                {
                    x = candidate;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                return Err(Error::Material("MCC local line search failed".into()));
            }
        }
        if !success {
            return Err(Error::Material("MCC local iteration limit".into()));
        }
        let p = x[0] * scale;
        let ratio = if qt > 0. { x[1] * scale / qt } else { 0. };
        let mut stress = trial.map(|s| s * ratio);
        for i in 0..3 {
            stress[i] = p + (trial[i] - pt) * ratio;
        }
        next.effective_stress_compression_positive_kpa = stress;
        next.preconsolidation_pressure_kpa = x[2].exp() * scale;
        for i in 0..3 {
            next.plastic_strain_compression_positive[i] +=
                ((trial[i] - pt) - (stress[i] - p)) / (2. * g) + (pt - p) / (3. * k);
        }
        for i in 3..6 {
            next.plastic_strain_compression_positive[i] += (trial[i] - stress[i]) / g;
        }
        Ok((next, true))
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialPath {
    pub schema: String,
    pub parameters: Parameters,
    pub initial_state: State,
    pub strain_increments_compression_positive: Vec<[f64; 6]>,
    pub maximum_substep: f64,
}

#[derive(Debug, Serialize)]
pub struct PathOutput {
    pub schema: &'static str,
    pub parameters: Parameters,
    pub states: Vec<State>,
    pub substeps: Vec<usize>,
    pub plastic_substeps: Vec<usize>,
}

pub fn run_path(path: MaterialPath) -> Result<PathOutput> {
    if path.schema != "soil-fem-mcc-path/1"
        || path.strain_increments_compression_positive.len() > 10000
    {
        return Err(Error::Input(
            "invalid MCC path schema or path length".into(),
        ));
    }
    path.parameters.validate_state(&path.initial_state)?;
    let mut output = PathOutput {
        schema: "soil-fem-mcc-path-result/1",
        parameters: path.parameters.clone(),
        states: vec![path.initial_state],
        substeps: vec![],
        plastic_substeps: vec![],
    };
    for increment in path.strain_increments_compression_positive {
        let old = output
            .states
            .last()
            .ok_or_else(|| Error::Material("MCC path has no state".into()))?;
        let update = path
            .parameters
            .integrate(old, increment, path.maximum_substep)?;
        output.states.push(update.state);
        output.substeps.push(update.substeps);
        output.plastic_substeps.push(update.plastic_substeps);
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    pub fn parameters() -> Parameters {
        Parameters {
            critical_state_m: 1.2,
            lambda: 0.2,
            kappa: 0.04,
            poisson: 0.3,
            reference_pressure_kpa: 1.,
            specific_volume_ncl_at_reference: 3.5,
        }
    }
    pub fn isotropic(p: f64, pc: f64) -> Result<State> {
        Ok(State {
            effective_stress_compression_positive_kpa: [p, p, p, 0., 0., 0.],
            void_ratio: parameters().reference_void_ratio(p, pc)?,
            preconsolidation_pressure_kpa: pc,
            plastic_strain_compression_positive: [0.; 6],
        })
    }
    #[test]
    fn isotropic_compression_unloading_follows_lambda_kappa() -> Result<()> {
        let model = parameters();
        let initial = isotropic(100., 100.)?;
        let load = model.integrate(&initial, [0.01, 0.01, 0.01, 0., 0., 0.], 1e-5)?;
        let p = invariants(&load.state.effective_stress_compression_positive_kpa).0;
        let expected_v = 1. + initial.void_ratio - model.lambda * (p / 100.).ln();
        assert!((1. + load.state.void_ratio - expected_v).abs() < 3e-5);
        assert!((p / load.state.preconsolidation_pressure_kpa - 1.).abs() < 1e-8);
        assert!(load.state.plastic_strain_compression_positive[0] > 0.);
        let unload = model.integrate(&load.state, [-0.002, -0.002, -0.002, 0., 0., 0.], 1e-5)?;
        let p2 = invariants(&unload.state.effective_stress_compression_positive_kpa).0;
        assert_eq!(unload.plastic_substeps, 0);
        assert_eq!(
            unload.state.preconsolidation_pressure_kpa,
            load.state.preconsolidation_pressure_kpa
        );
        let expected = 1. + load.state.void_ratio - model.kappa * (p2 / p).ln();
        assert!((1. + unload.state.void_ratio - expected).abs() < 3e-5);
        Ok(())
    }
    #[test]
    fn pure_shear_is_rotation_equivalent_and_refines() -> Result<()> {
        let model = parameters();
        let old = isotropic(100., 200.)?;
        let shear = model.integrate(&old, [0., 0., 0., 0.04, 0., 0.], 2e-5)?;
        let rotated = model.integrate(&old, [0.02, -0.02, 0., 0., 0., 0.], 2e-5)?;
        let (p, q) = invariants(&shear.state.effective_stress_compression_positive_kpa);
        let (pr, qr) = invariants(&rotated.state.effective_stress_compression_positive_kpa);
        assert!((p - pr).abs() < 1e-8 && (q - qr).abs() < 1e-8);
        assert!(shear.plastic_substeps > 0);
        let fine = model.integrate(&old, [0., 0., 0., 0.04, 0., 0.], 1e-5)?;
        assert!(
            (invariants(&fine.state.effective_stress_compression_positive_kpa).1 - q).abs() / q
                < 1e-3
        );
        Ok(())
    }
    #[test]
    fn invalid_state_and_failed_trial_are_not_committed() -> Result<()> {
        let model = parameters();
        let old = isotropic(100., 200.)?;
        let bytes = serde_json::to_string(&old)?;
        assert!(
            model
                .integrate(&old, [-10., -10., -10., 0., 0., 0.], 1e-4)
                .is_err()
        );
        assert_eq!(bytes, serde_json::to_string(&old)?);
        let mut bad = old;
        bad.preconsolidation_pressure_kpa = 50.;
        assert!(model.validate_state(&bad).is_err());
        Ok(())
    }

    #[test]
    fn restart_retains_all_hardening_and_plastic_history() -> Result<()> {
        let model = parameters();
        let old = isotropic(100., 100.)?;
        let first = model
            .integrate(&old, [0.005, 0.005, 0.005, 0.01, 0., 0.], 2e-5)?
            .state;
        let restored: State = serde_json::from_str(&serde_json::to_string(&first)?)?;
        let de = [-0.001, -0.001, -0.001, -0.005, 0., 0.];
        let direct = model.integrate(&first, de, 2e-5)?.state;
        let restart = model.integrate(&restored, de, 2e-5)?.state;
        assert_eq!(
            serde_json::to_string(&direct)?,
            serde_json::to_string(&restart)?
        );
        assert!(restart.preconsolidation_pressure_kpa > old.preconsolidation_pressure_kpa);
        Ok(())
    }

    #[test]
    fn arbitrary_basis_rotation_and_plastic_work() -> Result<()> {
        use nalgebra::{Matrix3, Rotation3};
        let model = parameters();
        let mut old = isotropic(100., 150.)?;
        old.effective_stress_compression_positive_kpa = [130., 80., 90., 10., -5., 3.];
        let rotation = Rotation3::from_euler_angles(0.4, -0.2, 0.7).into_inner();
        let rotate = |s: [f64; 6], engineering: bool| {
            let factor = if engineering { 0.5 } else { 1. };
            let tensor = Matrix3::new(
                s[0],
                s[3] * factor,
                s[5] * factor,
                s[3] * factor,
                s[1],
                s[4] * factor,
                s[5] * factor,
                s[4] * factor,
                s[2],
            );
            let t = rotation * tensor * rotation.transpose();
            [
                t[(0, 0)],
                t[(1, 1)],
                t[(2, 2)],
                t[(0, 1)] / factor,
                t[(1, 2)] / factor,
                t[(2, 0)] / factor,
            ]
        };
        let de = [2e-5, -1e-5, -1e-5, 5e-5, -2e-5, 3e-5];
        let mut rotated = old.clone();
        rotated.effective_stress_compression_positive_kpa =
            rotate(old.effective_stress_compression_positive_kpa, false);
        let mut count = 0;
        for _ in 0..800 {
            let (next, plastic) = model.substep(&old, de)?;
            let (next_rotated, _) = model.substep(&rotated, rotate(de, true))?;
            let expected = rotate(next.effective_stress_compression_positive_kpa, false);
            assert!(
                expected
                    .iter()
                    .zip(next_rotated.effective_stress_compression_positive_kpa)
                    .all(|(a, b)| (a - b).abs() < 1e-8)
            );
            if plastic {
                count += 1;
                let work: f64 = (0..6)
                    .map(|i| {
                        next.effective_stress_compression_positive_kpa[i]
                            * (next.plastic_strain_compression_positive[i]
                                - old.plastic_strain_compression_positive[i])
                    })
                    .sum();
                assert!(work >= -1e-12);
            }
            old = next;
            rotated = next_rotated;
        }
        assert!(count > 0);
        Ok(())
    }

    #[test]
    fn algorithmic_tangent_matches_half_step_and_independent_direction() -> Result<()> {
        let model = parameters();
        let old = isotropic(100., 180.)?;
        for strain in [
            [0.0001, -0.00005, 0.00002, 0.00003, 0., 0.],
            [0.008, -0.002, 0.001, 0.015, 0.002, 0.003],
        ] {
            let (update, d) = model.stress_tangent(&old, strain, 2e-4, 1e-7)?;
            let (_, half) = model.stress_tangent_partition(&old, strain, update.substeps, 5e-8)?;
            assert!((d - half).norm() / half.norm() < 1e-4);
            let direction =
                nalgebra::SVector::<f64, 6>::from_row_slice(&[0.2, -0.3, 0.4, 0.1, -0.5, 0.6]);
            let eps = 3e-8;
            let plus = model
                .integrate_partition(
                    &old,
                    std::array::from_fn(|i| strain[i] + eps * direction[i]),
                    update.substeps,
                )?
                .state;
            let minus = model
                .integrate_partition(
                    &old,
                    std::array::from_fn(|i| strain[i] - eps * direction[i]),
                    update.substeps,
                )?
                .state;
            let measured = nalgebra::SVector::<f64, 6>::from_fn(|i, _| {
                (plus.effective_stress_compression_positive_kpa[i]
                    - minus.effective_stress_compression_positive_kpa[i])
                    / (2. * eps)
            });
            assert!((d * direction - measured).norm() / measured.norm() < 1e-4);
        }
        Ok(())
    }
}
