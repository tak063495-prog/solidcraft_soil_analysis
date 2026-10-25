//! Invariant double hardening (IDH), a distinct compression-domain soil model.
//! Drucker-Prager shear + associated elliptic cap; NOT Hardening Soil/CYSoil.
//! Effective stresses and engineering strains are compression positive.
use crate::{Error, Result, elastic, mcc};
use nalgebra::{Matrix4, SVector, Vector4};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {
    pub unloading_young_ref_kpa: f64,
    pub poisson: f64,
    pub reference_pressure_kpa: f64,
    pub stiffness_exponent: f64,
    pub minimum_mean_effective_pressure_kpa: f64,
    pub shear_ratio_initial: f64,
    pub shear_ratio_ultimate: f64,
    pub shear_reference_plastic_strain: f64,
    pub dilation_ratio: f64,
    pub shear_pressure_shift_kpa: f64,
    pub cap_axis_ratio: f64,
    pub cap_hardening_modulus_ref_kpa: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Return {
    Elastic,
    Shear,
    Cap,
    ShearCap,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub effective_stress_compression_positive_kpa: [f64; 6],
    pub plastic_strain_compression_positive: [f64; 6],
    pub void_ratio: f64,
    pub shear_plastic_strain: f64,
    pub cap_plastic_volumetric_strain: f64,
    pub initial_cap_pressure_kpa: f64,
    pub cap_pressure_kpa: f64,
    pub plastic_work_kpa: f64,
    pub last_return: Return,
}

#[derive(Debug, Serialize)]
pub struct Update {
    pub state: State,
    pub substeps: usize,
    pub plastic_substeps: usize,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialPath {
    pub schema: String,
    pub parameters: Parameters,
    pub initial_state: State,
    pub maximum_substep: f64,
    pub increments: Vec<[f64; 6]>,
}

#[derive(Debug, Serialize)]
pub struct PathResult {
    pub schema: &'static str,
    pub parameters: Parameters,
    pub states: Vec<State>,
    pub plastic_substeps: usize,
    pub material_substeps: usize,
}

pub fn run_path(path: MaterialPath) -> Result<PathResult> {
    if path.schema != "soil-fem-hardening-path/1"
        || path.increments.len() > 100_000
        || !path.maximum_substep.is_finite()
        || !(1e-7..=1e-3).contains(&path.maximum_substep)
    {
        return Err(Error::Input("invalid IDH path schema or size".into()));
    }
    path.parameters.validate_state(&path.initial_state)?;
    let mut expected_substeps = 0.;
    for de in &path.increments {
        if !de.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("nonfinite IDH path strain".into()));
        }
        let count = (de.iter().map(|v| v.abs()).sum::<f64>() / path.maximum_substep)
            .ceil()
            .max(1.);
        expected_substeps += count;
        if count > 100_000. || expected_substeps > 2_000_000. {
            return Err(Error::Input(
                "IDH material path integration budget exceeded".into(),
            ));
        }
    }
    let mut state = path.initial_state;
    let mut states = vec![state.clone()];
    let mut plastic_substeps = 0;
    let mut material_substeps = 0;
    for de in path.increments {
        let update = path
            .parameters
            .integrate(&state, de, path.maximum_substep)?;
        state = update.state;
        plastic_substeps += update.plastic_substeps;
        material_substeps += update.substeps;
        states.push(state.clone());
    }
    Ok(PathResult {
        schema: "soil-fem-hardening-path-result/1",
        parameters: path.parameters,
        states,
        plastic_substeps,
        material_substeps,
    })
}

impl Parameters {
    pub fn validate(&self) -> Result<()> {
        let finite = [
            self.unloading_young_ref_kpa,
            self.poisson,
            self.reference_pressure_kpa,
            self.stiffness_exponent,
            self.minimum_mean_effective_pressure_kpa,
            self.shear_ratio_initial,
            self.shear_ratio_ultimate,
            self.shear_reference_plastic_strain,
            self.dilation_ratio,
            self.shear_pressure_shift_kpa,
            self.cap_axis_ratio,
            self.cap_hardening_modulus_ref_kpa,
        ]
        .iter()
        .all(|v| v.is_finite());
        if !finite
            || self.unloading_young_ref_kpa <= 0.
            || !(-1. ..0.5).contains(&self.poisson)
            || self.poisson <= -1.
            || self.reference_pressure_kpa <= 0.
            || !(0. ..=1.).contains(&self.stiffness_exponent)
            || self.minimum_mean_effective_pressure_kpa < 0.
            || self.shear_ratio_initial < 0.
            || self.shear_ratio_ultimate <= self.shear_ratio_initial
            || self.shear_ratio_ultimate >= 3.
            || self.shear_reference_plastic_strain <= 0.
            || self.dilation_ratio < 0.
            || self.dilation_ratio > self.shear_ratio_initial
            || self.shear_pressure_shift_kpa < 0.
            || self.cap_axis_ratio <= 0.
            || self.cap_hardening_modulus_ref_kpa <= 0.
        {
            return Err(Error::Input("invalid invariant double-hardening parameters; require positive stiffness/pressure/scales, 0<=m<=1, Mmax>M0>=dilation>=0".into()));
        }
        Ok(())
    }

    pub fn shear_ratio(&self, history: f64) -> f64 {
        self.shear_ratio_initial
            + (self.shear_ratio_ultimate - self.shear_ratio_initial) * history
                / (self.shear_reference_plastic_strain + history)
    }

    /// Exact integral of dpc/dkappa = Href (pc/pref)^m. No pressure clipping.
    fn harden_cap(&self, pc: f64, increment: f64) -> Result<f64> {
        let b = self.cap_hardening_modulus_ref_kpa / self.reference_pressure_kpa;
        let power = 1. - self.stiffness_exponent;
        let exponent = if power == 0. {
            b * increment
        } else {
            let x = power * b * increment / (pc / self.reference_pressure_kpa).powf(power);
            if x <= -1. {
                return Err(Error::Material("IDH cap hardening left its domain".into()));
            }
            x.ln_1p() / power
        };
        let result = pc * exponent.exp();
        if !result.is_finite() || result <= 0. {
            return Err(Error::Material(
                "nonfinite/nonpositive IDH cap pressure".into(),
            ));
        }
        Ok(result)
    }

    pub fn initial_state(
        &self,
        stress: [f64; 6],
        void_ratio: f64,
        cap_pressure: f64,
    ) -> Result<State> {
        let s = State {
            effective_stress_compression_positive_kpa: stress,
            plastic_strain_compression_positive: [0.; 6],
            void_ratio,
            shear_plastic_strain: 0.,
            cap_plastic_volumetric_strain: 0.,
            initial_cap_pressure_kpa: cap_pressure,
            cap_pressure_kpa: cap_pressure,
            plastic_work_kpa: 0.,
            last_return: Return::Elastic,
        };
        self.validate_state(&s)?;
        Ok(s)
    }

    pub fn yield_values(&self, s: &State) -> [f64; 2] {
        let (p, q) = mcc::invariants(&s.effective_stress_compression_positive_kpa);
        [
            q - self.shear_ratio(s.shear_plastic_strain) * (p + self.shear_pressure_shift_kpa),
            p.hypot(q / self.cap_axis_ratio) - s.cap_pressure_kpa,
        ]
    }

    pub fn normalized_yield(&self, s: &State) -> f64 {
        let f = self.yield_values(s);
        let (p, _) = mcc::invariants(&s.effective_stress_compression_positive_kpa);
        let shear_scale = (self.shear_ratio_ultimate * (p + self.shear_pressure_shift_kpa))
            .max(self.reference_pressure_kpa);
        (f[0] / shear_scale).max(f[1] / s.cap_pressure_kpa.max(self.reference_pressure_kpa))
    }

    pub fn validate_state(&self, s: &State) -> Result<()> {
        self.validate()?;
        let (p, _) = mcc::invariants(&s.effective_stress_compression_positive_kpa);
        if !s
            .effective_stress_compression_positive_kpa
            .iter()
            .chain(&s.plastic_strain_compression_positive)
            .all(|v| v.is_finite())
            || ![
                s.void_ratio,
                s.shear_plastic_strain,
                s.cap_plastic_volumetric_strain,
                s.initial_cap_pressure_kpa,
                s.cap_pressure_kpa,
                s.plastic_work_kpa,
            ]
            .iter()
            .all(|v| v.is_finite())
            || p <= self.minimum_mean_effective_pressure_kpa
            || s.void_ratio <= 0.
            || s.shear_plastic_strain < 0.
            || s.cap_plastic_volumetric_strain < 0.
            || s.initial_cap_pressure_kpa <= 0.
            || s.cap_pressure_kpa <= 0.
            || s.plastic_work_kpa < 0.
        {
            return Err(Error::Material(
                "invalid IDH history or mean pressure outside the declared compression domain"
                    .into(),
            ));
        }
        let expected =
            self.harden_cap(s.initial_cap_pressure_kpa, s.cap_plastic_volumetric_strain)?;
        let trace = s.plastic_strain_compression_positive[..3]
            .iter()
            .sum::<f64>();
        if (expected - s.cap_pressure_kpa).abs() > 5e-9 * expected.max(1.)
            || (trace - s.cap_plastic_volumetric_strain
                + self.dilation_ratio * s.shear_plastic_strain)
                .abs()
                > 1e-9 * (1. + trace.abs())
            || self.normalized_yield(s) > 2e-8
        {
            return Err(Error::Material(
                "IDH cap/shear history, plastic trace or yield consistency failed".into(),
            ));
        }
        Ok(())
    }

    pub fn moduli(&self, s: &State) -> Result<(f64, f64)> {
        let (p, _) = mcc::invariants(&s.effective_stress_compression_positive_kpa);
        if !p.is_finite() || p <= self.minimum_mean_effective_pressure_kpa {
            return Err(Error::Material(
                "IDH stiffness outside its declared pressure domain".into(),
            ));
        }
        // Cap pressure is a compression-history stress measure. Keeping it
        // fixed during elastic unloading/reloading makes closed elastic
        // strain cycles reversible; E(p_current) would be hypoelastic and
        // leave spurious residual shear stress on nonproportional cycles.
        let e = self.unloading_young_ref_kpa
            * (s.cap_pressure_kpa / self.reference_pressure_kpa).powf(self.stiffness_exponent);
        let k = e / (3. * (1. - 2. * self.poisson));
        let g = e / (2. * (1. + self.poisson));
        if !k.is_finite() || !g.is_finite() || k <= 0. || g <= 0. {
            return Err(Error::Material("nonfinite IDH elastic moduli".into()));
        }
        Ok((k, g))
    }

    fn equations(
        &self,
        old: &State,
        trial: [f64; 2],
        kg: [f64; 2],
        mask: u8,
        x: &Vector4<f64>,
    ) -> Result<Vector4<f64>> {
        let (p, q, shear, cap) = (x[0], x[1], x[2] / (3. * kg[1]), x[3] / kg[0]);
        let radius = p.hypot(q / self.cap_axis_ratio);
        if p <= self.minimum_mean_effective_pressure_kpa
            || radius <= 0.
            || old.shear_plastic_strain + shear <= -0.5 * self.shear_reference_plastic_strain
        {
            return Err(Error::Material(
                "IDH local Newton left its algebraic domain".into(),
            ));
        }
        let cap_increment = cap * p / radius;
        let pc = self.harden_cap(old.cap_pressure_kpa, cap_increment)?;
        let a = self.shear_ratio(old.shear_plastic_strain + shear);
        let result = Vector4::new(
            p - trial[0] - kg[0] * self.dilation_ratio * shear + x[3] * p / radius,
            q - trial[1]
                + x[2]
                + (3. * kg[1] / kg[0]) * x[3] * q / (self.cap_axis_ratio.powi(2) * radius),
            if mask & 1 != 0 {
                q - a * (p + self.shear_pressure_shift_kpa)
            } else {
                x[2]
            },
            if mask & 2 != 0 { radius - pc } else { x[3] },
        );
        if !result.iter().all(|v| v.is_finite()) {
            return Err(Error::Material("nonfinite IDH local equations".into()));
        }
        Ok(result)
    }

    fn return_candidate(
        &self,
        old: &State,
        trial: [f64; 2],
        kg: [f64; 2],
        mask: u8,
    ) -> Result<Vector4<f64>> {
        let scale = trial[0]
            .abs()
            .max(trial[1])
            .max(if mask & 2 != 0 {
                old.cap_pressure_kpa
            } else {
                0.
            })
            .max(self.reference_pressure_kpa);
        // A dilative shear correction can restore a positive mean stress even
        // when the elastic trial lies outside that domain. Keep the actual
        // trial stress in the equations; only the Newton starting iterate is
        // chosen inside the admissible domain. This is not pressure clipping.
        let initial_p = if trial[0] > self.minimum_mean_effective_pressure_kpa {
            trial[0]
        } else {
            mcc::invariants(&old.effective_stress_compression_positive_kpa).0
        };
        let mut x = Vector4::new(initial_p, trial[1], 0., 0.);
        for _ in 0..60 {
            let r = self.equations(old, trial, kg, mask, &x)?;
            if r.amax() <= 2e-11 * scale {
                if x[1] < -1e-10 * scale || x[2] < -1e-10 * scale || x[3] < -1e-10 * scale {
                    return Err(Error::Material(
                        "IDH negative plastic multiplier in active set".into(),
                    ));
                }
                return Ok(x);
            }
            let mut jac = Matrix4::zeros();
            for j in 0..4 {
                let h = 2e-6 * scale.max(x[j].abs());
                let mut plus = x;
                plus[j] += h;
                let mut minus = x;
                minus[j] -= h;
                let a = self.equations(old, trial, kg, mask, &plus);
                let b = self.equations(old, trial, kg, mask, &minus);
                let column = match (a, b) {
                    (Ok(a), Ok(b)) => (a - b) / (2. * h),
                    (Ok(a), Err(_)) => (a - r) / h,
                    (Err(_), Ok(b)) => (r - b) / h,
                    _ => return Err(Error::Material("IDH local derivative unavailable".into())),
                };
                jac.set_column(j, &column);
            }
            let step = jac
                .lu()
                .solve(&(-r))
                .ok_or_else(|| Error::Material("singular IDH active-set Jacobian".into()))?;
            let mut alpha = 1.;
            let mut accepted = false;
            for _ in 0..24 {
                let next = x + alpha * step;
                if let Ok(candidate) = self.equations(old, trial, kg, mask, &next)
                    && candidate.norm() < r.norm() * (1. - 1e-4 * alpha)
                {
                    x = next;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                return Err(Error::Material("IDH active-set line search failed".into()));
            }
        }
        Err(Error::Material("IDH active-set iteration limit".into()))
    }

    fn substep(&self, old: &State, strain: [f64; 6]) -> Result<(State, bool)> {
        if strain == [0.; 6] {
            return Ok((old.clone(), false));
        }
        let (k, g) = self.moduli(old)?;
        let e = 9. * k * g / (3. * k + g);
        let d = elastic::constitutive(e, self.poisson)?;
        let guess =
            SVector::<f64, 6>::from_row_slice(&old.effective_stress_compression_positive_kpa)
                + d * SVector::<f64, 6>::from_row_slice(&strain);
        let trial: [f64; 6] = guess.into();
        let (pt, qt) = mcc::invariants(&trial);
        let mut next = old.clone();
        next.effective_stress_compression_positive_kpa = trial;
        next.void_ratio = (1. + old.void_ratio) * (-strain[..3].iter().sum::<f64>()).exp() - 1.;
        next.last_return = Return::Elastic;
        if self.normalized_yield(&next) <= 1e-11 {
            self.validate_state(&next)?;
            return Ok((next, false));
        }
        let mut last_error = String::new();
        for mask in [1, 2, 3] {
            let x = match self.return_candidate(old, [pt, qt], [k, g], mask) {
                Ok(x) => x,
                Err(err) => {
                    last_error = err.to_string();
                    continue;
                }
            };
            let (p, q) = (x[0], x[1].max(0.));
            let ls = x[2].max(0.) / (3. * g);
            let lc = x[3].max(0.) / k;
            let radius = p.hypot(q / self.cap_axis_ratio);
            let cap_increment = lc * p / radius;
            let mut candidate = next.clone();
            candidate.shear_plastic_strain += ls;
            candidate.cap_plastic_volumetric_strain += cap_increment;
            candidate.cap_pressure_kpa = self.harden_cap(old.cap_pressure_kpa, cap_increment)?;
            candidate.last_return = match mask {
                1 => Return::Shear,
                2 => Return::Cap,
                _ => Return::ShearCap,
            };
            let mut flow = [0.; 6];
            for j in 0..6 {
                let dev = trial[j] - if j < 3 { pt } else { 0. };
                let nq = if qt > 0. {
                    (if j < 3 { 1.5 } else { 3. }) * dev / qt
                } else {
                    0.
                };
                candidate.effective_stress_compression_positive_kpa[j] =
                    if j < 3 { p } else { 0. } + if qt > 0. { q * dev / qt } else { 0. };
                flow[j] = ls * (nq - if j < 3 { self.dilation_ratio / 3. } else { 0. })
                    + lc * (if j < 3 { p / (3. * radius) } else { 0. }
                        + q * nq / (self.cap_axis_ratio.powi(2) * radius));
                candidate.plastic_strain_compression_positive[j] += flow[j];
            }
            let work: f64 = candidate
                .effective_stress_compression_positive_kpa
                .iter()
                .zip(flow)
                .map(|(s, e)| s * e)
                .sum();
            if work < -1e-10 || !work.is_finite() {
                continue;
            }
            candidate.plastic_work_kpa += work.max(0.);
            let elastic_strain = SVector::<f64, 6>::from_fn(|i, _| strain[i] - flow[i]);
            let closure = SVector::<f64, 6>::from_row_slice(
                &candidate.effective_stress_compression_positive_kpa,
            ) - SVector::<f64, 6>::from_row_slice(
                &old.effective_stress_compression_positive_kpa,
            ) - d * elastic_strain;
            if closure.amax() > 2e-8 * pt.max(qt).max(p).max(self.reference_pressure_kpa) {
                continue;
            }
            if self.validate_state(&candidate).is_ok() {
                return Ok((candidate, true));
            }
        }
        Err(Error::Material(format!(
            "no admissible IDH shear/cap active set; {last_error}"
        )))
    }

    pub fn integrate_partition(
        &self,
        old: &State,
        strain: [f64; 6],
        steps: usize,
    ) -> Result<Update> {
        self.validate_state(old)?;
        if steps == 0 || steps > 100_000 || !strain.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("invalid fixed IDH partition".into()));
        }
        let mut state = old.clone();
        let mut plastic = 0;
        for _ in 0..steps {
            let (next, yielded) = self.substep(&state, strain.map(|v| v / steps as f64))?;
            state = next;
            plastic += usize::from(yielded);
        }
        Ok(Update {
            state,
            substeps: steps,
            plastic_substeps: plastic,
        })
    }

    pub fn integrate(&self, old: &State, strain: [f64; 6], maximum: f64) -> Result<Update> {
        if !maximum.is_finite()
            || !(1e-7..=1e-3).contains(&maximum)
            || !strain.iter().all(|v| v.is_finite())
        {
            return Err(Error::Input("invalid IDH strain partition controls".into()));
        }
        let count = (strain.iter().map(|v| v.abs()).sum::<f64>() / maximum)
            .ceil()
            .max(1.);
        if count > 100_000. {
            return Err(Error::Material("IDH partition budget exceeded".into()));
        }
        self.integrate_partition(old, strain, count as usize)
    }

    pub fn stress_tangent_partition(
        &self,
        old: &State,
        strain: [f64; 6],
        steps: usize,
        h: f64,
    ) -> Result<(Update, elastic::Tensor)> {
        if !h.is_finite() || !(1e-9..=1e-5).contains(&h) {
            return Err(Error::Input("invalid IDH tangent step".into()));
        }
        let update = self.integrate_partition(old, strain, steps)?;
        let mut tangent = elastic::Tensor::zeros();
        for j in 0..6 {
            let mut a = strain;
            a[j] += h;
            let mut b = strain;
            b[j] -= h;
            let a = self.integrate_partition(old, a, steps)?.state;
            let b = self.integrate_partition(old, b, steps)?.state;
            for i in 0..6 {
                tangent[(i, j)] = (a.effective_stress_compression_positive_kpa[i]
                    - b.effective_stress_compression_positive_kpa[i])
                    / (2. * h);
            }
        }
        if !tangent.iter().all(|v| v.is_finite()) {
            return Err(Error::Material("nonfinite IDH tangent".into()));
        }
        Ok((update, tangent))
    }
}
