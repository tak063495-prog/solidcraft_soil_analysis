//! Second-order MCC with exact elastic evolution and midpoint plastic flow.
//! Local implicit sensitivities include hardening, all partition history and
//! the elastic first-yield event. No stress projection or tangent symmetrization.
use crate::{Error, Result, mcc};
use nalgebra::{SMatrix, SVector};
use std::ops::{Add, Div, Mul, Neg, Sub};
#[derive(Clone, Copy, Debug)]
struct D<const N: usize> {
    v: f64,
    g: [f64; N],
}
impl<const N: usize> D<N> {
    fn c(v: f64) -> Self {
        Self { v, g: [0.; N] }
    }
    fn seed(v: f64, i: usize) -> Self {
        let mut x = Self::c(v);
        x.g[i] = 1.;
        x
    }
    fn exp(self) -> Self {
        let v = self.v.exp();
        Self {
            v,
            g: self.g.map(|x| x * v),
        }
    }
    fn expm1(self) -> Self {
        Self {
            v: self.v.exp_m1(),
            g: self.g.map(|x| x * self.v.exp()),
        }
    }
    fn square(self) -> Self {
        self * self
    }
}
impl<const N: usize> Add for D<N> {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            v: self.v + b.v,
            g: std::array::from_fn(|i| self.g[i] + b.g[i]),
        }
    }
}
impl<const N: usize> Sub for D<N> {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            v: self.v - b.v,
            g: std::array::from_fn(|i| self.g[i] - b.g[i]),
        }
    }
}
// Product and quotient rules intentionally combine arithmetic operations.
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const N: usize> Mul for D<N> {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            g: std::array::from_fn(|i| self.g[i] * b.v + self.v * b.g[i]),
        }
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl<const N: usize> Div for D<N> {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        Self {
            v: self.v / b.v,
            g: std::array::from_fn(|i| (self.g[i] - self.v / b.v * b.g[i]) / b.v),
        }
    }
}
impl<const N: usize> Neg for D<N> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            g: self.g.map(|x| -x),
        }
    }
}
impl<const N: usize> Add<f64> for D<N> {
    type Output = Self;
    fn add(self, b: f64) -> Self {
        self + Self::c(b)
    }
}
impl<const N: usize> Sub<f64> for D<N> {
    type Output = Self;
    fn sub(self, b: f64) -> Self {
        self - Self::c(b)
    }
}
impl<const N: usize> Mul<f64> for D<N> {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        Self {
            v: self.v * b,
            g: self.g.map(|x| x * b),
        }
    }
}
impl<const N: usize> Div<f64> for D<N> {
    type Output = Self;
    fn div(self, b: f64) -> Self {
        self * (1. / b)
    }
}
fn sum<const N: usize>(a: impl IntoIterator<Item = D<N>>) -> D<N> {
    a.into_iter().fold(D::c(0.), |a, b| a + b)
}
fn exprel<const N: usize>(x: D<N>) -> D<N> {
    if x.v.abs() < 1e-5 {
        let mut result = D::c(1.);
        let mut term = D::c(1.);
        for i in 1..=9 {
            term = term * x / (i as f64 + 1.);
            result = result + term;
        }
        result
    } else {
        x.expm1() / x
    }
}
#[derive(Clone, Copy)]
struct Tracked<const N: usize> {
    pc_reference: f64,
    stress: [D<N>; 6],
    pc: D<N>,
    volume: D<N>,
    plastic: [D<N>; 6],
}
impl<const N: usize> Tracked<N> {
    fn from_state(s: &mcc::State) -> Self {
        Self {
            pc_reference: s.preconsolidation_pressure_kpa,
            stress: s.effective_stress_compression_positive_kpa.map(D::c),
            pc: D::c(0.),
            volume: D::c((1. + s.void_ratio).ln()),
            plastic: s.plastic_strain_compression_positive.map(D::c),
        }
    }
    fn values(self) -> mcc::State {
        mcc::State {
            effective_stress_compression_positive_kpa: self.stress.map(|x| x.v),
            preconsolidation_pressure_kpa: self.pc_reference * self.pc.v.exp(),
            void_ratio: self.volume.v.exp() - 1.,
            plastic_strain_compression_positive: self.plastic.map(|x| x.v),
        }
    }
}
fn yield_value<const N: usize>(p: &mcc::Parameters, s: &Tracked<N>) -> D<N> {
    let pressure = sum(s.stress[..3].iter().copied()) / 3.;
    let pc = s.pc.exp() * s.pc_reference;
    let m2 = p.critical_state_m.powi(2);
    let q2 = sum((0..6).map(|i| {
        if i < 3 {
            (s.stress[i] - pressure).square() * 1.5
        } else {
            s.stress[i].square() * 3.
        }
    }));
    (q2 + pressure * (pressure - pc) * m2) / (pc.square() * m2)
}
fn elastic<const N: usize>(
    p: &mcc::Parameters,
    old: &Tracked<N>,
    de: [D<N>; 6],
) -> Result<Tracked<N>> {
    let ev = sum(de[..3].iter().copied());
    let v = old.volume.exp();
    let pressure = sum(old.stress[..3].iter().copied()) / 3.;
    let a = -v / p.kappa * (-ev).expm1();
    let integrated = v * pressure / p.kappa * exprel(-ev) * exprel(a);
    let ratio = 3. * (1. - 2. * p.poisson) / (2. * (1. + p.poisson));
    let mut next = *old;
    next.volume = old.volume - ev;
    next.stress = std::array::from_fn(|i| {
        old.stress[i]
            + if i < 3 {
                integrated * (de[i] - ev / 3.) * (2. * ratio) + pressure * a.expm1()
            } else {
                integrated * de[i] * ratio
            }
    });
    let values = next.values();
    if values.void_ratio <= 0.
        || !values.void_ratio.is_finite()
        || mcc::invariants(&values.effective_stress_compression_positive_kpa).0 <= 0.
        || !values
            .effective_stress_compression_positive_kpa
            .iter()
            .all(|x| x.is_finite())
    {
        return Err(Error::Material(
            "midpoint elastic path left p>0 and v>1".into(),
        ));
    }
    Ok(next)
}
fn equations<const N: usize>(
    p: &mcc::Parameters,
    old: &Tracked<N>,
    de: [D<N>; 6],
    x: [D<N>; 8],
) -> ([D<N>; 8], [D<N>; 6]) {
    let scale = old.pc.exp() * old.pc_reference;
    let s0 = old.stress.map(|x| x / scale);
    let sm: [D<N>; 6] = std::array::from_fn(|i| (s0[i] + x[i]) / 2.);
    let pm = sum(sm[..3].iter().copied()) / 3.;
    let pcm = (x[6] / 2.).exp();
    let m2 = p.critical_state_m.powi(2);
    let trace = (pm * 2. - pcm) * m2;
    let normal: [D<N>; 6] = std::array::from_fn(|i| {
        if i < 3 {
            (sm[i] - pm) * 3. + trace / 3.
        } else {
            sm[i] * 6.
        }
    });
    let vm = (old.volume - sum(de[..3].iter().copied()) / 2.).exp();
    let bulk = vm * pm / p.kappa;
    let shear = bulk * (3. * (1. - 2. * p.poisson) / (2. * (1. + p.poisson)));
    let flow = normal.map(|n| n * x[7]);
    let remainder: [D<N>; 6] = std::array::from_fn(|i| de[i] - flow[i]);
    let ev = sum(remainder[..3].iter().copied());
    let mut r = std::array::from_fn(|i| {
        if i < 6 {
            x[i] - s0[i]
                - if i < 3 {
                    bulk * ev + shear * (remainder[i] - ev / 3.) * 2.
                } else {
                    shear * remainder[i]
                }
        } else {
            D::c(0.)
        }
    });
    r[6] = x[6] - vm / (p.lambda - p.kappa) * sum(flow[..3].iter().copied());
    let p1 = (x[0] + x[1] + x[2]) / 3.;
    let pc1 = x[6].exp();
    let q2 = sum((0..6).map(|i| {
        if i < 3 {
            (x[i] - p1).square() * 1.5
        } else {
            x[i].square() * 3.
        }
    }));
    r[7] = (q2 + p1 * (p1 - pc1) * m2) / (pc1.square() * m2);
    (r, flow)
}
fn plastic<const N: usize>(
    p: &mcc::Parameters,
    old: &Tracked<N>,
    de: [D<N>; 6],
    tangent: bool,
) -> Result<Tracked<N>> {
    let old8 = Tracked::<8>::from_state(&old.values());
    let strain8 = de.map(|x| D::<8>::c(x.v));
    let trial = elastic(p, &old8, strain8)?;
    let scale = old.pc_reference * old.pc.v.exp();
    let mut x = SVector::<f64, 8>::zeros();
    for i in 0..6 {
        x[i] = trial.stress[i].v / scale;
    }
    let mut solution = None;
    for _ in 0..40 {
        let variables = std::array::from_fn(|i| D::<8>::seed(x[i], i));
        let r = equations(p, &old8, strain8, variables).0;
        let value = SVector::<f64, 8>::from_fn(|i, _| r[i].v);
        let jac = SMatrix::<f64, 8, 8>::from_fn(|i, j| r[i].g[j]);
        if value.norm() < 1e-12 {
            solution = Some((x, jac));
            break;
        }
        let correction = jac
            .lu()
            .solve(&(-value))
            .ok_or_else(|| Error::Material("singular midpoint local Jacobian".into()))?;
        let mut accepted = false;
        let mut alpha = 1.;
        for _ in 0..30 {
            let candidate = x + correction * alpha;
            if candidate.iter().all(|x| x.is_finite())
                && candidate[6].abs() < 100.
                && candidate[7] >= 0.
                && candidate[0] + candidate[1] + candidate[2] > 0.
            {
                let variables = std::array::from_fn(|i| D::<8>::c(candidate[i]));
                let residual = equations(p, &old8, strain8, variables).0;
                if sum(residual.map(D::square)).v.sqrt() < value.norm() {
                    x = candidate;
                    accepted = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if !accepted {
            return Err(Error::Material("midpoint local line search failed".into()));
        }
    }
    let (x, jac) =
        solution.ok_or_else(|| Error::Material("midpoint local iteration limit".into()))?;
    let mut variables = std::array::from_fn(|i| D::<N>::c(x[i]));
    if tangent {
        let residual = equations(p, old, de, variables).0;
        let rhs = SMatrix::<f64, 8, N>::from_fn(|i, j| -residual[i].g[j]);
        let derivative = jac
            .lu()
            .solve(&rhs)
            .ok_or_else(|| Error::Material("singular midpoint sensitivity system".into()))?;
        for i in 0..8 {
            variables[i].g = std::array::from_fn(|j| derivative[(i, j)]);
        }
    }
    let mut next = *old;
    let scale = old.pc.exp() * old.pc_reference;
    next.stress = std::array::from_fn(|i| variables[i] * scale);
    next.pc = old.pc + variables[6];
    next.volume = old.volume - sum(de[..3].iter().copied());
    let flow = equations(p, old, de, variables).1;
    let work = (0..6)
        .map(|i| (old.stress[i].v + next.stress[i].v) * flow[i].v / 2.)
        .sum::<f64>();
    let scale = (0..6)
        .map(|i| ((old.stress[i].v + next.stress[i].v) * flow[i].v / 2.).abs())
        .sum::<f64>();
    if work < -1e-12 * scale.max(1e-15) {
        return Err(Error::Material(
            "negative midpoint plastic work; reduce the material increment".into(),
        ));
    }
    next.plastic = std::array::from_fn(|i| old.plastic[i] + flow[i]);
    p.validate_state(&next.values())?;
    Ok(next)
}
fn lift<const N: usize>(x: D<N>) -> D<7> {
    D {
        v: x.v,
        g: std::array::from_fn(|i| {
            if i == 0 {
                0.
            } else {
                x.g.get(i - 1).copied().unwrap_or(0.)
            }
        }),
    }
}
fn step<const N: usize>(
    p: &mcc::Parameters,
    old: &Tracked<N>,
    de: [D<N>; 6],
    tangent: bool,
) -> Result<(Tracked<N>, bool)> {
    let trial = elastic(p, old, de)?;
    let f0 = yield_value(p, old).v;
    let f1 = yield_value(p, &trial).v;
    let s = old.values();
    let (pressure, _) = mcc::invariants(&s.effective_stress_compression_positive_kpa);
    let pc = s.preconsolidation_pressure_kpa;
    let m2 = p.critical_state_m.powi(2);
    let delta: [f64; 6] = std::array::from_fn(|i| trial.stress[i].v - old.stress[i].v);
    let normal: [f64; 6] = std::array::from_fn(|i| {
        if i < 3 {
            3. * (old.stress[i].v - pressure) + m2 * (2. * pressure - pc) / 3.
        } else {
            6. * old.stress[i].v
        }
    });
    let loading = (0..6).map(|i| normal[i] * delta[i]).sum::<f64>() / (m2 * pc * pc);
    if f1 <= 1e-11 || (loading <= 0. && f1 <= f0) {
        return Ok((trial, false));
    }
    if f0 < -2e-9 || loading < 0. {
        let a = f1 - f0 - loading;
        let mut lo = if loading < 0. {
            if a <= 0. {
                return Err(Error::Material("invalid unloading yield bracket".into()));
            }
            -loading / (2. * a)
        } else {
            0.
        };
        let mut hi = 1.;
        if !lo.is_finite() || !(0. ..1.).contains(&lo) {
            return Err(Error::Material(
                "yield minimum is outside the elastic interval".into(),
            ));
        }
        for _ in 0..60 {
            let eta = (lo + hi) / 2.;
            let mut state = s.clone();
            state.effective_stress_compression_positive_kpa =
                std::array::from_fn(|i| old.stress[i].v + eta * delta[i]);
            if p.normalized_yield(&state) > 0. {
                hi = eta;
            } else {
                lo = eta;
            }
        }
        let eta = (lo + hi) / 2.;
        let ev = sum(de[..3].iter().copied()).v;
        let volume = old.volume.v.exp();
        let a = -volume / p.kappa * (-ev).exp_m1();
        let fraction = if ev.abs() < 1e-14 {
            eta
        } else {
            -(-p.kappa / volume * (eta * a.exp_m1()).ln_1p()).ln_1p() / ev
        };
        if !(0. ..=1.).contains(&fraction) {
            return Err(Error::Material("invalid elastic yield fraction".into()));
        }
        let mut t = D::<N>::c(fraction);
        if tangent {
            let old7 = Tracked {
                pc_reference: old.pc_reference,
                stress: old.stress.map(lift),
                pc: lift(old.pc),
                volume: lift(old.volume),
                plastic: old.plastic.map(lift),
            };
            let scalar = D::<7>::seed(fraction, 0);
            let prefix = elastic(p, &old7, de.map(|x| lift(x) * scalar))?;
            let equation = yield_value(p, &prefix);
            if !equation.g[0].is_finite() || equation.g[0] <= 0. {
                return Err(Error::Material(
                    "singular elastic yield event derivative".into(),
                ));
            }
            t.g = std::array::from_fn(|i| -equation.g[i + 1] / equation.g[0]);
        }
        let prefix = elastic(p, old, de.map(|x| x * t))?;
        if yield_value(p, &prefix).v.abs() > 1e-10 {
            return Err(Error::Material(
                "elastic first yield residual exceeds tolerance".into(),
            ));
        }
        return Ok((
            plastic(p, &prefix, de.map(|x| x * (D::c(1.) - t)), tangent)?,
            true,
        ));
    }
    Ok((plastic(p, old, de, tangent)?, true))
}
fn integrate_tracked<const N: usize>(
    p: &mcc::Parameters,
    committed: &mcc::State,
    strain: [f64; 6],
    steps: usize,
    tangent: bool,
) -> Result<(mcc::Update, crate::elastic::Tensor)> {
    p.validate_state(committed)?;
    if steps == 0 || steps > 100_000 || !strain.iter().all(|x| x.is_finite()) {
        return Err(Error::Input("invalid midpoint fixed partition".into()));
    }
    let mut state = Tracked::<N>::from_state(committed);
    let de = std::array::from_fn(|i| {
        let mut x = D::c(strain[i] / steps as f64);
        if tangent && let Some(gradient) = x.g.get_mut(i) {
            *gradient = 1. / steps as f64;
        }
        x
    });
    let mut plastic_count = 0;
    for _ in 0..steps {
        let (next, yielded) = step(p, &state, de, tangent)?;
        state = next;
        plastic_count += usize::from(yielded);
    }
    let final_state = state.values();
    p.validate_state(&final_state)?;
    let tensor =
        crate::elastic::Tensor::from_fn(|i, j| state.stress[i].g.get(j).copied().unwrap_or(0.));
    if !tensor.iter().all(|x| x.is_finite()) {
        return Err(Error::Material(
            "nonfinite midpoint fixed partition tangent".into(),
        ));
    }
    Ok((
        mcc::Update {
            state: final_state,
            substeps: steps,
            plastic_substeps: plastic_count,
        },
        tensor,
    ))
}

/// State-only trials retain the same primal map without six unused sensitivities.
pub fn integrate(
    p: &mcc::Parameters,
    committed: &mcc::State,
    strain: [f64; 6],
    steps: usize,
    tangent: bool,
) -> Result<(mcc::Update, crate::elastic::Tensor)> {
    if tangent {
        integrate_tracked::<6>(p, committed, strain, steps, true)
    } else {
        integrate_tracked::<0>(p, committed, strain, steps, false)
    }
}
