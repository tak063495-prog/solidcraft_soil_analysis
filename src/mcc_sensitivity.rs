//! Fixed-partition MCC sensitivities, using the existing primal return.
//! No stress, history, multiplier or volumetric tangent is clipped or shifted.
use super::*;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy)]
struct Dual {
    value: f64,
    derivative: [f64; 6],
}
impl Dual {
    fn constant(value: f64) -> Self {
        Self {
            value,
            derivative: [0.; 6],
        }
    }
    fn square(self) -> Self {
        self * self
    }
    fn ln(self) -> Self {
        Self {
            value: self.value.ln(),
            derivative: self.derivative.map(|x| x / self.value),
        }
    }
    fn exp(self) -> Self {
        let value = self.value.exp();
        Self {
            value,
            derivative: self.derivative.map(|x| value * x),
        }
    }
    fn sqrt_norm(self) -> Self {
        let value = self.value.sqrt();
        // The scalar q norm is nondifferentiable at zero. The stress map has
        // a unique isotropic tensor limit: its deviatoric factor is evaluated
        // without dividing by q, and p/pc do not couple to dq at q=0.
        Self {
            value,
            derivative: if value == 0. {
                [0.; 6]
            } else {
                self.derivative.map(|x| x / (2. * value))
            },
        }
    }
}
impl Add for Dual {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            value: self.value + b.value,
            derivative: std::array::from_fn(|i| self.derivative[i] + b.derivative[i]),
        }
    }
}
impl Sub for Dual {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            value: self.value - b.value,
            derivative: std::array::from_fn(|i| self.derivative[i] - b.derivative[i]),
        }
    }
}
impl Mul for Dual {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, b: Self) -> Self {
        Self {
            value: self.value * b.value,
            derivative: std::array::from_fn(|i| {
                self.derivative[i] * b.value + self.value * b.derivative[i]
            }),
        }
    }
}
impl Div for Dual {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, b: Self) -> Self {
        Self {
            value: self.value / b.value,
            derivative: std::array::from_fn(|i| {
                (self.derivative[i] - self.value / b.value * b.derivative[i]) / b.value
            }),
        }
    }
}
impl Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            derivative: self.derivative.map(|x| -x),
        }
    }
}
macro_rules! scalar {
    ($trait:ident,$method:ident,$op:tt) => {
        impl $trait<f64> for Dual {
            type Output=Self;
            fn $method(self,b:f64)->Self {self $op Dual::constant(b)}
        }
        impl $trait<Dual> for f64 {
            type Output=Dual;
            fn $method(self,b:Dual)->Dual {Dual::constant(self) $op b}
        }
    }
}
scalar!(Add,add,+);
scalar!(Sub,sub,-);
scalar!(Mul,mul,*);
scalar!(Div,div,/);

fn pressure(s: &[Dual; 6]) -> Dual {
    (s[0] + s[1] + s[2]) / 3.
}
fn q(s: &[Dual; 6]) -> Dual {
    let p = pressure(s);
    (1.5 * ((s[0] - p).square()
        + (s[1] - p).square()
        + (s[2] - p).square()
        + 2. * (s[3].square() + s[4].square() + s[5].square())))
    .sqrt_norm()
}

struct ReturnInput {
    trial: [Dual; 6],
    k: Dual,
    g: Dual,
    volume: Dual,
    pc: Dual,
}

fn plastic(
    model: &Parameters,
    old: &State,
    next: &State,
    input: ReturnInput,
) -> Result<([Dual; 6], Dual)> {
    let ReturnInput {
        trial,
        k,
        g,
        volume,
        pc,
    } = input;
    let pt = pressure(&trial);
    let qt = q(&trial);
    let (p, q) = invariants(&next.effective_stress_compression_positive_kpa);
    let scale = old
        .preconsolidation_pressure_kpa
        .max(pt.value)
        .max(qt.value);
    let (pn, qn, cn) = (
        p / scale,
        q / scale,
        next.preconsolidation_pressure_kpa / scale,
    );
    let (kn, gn) = (k / scale, g / scale);
    let a = volume / (model.lambda - model.kappa);
    let m2 = model.critical_state_m.powi(2);
    let h = m2 * (2. * pn - cn);
    // Reconstruct the multiplier from all three converged return equations.
    // Normalization is frozen as a coordinate scaling, not a material input.
    // This avoids division by q or by the volumetric flow at critical state.
    let coefficients = [kn.value * h, 6. * gn.value * qn, a.value * h];
    let rhs = [
        pt.value / scale - pn,
        qt.value / scale - qn,
        cn.ln() - (old.preconsolidation_pressure_kpa / scale).ln(),
    ];
    let denominator = coefficients.iter().map(|x| x * x).sum::<f64>();
    if !denominator.is_finite() || denominator <= 0. {
        return Err(Error::Material(
            "degenerate MCC sensitivity multiplier equations".into(),
        ));
    }
    let multiplier = coefficients
        .iter()
        .zip(rhs)
        .map(|(a, b)| a * b)
        .sum::<f64>()
        / denominator;
    if !multiplier.is_finite() || multiplier < -1e-10 {
        return Err(Error::Material(
            "inadmissible MCC sensitivity multiplier".into(),
        ));
    }
    let residual = Vector4::new(
        pn - pt.value / scale + kn.value * multiplier * h,
        qn - qt.value / scale + 6. * gn.value * multiplier * qn,
        cn.ln() - (old.preconsolidation_pressure_kpa / scale).ln() - a.value * multiplier * h,
        qn * qn + m2 * pn * (pn - cn),
    );
    if !residual.norm().is_finite() || residual.norm() > 1e-8 {
        return Err(Error::Material(
            "MCC sensitivity return reconstruction residual".into(),
        ));
    }
    let jacobian = Matrix4::from_row_slice(&[
        1. + 2. * kn.value * multiplier * m2,
        0.,
        -kn.value * multiplier * m2 * cn,
        kn.value * h,
        0.,
        1. + 6. * gn.value * multiplier,
        0.,
        6. * gn.value * qn,
        -2. * a.value * multiplier * m2,
        0.,
        1. + a.value * multiplier * m2 * cn,
        -a.value * h,
        m2 * (2. * pn - cn),
        2. * qn,
        -m2 * pn * cn,
        0.,
    ]);
    // x is held fixed here. Only trial stress, old moduli, old v and old pc
    // vary. Solve J dx = -dF for all six complete-increment directions at once.
    let fixed = [
        Dual::constant(pn) - pt / scale + kn * multiplier * h,
        Dual::constant(qn) - qt / scale + 6. * gn * multiplier * qn,
        Dual::constant(cn.ln()) - (pc / scale).ln() - a * multiplier * h,
        Dual::constant(0.),
    ];
    let rhs = nalgebra::SMatrix::<f64, 4, 6>::from_fn(|r, c| -fixed[r].derivative[c]);
    let sensitivity = jacobian
        .lu()
        .solve(&rhs)
        .ok_or_else(|| Error::Material("singular MCC sensitivity Jacobian".into()))?;
    let pressure = Dual {
        value: p,
        derivative: std::array::from_fn(|i| scale * sensitivity[(0, i)]),
    };
    let multiplier = Dual {
        value: multiplier,
        derivative: std::array::from_fn(|i| sensitivity[(3, i)]),
    };
    let factor = 1. / (1. + 6. * gn * multiplier);
    let stress = std::array::from_fn(|i| {
        if i < 3 {
            pressure + factor * (trial[i] - pt)
        } else {
            factor * trial[i]
        }
    });
    let pc = Dual {
        value: next.preconsolidation_pressure_kpa,
        derivative: std::array::from_fn(|i| {
            next.preconsolidation_pressure_kpa * sensitivity[(2, i)]
        }),
    };
    Ok((stress, pc))
}

pub(super) fn integrate(
    model: &Parameters,
    committed: &State,
    strain: [f64; 6],
    steps: usize,
) -> Result<(Update, crate::elastic::Tensor)> {
    model.validate_state(committed)?;
    if steps == 0 || steps > 100_000 || !strain.iter().all(|x| x.is_finite()) {
        return Err(Error::Input(
            "invalid fixed MCC sensitivity partition".into(),
        ));
    }
    let increment = strain.map(|x| x / steps as f64);
    let de: [Dual; 6] = std::array::from_fn(|i| Dual {
        value: increment[i],
        derivative: std::array::from_fn(|j| if i == j { 1. / steps as f64 } else { 0. }),
    });
    let ev = de[0] + de[1] + de[2];
    let mut state = committed.clone();
    let mut stress = state
        .effective_stress_compression_positive_kpa
        .map(Dual::constant);
    let mut void = Dual::constant(state.void_ratio);
    let mut pc = Dual::constant(state.preconsolidation_pressure_kpa);
    let mut count = 0;
    for _ in 0..steps {
        let k = (1. + void) * pressure(&stress) / model.kappa;
        let g = 3. * k * (1. - 2. * model.poisson) / (2. * (1. + model.poisson));
        let mut trial = stress;
        for i in 0..3 {
            trial[i] = trial[i] + (k * ev + 2. * g * (de[i] - ev / 3.));
        }
        for i in 3..6 {
            trial[i] = trial[i] + g * de[i];
        }
        // Every primal value and branch comes from the already verified return.
        let (next, yielded) = model.substep(&state, increment)?;
        let next_void = (1. + void) * (-ev).exp() - 1.;
        if yielded {
            (stress, pc) = plastic(
                model,
                &state,
                &next,
                ReturnInput {
                    trial,
                    k,
                    g,
                    volume: 1. + void,
                    pc,
                },
            )?;
            count += 1;
        } else {
            stress = trial;
        }
        void = next_void;
        for (i, s) in stress.iter_mut().enumerate() {
            s.value = next.effective_stress_compression_positive_kpa[i];
        }
        void.value = next.void_ratio;
        pc.value = next.preconsolidation_pressure_kpa;
        state = next;
    }
    model.validate_state(&state)?;
    let tangent = crate::elastic::Tensor::from_fn(|i, j| stress[i].derivative[j]);
    if !tangent.iter().all(|x| x.is_finite()) {
        return Err(Error::Material(
            "nonfinite MCC fixed-partition sensitivity".into(),
        ));
    }
    Ok((
        Update {
            state,
            substeps: steps,
            plastic_substeps: count,
        },
        tangent,
    ))
}
