//! Compression-positive 3D Mohr-Coulomb perfect plasticity. Exact polyhedral
//! multisurface return in trial principal axes; no smoothing or pressure clip.
use crate::{Error, Result, elastic};
use nalgebra::{DMatrix, DVector, Matrix3, SVector, Vector3};
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;

fn candidate_sets(tension: bool) -> &'static [usize] {
    static PURE: OnceLock<Vec<usize>> = OnceLock::new();
    static CAPPED: OnceLock<Vec<usize>> = OnceLock::new();
    let cache = if tension { &CAPPED } else { &PURE };
    cache.get_or_init(|| {
        // The smooth face and the two adjacent edges in sorted trial axes,
        // followed by Rankine min, min/mid, apex and the MC/Rankine edge.
        // Every remaining nonsingular cone representation is still enumerated.
        let mut candidates = vec![2, 3, 10];
        if tension {
            candidates.extend([256, 384, 448, 258, 259, 266]);
        }
        for mask in 1usize..(1usize << if tension { 9 } else { 6 }) {
            if mask.count_ones() <= 3 && !candidates.contains(&mask) {
                candidates.push(mask);
            }
        }
        candidates
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {
    pub young_kpa: f64,
    pub poisson: f64,
    pub cohesion_kpa: f64,
    pub friction_deg: f64,
    pub dilation_deg: f64,
    /// Explicit associated Rankine surface; None means pure MC, including its
    /// genuine no-solution tension paths at zero dilation.
    pub tension_cutoff_kpa: Option<f64>,
    pub flow: Flow,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Flow {
    NonAssociated,
    Associated,
    Davis,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Return {
    Elastic,
    ShearFace,
    CompressionEdge,
    ExtensionEdge,
    ShearApex,
    TensionFace,
    TensionEdge,
    TensionApex,
    Mixed,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct State {
    pub effective_stress_compression_positive_kpa: [f64; 6],
    pub plastic_strain_compression_positive: [f64; 6],
    pub void_ratio: f64,
    pub equivalent_plastic_strain: f64,
    pub plastic_work_kpa: f64,
    pub last_return: Return,
}
#[derive(Debug, Serialize)]
pub struct Update {
    pub state: State,
    pub substeps: usize,
    pub plastic_substeps: usize,
}
#[derive(Clone)]
struct Surface {
    n: Vector3<f64>,
    m: Vector3<f64>,
    c: f64,
}

pub fn tensor(s: &[f64; 6]) -> Matrix3<f64> {
    Matrix3::new(s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2])
}
pub fn voigt(s: &Matrix3<f64>) -> [f64; 6] {
    [
        s[(0, 0)],
        s[(1, 1)],
        s[(2, 2)],
        s[(0, 1)],
        s[(1, 2)],
        s[(2, 0)],
    ]
}
pub fn principal(s: &[f64; 6]) -> (Vector3<f64>, Matrix3<f64>) {
    let e = tensor(s).symmetric_eigen();
    let mut order = [0, 1, 2];
    order.sort_by(|&a, &b| e.eigenvalues[b].total_cmp(&e.eigenvalues[a]));
    (
        Vector3::from_fn(|i, _| e.eigenvalues[order[i]]),
        Matrix3::from_fn(|i, j| e.eigenvectors[(i, order[j])]),
    )
}

impl Parameters {
    pub fn validate(&self) -> Result<()> {
        elastic::constitutive(self.young_kpa, self.poisson)?;
        if ![self.cohesion_kpa, self.friction_deg, self.dilation_deg]
            .iter()
            .all(|x| x.is_finite())
            || self.cohesion_kpa < 0.
            || !(0. ..89.).contains(&self.friction_deg)
            || !(0. ..=self.friction_deg).contains(&self.dilation_deg)
            || (self.flow == Flow::Associated
                && (self.dilation_deg - self.friction_deg).abs() > 1e-12)
            || self
                .tension_cutoff_kpa
                .is_some_and(|v| !v.is_finite() || v < 0.)
        {
            return Err(Error::Input("MC requires c>=0, 0<=psi<=phi<89 degrees; associated requires psi=phi; finite nonnegative tension cutoff".into()));
        }
        let (c, phi, _) = self.effective_strength();
        if self
            .tension_cutoff_kpa
            .is_some_and(|v| phi > 0. && v > c / phi.tan() + 1e-10)
        {
            return Err(Error::Input("MC tension cutoff exceeds effective c*cot(phi) apex; choose an explicit lower value".into()));
        }
        Ok(())
    }
    /// Davis changes the yield strength AND uses associated flow. Input values
    /// remain recorded; it is never silently enabled for nonassociated MC.
    pub fn effective_strength(&self) -> (f64, f64, f64) {
        let phi = self.friction_deg.to_radians();
        let psi = self.dilation_deg.to_radians();
        if self.flow == Flow::Davis {
            let beta = phi.cos() * psi.cos() / (1. - phi.sin() * psi.sin());
            let p = (beta * phi.tan()).atan();
            (beta * self.cohesion_kpa, p, p)
        } else {
            (self.cohesion_kpa, phi, psi)
        }
    }
    fn surfaces(&self) -> Vec<Surface> {
        let (c, p, g) = self.effective_strength();
        let mut out = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    let mut n = Vector3::zeros();
                    let mut m = n;
                    n[i] = 1. - p.sin();
                    n[j] = -1. - p.sin();
                    m[i] = 1. - g.sin();
                    m[j] = -1. - g.sin();
                    out.push(Surface {
                        n,
                        m,
                        c: 2. * c * p.cos(),
                    });
                }
            }
        }
        if let Some(t) = self.tension_cutoff_kpa {
            for i in 0..3 {
                let mut n = Vector3::zeros();
                n[i] = -1.;
                out.push(Surface { n, m: n, c: t });
            }
        }
        out
    }
    pub fn normalized_yield(&self, state: &State) -> f64 {
        let (s, _) = principal(&state.effective_stress_compression_positive_kpa);
        let scale = s.amax().max(self.cohesion_kpa).max(1.);
        self.surfaces()
            .iter()
            .map(|f| (f.n.dot(&s) - f.c) / scale)
            .fold(f64::NEG_INFINITY, f64::max)
    }
    pub fn validate_state(&self, s: &State) -> Result<()> {
        self.validate()?;
        if !s
            .effective_stress_compression_positive_kpa
            .iter()
            .chain(&s.plastic_strain_compression_positive)
            .chain(
                [
                    s.void_ratio,
                    s.equivalent_plastic_strain,
                    s.plastic_work_kpa,
                ]
                .iter(),
            )
            .all(|v| v.is_finite())
            || s.void_ratio <= 0.
            || s.equivalent_plastic_strain < 0.
            || s.plastic_work_kpa < -1e-9
            || self.normalized_yield(s) > 1e-8
        {
            return Err(Error::Material(
                "inadmissible/nonfinite MC committed state".into(),
            ));
        }
        Ok(())
    }
    pub fn moduli(&self) -> (f64, f64) {
        (
            self.young_kpa / (3. * (1. - 2. * self.poisson)),
            self.young_kpa / (2. * (1. + self.poisson)),
        )
    }
    fn substep(&self, old: &State, de: [f64; 6]) -> Result<State> {
        let d = elastic::constitutive(self.young_kpa, self.poisson)?;
        let trial = SVector::from_row_slice(&old.effective_stress_compression_positive_kpa)
            + d * SVector::from_row_slice(&de);
        let (s, v) = principal(&std::array::from_fn(|i| trial[i]));
        let surfaces = self.surfaces();
        let scale = s.amax().max(self.cohesion_kpa).max(1.);
        let tol = 1e-10 * scale;
        let f: Vec<_> = surfaces.iter().map(|a| a.n.dot(&s) - a.c).collect();
        let mut out = old.clone();
        out.last_return = Return::Elastic;
        let total_vol = de[..3].iter().sum::<f64>();
        out.void_ratio = (1. + old.void_ratio) * (-total_vol).exp() - 1.;
        if f.iter().all(|&x| x <= tol) {
            out.effective_stress_compression_positive_kpa = std::array::from_fn(|i| trial[i]);
            return Ok(out);
        }
        let normal = Matrix3::from_fn(|i, j| d[(i, j)]);
        // Conic Caratheodory: a feasible 3D correction has a representation
        // with at most three generators, even when six planes meet at apex.
        for &mask in candidate_sets(self.tension_cutoff_kpa.is_some()) {
            let size = mask.count_ones() as usize;
            if size > 3 {
                continue;
            }
            let ids: Vec<_> = (0..surfaces.len())
                .filter(|&i| mask & (1 << i) != 0)
                .collect();
            let a = DMatrix::from_fn(size, size, |i, j| {
                surfaces[ids[i]].n.dot(&(normal * surfaces[ids[j]].m))
            });
            let rhs = DVector::from_fn(size, |i, _| f[ids[i]]);
            let Some(mult) = a.clone().lu().solve(&rhs) else {
                continue;
            };
            if (&a * &mult - &rhs).norm() > 1e-9 * scale
                || mult.iter().any(|x| !x.is_finite() || *x < -1e-12)
            {
                continue;
            }
            let dp = ids
                .iter()
                .enumerate()
                .fold(Vector3::zeros(), |sum, (i, &j)| {
                    sum + mult[i].max(0.) * surfaces[j].m
                });
            let corrected = s - normal * dp;
            if corrected[0] + tol < corrected[1]
                || corrected[1] + tol < corrected[2]
                || surfaces.iter().any(|g| g.n.dot(&corrected) - g.c > tol)
                || ids
                    .iter()
                    .any(|&i| (surfaces[i].n.dot(&corrected) - surfaces[i].c).abs() > tol)
            {
                continue;
            }
            let stress = v * Matrix3::from_diagonal(&corrected) * v.transpose();
            let plastic = v * Matrix3::from_diagonal(&dp) * v.transpose();
            let pv = voigt(&plastic);
            for (i, value) in pv.iter().enumerate() {
                out.plastic_strain_compression_positive[i] += value * if i < 3 { 1. } else { 2. };
            }
            out.effective_stress_compression_positive_kpa = voigt(&stress);
            let dev = plastic - Matrix3::identity() * plastic.trace() / 3.;
            out.equivalent_plastic_strain += (2. / 3. * dev.norm_squared()).sqrt();
            let work = corrected.dot(&dp);
            if work < -1e-9 * scale * dp.norm().max(1e-12) {
                return Err(Error::Material("negative MC plastic dissipation".into()));
            }
            out.plastic_work_kpa += work.max(0.);
            let shear = ids.iter().any(|&i| i < 6);
            let tension = ids.iter().any(|&i| i >= 6);
            out.last_return = if shear && tension {
                Return::Mixed
            } else if tension {
                match corrected
                    .iter()
                    .filter(|&&s| {
                        self.tension_cutoff_kpa
                            .is_some_and(|t| (s + t).abs() <= tol)
                    })
                    .count()
                {
                    1 => Return::TensionFace,
                    2 => Return::TensionEdge,
                    _ => Return::TensionApex,
                }
            } else if (corrected[0] - corrected[2]).abs() <= tol {
                Return::ShearApex
            } else if (corrected[0] - corrected[1]).abs() <= tol {
                Return::ExtensionEdge
            } else if (corrected[1] - corrected[2]).abs() <= tol {
                Return::CompressionEdge
            } else {
                Return::ShearFace
            };
            return Ok(out);
        }
        Err(Error::Material("MC return has no admissible nonnegative active-set solution (including zero-dilation pure-MC tension/apex paths); no stress clipping was applied".into()))
    }
    /// A parameter jump followed by the prescribed strain subdivision. The
    /// committed state must be admissible under `from`; only this explicit
    /// transition may begin outside the new yield surface. Elastic constants
    /// and flow family are held fixed; no history is cleared or relabelled.
    pub fn integrate_transition_partition(
        &self,
        from: &Self,
        old: &State,
        de: [f64; 6],
        steps: usize,
    ) -> Result<Update> {
        from.validate_state(old)?;
        self.validate()?;
        if self.young_kpa != from.young_kpa
            || self.poisson != from.poisson
            || self.flow != from.flow
        {
            return Err(Error::Input(
                "MC parameter transition requires unchanged elastic constants and flow family"
                    .into(),
            ));
        }
        if steps == 0 || steps > 100_000 || !de.iter().all(|v| v.is_finite()) {
            return Err(Error::Input(
                "invalid MC transition partition/strain".into(),
            ));
        }
        let mut state = old.clone();
        let mut plastic = 0;
        for _ in 0..steps {
            state = self.substep(&state, de.map(|v| v / steps as f64))?;
            self.validate_state(&state)?;
            plastic += usize::from(state.last_return != Return::Elastic);
        }
        Ok(Update {
            state,
            substeps: steps,
            plastic_substeps: plastic,
        })
    }
    /// Tangent of the same old-state/new-parameter map used for every trial.
    /// It never reuses a previous Newton perturbation's plastic history.
    pub fn stress_tangent_transition_partition(
        &self,
        from: &Self,
        old: &State,
        de: [f64; 6],
        steps: usize,
        h: f64,
    ) -> Result<(Update, elastic::Tensor)> {
        if !h.is_finite() || !(1e-9..=1e-5).contains(&h) {
            return Err(Error::Input(
                "invalid MC transition tangent perturbation".into(),
            ));
        }
        let update = self.integrate_transition_partition(from, old, de, steps)?;
        let mut d = elastic::Tensor::zeros();
        for j in 0..6 {
            let mut a = de;
            let mut b = de;
            a[j] += h;
            b[j] -= h;
            let a = self
                .integrate_transition_partition(from, old, a, steps)?
                .state;
            let b = self
                .integrate_transition_partition(from, old, b, steps)?
                .state;
            for i in 0..6 {
                d[(i, j)] = (a.effective_stress_compression_positive_kpa[i]
                    - b.effective_stress_compression_positive_kpa[i])
                    / (2. * h);
            }
        }
        if !d.iter().all(|v| v.is_finite()) {
            return Err(Error::Material(
                "nonfinite MC parameter transition tangent".into(),
            ));
        }
        Ok((update, d))
    }
    pub fn integrate_partition(&self, old: &State, de: [f64; 6], steps: usize) -> Result<Update> {
        self.validate_state(old)?;
        if steps == 0 || steps > 100_000 || !de.iter().all(|v| v.is_finite()) {
            return Err(Error::Input("invalid MC partition/strain".into()));
        }
        let mut state = old.clone();
        let mut plastic = 0;
        for _ in 0..steps {
            state = self.substep(&state, de.map(|v| v / steps as f64))?;
            plastic += usize::from(state.last_return != Return::Elastic);
        }
        self.validate_state(&state)?;
        Ok(Update {
            state,
            substeps: steps,
            plastic_substeps: plastic,
        })
    }
    pub fn integrate(&self, old: &State, de: [f64; 6], maximum: f64) -> Result<Update> {
        if !maximum.is_finite() || !(1e-7..=1e-3).contains(&maximum) {
            return Err(Error::Input("invalid MC substep size".into()));
        }
        let steps = (de.iter().map(|x| x.abs()).fold(0., f64::max) / maximum)
            .ceil()
            .max(1.) as usize;
        self.integrate_partition(old, de, steps)
    }
    pub fn stress_tangent_partition(
        &self,
        old: &State,
        de: [f64; 6],
        steps: usize,
        h: f64,
    ) -> Result<(Update, elastic::Tensor)> {
        if !h.is_finite() || !(1e-9..=1e-5).contains(&h) {
            return Err(Error::Input("invalid MC tangent perturbation".into()));
        }
        let update = self.integrate_partition(old, de, steps)?;
        let mut d = elastic::Tensor::zeros();
        for j in 0..6 {
            let mut a = de;
            let mut b = de;
            a[j] += h;
            b[j] -= h;
            let a = self.integrate_partition(old, a, steps)?.state;
            let b = self.integrate_partition(old, b, steps)?.state;
            for i in 0..6 {
                d[(i, j)] = (a.effective_stress_compression_positive_kpa[i]
                    - b.effective_stress_compression_positive_kpa[i])
                    / (2. * h);
            }
        }
        if !d.iter().all(|v| v.is_finite()) {
            return Err(Error::Material("nonfinite MC tangent".into()));
        }
        Ok((update, d))
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Path {
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
pub fn run_path(p: Path) -> Result<PathOutput> {
    if p.schema != "soil-fem-mc-path/1" || p.strain_increments_compression_positive.len() > 10000 {
        return Err(Error::Input("invalid MC path schema/length".into()));
    }
    p.parameters.validate_state(&p.initial_state)?;
    let mut out = PathOutput {
        schema: "soil-fem-mc-path-result/1",
        parameters: p.parameters.clone(),
        states: vec![p.initial_state],
        substeps: vec![],
        plastic_substeps: vec![],
    };
    for de in p.strain_increments_compression_positive {
        let u = p.parameters.integrate(
            out.states
                .last()
                .ok_or_else(|| Error::Material("empty MC path".into()))?,
            de,
            p.maximum_substep,
        )?;
        out.states.push(u.state);
        out.substeps.push(u.substeps);
        out.plastic_substeps.push(u.plastic_substeps);
    }
    Ok(out)
}
