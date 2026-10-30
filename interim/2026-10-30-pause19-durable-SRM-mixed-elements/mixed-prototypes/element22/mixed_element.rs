//! Isolated Tet10 u / P1 mean-p / independent or P1 volume element.
//! No global construction solver or checkpoint/SRM claim is provided here.
use crate::soil_split;
use soil_fem_3d::{Error, Result, elements, soil};
#[derive(Clone, Copy)]
pub enum Volume {
    Independent,
    Projected,
}
pub struct Element {
    pub points: Vec<elements::VolumePoint>,
    pub nodes: Vec<[f64; 3]>,
    pub volume: Volume,
    pub nt: usize,
}
pub struct Evaluation {
    pub residual: Vec<f64>,
    pub tangent: Option<Vec<Vec<f64>>>,
    pub actual_internal: Vec<f64>,
    pub mean: Vec<f64>,
    pub maximum_pressure_gap: f64,
    pub states: Vec<soil::State>,
}
impl Element {
    pub fn new(nodes: Vec<[f64; 3]>, volume: Volume) -> Result<Self> {
        if nodes.len() != 10 {
            return Err(Error::Input("mixed element needs Tet10 nodes".into()));
        }
        let points = elements::volume_points(&nodes, 3)?;
        let nt = match volume {
            Volume::Independent => points.len(),
            Volume::Projected => 4,
        };
        Ok(Self {
            points,
            nodes,
            volume,
            nt,
        })
    }
    pub fn count(&self) -> usize {
        self.nt + 34
    }
    pub fn theta_rows(&self, qi: usize, n: [f64; 4]) -> Vec<(usize, f64)> {
        match self.volume {
            Volume::Independent => vec![(qi, 1.)],
            Volume::Projected => n.into_iter().enumerate().collect(),
        }
    }
    pub fn evaluate(
        &self,
        model: &soil::Model,
        old: &[soil::State],
        x: &[f64],
        tangent: bool,
    ) -> Result<Evaluation> {
        if old.len() != self.points.len()
            || x.len() != self.count()
            || !x.iter().all(|v| v.is_finite())
        {
            return Err(Error::Input(
                "mixed element input inventory/value differs".into(),
            ));
        }
        let mut out = Evaluation {
            residual: vec![0.; self.count()],
            tangent: tangent.then(|| vec![vec![0.; self.count()]; self.count()]),
            actual_internal: vec![0.; 30],
            mean: vec![],
            maximum_pressure_gap: 0.,
            states: Vec::new(),
        };
        for (qi, q) in self.points.iter().enumerate() {
            let n = [
                1. - q.natural.iter().sum::<f64>(),
                q.natural[0],
                q.natural[1],
                q.natural[2],
            ];
            let rows = self.theta_rows(qi, n);
            let theta = rows.iter().map(|&(i, w)| w * x[i]).sum::<f64>();
            let de: [f64; 6] = std::array::from_fn(|a| {
                -(0..30).map(|i| q.b[(a, i)] * x[self.nt + i]).sum::<f64>()
            });
            let div = de[..3].iter().sum::<f64>();
            let mut dev = std::array::from_fn(|i| de[i] - if i < 3 { div / 3. } else { 0. });
            dev[2] = -(dev[0] + dev[1]);
            let (values, t) = if tangent {
                let method = if matches!(model, soil::Model::Mcc(_)) {
                    soil_split::Method::ImplicitMcc
                } else {
                    soil_split::Method::FiniteDifference(1e-7)
                };
                let r = soil_split::evaluate(model, &old[qi], dev, theta, 1, method)?;
                (r.values, Some(r.tangent))
            } else {
                (soil_split::integrate(model, &old[qi], dev, theta, 1)?, None)
            };
            let mean = values.mean_pressure_compression_positive_kpa;
            let pressure = (0..4).map(|a| n[a] * x[self.nt + 30 + a]).sum::<f64>();
            let stress: [f64; 6] = std::array::from_fn(|a| {
                values.deviatoric_stress_compression_positive_kpa[a]
                    + if a < 3 { pressure } else { 0. }
            });
            let actual = *values.update.state.stress();
            let w = q.weight_m3;
            out.mean.push(mean);
            out.maximum_pressure_gap = out.maximum_pressure_gap.max((mean - pressure).abs());
            for i in 0..30 {
                out.residual[self.nt + i] -=
                    w * (0..6).map(|a| q.b[(a, i)] * stress[a]).sum::<f64>();
                out.actual_internal[i] -= w * (0..6).map(|a| q.b[(a, i)] * actual[a]).sum::<f64>();
            }
            for &(ti, wi) in &rows {
                out.residual[ti] += w * wi * (mean - pressure);
            }
            for a in 0..4 {
                out.residual[self.nt + 30 + a] -= w * n[a] * (-div + theta);
            }
            if let (Some(matrix), Some(t)) = (&mut out.tangent, t) {
                for &(ti, wi) in &rows {
                    for &(tj, wj) in &rows {
                        matrix[ti][tj] += w * wi * wj * t.mean_by_volume;
                    }
                }
                for i in 0..30 {
                    let bi: [f64; 6] = std::array::from_fn(|a| q.b[(a, i)]);
                    let vi = bi[..3].iter().sum::<f64>();
                    let ui = self.nt + i;
                    for &(ti, wi) in &rows {
                        matrix[ui][ti] -= w
                            * wi
                            * (0..6)
                                .map(|a| bi[a] * t.deviatoric_by_volume[a])
                                .sum::<f64>();
                        matrix[ti][ui] -=
                            w * wi * (0..6).map(|a| t.mean_by_deviatoric[a] * bi[a]).sum::<f64>();
                    }
                    for k in 0..30 {
                        let mut value = 0.;
                        for a in 0..6 {
                            for b in 0..6 {
                                value += bi[a] * t.deviatoric_by_deviatoric[a][b] * q.b[(b, k)];
                            }
                        }
                        matrix[ui][self.nt + k] += w * value;
                    }
                    for a in 0..4 {
                        let pr = self.nt + 30 + a;
                        matrix[ui][pr] -= w * vi * n[a];
                        matrix[pr][ui] -= w * n[a] * vi;
                    }
                }
                for a in 0..4 {
                    let pr = self.nt + 30 + a;
                    for &(ti, wi) in &rows {
                        matrix[ti][pr] -= w * wi * n[a];
                        matrix[pr][ti] -= w * n[a] * wi;
                    }
                }
            }
            out.states.push(values.update.state);
        }
        Ok(out)
    }
    pub fn mass(&self) -> Vec<Vec<f64>> {
        let mut matrix = vec![vec![0.; 4]; 4];
        for q in &self.points {
            let n = [
                1. - q.natural.iter().sum::<f64>(),
                q.natural[0],
                q.natural[1],
                q.natural[2],
            ];
            for i in 0..4 {
                for j in 0..4 {
                    matrix[i][j] += q.weight_m3 * n[i] * n[j];
                }
            }
        }
        matrix
    }
    pub fn project(&self, values: &[f64]) -> Vec<f64> {
        assert_eq!(values.len(), self.points.len());
        let mut rhs = vec![0.; 4];
        for (q, &value) in self.points.iter().zip(values) {
            let n = [
                1. - q.natural.iter().sum::<f64>(),
                q.natural[0],
                q.natural[1],
                q.natural[2],
            ];
            for i in 0..4 {
                rhs[i] += q.weight_m3 * n[i] * value;
            }
        }
        solve(self.mass(), rhs)
    }
}
pub fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    for i in 0..b.len() {
        let pivot = (i..b.len())
            .max_by(|&j, &k| a[j][i].abs().total_cmp(&a[k][i].abs()))
            .unwrap();
        a.swap(i, pivot);
        b.swap(i, pivot);
        assert!(a[i][i].abs() > 1e-16);
        let d = a[i][i];
        for j in i..b.len() {
            a[i][j] /= d;
        }
        b[i] /= d;
        for k in 0..b.len() {
            if k != i {
                let d = a[k][i];
                for j in i..b.len() {
                    a[k][j] -= d * a[i][j];
                }
                b[k] -= d * b[i];
            }
        }
    }
    b
}
