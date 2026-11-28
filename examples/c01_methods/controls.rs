use serde::Deserialize;
use soil_fem_3d::{Error, Result, stage_control::ArcSample};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

#[derive(Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub seed: Option<PathBuf>,
    pub full_field_grid: Option<std::num::NonZeroU32>,
    pub dissipation_history: Option<PathBuf>,
    pub arc_scale_m: Option<f64>,
}
type Snapshot = (f64, Vec<[f64; 3]>);
pub struct Seed(Vec<Vec<Snapshot>>);

impl Seed {
    pub fn read(path: &Path, nodes: usize) -> Result<Self> {
        let mut stages = vec![Vec::new(); 7];
        stages[0].push((0., vec![[0.; 3]; nodes]));
        let mut phase = 0;
        for line in BufReader::new(File::open(path)?).lines() {
            let (s, t, u): (usize, f64, Vec<[f64; 3]>) = serde_json::from_str(&line?)?;
            if s >= 7
                || !(0. < t && t <= 1.)
                || u.len() != nodes
                || !u.iter().flatten().all(|v| v.is_finite())
            {
                return Err(Error::Input("invalid C-01 coarse history".into()));
            }
            if s != phase {
                let (time, last) = stages[phase].last().unwrap();
                if s != phase + 1 || *time != 1. {
                    return Err(Error::Input("missing coarse stage endpoint".into()));
                }
                let origin = last.clone();
                stages[s].push((0., origin));
                phase = s;
            }
            if t <= stages[s].last().unwrap().0 {
                return Err(Error::Input("unordered coarse history".into()));
            }
            stages[s].push((t, u));
        }
        if phase != 6 || stages[phase].last().unwrap().0 != 1. {
            return Err(Error::Input("incomplete coarse history".into()));
        }
        Ok(Self(stages))
    }
    pub fn interpolate(&self, phase: usize, t: f64) -> Result<Vec<[f64; 3]>> {
        if phase >= 7 || !(0. ..=1.).contains(&t) {
            return Err(Error::Input("invalid seed coordinate".into()));
        }
        let rows = &self.0[phase];
        let i = rows.partition_point(|x| x.0 < t);
        if rows[i].0 == t {
            return Ok(rows[i].1.clone());
        }
        let (a, left) = &rows[i - 1];
        let (b, right) = &rows[i];
        let f = (t - a) / (b - a);
        Ok(left
            .iter()
            .zip(right)
            .map(|(l, r)| std::array::from_fn(|k| l[k] + f * (r[k] - l[k])))
            .collect())
    }
}

#[derive(Clone, Copy)]
struct Bracket {
    stage: usize,
    start: f64,
    expected: f64,
    lo: f64,
    low: f64,
    hi: f64,
    high: f64,
}

pub struct Sphere {
    scale: f64,
    pending: Option<Bracket>,
    pub reviews: usize,
}
impl Sphere {
    pub fn new(scale: f64) -> Result<Self> {
        if !scale.is_finite() || scale <= 0. {
            return Err(Error::Input("positive displacement scale required".into()));
        }
        Ok(Self {
            scale,
            pending: None,
            reviews: 0,
        })
    }
    pub fn retry(&mut self, s: ArcSample<'_>) -> Option<f64> {
        if s.stage == 0 || s.path_rejected {
            self.pending = None;
            return None;
        }
        self.reviews += 1;
        let sum: f64 = s
            .trial_u
            .iter()
            .flatten()
            .zip(s.old_u.iter().flatten())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        let rms = (sum / (3 * s.old_u.len()) as f64).sqrt() / self.scale;
        let length = (s.end - s.start).hypot(rms);
        let mut b = match self.pending.take() {
            Some(b) if (b.stage, b.start, b.expected) == (s.stage, s.start, s.end) => b,
            _ if length <= s.base => return None, // original error/boundary limit lies inside the sphere
            _ => Bracket {
                stage: s.stage,
                start: s.start,
                expected: s.end,
                lo: s.start,
                low: 0.,
                hi: s.end,
                high: length,
            },
        };
        if length == s.base {
            return None;
        }
        if length < s.base {
            b.lo = s.end;
            b.low = length;
        } else {
            b.hi = s.end;
            b.high = length;
        }
        let estimate = b.lo + (b.hi - b.lo) * (s.base - b.low) / (b.high - b.low);
        if length < s.base && (estimate <= s.end || b.lo.next_up() >= b.hi) {
            return None;
        }
        b.expected = if b.lo.next_up() >= b.hi {
            b.lo
        } else {
            estimate.clamp(b.lo.next_up(), b.hi.next_down())
        };
        self.pending = Some(b);
        Some(b.expected - s.start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sphere_counts_an_outside_trial_before_accepting_the_root() {
        let old = [[0.; 3]];
        let mut sphere = Sphere::new(1.).unwrap();
        let trial = [[0.75; 3]];
        let next = sphere
            .retry(ArcSample {
                stage: 1,
                start: 0.,
                end: 1.,
                base: 1.,
                old_u: &old,
                trial_u: &trial,
                path_rejected: false,
            })
            .unwrap();
        assert_eq!(next, 0.8);
        let trial = [[0.6; 3]];
        assert!(
            sphere
                .retry(ArcSample {
                    stage: 1,
                    start: 0.,
                    end: next,
                    base: 1.,
                    old_u: &old,
                    trial_u: &trial,
                    path_rejected: false
                })
                .is_none()
        );
    }
}
