//! Coupled high-order wet-column study. Input contains the elastic Schur
//! complement and the two load vectors. All active interface history,
//! tractions, tangents and equilibrium increments are solved here.
use crate::{
    Error, Result,
    joint_law::{Material, PointState},
};
use nalgebra::{DMatrix, DVector, Vector3};
use serde::Deserialize;
use serde_json::json;
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    elements: usize,
    stiffness: Vec<Vec<f64>>,
    rhs_load: Vec<f64>,
    rhs_boundary: Vec<f64>,
    load: Vec<f64>,
    initial: Vec<f64>,
    trace: Vec<Vec<f64>>,
    observation: Vec<Vec<f64>>,
    area: Vec<f64>,
    normal: Vec<[f64; 3]>,
    observation_normal: Vec<[f64; 3]>,
    material: Material,
    targets: Vec<f64>,
    grids: Vec<Vec<f64>>,
}

type Trace = Vec<Vec<(usize, f64)>>;

fn trace(rows: &[Vec<f64>]) -> Trace {
    rows.iter()
        .map(|row| {
            row.iter()
                .copied()
                .enumerate()
                .filter(|(_, w)| *w != 0.)
                .collect()
        })
        .collect()
}

fn jumps(rows: &Trace, x: &DVector<f64>) -> Vec<[f64; 3]> {
    rows.iter()
        .map(|row| std::array::from_fn(|d| row.iter().map(|(i, w)| w * x[3 * i + d]).sum()))
        .collect()
}

fn states(normals: &[[f64; 3]]) -> Result<Vec<PointState>> {
    normals.iter().map(|n| PointState::initial(*n)).collect()
}

/// Run the seven-stage high-order column on supplied, statically condensed elastic operators.
pub fn run(operators: &Path, output: &Path) -> Result<()> {
    let data: Input = serde_json::from_slice(&fs::read(operators)?)?;
    let nd = data.rhs_load.len();
    let limit = data
        .elements
        .checked_mul(2065)
        .ok_or_else(|| Error::Input("invalid element count".into()))?
        .div_ceil(12);
    if nd == 0
        || !nd.is_multiple_of(3)
        || data.stiffness.len() != nd
        || data.stiffness.iter().any(|r| r.len() != nd)
        || data.rhs_boundary.len() != nd
        || data.load.len() != data.initial.len()
        || data.load.len() < nd
        || data.area.is_empty()
        || data.area.len() != data.trace.len()
        || data.normal.len() != data.trace.len()
        || data.observation.is_empty()
        || data.observation.len() != data.observation_normal.len()
        || data.area.iter().any(|&w| w <= 0.)
        || data
            .trace
            .iter()
            .chain(&data.observation)
            .any(|r| r.len() != nd / 3)
        || data.targets.len() != 7
        || data.grids.len() != 7
        || data.targets[0] != 0.
        || data.grids.iter().any(|g| {
            g.last() != Some(&1.)
                || g[0] <= 0.
                || g.windows(2).any(|v| v[1] <= v[0])
                || g.len() > limit
        })
    {
        return Err(Error::Input(
            "inconsistent wet-column operators, area weights or stage grids".into(),
        ));
    }
    data.material.validate()?;
    let mut output = BufWriter::new(
        fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(output)?,
    );
    let stiff = DMatrix::from_fn(nd, nd, |i, j| data.stiffness[i][j]);
    let load = DVector::from_vec(data.rhs_load);
    let bc = DVector::from_vec(data.rhs_boundary);
    let active = trace(&data.trace);
    let passive = trace(&data.observation);
    let mut old = states(&data.normal)?;
    let mut probe = states(&data.observation_normal)?;
    let mut x = DVector::<f64>::zeros(nd);
    let mut origin = 0.;
    let law = &data.material;
    for (phase, (&target, grid)) in data.targets.iter().zip(&data.grids).enumerate() {
        for &t in grid {
            let factor = origin + (target - origin) * t;
            let fraction = if phase == 0 { t } else { 1. };
            let rhs = fraction * &load + factor * &bc;
            let scale = data
                .load
                .iter()
                .zip(&data.initial)
                .map(|(l, i)| (fraction * l + (1. - fraction) * i).powi(2))
                .sum::<f64>()
                .sqrt()
                .max(1.);
            let mut counts = vec![1; old.len()];
            let evaluate = |x: &DVector<f64>, counts: &[usize]| -> Result<_> {
                let mut residual = &stiff * x - &rhs;
                let mut new = Vec::with_capacity(old.len());
                for (q, jump) in jumps(&active, x).into_iter().enumerate() {
                    let state = law.integrate_partition(&old[q], jump, counts[q])?;
                    for &(i, w) in &active[q] {
                        for d in 0..3 {
                            residual[3 * i + d] += data.area[q] * w * state.traction_kpa[d];
                        }
                    }
                    new.push(state);
                }
                Ok((residual, new))
            };
            let mut max_linear = 0.0_f64;
            let (new, rel, nit) = 'newton: {
                for it in 0..60 {
                    let (res, new) = evaluate(&x, &counts)?;
                    let norm = res.norm();
                    if norm < 1e-10 * scale {
                        let mut changed = false;
                        for (q, s) in new.iter().enumerate() {
                            let needed = law.partition(&old[q], s.jump_m, 5e-5)?;
                            if needed > counts[q] {
                                counts[q] = needed;
                                changed = true;
                            }
                        }
                        if changed {
                            continue;
                        }
                        break 'newton (new, norm / scale, it);
                    }
                    let mut jac = stiff.clone();
                    for (q, row) in active.iter().enumerate() {
                        let d =
                            law.tangent(&old[q], new[q].jump_m, counts[q], 2e-9)? * data.area[q];
                        for &(j, b) in row {
                            let block = b * d;
                            for &(i, a) in row {
                                for col in 0..3 {
                                    for r in 0..3 {
                                        jac[(3 * i + r, 3 * j + col)] += a * block[(r, col)];
                                    }
                                }
                            }
                        }
                    }
                    let dx = jac
                        .clone()
                        .lu()
                        .solve(&(-&res))
                        .ok_or_else(|| Error::Solve("singular JOINT tangent".into()))?;
                    let lr = (&jac * &dx + &res).norm() / norm;
                    if !lr.is_finite() || lr > 1e-10 {
                        return Err(Error::Solve(format!("linear residual {lr} at {phase} {t}")));
                    }
                    max_linear = max_linear.max(lr);
                    let mut improved = false;
                    for power in 0..16 {
                        let trial = &x + &dx * 0.5_f64.powi(power);
                        if evaluate(&trial, &counts)?.0.norm() < norm {
                            x = trial;
                            improved = true;
                            break;
                        }
                    }
                    if !improved {
                        return Err(Error::Solve(format!("line search at {phase} {t}")));
                    }
                }
                return Err(Error::Solve(format!("Newton cap at {phase} {t}")));
            };
            old = new;
            for (s, jump) in probe.iter_mut().zip(jumps(&passive, &x)) {
                *s = law.integrate(s, jump, 5e-5)?;
            }
            let mut row = json!({"phase":phase,"t":t,"x":x.as_slice(),"relative_force":rel,
                                 "linear":max_linear,"iterations":nit,"traction":old.iter().map(|s| s.traction_kpa).collect::<Vec<_>>()});
            if (t * 64.).fract() == 0. {
                let integral = old
                    .iter()
                    .zip(&data.area)
                    .fold(Vector3::zeros(), |sum, (s, &w)| {
                        sum + w * Vector3::from(s.plastic_slip_m)
                    });
                row["integral"] = json!(integral.as_slice());
                row["slip"] = json!(probe.iter().map(|s| s.plastic_slip_m).collect::<Vec<_>>());
            }
            serde_json::to_writer(&mut output, &row)?;
            writeln!(output)?;
        }
        origin = target;
    }
    output.flush()?;
    Ok(())
}
