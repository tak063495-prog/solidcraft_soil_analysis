//! Opt-in replay of committed MCC material substeps from the C-03 native probe.
use serde::{Deserialize, Serialize};
use soil_fem_3d::{elements::PointInfo, geostatic::Seed, ground::Ground, mcc, mesh::Mesh, soil};
use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

const PHASES: [&str; 7] = [
    "geostatic",
    "shear",
    "open",
    "recontact",
    "reverse",
    "reload",
    "return",
];

#[derive(Deserialize)]
struct Run {
    ground: Ground,
    complete: bool,
    integration_points: Vec<PointInfo>,
}

#[derive(Deserialize)]
struct Increment {
    stage: usize,
    progress: f64,
    strain_increment: Vec<[f64; 6]>,
    partitions: Vec<usize>,
    relative_force_residual: f64,
}

#[derive(Deserialize)]
struct Accepted {
    stage_index: usize,
    stage_progress: f64,
    material_states: Vec<mcc::State>,
}

struct Point<'a> {
    model: &'a mcc::Parameters,
    state: mcc::State,
    weight: f64,
}

#[derive(Default, Serialize)]
struct WorkRow {
    stage: &'static str,
    progress: f64,
    relative_force_residual: f64,
    work_kn_m: [f64; 3],
    substeps: usize,
    plastic_substeps: usize,
    max_local_residual_kpa: f64,
    sum_abs_residual_kn_m: f64,
}

fn work(stress: &[f64; 6], dep: [f64; 6]) -> [f64; 3] {
    let tensor = stress.iter().zip(dep).map(|(s, e)| s * e).sum();
    let (p, q) = mcc::invariants(stress);
    let ev: f64 = dep[..3].iter().sum();
    let normal_dev: f64 = dep[..3].iter().map(|e| (e - ev / 3.).powi(2)).sum();
    let engineering_shear: f64 = dep[3..].iter().map(|g| g * g).sum();
    let es = ((2. / 3.) * (normal_dev + 0.5 * engineering_shear)).sqrt();
    [tensor, p * ev, q * es]
}

impl Point<'_> {
    fn advance(
        &mut self,
        strain: [f64; 6],
        count: usize,
        row: &mut WorkRow,
    ) -> Result<(), Box<dyn Error>> {
        if !(1..=100_000).contains(&count) {
            return Err("invalid recorded material partition".into());
        }
        let de = strain.map(|v| v / count as f64);
        for _ in 0..count {
            let next = self.model.integrate_partition(&self.state, de, 1)?.state;
            let dep = std::array::from_fn(|i| {
                next.plastic_strain_compression_positive[i]
                    - self.state.plastic_strain_compression_positive[i]
            });
            let terms = work(&next.effective_stress_compression_positive_kpa, dep);
            let discrepancy = (terms[0] - (terms[1] + terms[2])).abs();
            row.max_local_residual_kpa = row.max_local_residual_kpa.max(discrepancy);
            row.sum_abs_residual_kn_m += self.weight * discrepancy;
            for (total, term) in row.work_kn_m.iter_mut().zip(terms) {
                *total += self.weight * term;
            }
            row.substeps += 1;
            row.plastic_substeps += usize::from(dep.iter().any(|x| *x != 0.));
            self.state = next;
        }
        Ok(())
    }
}

fn point_layout<'a>(run: &'a Run, mesh: &Mesh) -> Result<Vec<Point<'a>>, Box<dyn Error>> {
    let config = run
        .ground
        .initialization
        .geostatic
        .as_ref()
        .ok_or("missing geostatic input")?;
    let Seed::BodyStates { states } = &config.seed else {
        return Err("body seeds required".into());
    };
    let mut points = Vec::new();
    for (i, info) in run.integration_points.iter().enumerate() {
        let e = info.element_index;
        if e != i / 27 || info.local_index != i % 27 || info.weight_m3 <= 0. {
            return Err("unexpected Tet10 integration point layout".into());
        }
        let region = mesh
            .regions
            .iter()
            .find(|r| r.id == mesh.region_ids[e])
            .ok_or("unknown region")?;
        let mat = run
            .ground
            .materials
            .iter()
            .find(|m| m.id == mesh.material_ids[e])
            .ok_or("unknown material")?;
        let (soil::Model::Mcc(model), Some(soil::State::Mcc(seed))) =
            (&mat.model, states.get(&region.body))
        else {
            return Err("MCC material and seed required".into());
        };
        if model.integration != mcc::Integration::BackwardEuler {
            return Err("backward Euler required".into());
        }
        points.push(Point {
            model,
            state: seed.clone(),
            weight: info.weight_m3,
        });
    }
    Ok(points)
}

fn same_state(a: &mcc::State, b: &mcc::State) -> bool {
    a.effective_stress_compression_positive_kpa == b.effective_stress_compression_positive_kpa
        && a.plastic_strain_compression_positive == b.plastic_strain_compression_positive
        && a.void_ratio == b.void_ratio
        && a.preconsolidation_pressure_kpa == b.preconsolidation_pressure_kpa
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let [mesh, raw, output] = args.as_slice() else {
        return Err("usage: MESH RAW NEW_WORK_JSONL".into());
    };
    let mesh: Mesh = serde_json::from_reader(BufReader::new(File::open(mesh)?))?;
    let run: Run = serde_json::from_reader(BufReader::new(File::open(raw.join("run.json"))?))?;
    if !run.complete || mesh.tetrahedra10.len() != 12 || run.integration_points.len() != 324 {
        return Err("complete 12 Tet10 history required".into());
    }
    let mut points = point_layout(&run, &mesh)?;
    let mut diagnostics = BufReader::new(File::open(raw.join("constitutive.jsonl"))?).lines();
    let mut accepted = BufReader::new(File::open(raw.join("accepted.jsonl"))?).lines();
    let mut sink = BufWriter::new(File::create_new(output)?);
    for (phase, &stage) in PHASES.iter().enumerate() {
        for step in 1..=32 {
            let progress = f64::from(step) / 32.;
            let d: Increment =
                serde_json::from_str(&diagnostics.next().ok_or("missing material partitions")??)?;
            let a: Accepted =
                serde_json::from_str(&accepted.next().ok_or("missing accepted state")??)?;
            // The committed geostatic end is stored at the start coordinate of the following stage.
            let stored = if phase == 0 && step == 32 {
                0.
            } else {
                progress
            };
            if (d.stage, d.progress) != (phase, progress)
                || (a.stage_index, a.stage_progress) != (phase, stored)
                || d.strain_increment.len() != points.len()
                || a.material_states.len() != points.len()
                || d.partitions.len() != points.len() + 32
            {
                return Err("accepted coordinate or GP inventory differs".into());
            }
            let mut row = WorkRow {
                stage,
                progress,
                relative_force_residual: d.relative_force_residual,
                ..Default::default()
            };
            for (i, point) in points.iter_mut().enumerate() {
                point.advance(d.strain_increment[i], d.partitions[i], &mut row)?;
                if !same_state(&point.state, &a.material_states[i]) {
                    return Err(format!("substep replay differs from the committed state at {stage} {progress} GP {i}").into());
                }
            }
            serde_json::to_writer(&mut sink, &row)?;
            writeln!(sink)?;
        }
    }
    if diagnostics.next().is_some() || accepted.next().is_some() {
        return Err("unexpected trailing history".into());
    }
    sink.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::work;
    #[test]
    fn hydrostatic_and_deviatoric_work_use_engineering_shear() {
        assert_eq!(
            work(&[8., 8., 8., 0., 0., 0.], [0.25, 0.25, 0.25, 0., 0., 0.]),
            [6., 6., 0.]
        );
        let result = work(&[2., -1., -1., 0., 0., 0.], [2., -1., -1., 0., 0., 0.]);
        assert_eq!(result, [6., 0., 6.]);
        let shear = work(&[0., 0., 0., 4., 0., 0.], [0., 0., 0., 0.5, 0., 0.]);
        assert_eq!(shear[0], 2.);
        let rotated = work(&[0., 0., 0., 4., 0., 0.], [0., 0., 0., 0., 0.5, 0.]);
        assert_eq!(rotated[0], 0.);
        assert!(rotated[2] > 0.);
    }
}
