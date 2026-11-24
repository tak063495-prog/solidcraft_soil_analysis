//! Drained homogeneous triaxial specimen, using the production material update.
use nalgebra::{SVector, Vector3};
use serde_json::{Value, json};
use soil_fem_3d::{
    elastic,
    mcc::{Integration, Parameters, State, Update, invariants},
};
use std::collections::BTreeMap;
use std::error::Error;

type Trial<T> = Result<T, Box<dyn Error>>;

const SPEC: &str = include_str!("../../interim/2026-11-24-t9d-benchmarks/criteria.json");

#[derive(Clone, Copy, PartialEq)]
enum ElasticPath {
    Default,
    Analytic,
}

struct MaterialStep<'a> {
    model: &'a Parameters,
    h: f64,
    shear: Option<f64>,
    path: ElasticPath,
}

impl MaterialStep<'_> {
    fn integrate(&self, state: &State, strain: [f64; 6]) -> Trial<Update> {
        Ok(match self.path {
            ElasticPath::Default => self
                .model
                .integrate_with_shear(state, strain, self.h, self.shear)?,
            ElasticPath::Analytic => self.model.integrate_with_elastic_split(
                state,
                strain,
                self.h,
                self.shear.ok_or("analytic path requires constant G")?,
            )?,
        })
    }
}

struct Cell {
    ids: [usize; 4],
    b: elastic::BMatrix,
    vol: f64,
    state: State,
}

struct Cube {
    nodes: Vec<[f64; 3]>,
    cells: Vec<Cell>,
    faces: Vec<([usize; 3], Vector3<f64>)>,
}

impl Cube {
    fn new(n: usize, initial: &State) -> Trial<Self> {
        // Freudenthal triangulation on 1x1x2 or 2x2x4 brick grids.
        let nz = 2 * n;
        let id = |x, y, z| (z * (n + 1) + y) * (n + 1) + x;
        let mut nodes = Vec::new();
        for z in 0..=nz {
            for y in 0..=n {
                for x in 0..=n {
                    nodes.push([
                        x as f64 / n as f64,
                        y as f64 / n as f64,
                        z as f64 / nz as f64,
                    ]);
                }
            }
        }
        let mut cells = Vec::new();
        let mut exposed = BTreeMap::new();
        for z in 0..nz {
            for y in 0..n {
                for x in 0..n {
                    for axes in [
                        [0, 1, 2],
                        [0, 2, 1],
                        [1, 0, 2],
                        [1, 2, 0],
                        [2, 0, 1],
                        [2, 1, 0],
                    ] {
                        let mut corner = [x, y, z];
                        let mut ids = [id(x, y, z); 4];
                        for (i, axis) in axes.into_iter().enumerate() {
                            corner[axis] += 1;
                            ids[i + 1] = id(corner[0], corner[1], corner[2]);
                        }
                        let point = |i: usize| Vector3::from(nodes[i]);
                        let e = [1, 2, 3].map(|k| point(ids[k]) - point(ids[0]));
                        if e[0].cross(&e[1]).dot(&e[2]) < 0.0 {
                            ids.swap(1, 2);
                        }
                        let (vol, b, _) =
                            elastic::element(&ids.map(|i| nodes[i]), &elastic::Tensor::identity())?;
                        for face in soil_fem_3d::elements::FACES {
                            let f = face.map(|i| ids[i]);
                            let mut key = f;
                            key.sort_unstable();
                            if exposed.remove(&key).is_none() {
                                let area = (point(f[1]) - point(f[0]))
                                    .cross(&(point(f[2]) - point(f[0])))
                                    / 2.;
                                exposed.insert(key, (f, area));
                            }
                        }
                        cells.push(Cell {
                            ids,
                            b,
                            vol,
                            state: initial.clone(),
                        });
                    }
                }
            }
        }
        Ok(Self {
            nodes,
            cells,
            faces: exposed.into_values().collect(),
        })
    }

    fn volume_mean(&self) -> [f64; 4] {
        let mut total = [0.; 4];
        let mut volume = 0.;
        for c in &self.cells {
            let (p, q) = invariants(&c.state.effective_stress_compression_positive_kpa);
            let values = [
                p,
                q,
                1. + c.state.void_ratio,
                c.state.preconsolidation_pressure_kpa,
            ];
            for i in 0..4 {
                total[i] += c.vol * values[i];
            }
            volume += c.vol;
        }
        total.map(|x| x / volume)
    }

    /// Reconstruct from nodal increments, integrate each cell's own parent, and
    /// test all free nodal equations against the prescribed pressure traction.
    fn step(
        &mut self,
        method: &MaterialStep,
        strain: [f64; 6],
        cell_pressure: f64,
        platen: bool,
    ) -> Trial<(f64, f64, f64, usize)> {
        let mut residual = vec![Vector3::zeros(); self.nodes.len()];
        let mut external = residual.clone();
        let mut p_range = [f64::INFINITY, f64::NEG_INFINITY];
        let mut substeps = 0;
        for c in &mut self.cells {
            let du = SVector::<f64, 12>::from_fn(|k, _| {
                -self.nodes[c.ids[k / 3]][k % 3] * strain[k % 3]
            });
            let de = -(c.b * du);
            let update = method.integrate(&c.state, de.into())?;
            substeps += update.substeps;
            c.state = update.state;
            let stress = SVector::<f64, 6>::from(c.state.effective_stress_compression_positive_kpa);
            let local = -c.vol * c.b.transpose() * stress;
            for k in 0..12 {
                residual[c.ids[k / 3]][k % 3] += local[k];
            }
            let p = invariants(&c.state.effective_stress_compression_positive_kpa).0;
            p_range = [p_range[0].min(p), p_range[1].max(p)];
        }
        for &(face, area) in &self.faces {
            if platen && area.z != 0.0 {
                continue;
            }
            for i in face {
                external[i] -= area * (cell_pressure / 3.);
            }
        }
        let mut free2 = 0.;
        let mut norm2 = 0.;
        let mut axial_reaction = 0.;
        for (i, xyz) in self.nodes.iter().enumerate() {
            for a in 0..3 {
                let r = residual[i][a] - external[i][a];
                if xyz[a] != 0. && !(platen && a == 2 && xyz[2] == 1.) {
                    free2 += r * r;
                }
                if platen && a == 2 && xyz[2] == 1. {
                    axial_reaction -= r;
                }
            }
            norm2 += external[i].norm_squared();
        }
        Ok((
            free2.sqrt() / norm2.sqrt(),
            p_range[1] - p_range[0],
            axial_reaction,
            substeps,
        ))
    }
}

fn specimen(
    ocr_target: f64,
    n: usize,
    h: f64,
    cfg: &Value,
    constant_g: bool,
    observe_strain: bool,
    comparison: Option<ElasticPath>,
) -> Trial<Value> {
    let start = std::time::Instant::now();
    let get = |key: &str| cfg[key].as_f64().unwrap();
    let p0 = get("p0_kpa");
    // Match Table 1's G at the common normally consolidated starting state.
    let v0 = get("v_lambda") - get("lambda") * (p0 / get("p1_kpa")).ln();
    let bulk0 = v0 * p0 / get("kappa");
    let shear0 = 250.0 * get("p1_kpa");
    let model = Parameters {
        critical_state_m: get("M"),
        lambda: get("lambda"),
        kappa: get("kappa"),
        poisson: (3. * bulk0 - 2. * shear0) / (2. * (3. * bulk0 + shear0)),
        reference_pressure_kpa: get("p1_kpa"),
        specific_volume_ncl_at_reference: get("v_lambda"),
        integration: Integration::BackwardEuler,
    };
    let initial_pc = if constant_g { ocr_target } else { p0 };
    let shear = constant_g.then_some(shear0);
    let method = MaterialStep {
        model: &model,
        h,
        shear,
        path: comparison.unwrap_or(ElasticPath::Default),
    };
    let (mut control_steps, mut cell_steps) = (0, 0);
    let mut state = State {
        effective_stress_compression_positive_kpa: [p0, p0, p0, 0., 0., 0.],
        preconsolidation_pressure_kpa: initial_pc,
        void_ratio: model.reference_void_ratio(p0, initial_pc)?,
        plastic_strain_compression_positive: [0.; 6],
    };
    let mut mesh = Cube::new(n, &state)?;
    let mut history = Vec::new();
    let mut axial_strain = 0.;
    let mut lateral_strain = 0.;
    let mut max_free: f64 = 0.;
    let mut max_spread: f64 = 0.;
    let mut max_side_error: f64 = 0.;
    let mut max_pc_error: f64 = 0.;
    let mut max_v_error: f64 = 0.;
    let mut max_q: f64 = 0.;
    let axial_end = get("terminal_axial_strain");
    for (stage, count) in [("consolidate", 128), ("unload", 128), ("shear", 512)] {
        let platen = stage == "shear";
        if constant_g && !platen {
            continue;
        }
        for i in 1..=count {
            let f = f64::from(i) / f64::from(count);
            let pressure = match stage {
                "consolidate" => p0 * (ocr_target / p0).powf(f),
                "unload" => ocr_target * (p0 / ocr_target).powf(f),
                _ => p0,
            };
            let dz = axial_end / f64::from(count);
            let strain = |lateral| {
                [
                    lateral,
                    lateral,
                    if platen { dz } else { lateral },
                    0.,
                    0.,
                    0.,
                ]
            };
            // The root is the physical strain giving the prescribed cell pressure.
            // Every evaluation starts from the same committed material state.
            let mut evaluate = |x| -> Trial<(f64, State)> {
                let update = method.integrate(&state, strain(x))?;
                control_steps += update.substeps;
                let s = update.state;
                Ok((s.effective_stress_compression_positive_kpa[0] - pressure, s))
            };
            let (mut lo, mut hi) = (-0.002, 0.002);
            if evaluate(lo)?.0 > 0. || evaluate(hi)?.0 < 0. {
                return Err("unbracketed cell pressure".into());
            }
            let mut success = None;
            for _ in 0..60 {
                let x = (lo + hi) / 2.;
                let (err, s) = evaluate(x)?;
                if err.abs() <= 1e-12 * pressure {
                    success = Some((x, s));
                    break;
                }
                if err < 0. {
                    lo = x;
                } else {
                    hi = x;
                }
            }
            let (lateral, next) = success.ok_or("pressure control did not converge")?;
            state = next;
            let (free, spread, reaction, substeps) =
                mesh.step(&method, strain(lateral), pressure, platen)?;
            cell_steps += substeps;
            max_free = max_free.max(free);
            max_spread = max_spread.max(spread);
            max_side_error = max_side_error
                .max((state.effective_stress_compression_positive_kpa[0] - pressure).abs());
            let [p, q, v, pc_history] = mesh.volume_mean();
            if platen {
                max_q = max_q.max(q);
            } else {
                let exact_pc = if stage == "consolidate" {
                    pressure
                } else {
                    ocr_target
                };
                max_pc_error = max_pc_error.max((pc_history - exact_pc).abs());
                max_v_error = max_v_error
                    .max((v - 1. - model.reference_void_ratio(pressure, exact_pc)?).abs());
            }
            if constant_g || i % 16 == 0 || i == 1 {
                let mut row = json!({"stage":stage,"fraction":f,"p_kpa":p,"q_kpa":q,
                    "pc_kpa":pc_history,"v":v,
                    "free_relative":free,"top_compression_reaction_kn":reaction});
                if observe_strain {
                    let increment = strain(lateral);
                    axial_strain += increment[2];
                    lateral_strain += increment[0];
                    row["axial_strain"] = json!(axial_strain);
                    row["lateral_strain"] = json!(lateral_strain);
                    row["volumetric_strain"] = json!(axial_strain + 2. * lateral_strain);
                }
                if comparison.is_some() {
                    row["plastic_volumetric_strain"] = json!(
                        state.plastic_strain_compression_positive[..3]
                            .iter()
                            .sum::<f64>()
                    );
                }
                history.push(row);
            }
        }
    }
    let [p, q, v, _] = mesh.volume_mean();
    let pc = p0 / (1. - model.critical_state_m / 3.);
    let predicted = [
        pc,
        model.critical_state_m * pc,
        1. + model.reference_void_ratio(pc, 2. * pc)?,
    ];
    let values = [p, q, v];
    let errors: [f64; 3] = std::array::from_fn(|i| (values[i] / predicted[i] - 1.).abs());
    let mut output = json!({"elements":mesh.cells.len(),"poisson":model.poisson,"material_substep":h,"ocr":ocr_target/p0,"history":history,
        "endpoint_p_q_v":values,"critical_p_q_v":predicted,"relative_errors":errors,
        "endpoint_passed":errors.iter().all(|e| *e < get("endpoint_relative_max_exclusive")),
        "max_free_relative":max_free,"max_cell_p_spread_kpa":max_spread,"max_lateral_error_kpa":max_side_error,
        "history_max_pc_error_kpa":max_pc_error,"history_max_v_error":max_v_error,"peak_q_kpa":max_q});
    if constant_g {
        let fields = output.as_object_mut().unwrap();
        fields.remove("poisson");
        fields.remove("history_max_pc_error_kpa");
        fields.remove("history_max_v_error");
        fields.insert("constant_shear_kpa".into(), json!(shear0));
        fields.insert("initial_pc_kpa".into(), json!(initial_pc));
    }
    if let Some(path) = comparison {
        output["analytic_elastic"] = json!(path == ElasticPath::Analytic);
        output["control_substeps"] = json!(control_steps);
        output["cell_substeps"] = json!(cell_steps);
        output["seconds"] = json!(start.elapsed().as_secs_f64());
    }
    Ok(output)
}

fn main() -> Trial<()> {
    let cfg: Value = serde_json::from_str(SPEC)?;
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args == ["--compare-ocr8"] {
        let mut records = Vec::new();
        for n in [1, 2] {
            for h in [1e-4, 5e-5] {
                for path in [ElasticPath::Default, ElasticPath::Analytic] {
                    records.push(specimen(40., n, h, &cfg["mcc"], true, true, Some(path))?);
                }
            }
        }
        println!("{}", serde_json::to_string_pretty(&records)?);
        return Ok(());
    }
    let (constant_g, observe_strain) = match args.as_slice() {
        [] => (false, false),
        [arg] if arg == "--constant-g" => (true, false),
        [arg, observation] if arg == "--constant-g" && observation == "--strain-history" => {
            (true, true)
        }
        _ => return Err("usage: t9d_mcc [--constant-g]".into()),
    };
    let mut records = Vec::new();
    for pressure in cfg["mcc"]["pc_targets_kpa"].as_array().unwrap() {
        for n in [1, 2] {
            for h in [1e-4, 5e-5] {
                records.push(specimen(
                    pressure.as_f64().unwrap(),
                    n,
                    h,
                    &cfg["mcc"],
                    constant_g,
                    observe_strain,
                    None,
                )?);
            }
        }
    }
    println!("{}", serde_json::to_string_pretty(&records)?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drained_history_reaches_the_published_critical_state() -> Trial<()> {
        let spec: Value = serde_json::from_str(SPEC)?;
        for pc in [8., 40.] {
            let result = specimen(pc, 1, 1e-4, &spec["mcc"], false, false, None)?;
            assert_eq!(result["endpoint_passed"], true);
            assert!(
                result["max_free_relative"].as_f64().unwrap()
                    <= spec["linear_relative_max"].as_f64().unwrap()
            );
            let history = result["history"].as_array().unwrap();
            let handoff = history
                .iter()
                .rev()
                .find(|r| r["stage"] == "unload")
                .unwrap();
            assert!(
                (handoff["pc_kpa"].as_f64().unwrap() - pc).abs()
                    < (handoff["pc_kpa"].as_f64().unwrap() - 5.).abs()
            );
        }
        Ok(())
    }
}
