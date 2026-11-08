//! Read-only material-map measurements; not a reference solution or a solver.
use serde::{Deserialize, Serialize};
use soil_fem_3d::{Result, elastic::Tensor, mcc};
use std::{fs, path::Path};

#[derive(Deserialize)]
struct Request {
    parameters: mcc::Parameters,
    state: mcc::State,
    strain: [f64; 6],
    partitions: usize,
}
#[derive(Serialize)]
struct Difference {
    width: f64,
    tangent: Vec<Vec<f64>>,
    relative_to_implicit: f64,
    plus_minus_plastic_substeps: Vec<[usize; 2]>,
}
fn rows(t: &Tensor) -> Vec<Vec<f64>> {
    (0..6)
        .map(|i| (0..6).map(|j| t[(i, j)]).collect())
        .collect()
}
fn measure(r: &Request) -> Result<serde_json::Value> {
    let before = serde_json::to_value(&r.state)?;
    let (implicit, tangent) =
        r.parameters
            .implicit_stress_tangent_partition(&r.state, r.strain, r.partitions)?;
    let primal = r
        .parameters
        .integrate_partition(&r.state, r.strain, r.partitions)?;
    let same = serde_json::to_value(&implicit.state)? == serde_json::to_value(&primal.state)?;
    let mut differences = Vec::new();
    for width in [1e-7, 1e-8, 1e-9] {
        let (_, fd) =
            r.parameters
                .stress_tangent_partition(&r.state, r.strain, r.partitions, width)?;
        let mut branches = Vec::new();
        for j in 0..6 {
            let mut plus = r.strain;
            let mut minus = r.strain;
            plus[j] += width;
            minus[j] -= width;
            branches.push([
                r.parameters
                    .integrate_partition(&r.state, plus, r.partitions)?
                    .plastic_substeps,
                r.parameters
                    .integrate_partition(&r.state, minus, r.partitions)?
                    .plastic_substeps,
            ]);
        }
        differences.push(Difference {
            width,
            tangent: rows(&fd),
            relative_to_implicit: (fd - tangent).norm() / tangent.norm(),
            plus_minus_plastic_substeps: branches,
        });
    }
    Ok(serde_json::json!({
        "implicit_tangent": rows(&tangent), "finite_differences": differences,
        "primal_state_exact": same,
        "input_state_unchanged": before == serde_json::to_value(&r.state)?,
        "plastic_substeps": primal.plastic_substeps,
        "normalized_yield": r.parameters.normalized_yield(&primal.state),
        "updated_state": primal.state,
    }))
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 || Path::new(&args[2]).exists() {
        return Err(soil_fem_3d::Error::Input(
            "nc_material_probe INPUT NEW_OUTPUT".into(),
        ));
    }
    let requests: Vec<Request> = serde_json::from_slice(&fs::read(&args[1])?)?;
    let result: Vec<_> = requests
        .iter()
        .map(|r| match measure(r) {
            Ok(v) => v,
            Err(e) => serde_json::json!({"error": e.to_string()}),
        })
        .collect();
    fs::write(
        &args[2],
        serde_json::to_vec(&serde_json::json!({
            "passed": false, "independent_reference": false, "points": result,
        }))?,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagnostic_keeps_primal_return_and_committed_state() -> Result<()> {
        let parameters: mcc::Parameters = serde_json::from_value(
            serde_json::from_str::<serde_json::Value>(include_str!("ground-joint-mcc.json"))?
                ["materials"][0]["mcc"].clone(),
        )?;
        let state = mcc::State {
            effective_stress_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
            preconsolidation_pressure_kpa: 100.,
            void_ratio: parameters.reference_void_ratio(100., 100.)?,
            plastic_strain_compression_positive: [0.; 6],
        };
        let result = measure(&Request {
            parameters,
            state,
            strain: [1e-5, 1e-5, 1e-5, 0., 0., 0.],
            partitions: 2,
        })?;
        assert_eq!(result["primal_state_exact"], true);
        assert_eq!(result["input_state_unchanged"], true);
        assert_eq!(
            result["finite_differences"].as_array().map(Vec::len),
            Some(3)
        );
        assert!(result["plastic_substeps"].as_u64().is_some_and(|n| n > 0));
        Ok(())
    }
}
