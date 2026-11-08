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
// Reconstruct the existing backward-Euler local equations from each exact
// substep return. This is a residual audit, not a second constitutive model.
fn local_residuals(r: &Request) -> Result<serde_json::Value> {
    if r.parameters.integration != mcc::Integration::BackwardEuler {
        return Err(soil_fem_3d::Error::Input(
            "local residual audit requires backward Euler".into(),
        ));
    }
    let mut old = r.state.clone();
    let de = r.strain.map(|v| v / r.partitions as f64);
    let mut records = Vec::new();
    for _ in 0..r.partitions {
        let update = r.parameters.integrate_partition(&old, de, 1)?;
        let (k, g) = r.parameters.moduli(&old)?;
        let ev = de[..3].iter().sum::<f64>();
        let mut trial = old.effective_stress_compression_positive_kpa;
        for i in 0..3 {
            trial[i] += k * ev + 2. * g * (de[i] - ev / 3.);
        }
        for i in 3..6 {
            trial[i] += g * de[i];
        }
        let (pt, qt) = mcc::invariants(&trial);
        let (p, q) = mcc::invariants(&update.state.effective_stress_compression_positive_kpa);
        let c = update.state.preconsolidation_pressure_kpa;
        let scale = old.preconsolidation_pressure_kpa.max(pt).max(qt);
        let m2 = r.parameters.critical_state_m.powi(2);
        let pn = p / scale;
        let qn = q / scale;
        let cn = c / scale;
        let h = m2 * (2. * pn - cn);
        let kn = k / scale;
        let gn = g / scale;
        // Least-squares reconstruction avoids division by h at critical state.
        let a = kn * h;
        let b = 6. * gn * qn;
        let l = if update.plastic_substeps > 0 {
            (a * (pt - p) / scale + b * (qt - q) / scale) / (a * a + b * b)
        } else {
            0.
        };
        let hardening = (1. + old.void_ratio) / (r.parameters.lambda - r.parameters.kappa);
        let residual = if update.plastic_substeps > 0 {
            vec![
                pn - pt / scale + kn * l * h,
                qn - qt / scale + 6. * gn * l * qn,
                (cn / (old.preconsolidation_pressure_kpa / scale)).ln() - hardening * l * h,
                qn * qn + m2 * pn * (pn - cn),
            ]
        } else {
            vec![0.; 4]
        };
        records.push(serde_json::json!({"plastic": update.plastic_substeps > 0,
            "normalized_yield": r.parameters.normalized_yield(&update.state),
            "reconstructed_local_residual_norm": residual.iter().map(|x| x*x).sum::<f64>().sqrt(),
            "local_equations": residual,
            "specific_volume_error": (1.+update.state.void_ratio)-(1.+old.void_ratio)*(-ev).exp(),
        }));
        old = update.state;
    }
    Ok(
        serde_json::json!({"substeps": records, "replayed_state": old,
        "residual_kind": "four existing BE equations; multiplier reconstructed by least squares"}),
    )
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
        "local_residuals": local_residuals(r)?,
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
        assert_eq!(
            result["local_residuals"]["replayed_state"],
            result["updated_state"]
        );
        for step in result["local_residuals"]["substeps"].as_array().unwrap() {
            assert!(step["reconstructed_local_residual_norm"].as_f64().unwrap() < 1e-11);
            assert_eq!(step["specific_volume_error"], 0.);
        }
        Ok(())
    }
}
