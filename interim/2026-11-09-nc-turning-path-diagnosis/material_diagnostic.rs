//! Read-only fixed-history MCC map and derivative measurements. No equilibrium solve.
use nalgebra::SMatrix;
use serde::Deserialize;
use soil_fem_3d::{Error, Result, mcc};
use std::{fs, path::Path};

#[derive(Deserialize)]
struct Request {
    label: String,
    point: usize,
    parameters: mcc::Parameters,
    state: mcc::State,
    strain: [f64; 6],
    partitions: usize,
    widths: Vec<f64>,
}
fn rows(t: &SMatrix<f64, 6, 6>) -> Vec<Vec<f64>> {
    (0..6)
        .map(|i| (0..6).map(|j| t[(i, j)]).collect())
        .collect()
}
fn measure(r: &Request) -> Result<serde_json::Value> {
    if r.parameters.integration != mcc::Integration::Midpoint {
        return Err(Error::Input(
            "this diagnostic requires midpoint integration".into(),
        ));
    }
    let (ad, k) =
        r.parameters
            .implicit_stress_tangent_partition(&r.state, r.strain, r.partitions)?;
    let primal = r
        .parameters
        .integrate_partition(&r.state, r.strain, r.partitions)?;
    let mut finite = Vec::new();
    for &h in &r.widths {
        if !(1e-14..=1e-5).contains(&h) {
            return Err(Error::Input("diagnostic width out of bounds".into()));
        }
        let mut central = SMatrix::<f64, 6, 6>::zeros();
        let mut forward = central;
        let mut backward = central;
        let mut branches = Vec::new();
        for j in 0..6 {
            let mut plus = r.strain;
            plus[j] += h;
            let mut minus = r.strain;
            minus[j] -= h;
            let a = r
                .parameters
                .integrate_partition(&r.state, plus, r.partitions)?;
            let b = r
                .parameters
                .integrate_partition(&r.state, minus, r.partitions)?;
            branches.push([a.plastic_substeps, b.plastic_substeps]);
            for i in 0..6 {
                let base = primal.state.effective_stress_compression_positive_kpa[i];
                central[(i, j)] = (a.state.effective_stress_compression_positive_kpa[i]
                    - b.state.effective_stress_compression_positive_kpa[i])
                    / (2. * h);
                forward[(i, j)] = (a.state.effective_stress_compression_positive_kpa[i] - base) / h;
                backward[(i, j)] =
                    (base - b.state.effective_stress_compression_positive_kpa[i]) / h;
            }
        }
        finite.push(serde_json::json!({"width":h,"central":rows(&central),
            "forward":rows(&forward),"backward":rows(&backward),
            "central_relative_to_ad":(central-k).norm()/k.norm(),
            "forward_relative_to_ad":(forward-k).norm()/k.norm(),
            "backward_relative_to_ad":(backward-k).norm()/k.norm(),
            "plastic_substeps_plus_minus":branches}));
    }
    let path = mcc::diagnostic_midpoint_path(&r.parameters, &r.state, r.strain, r.partitions)?;
    Ok(
        serde_json::json!({"label":r.label,"point":r.point,"implicit_tangent":rows(&k),
        "updated_state":primal.state,"plastic_substeps":primal.plastic_substeps,
        "ad_primal_exact":serde_json::to_value(&ad.state)?==serde_json::to_value(&primal.state)?,
        "substep_replay_exact":path["state"]==serde_json::to_value(&primal.state)?,
        "substeps":path["substeps"],"finite_differences":finite}),
    )
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    if a.len() != 3 || Path::new(&a[2]).exists() {
        return Err(Error::Input(
            "material_diagnostic REQUESTS NEW_OUTPUT".into(),
        ));
    }
    let requests: Vec<Request> = serde_json::from_slice(&fs::read(&a[1])?)?;
    if requests.len() > 6000 {
        return Err(Error::Input("small-case diagnostic limit".into()));
    }
    let output: Vec<_> = requests
        .iter()
        .map(|r| match measure(r) {
            Ok(v) => v,
            Err(e) => serde_json::json!({"label":r.label,"point":r.point,"error":e.to_string()}),
        })
        .collect();
    fs::write(
        &a[2],
        serde_json::to_vec(
            &serde_json::json!({"passed":false,"independent_reference":false,"points":output}),
        )?,
    )?;
    Ok(())
}
