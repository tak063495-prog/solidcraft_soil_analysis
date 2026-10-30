//! Prototype only: independent deviatoric/volume inputs to the unchanged MC map.
//! A mixed FE formulation must separately enforce volume compatibility and
//! pressure consistency. This constitutive split proves no mixed equilibrium.
use soil_fem_3d::{Error, Result, mc};

pub struct Tangent {
    pub deviatoric_by_deviatoric: [[f64; 6]; 6],
    pub deviatoric_by_volume: [f64; 6],
    pub mean_by_deviatoric: [f64; 6],
    pub mean_by_volume: f64,
}
pub struct Values {
    pub update: mc::Update,
    pub mean_pressure_compression_positive_kpa: f64,
    pub deviatoric_stress_compression_positive_kpa: [f64; 6],
}
pub struct Response {
    pub values: Values,
    pub tangent: Tangent,
}
fn strain(deviatoric: [f64; 6], volume: f64) -> Result<[f64; 6]> {
    let magnitude = deviatoric.iter().map(|v| v.abs()).sum::<f64>();
    let trace = deviatoric[..3].iter().sum::<f64>();
    if !volume.is_finite()
        || !magnitude.is_finite()
        || !deviatoric.iter().all(|v| v.is_finite())
        || trace.abs() > 64. * f64::EPSILON * magnitude.max(1e-12)
    {
        return Err(Error::Input(
            "MC split requires finite volume and trace-free deviatoric engineering strain".into(),
        ));
    }
    let mut out = deviatoric;
    for v in &mut out[..3] {
        *v += volume / 3.;
    }
    Ok(out)
}
fn values(update: mc::Update) -> Values {
    let stress = update.state.effective_stress_compression_positive_kpa;
    let mean = stress[..3].iter().sum::<f64>() / 3.;
    Values {
        update,
        mean_pressure_compression_positive_kpa: mean,
        deviatoric_stress_compression_positive_kpa: std::array::from_fn(|i| {
            stress[i] - if i < 3 { mean } else { 0. }
        }),
    }
}
pub fn integrate(
    parameters: &mc::Parameters,
    committed: &mc::State,
    deviatoric: [f64; 6],
    volume: f64,
    partitions: usize,
) -> Result<Values> {
    Ok(values(parameters.integrate_partition(
        committed,
        strain(deviatoric, volume)?,
        partitions,
    )?))
}
pub fn evaluate(
    parameters: &mc::Parameters,
    committed: &mc::State,
    deviatoric: [f64; 6],
    volume: f64,
    partitions: usize,
    perturbation: f64,
) -> Result<Response> {
    let (update, full) = parameters.stress_tangent_partition(
        committed,
        strain(deviatoric, volume)?,
        partitions,
        perturbation,
    )?;
    Ok(response(update, full))
}
fn response(update: mc::Update, full: soil_fem_3d::elastic::Tensor) -> Response {
    let normal: [f64; 6] = [1., 1., 1., 0., 0., 0.];
    let projection: [[f64; 6]; 6] = std::array::from_fn(|i| {
        std::array::from_fn(|j| if i == j { 1. } else { 0. } - normal[i] * normal[j] / 3.)
    });
    let deviatoric_by_deviatoric = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let mut value = 0.;
            for a in 0..6 {
                for b in 0..6 {
                    value += projection[i][a] * full[(a, b)] * projection[b][j];
                }
            }
            value
        })
    });
    let deviatoric_by_volume = std::array::from_fn(|i| {
        let mut value = 0.;
        for a in 0..6 {
            for b in 0..3 {
                value += projection[i][a] * full[(a, b)] / 3.;
            }
        }
        value
    });
    let mean_by_deviatoric = std::array::from_fn(|j| {
        let mut value = 0.;
        for a in 0..3 {
            for b in 0..6 {
                value += full[(a, b)] * projection[b][j] / 3.;
            }
        }
        value
    });
    let mut mean_by_volume = 0.;
    for a in 0..3 {
        for b in 0..3 {
            mean_by_volume += full[(a, b)] / 9.;
        }
    }
    Response {
        values: values(update),
        tangent: Tangent {
            deviatoric_by_deviatoric,
            deviatoric_by_volume,
            mean_by_deviatoric,
            mean_by_volume,
        },
    }
}

/// Explicit old-parameter to new-parameter map. The original MC transition
/// validates the old state under `from` without clearing or relabelling it.
pub fn integrate_transition(
    parameters: &mc::Parameters,
    from: &mc::Parameters,
    committed: &mc::State,
    deviatoric: [f64; 6],
    volume: f64,
    partitions: usize,
) -> Result<Values> {
    Ok(values(parameters.integrate_transition_partition(
        from,
        committed,
        strain(deviatoric, volume)?,
        partitions,
    )?))
}
pub fn evaluate_transition(
    parameters: &mc::Parameters,
    from: &mc::Parameters,
    committed: &mc::State,
    deviatoric: [f64; 6],
    volume: f64,
    partitions: usize,
    perturbation: f64,
) -> Result<Response> {
    let (update, full) = parameters.stress_tangent_transition_partition(
        from,
        committed,
        strain(deviatoric, volume)?,
        partitions,
        perturbation,
    )?;
    Ok(response(update, full))
}
