//! MCC constitutive response for displacement/mean-pressure/volume formulations.
//!
//! The independent scalar is compression-positive volumetric strain, not a
//! replacement pressure added after a conventional return. The full MCC return
//! retains its own admissible mean stress and all history. A mixed formulation
//! must enforce pressure consistency and weak volume compatibility separately.
use crate::{Error, Result, elastic::Tensor, mcc};
use nalgebra::{SMatrix, SVector};

/// Derivatives of the complete fixed-partition increment map. Engineering
/// shear components use the same convention as `mcc`; no symmetrization.
pub struct Tangent {
    pub deviatoric_by_deviatoric: Tensor,
    pub deviatoric_by_volume: [f64; 6],
    pub mean_by_deviatoric: [f64; 6],
    pub mean_by_volume: f64,
}

pub struct Response {
    pub update: mcc::Update,
    pub mean_pressure_compression_positive_kpa: f64,
    pub deviatoric_stress_compression_positive_kpa: [f64; 6],
    pub tangent: Tangent,
}

pub struct Values {
    pub update: mcc::Update,
    pub mean_pressure_compression_positive_kpa: f64,
    pub deviatoric_stress_compression_positive_kpa: [f64; 6],
}

fn strain(deviatoric_strain: [f64; 6], volumetric_strain: f64) -> Result<[f64; 6]> {
    let magnitude = deviatoric_strain.iter().map(|x| x.abs()).sum::<f64>();
    let trace = deviatoric_strain[..3].iter().sum::<f64>();
    if !volumetric_strain.is_finite()
        || !magnitude.is_finite()
        || !deviatoric_strain.iter().all(|x| x.is_finite())
        || trace.abs() > 64. * f64::EPSILON * magnitude.max(1e-12)
    {
        return Err(Error::Input(
            "MCC split response requires finite volume and trace-free deviatoric strain".into(),
        ));
    }
    let mut out = deviatoric_strain;
    for value in &mut out[..3] {
        *value += volumetric_strain / 3.;
    }
    Ok(out)
}

pub fn integrate(
    parameters: &mcc::Parameters,
    committed: &mcc::State,
    deviatoric_strain: [f64; 6],
    volumetric_strain: f64,
    partition: usize,
) -> Result<Values> {
    let update = parameters.integrate_partition(
        committed,
        strain(deviatoric_strain, volumetric_strain)?,
        partition,
    )?;
    let stress = update.state.effective_stress_compression_positive_kpa;
    let mean = stress[..3].iter().sum::<f64>() / 3.;
    Ok(Values {
        update,
        mean_pressure_compression_positive_kpa: mean,
        deviatoric_stress_compression_positive_kpa: std::array::from_fn(|i| {
            stress[i] - if i < 3 { mean } else { 0. }
        }),
    })
}

/// Evaluate an immutable trial from committed history with a fixed material
/// partition. Failure must be handled by the caller's full state rollback.
/// A nonpositive volumetric tangent is retained: it is not clipped or inverted.
pub fn evaluate(
    parameters: &mcc::Parameters,
    committed: &mcc::State,
    deviatoric_strain: [f64; 6],
    volumetric_strain: f64,
    partition: usize,
    perturbation: f64,
) -> Result<Response> {
    let strain = strain(deviatoric_strain, volumetric_strain)?;
    let (update, full) =
        parameters.stress_tangent_partition(committed, strain, partition, perturbation)?;
    Ok(response(update, full))
}

/// Branch-consistent implicit derivative of the unchanged primal return.
pub fn evaluate_implicit(
    parameters: &mcc::Parameters,
    committed: &mcc::State,
    deviatoric_strain: [f64; 6],
    volumetric_strain: f64,
    partition: usize,
) -> Result<Response> {
    let (update, full) = parameters.implicit_stress_tangent_partition(
        committed,
        strain(deviatoric_strain, volumetric_strain)?,
        partition,
    )?;
    Ok(response(update, full))
}

fn response(update: mcc::Update, full: Tensor) -> Response {
    let normal = SVector::<f64, 6>::from_row_slice(&[1., 1., 1., 0., 0., 0.]);
    let projection = SMatrix::<f64, 6, 6>::identity() - normal * normal.transpose() / 3.;
    let stress = SVector::from_row_slice(&update.state.effective_stress_compression_positive_kpa);
    let mean = normal.dot(&stress) / 3.;
    let dev = projection * stress;
    let dev_volume = projection * full * normal / 3.;
    let mean_dev = normal.transpose() * full * projection / 3.;
    let mean_volume = normal.dot(&(full * normal)) / 9.;
    Response {
        update,
        mean_pressure_compression_positive_kpa: mean,
        deviatoric_stress_compression_positive_kpa: std::array::from_fn(|i| dev[i]),
        tangent: Tangent {
            deviatoric_by_deviatoric: projection * full * projection,
            deviatoric_by_volume: std::array::from_fn(|i| dev_volume[i]),
            mean_by_deviatoric: std::array::from_fn(|i| mean_dev[i]),
            mean_by_volume: mean_volume,
        },
    }
}
