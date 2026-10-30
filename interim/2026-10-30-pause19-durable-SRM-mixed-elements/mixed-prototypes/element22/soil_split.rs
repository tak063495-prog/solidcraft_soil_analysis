//! Prototype typed split dispatch. MC and MCC keep their own complete states.
use crate::mc_split;
use soil_fem_3d::{Error, Result, mcc_split, soil};
#[derive(Clone, Copy)]
pub enum Method {
    FiniteDifference(f64),
    ImplicitMcc,
}
pub struct Values {
    pub update: soil::Update,
    pub mean_pressure_compression_positive_kpa: f64,
    pub deviatoric_stress_compression_positive_kpa: [f64; 6],
}
pub struct Response {
    pub values: Values,
    pub tangent: mc_split::Tangent,
}
fn mc_values(v: mc_split::Values) -> Values {
    Values {
        update: soil::Update {
            state: soil::State::MohrCoulomb(v.update.state),
            substeps: v.update.substeps,
            plastic_substeps: v.update.plastic_substeps,
        },
        mean_pressure_compression_positive_kpa: v.mean_pressure_compression_positive_kpa,
        deviatoric_stress_compression_positive_kpa: v.deviatoric_stress_compression_positive_kpa,
    }
}
fn mcc_values(v: mcc_split::Values) -> Values {
    Values {
        update: soil::Update {
            state: soil::State::Mcc(v.update.state),
            substeps: v.update.substeps,
            plastic_substeps: v.update.plastic_substeps,
        },
        mean_pressure_compression_positive_kpa: v.mean_pressure_compression_positive_kpa,
        deviatoric_stress_compression_positive_kpa: v.deviatoric_stress_compression_positive_kpa,
    }
}
pub fn integrate(
    model: &soil::Model,
    committed: &soil::State,
    dev: [f64; 6],
    volume: f64,
    partitions: usize,
) -> Result<Values> {
    match (model, committed) {
        (soil::Model::Mcc(p), soil::State::Mcc(s)) => Ok(mcc_values(mcc_split::integrate(
            p, s, dev, volume, partitions,
        )?)),
        (soil::Model::MohrCoulomb(p), soil::State::MohrCoulomb(s)) => Ok(mc_values(
            mc_split::integrate(p, s, dev, volume, partitions)?,
        )),
        _ => Err(Error::Input(
            "split material/state type mismatch or unsupported split model".into(),
        )),
    }
}
pub fn evaluate(
    model: &soil::Model,
    committed: &soil::State,
    dev: [f64; 6],
    volume: f64,
    partitions: usize,
    method: Method,
) -> Result<Response> {
    match (model, committed, method) {
        (soil::Model::Mcc(p), soil::State::Mcc(s), method) => {
            let r = match method {
                Method::ImplicitMcc => mcc_split::evaluate_implicit(p, s, dev, volume, partitions)?,
                Method::FiniteDifference(h) => {
                    mcc_split::evaluate(p, s, dev, volume, partitions, h)?
                }
            };
            let tangent = mc_split::Tangent {
                deviatoric_by_deviatoric: std::array::from_fn(|i| {
                    std::array::from_fn(|j| r.tangent.deviatoric_by_deviatoric[(i, j)])
                }),
                deviatoric_by_volume: r.tangent.deviatoric_by_volume,
                mean_by_deviatoric: r.tangent.mean_by_deviatoric,
                mean_by_volume: r.tangent.mean_by_volume,
            };
            let values = mcc_values(mcc_split::Values {
                update: r.update,
                mean_pressure_compression_positive_kpa: r.mean_pressure_compression_positive_kpa,
                deviatoric_stress_compression_positive_kpa: r
                    .deviatoric_stress_compression_positive_kpa,
            });
            Ok(Response { values, tangent })
        }
        (soil::Model::MohrCoulomb(p), soil::State::MohrCoulomb(s), Method::FiniteDifference(h)) => {
            let r = mc_split::evaluate(p, s, dev, volume, partitions, h)?;
            Ok(Response {
                values: mc_values(r.values),
                tangent: r.tangent,
            })
        }
        _ => Err(Error::Input(
            "split material/state/tangent method mismatch or unsupported split model".into(),
        )),
    }
}
pub fn integrate_transition(
    to: &soil::Model,
    from: &soil::Model,
    committed: &soil::State,
    dev: [f64; 6],
    volume: f64,
    partitions: usize,
) -> Result<Values> {
    match (to, from, committed) {
        (
            soil::Model::MohrCoulomb(p),
            soil::Model::MohrCoulomb(old),
            soil::State::MohrCoulomb(s),
        ) => Ok(mc_values(mc_split::integrate_transition(
            p, old, s, dev, volume, partitions,
        )?)),
        _ => Err(Error::Input(
            "explicit split transition requires matching old/new MC models and an MC state".into(),
        )),
    }
}

/// An explicitly typed fixed-partition parameter jump. No implicit MC
/// approximation or cross-model history conversion is offered.
pub fn evaluate_transition(
    to: &soil::Model,
    from: &soil::Model,
    committed: &soil::State,
    dev: [f64; 6],
    volume: f64,
    partitions: usize,
    method: Method,
) -> Result<Response> {
    match (to,from,committed,method) {
        (soil::Model::MohrCoulomb(p),soil::Model::MohrCoulomb(old),soil::State::MohrCoulomb(s),Method::FiniteDifference(h)) => {
            let r=mc_split::evaluate_transition(p,old,s,dev,volume,partitions,h)?;
            Ok(Response {values:mc_values(r.values),tangent:r.tangent})
        },
        _ => Err(Error::Input("split transition tangent requires matching MC types and explicit fixed-partition finite difference".into()))
    }
}
