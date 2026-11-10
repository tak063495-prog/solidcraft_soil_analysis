//! Small constitutive evidence driver; independent integration lives in Python.
use serde::{Deserialize, Serialize};
use soil_fem_3d::{
    Result,
    joint_law::{Material, PointState},
};
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    material: Material,
    normal: [f64; 3],
    reference: [f64; 3],
    traction: [f64; 3],
    increments: Vec<[f64; 3]>,
    maximum_increment_m: f64,
}
#[derive(Serialize)]
struct Output {
    states: Vec<PointState>,
}
fn main() -> Result<()> {
    let input: Input = serde_json::from_reader(std::io::stdin())?;
    let mut s = input
        .material
        .birth(input.normal, input.reference, input.traction)?;
    let mut states = vec![s.clone()];
    for increment in input.increments {
        s = input.material.integrate(
            &s,
            std::array::from_fn(|a| input.reference[a] + increment[a]),
            input.maximum_increment_m,
        )?;
        states.push(s.clone());
    }
    serde_json::to_writer_pretty(std::io::stdout(), &Output { states })?;
    Ok(())
}
