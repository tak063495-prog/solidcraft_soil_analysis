//! Integration-only probe for independent affine cut-volume verification.
use serde::Deserialize;
use soil_fem_3d::{Error, Result, hydrostatic_volume};
use std::{fs, path::Path};

#[derive(Deserialize)]
struct Request {
    points: Vec<[f64; 3]>,
    level: f64,
    gamma: f64,
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 || Path::new(&args[2]).exists() {
        return Err(Error::Input(
            "hydrostatic_volume_probe INPUT NEW_OUTPUT".into(),
        ));
    }
    let requests: Vec<Request> = serde_json::from_slice(&fs::read(&args[1])?)?;
    let values: Vec<_> = requests
        .iter()
        .map(|r| hydrostatic_volume::affine_load(&r.points, r.level, r.gamma))
        .collect::<Result<_>>()?;
    fs::write(&args[2], serde_json::to_vec(&values)?)?;
    Ok(())
}
