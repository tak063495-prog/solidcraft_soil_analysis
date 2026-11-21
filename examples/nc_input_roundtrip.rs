//! Print the native schema representation of diagnostic inputs, without solving.
use soil_fem_3d::{ground::Ground, mesh::Mesh, stages::StageCase};
use std::io;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inputs: (Mesh, Ground, StageCase) = serde_json::from_reader(io::stdin().lock())?;
    serde_json::to_writer(io::stdout().lock(), &inputs)?;
    Ok(())
}
