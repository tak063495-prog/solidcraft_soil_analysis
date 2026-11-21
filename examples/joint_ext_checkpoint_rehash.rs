//! Test-only adversarial re-sealing; resume must still authenticate the path.
use soil_fem_3d::{
    Result,
    checkpoint::{Checkpoint, hash_json},
};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut saved: Checkpoint = serde_json::from_slice(&std::fs::read(&args[1])?)?;
    saved.run_sha256 = hash_json(&saved.run)?;
    std::fs::write(&args[2], serde_json::to_vec_pretty(&saved)?)?;
    Ok(())
}
