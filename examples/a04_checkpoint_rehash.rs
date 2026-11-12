//! A-04 diagnostic only: rehash a typed construction checkpoint for tamper tests.
use soil_fem_3d::{Result, checkpoint::hash_json, one_way_construction::Checkpoint};
use std::{fs, fs::OpenOptions, io::Write};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "INPUT NEW_OUTPUT");
    let mut saved: Checkpoint = serde_json::from_slice(&fs::read(&args[0])?)?;
    saved.run_sha256 = hash_json(&saved.run)?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[1])?;
    file.write_all(&serde_json::to_vec(&saved)?)?;
    Ok(())
}
