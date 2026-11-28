//! JSON-lines evaluation of a trial against a caller-owned accepted checkpoint.
use soil_fem_3d::{Error, Result, stages};
use std::{
    env, fs,
    io::{self, BufRead, Write},
};

fn main() -> Result<()> {
    let args: Vec<_> = env::args().skip(1).collect();
    let [mesh, ground, case] = args.as_slice() else {
        return Err(Error::Input("t7_boundary MESH GROUND CASE".into()));
    };
    let mesh = serde_json::from_slice(&fs::read(mesh)?)?;
    let ground = serde_json::from_slice(&fs::read(ground)?)?;
    let case = serde_json::from_slice(&fs::read(case)?)?;
    let mut output = io::stdout().lock();
    for line in io::stdin().lock().lines() {
        let request = serde_json::from_str(&line?)?;
        let result = stages::t7::evaluate(&mesh, &ground, &case, request)
            .unwrap_or_else(|e| serde_json::json!({"error":e.to_string()}));
        serde_json::to_writer(&mut output, &result)?;
        writeln!(output)?;
        output.flush()?;
    }
    Ok(())
}
