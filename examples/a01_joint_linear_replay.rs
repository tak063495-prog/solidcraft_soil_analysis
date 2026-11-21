//! Replay captured A,b only; no nonlinear continuation or acceptance fallback.
use serde_json::{Value, json};
use soil_fem_3d::{
    Result,
    linalg::{self, CsrMatrix, CsrPattern, GmresOptions, IlutOptions},
};
use std::{fs, io::Write};
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "SYSTEMS NEW_OUTPUT");
    let systems: Vec<Value> = serde_json::from_slice(&fs::read(&args[0])?)?;
    let mut output = Vec::new();
    for system in systems {
        let rhs: Vec<f64> = serde_json::from_value(system["rhs"].clone())?;
        let n = rhs.len();
        assert!(n > 0 && n <= 256);
        let row: Vec<usize> = serde_json::from_value(system["storage"]["row_ptr"].clone())?;
        let cols: Vec<usize> = serde_json::from_value(system["storage"]["columns"].clone())?;
        let values: Vec<f64> = serde_json::from_value(system["storage"]["values"].clone())?;
        assert_eq!(row.len(), n + 1);
        assert_eq!(cols.len(), values.len());
        let connections =
            (0..n).flat_map(|i| cols[row[i]..row[i + 1]].iter().map(move |&j| vec![i, j]));
        let mut a = CsrMatrix::zeros(CsrPattern::from_connectivity(n, connections)?);
        assert_eq!(a.row_ptr(), row);
        assert_eq!(a.columns(), cols);
        for i in 0..n {
            for ij in row[i]..row[i + 1] {
                a.add(i, cols[ij], values[ij])?;
            }
        }
        let original: GmresOptions =
            serde_json::from_value(system["settings"]["linear_solver"]["options"].clone())?;
        let tol = system["settings"]["relative_linear_tolerance"]
            .as_f64()
            .unwrap();
        let mut longer = original.clone();
        longer.maximum_iterations = (original.maximum_iterations * 4).min(1_000_000);
        let mut full = longer.clone();
        full.restart = n;
        let mut ilut = full.clone();
        ilut.ilut = Some(IlutOptions {
            drop_tolerance: 0.,
            fill_limit: n,
        });
        let mut variants = Vec::new();
        for (name, options) in [
            ("original", original),
            ("iteration_cap_x4", longer),
            ("full_restart", full),
            ("existing_ilut_full_fill", ilut),
        ] {
            let mut samples = Vec::new();
            let result = linalg::gmres_ilu0_a01_residual_observed(
                &a,
                &rhs,
                tol,
                &options,
                || false,
                &mut |iteration, residual| {
                    samples.push(json!({"iteration":iteration,"true_relative_residual":residual}))
                },
            );
            let outcome = match result {
                Ok((x, report)) => json!({"solution":x,"report":report}),
                Err(e) => json!({"error":e.to_string(),"structured_error":e}),
            };
            variants
                .push(json!({"name":name,"options":options,"samples":samples,"outcome":outcome}));
        }
        output.push(json!({"id":system["id"],"variants":variants}));
    }
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[1])?;
    serde_json::to_writer(&mut file, &output)?;
    writeln!(file)?;
    Ok(())
}
