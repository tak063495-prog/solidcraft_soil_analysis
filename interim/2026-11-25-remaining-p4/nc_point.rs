// Replay the rejected material call at two uniform stress units.
use serde_json::{Value, json};
use soil_fem_3d::mcc::{Parameters, State};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let row: Value = serde_json::from_str(&line?)?;
        for exponent in [448, 480] {
            let factor = 2.0_f64.powi(exponent);
            let mut state: State = serde_json::from_value(row["state"].clone())?;
            let mut law: Parameters = serde_json::from_value(row["law"].clone())?;
            let strain = serde_json::from_value(row["strain"].clone())?;
            let count = serde_json::from_value(row["count"].clone())?;
            law.reference_pressure_kpa *= factor;
            state.preconsolidation_pressure_kpa *= factor;
            for x in &mut state.effective_stress_compression_positive_kpa {
                *x *= factor;
            }
            let square = state.preconsolidation_pressure_kpa.powi(2);
            let result = law.implicit_stress_tangent_partition(&state, strain, count);
            println!(
                "{}",
                json!({"point": row["point"], "exponent": exponent,
                "scaled_pc_squared": square, "subnormal_square": square.is_subnormal(),
                "reason": result.err().map(|e|e.to_string())})
            );
        }
    }
    Ok(())
}
