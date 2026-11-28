// Input is an immutable, captured local update. No return value feeds a global solve.
use serde_json::{json, Value};
use soil_fem_3d::mcc::{Parameters, State};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let row: Value = serde_json::from_str(&line?)?;
        let mut record = row["state"].clone();
        record.as_object_mut().unwrap().remove("stress_exponent");
        let mut state: State = serde_json::from_value(record)?;
        let mut law: Parameters = serde_json::from_value(row["law"].clone())?;
        let shift = match row["mode"].as_str().unwrap() {
            "homogeneous" => -(state.preconsolidation_pressure_kpa.log2().floor() as i32),
            "units448" => 448,
            "comparison480" => 480,
            _ => panic!("unknown diagnostic"),
        };
        let factor = 2_f64.powi(shift);
        if row["mode"] != "homogeneous" {
            law.reference_pressure_kpa *= factor;
        }
        state.preconsolidation_pressure_kpa *= factor;
        for s in &mut state.effective_stress_compression_positive_kpa {
            *s *= factor;
        }
        let pc = state.preconsolidation_pressure_kpa;
        let strain = serde_json::from_value(row["strain"].clone())?;
        let count = serde_json::from_value(row["count"].clone())?;
        let reason = law.implicit_stress_tangent_partition(&state, strain, count)
            .err().map(|e| e.to_string());
        println!("{}",json!({"mode":row["mode"],"point":row["point"],"shift":shift,
            "pc":pc,"pc_squared":pc*pc,"subnormal":(pc*pc).is_subnormal(),"reason":reason}));
    }
    Ok(())
}
