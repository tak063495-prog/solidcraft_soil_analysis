// Homogeneity of the committed MCC stress update. Only stress, pc and
// the stress tangent have degree one. Histories of strain and volume have degree zero.
use serde_json::{Value, json};
use soil_fem_3d::mcc::{Parameters, State};
use std::io::{self, BufReader, BufWriter};

type Request = (
    Vec<Parameters>,
    Vec<usize>,
    Vec<Value>,
    Vec<[f64; 6]>,
    Vec<usize>,
);

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (laws, ids, old, strains, counts): Request =
        serde_json::from_reader(BufReader::new(io::stdin()))?;
    assert_eq!(
        [ids.len(), old.len(), strains.len(), counts.len()],
        [1296; 4]
    );
    let mut output = Vec::with_capacity(old.len());
    for (i, mut record) in old.into_iter().enumerate() {
        let exponent = record
            .as_object_mut()
            .unwrap()
            .remove("stress_exponent")
            .and_then(|e| e.as_i64())
            .unwrap_or(0);
        let mut state: State = serde_json::from_value(record)?;
        let shift = state.preconsolidation_pressure_kpa.log2().floor() as i32;
        let factor = 2_f64.powi(-shift);
        assert!(factor.is_normal());
        for value in &mut state.effective_stress_compression_positive_kpa {
            *value *= factor;
        }
        state.preconsolidation_pressure_kpa *= factor;
        let law = &laws[ids[i]];
        let (update, tangent) = law
            .implicit_stress_tangent_partition(&state, strains[i], counts[i])
            .map_err(|e| format!("active point {i}, tangent: {e}"))?;
        let primal = law.integrate_partition(&state, strains[i], counts[i])?;
        let mut next = serde_json::to_value(&primal.state)?;
        assert_eq!(serde_json::to_value(&update.state)?, next);
        law.validate_state(&primal.state)?;
        next["stress_exponent"] = json!(exponent + i64::from(shift));
        output.push((next, tangent.as_slice().to_vec(), primal.plastic_substeps));
    }
    serde_json::to_writer(BufWriter::new(io::stdout()), &output)?;
    Ok(())
}
