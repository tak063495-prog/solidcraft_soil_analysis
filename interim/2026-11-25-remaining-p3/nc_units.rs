//! Same midpoint primal and tangent as P2, evaluated in one common stress unit.
//! The power-of-two change is exact for normal binary64 operands. External
//! loads, displacements, strain, lambda and physical residuals remain in P2 units.
use soil_fem_3d::mcc::{Parameters, State};
use std::io::{self, BufReader, BufWriter};

type Request = (
    Vec<Parameters>,
    Vec<usize>,
    Vec<State>,
    Vec<[f64; 6]>,
    Vec<usize>,
);

fn scale(state: &mut State, factor: f64) {
    state.preconsolidation_pressure_kpa *= factor;
    for s in &mut state.effective_stress_compression_positive_kpa {
        *s *= factor;
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (laws, ids, mut old, strains, counts): Request =
        serde_json::from_reader(BufReader::new(io::stdin()))?;
    assert_eq!(
        [ids.len(), old.len(), strains.len(), counts.len()],
        [1296; 4]
    );
    let factor = 2.0_f64.powi(448);
    let mut result = Vec::with_capacity(old.len());
    for (i, parent) in old.iter_mut().enumerate() {
        let original = &laws[ids[i]];
        let mut law = original.clone();
        law.reference_pressure_kpa *= factor;
        scale(parent, factor);
        let (update, tangent) = law
            .implicit_stress_tangent_partition(parent, strains[i], counts[i])
            .map_err(|e| format!("active point {i}, tangent: {e}"))?;
        let primal = law.integrate_partition(parent, strains[i], counts[i])?;
        assert_eq!(
            serde_json::to_value(&update.state)?,
            serde_json::to_value(&primal.state)?
        );
        let mut state = primal.state;
        scale(&mut state, factor.recip());
        if cfg!(roundtrip_units) {
            // Repeat the unmodified admissibility check on the stored kPa
            // values after a round trip to the calculation's common unit.
            // A zero or nonfinite physical pressure cannot survive this check.
            let mut checked = state.clone();
            scale(&mut checked, factor);
            law.validate_state(&checked)
                .map_err(|e| format!("active point {i}, round-trip admissibility: {e}"))?;
        } else {
            original.validate_state(&state)?;
        }
        let c: Vec<_> = tangent.as_slice().iter().map(|v| v / factor).collect();
        result.push((state, c, primal.plastic_substeps));
    }
    serde_json::to_writer(BufWriter::new(io::stdout()), &result)?;
    Ok(())
}
