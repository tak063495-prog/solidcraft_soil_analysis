//! Evaluate the unchanged MCC increment and its consistent tangent on a fixed parent.
use soil_fem_3d::mcc::{Parameters, State};
use std::io::{self, BufReader, BufWriter};

type Request = (
    Vec<Parameters>,
    Vec<usize>,
    Vec<State>,
    Vec<[f64; 6]>,
    Vec<usize>,
);

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (laws, ids, old, strains, counts): Request =
        serde_json::from_reader(BufReader::new(io::stdin()))?;
    assert_eq!(old.len(), 1296);
    assert_eq!([ids.len(), strains.len(), counts.len()], [old.len(); 3]);
    let mut result = Vec::with_capacity(old.len());
    for i in 0..old.len() {
        let law = &laws[ids[i]];
        let (update, tangent) =
            law.implicit_stress_tangent_partition(&old[i], strains[i], counts[i])?;
        let primal = law.integrate_partition(&old[i], strains[i], counts[i])?;
        assert_eq!(
            serde_json::to_value(&update.state)?,
            serde_json::to_value(&primal.state)?
        );
        result.push((
            primal.state,
            tangent.as_slice().to_vec(),
            primal.plastic_substeps,
        ));
    }
    serde_json::to_writer(BufWriter::new(io::stdout()), &result)?;
    Ok(())
}
