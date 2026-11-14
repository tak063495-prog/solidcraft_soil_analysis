
/// Observe the identical N=0 fixed-partition path, retaining its pc reference.
pub fn diagnostic_path(
    p: &mcc::Parameters,
    s: &mcc::State,
    strain: [f64; 6],
    steps: usize,
) -> Result<serde_json::Value> {
    p.validate_state(s)?;
    if steps == 0 || steps > 100_000 || !strain.iter().all(|x| x.is_finite()) {
        return Err(Error::Input("diagnostic partition bounds".into()));
    }
    let mut state = Tracked::<0>::from_state(s);
    let de = strain.map(|v| D::<0>::c(v / steps as f64));
    let mut records = Vec::new();
    for _ in 0..steps {
        let trial = elastic(p, &state, de)?;
        let before = yield_value(p, &state).v;
        let elastic_yield = yield_value(p, &trial).v;
        let (next, plastic) = step(p, &state, de, false)?;
        let values = next.values();
        let (pressure, q) = mcc::invariants(&values.effective_stress_compression_positive_kpa);
        records.push(serde_json::json!({"plastic":plastic,"normalized_yield_before":before,
            "elastic_trial_normalized_yield":elastic_yield,"normalized_yield_after":yield_value(p,&next).v,
            "p_kpa":pressure,"q_kpa":q,"pc_kpa":values.preconsolidation_pressure_kpa}));
        state = next;
    }
    Ok(serde_json::json!({"substeps":records,"state":state.values()}))
}
