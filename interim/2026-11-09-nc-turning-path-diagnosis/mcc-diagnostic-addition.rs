
/// Read-only isolated diagnostic; absent from production.
pub fn diagnostic_midpoint_path(
    p: &Parameters,
    s: &State,
    strain: [f64; 6],
    steps: usize,
) -> Result<serde_json::Value> {
    midpoint::diagnostic_path(p, s, strain, steps)
}
