//! Diagnostic only: use the production point law on prescribed local paths.
use serde::{Deserialize, Serialize};
use soil_fem_3d::{
    Result,
    joint_law::{Material, PointState},
};
use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    material: Material,
    normal: [f64; 3],
    maximum_increment_m: f64,
    runs: Vec<Run>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Run {
    id: String,
    policy: String,
    fixed_count: usize,
    steps: Vec<Step>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Step {
    jump_m: [f64; 3],
    trial_jumps_m: Vec<[f64; 3]>,
}
#[derive(Serialize)]
struct Sample {
    endpoint_lower_bound: usize,
    actual_count: usize,
    plastic_dissipation_increment_kpa_m: f64,
    elastic_energy_kpa_m: f64,
    state: PointState,
}
#[derive(Serialize)]
struct OutputRun {
    id: String,
    policy: String,
    samples: Vec<Sample>,
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "REQUEST NEW_OUTPUT");
    let request: Request = serde_json::from_slice(&fs::read(&args[0])?)?;
    request.material.validate()?;
    let mut output = Vec::new();
    for run in request.runs {
        let mut state = PointState::initial(request.normal)?;
        let mut samples = Vec::new();
        for step in run.steps {
            let lower =
                request
                    .material
                    .partition(&state, step.jump_m, request.maximum_increment_m)?;
            let count = match run.policy.as_str() {
                "fixed" => run.fixed_count,
                "endpoint_lower" => lower,
                "frozen_rule_replay" => {
                    // stages.rs uses this rule during trial expansion. This is a
                    // prescribed trial sequence, not an equilibrium/Newton solve.
                    let mut frozen = run.fixed_count;
                    for candidate in step.trial_jumps_m.iter().chain([&step.jump_m]) {
                        let needed = request.material.partition(
                            &state,
                            *candidate,
                            request.maximum_increment_m,
                        )?;
                        if needed > frozen {
                            frozen = (2 * needed).min(50_000);
                        }
                    }
                    frozen
                }
                _ => panic!("unknown diagnostic partition policy"),
            };
            let previous_d = state.plastic_dissipation_kpa_m;
            state = request
                .material
                .integrate_partition(&state, step.jump_m, count)?;
            request.material.validate_state(&state)?;
            samples.push(Sample {
                endpoint_lower_bound: lower,
                actual_count: count,
                plastic_dissipation_increment_kpa_m: state.plastic_dissipation_kpa_m - previous_d,
                elastic_energy_kpa_m: request.material.energy(&state),
                state: state.clone(),
            });
        }
        output.push(OutputRun {
            id: run.id,
            policy: run.policy,
            samples,
        });
    }
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(Path::new(&args[1]))?;
    let mut stream = BufWriter::new(file);
    serde_json::to_writer(&mut stream, &output)?;
    writeln!(stream)?;
    Ok(())
}
