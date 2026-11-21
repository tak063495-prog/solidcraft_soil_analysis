//! Local observation only: expose every MC return on the saved fixed partition.
use serde::{Deserialize, Serialize};
use soil_fem_3d::{Result, mc};
use std::{collections::BTreeMap, fs, io::Write};

#[derive(Deserialize)]
struct Request {
    route: String,
    increment: usize,
    gp: usize,
    parameters: mc::Parameters,
    old: mc::State,
    strain: [f64; 6],
    count: usize,
}
#[derive(Serialize)]
struct Record {
    route: String,
    increment: usize,
    gp: usize,
    count: usize,
    returns: BTreeMap<String, usize>,
    state: mc::State,
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "REQUEST_JSON NEW_JSON");
    let requests: Vec<Request> = serde_json::from_slice(&fs::read(&args[0])?)?;
    let mut records = Vec::new();
    for r in requests {
        assert!(r.count > 0 && r.count <= 50_000);
        let mut state = r.old;
        let strain = r.strain.map(|x| x / r.count as f64);
        let mut returns = BTreeMap::new();
        for _ in 0..r.count {
            state = r.parameters.integrate_partition(&state, strain, 1)?.state;
            let category = serde_json::to_value(&state)?["last_return"]
                .as_str()
                .expect("MC category")
                .to_owned();
            *returns.entry(category).or_insert(0) += 1;
        }
        records.push(Record {
            route: r.route,
            increment: r.increment,
            gp: r.gp,
            count: r.count,
            returns,
            state,
        });
    }
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[1])?;
    serde_json::to_writer(&mut out, &records)?;
    writeln!(out)?;
    Ok(())
}
