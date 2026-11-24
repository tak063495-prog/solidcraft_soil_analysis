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
    #[serde(default)]
    save_substeps: bool,
    #[serde(default)]
    locate_on: Option<usize>,
}
#[derive(Serialize)]
struct Substep {
    stress_kpa: [f64; 6],
    plastic_increment: [f64; 6],
    label: mc::Return,
}

#[derive(Serialize)]
struct Record {
    route: String,
    increment: usize,
    gp: usize,
    count: usize,
    returns: BTreeMap<String, usize>,
    state: mc::State,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    substeps: Vec<Substep>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    boundaries: Vec<Boundary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    search_steps: Option<usize>,
}

#[derive(Serialize)]
struct Boundary {
    partition: usize,
    from: mc::Return,
    to: mc::Return,
    fraction: f64,
}

// Trace the specified uniform partition from the actual incoming state. Bisect
// within a switching substep from its unchanged parent, using the native return.
fn boundaries(
    r: &Request,
    start: &mc::State,
    count: usize,
    calls: &mut usize,
) -> Result<Vec<Boundary>> {
    let mut found = Vec::new();
    let mut state = start.clone();
    let strain = r.strain.map(|x| x / count as f64);
    for i in 0..count {
        let next = r.parameters.integrate_partition(&state, strain, 1)?.state;
        *calls += 1;
        use mc::Return::{ExtensionEdge, ShearFace};
        if matches!(
            (state.last_return, next.last_return),
            (ShearFace, ExtensionEdge) | (ExtensionEdge, ShearFace)
        ) {
            let (mut lo, mut hi) = (0., 1.);
            for _ in 0..48 {
                let mid = (lo + hi) / 2.;
                let probe = r
                    .parameters
                    .integrate_partition(&state, strain.map(|x| mid * x), 1)?;
                *calls += 1;
                if probe.state.last_return == next.last_return {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            found.push(Boundary {
                partition: count,
                from: state.last_return,
                to: next.last_return,
                fraction: (i as f64 + (lo + hi) / 2.) / count as f64,
            });
        }
        state = next;
    }
    Ok(found)
}
fn main() -> Result<()> {
    let mut args: Vec<_> = std::env::args().skip(1).collect();
    let history = args.last().is_some_and(|x| x == "--history");
    if history {
        args.pop();
    }
    assert_eq!(args.len(), 2, "REQUEST_JSON NEW_JSON");
    let requests: Vec<Request> = serde_json::from_slice(&fs::read(&args[0])?)?;
    let mut records = Vec::new();
    let mut parents = BTreeMap::<(String, usize), (usize, mc::State)>::new();
    for r in requests {
        assert!(r.count > 0 && r.count <= 50_000);
        let key = (r.route.clone(), r.gp);
        let mut state = if history && r.increment != 0 {
            let (index, parent) = parents.remove(&key).expect("missing preceding frame");
            assert_eq!(index + 1, r.increment);
            parent
        } else {
            r.old.clone()
        };
        let mut calls = 0;
        let mut found = Vec::new();
        if let Some(native) = r.locate_on {
            assert!(history && native > 0 && native <= r.count);
            found = boundaries(&r, &state, native, &mut calls)?;
            if native != r.count {
                found.extend(boundaries(&r, &state, r.count, &mut calls)?);
            }
        }
        let mut cuts: Vec<_> = found.iter().map(|b| b.fraction).collect();
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let mut cuts = cuts.into_iter().peekable();
        let strain = r.strain.map(|x| x / r.count as f64);
        let mut returns = BTreeMap::new();
        let mut substeps = Vec::new();
        let mut accepted = 0;
        for i in 0..r.count {
            let end = (i + 1) as f64 / r.count as f64;
            let mut start = i as f64 / r.count as f64;
            let mut split = false;
            loop {
                let stop = match cuts.next_if(|&x| x <= end) {
                    Some(x) if x <= start => continue,
                    Some(x) => x,
                    None => end,
                };
                let delta = if split || stop < end {
                    r.strain.map(|x| x * (stop - start))
                } else {
                    strain
                };
                accepted += 1;
                let previous = state.plastic_strain_compression_positive;
                state = r.parameters.integrate_partition(&state, delta, 1)?.state;
                let category = serde_json::to_value(&state)?["last_return"]
                    .as_str()
                    .expect("MC category")
                    .to_owned();
                *returns.entry(category).or_insert(0) += 1;
                if r.save_substeps {
                    substeps.push(Substep {
                        stress_kpa: state.effective_stress_compression_positive_kpa,
                        plastic_increment: std::array::from_fn(|i| {
                            state.plastic_strain_compression_positive[i] - previous[i]
                        }),
                        label: state.last_return,
                    });
                }
                if stop == end {
                    break;
                }
                start = stop;
                split = true;
            }
        }
        if history {
            parents.insert(key, (r.increment, state.clone()));
        }
        records.push(Record {
            route: r.route,
            increment: r.increment,
            gp: r.gp,
            count: accepted,
            boundaries: found,
            search_steps: r.locate_on.map(|_| calls),
            returns,
            state,
            substeps,
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
