// Diagnostic executable: only calls existing public production APIs.
use serde::de::DeserializeOwned;
use soil_fem_3d::{
    Result, beams,
    checkpoint::{Checkpoint, hash_json},
    ground::Ground,
    mesh::Mesh,
    stages::{self, Run, StageCase},
};
use std::{
    fs::{self, File, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
};

fn read<T: DeserializeOwned>(p: &str) -> Result<T> {
    Ok(serde_json::from_slice(&fs::read(p)?)?)
}
fn write(p: &Path, v: &impl serde::Serialize) -> Result<()> {
    let f = OpenOptions::new().write(true).create_new(true).open(p)?;
    serde_json::to_writer(BufWriter::new(f), v)?;
    Ok(())
}
fn count(r: &Run) -> usize {
    r.attempts.iter().filter(|a| a.accepted).count()
        + r.geostatic
            .as_ref()
            .map_or(0, |g| g.attempts.iter().filter(|a| a.accepted).count())
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a[1] == "rehash" {
        let mut c: Checkpoint = read(&a[2])?;
        c.run_sha256 = hash_json(&c.run)?;
        return write(Path::new(&a[3]), &c);
    }
    let m: Mesh = read(&a[2])?;
    let g: Ground = read(&a[3])?;
    let c: StageCase = read(&a[4])?;
    let out = Path::new(&a[5]);
    fs::create_dir(out)?;
    if a[1] == "plain" {
        let r = stages::solve(&m, g, c)?;
        assert!(r.complete, "{:?}", r.failure_reason);
        return write(&out.join("run.json"), &r);
    }
    let model = beams::Model::new(&m, c.beam_system.clone().unwrap())?;
    let quadrature: Vec<_> = (0..m.tetrahedra.len())
        .map(|i| m.integration(i))
        .collect::<Result<_>>()?;
    let mut history = BufWriter::new(File::create(out.join("accepted.jsonl"))?);
    let mut prior = 0;
    let mut r = if a[1] == "resume" {
        let saved: Checkpoint = read(&a[6])?;
        prior = count(&saved.run);
        saved.resume(&m, g.clone(), c.clone(), Some(1))?
    } else {
        stages::solve_until(&m, g.clone(), c.clone(), 1)?
    };
    loop {
        let n = count(&r);
        if n > prior {
            assert_eq!(n, prior + 1);
            let geo = r.attempts.iter().all(|x| !x.accepted);
            let s = if geo {
                &r.geostatic.as_ref().unwrap().last_committed
            } else {
                &r.last_committed
            };
            let progress = if geo {
                r.geostatic.as_ref().unwrap().progress
            } else {
                s.stage_progress
            };
            let e = model.evaluate(
                s.beams.as_ref().unwrap(),
                &s.displacement_m,
                &vec![true; m.nodes.len()],
                s.boundary_factor,
            )?;
            // Production quadrature + actual saved production stress. This is
            // observation/reassembly, never a new independent solution.
            let mut soil_internal = vec![0.; 3 * m.nodes.len()];
            let mut point = 0;
            for (ei, points) in quadrature.iter().enumerate() {
                let nodes = m.cell_nodes(ei);
                for q in points {
                    let stress = s.material_states[point].stress();
                    for (i, node) in nodes.iter().enumerate() {
                        for a in 0..3 {
                            soil_internal[3 * node + a] -= q.weight_m3
                                * (0..6).map(|j| q.b[(j, 3 * i + a)] * stress[j]).sum::<f64>();
                        }
                    }
                    point += 1;
                }
            }
            serde_json::to_writer(
                &mut history,
                &serde_json::json!({
                    "accepted_index": n, "stage_index": if geo {0} else {s.stage_index},
                    "progress": progress, "state": s, "beam_results": e.report,
                    "soil_internal_kn": soil_internal, "interface_dofs": e.interface_dofs,
                    "interface_internal_kn": e.interface_internal_kn.as_slice()
                }),
            )?;
            history.write_all(b"\n")?;
            history.flush()?;
            prior = n;
        }
        if r.complete {
            write(&out.join("run.json"), &r)?;
            break;
        }
        assert_eq!(r.failure_reason.as_deref(), Some("checkpoint_requested"));
        let saved = Checkpoint::seal(&m, r)?;
        if a[1] == "solve"
            && a.get(6)
                .is_some_and(|x| x.split(',').any(|v| v.parse::<usize>().ok() == Some(n)))
        {
            write(&out.join(format!("checkpoint-{n}.json")), &saved)?;
        }
        r = saved.resume(&m, g.clone(), c.clone(), Some(1))?;
    }
    Ok(())
}
