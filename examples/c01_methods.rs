//! Explicit small C-01 route driver. A saved displacement never carries material/contact history.
#[path = "c01_methods/controls.rs"]
mod controls;
use controls::{Options, Seed, Sphere};
use serde_json::json;
use soil_fem_3d::{
    Result,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Observer, Progress},
    stages,
};
use std::{
    cell::{Cell, RefCell},
    fs::{self, File},
    io::{BufWriter, Write},
    path::PathBuf,
};

struct Study {
    grid: Option<u32>,
    seed: Option<Seed>,
    guesses: Cell<usize>,
    arc: Option<RefCell<Sphere>>,
    sink: Option<RefCell<BufWriter<File>>>,
    history: Option<RefCell<BufWriter<File>>>,
}
impl Observer for Study {
    fn on_progress(&self, _: &Progress) -> Result<()> {
        Ok(())
    }
    fn c01_initial_displacement(
        &self,
        stage: usize,
        progress: f64,
    ) -> Result<Option<Vec<[f64; 3]>>> {
        if let Some(seed) = &self.seed {
            self.guesses.set(self.guesses.get() + 1);
            Ok(Some(seed.interpolate(stage, progress)?))
        } else {
            Ok(None)
        }
    }
    fn c01_arc_retry(&self, sample: soil_fem_3d::stage_control::ArcSample<'_>) -> Option<f64> {
        self.arc
            .as_ref()
            .and_then(|arc| arc.borrow_mut().retry(sample))
    }
    fn on_accepted_state(
        &self,
        stage: usize,
        progress: f64,
        state: &stages::Committed,
    ) -> Result<()> {
        if self
            .grid
            .is_some_and(|n| (progress * f64::from(n)).fract() == 0.)
        {
            let mut out = std::io::stdout().lock();
            serde_json::to_writer(&mut out, &(stage, progress, state))?;
            writeln!(out)?;
        }
        if let Some(history) = &self.history {
            let joint = state.joints.as_ref().expect("C-01 joint state");
            let field: Vec<_> = joint
                .points
                .iter()
                .map(|q| q.plastic_dissipation_kpa_m)
                .collect();
            let mut sink = history.borrow_mut();
            serde_json::to_writer(&mut *sink, &(stage, progress, field))?;
            writeln!(sink)?;
        }
        self.on_accepted_displacement(stage, progress, &state.displacement_m)
    }
    fn on_accepted_displacement(&self, stage: usize, progress: f64, u: &[[f64; 3]]) -> Result<()> {
        if let Some(sink) = &self.sink {
            let mut sink = sink.borrow_mut();
            serde_json::to_writer(&mut *sink, &(stage, progress, u))?;
            writeln!(sink)?;
        }
        Ok(())
    }
}

fn main() -> Result<()> {
    let mut args: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let options = if args.len() >= 6 && args[args.len() - 2].as_os_str() == "--methods" {
        let path = args.pop().unwrap();
        args.pop();
        serde_json::from_slice(&fs::read(path)?)?
    } else {
        Options::default()
    };
    assert!(args.len() == 4 || args.len() == 5);
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    let ground: Ground = serde_json::from_slice(&fs::read(&args[1])?)?;
    let case: stages::StageCase = serde_json::from_slice(&fs::read(&args[2])?)?;
    assert_eq!(mesh.tetrahedra10.len(), 12);
    let sink = args
        .get(4)
        .map(File::create_new)
        .transpose()?
        .map(BufWriter::new)
        .map(RefCell::new);
    let attempt_budget = case.settings.attempt_budget_policy.stage_limit(&mesh)?;
    let study = Study {
        grid: options.full_field_grid.map(|n| n.get()),
        history: options
            .dissipation_history
            .map(File::create_new)
            .transpose()?
            .map(BufWriter::new)
            .map(RefCell::new),
        sink,
        seed: options
            .seed
            .as_ref()
            .map(|p| Seed::read(p, mesh.nodes.len()))
            .transpose()?,
        guesses: Cell::new(0),
        arc: options
            .arc_scale_m
            .map(Sphere::new)
            .transpose()?
            .map(RefCell::new),
    };
    let run = stages::solve_observed(&mesh, ground, case, &study)?;
    let mut groups = Vec::new();
    let geo = &run.geostatic.as_ref().unwrap().attempts;
    let mut groups_of_attempts = vec![("geostatic", geo.iter().collect::<Vec<_>>())];
    for stage in ["shear", "open", "recontact", "reverse", "reload", "return"] {
        groups_of_attempts.push((
            stage,
            run.attempts.iter().filter(|a| a.stage == stage).collect(),
        ));
    }
    for (name, attempts) in groups_of_attempts {
        assert!(attempts.len() <= attempt_budget);
        let accepted: Vec<_> = attempts.iter().filter(|a| a.accepted).collect();
        let mut reasons = std::collections::BTreeMap::new();
        for a in attempts.iter().filter(|a| !a.accepted) {
            assert!(a.committed_state_unchanged_on_rejection);
            *reasons.entry(&a.reason).or_insert(0) += 1;
        }
        let force_max = accepted
            .iter()
            .map(|a| a.relative_force_residual.expect("accepted residual"))
            .fold(0., f64::max);
        for a in &accepted {
            let e = a.joint_path_error.as_ref().unwrap();
            assert!(
                a.relative_force_residual.unwrap() <= run.case.settings.relative_force_tolerance
                    && e.maximum_slip_error_m <= e.slip_budget_m
                    && e.maximum_dissipation_error_kpa_m <= e.dissipation_budget_kpa_m
            );
        }
        let rebalance = attempts
            .iter()
            .filter_map(|a| a.joint_path_error.as_ref())
            .filter_map(|e| e.endpoint_rebalance_from)
            .collect::<Vec<_>>();
        let mut summary = json!({"stage":name,"attempts":attempts.len(),"accepted":accepted.len(),
        "last_accepted":accepted.last().map(|a|a.progress_to),"rejections":reasons,"max_force":force_max});
        if !rebalance.is_empty() {
            let counts: [usize; 3] =
                std::array::from_fn(|i| rebalance.iter().filter(|r| r[i].is_some()).count());
            summary["reequilibrated_coarse_mid_fine"] = json!(counts);
        }
        groups.push(summary);
    }
    if let Some(s) = study.history {
        s.into_inner().flush()?;
    }
    if let Some(s) = study.sink {
        s.into_inner().flush()?;
    }
    let mut out = BufWriter::new(File::create_new(&args[3])?);
    serde_json::to_writer(
        &mut out,
        &json!({"complete":run.complete,"stop":run.failure_reason,"stages":groups,"seed_queries":study.guesses.get(),"arc_reviews":study.arc.as_ref().map_or(0,|a|a.borrow().reviews)}),
    )?;
    out.flush()?;
    Ok(())
}
