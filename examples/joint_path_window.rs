//! Observe a 12-element C-01 prefix without resetting the per-stage attempt count.
use soil_fem_3d::{
    Result,
    ground::Ground,
    mesh::Mesh,
    stage_control::{Kind, Observer, Progress},
    stages::{StageCase, solve_observed},
};
use std::{cell::Cell, fs, io::BufWriter};

struct Window {
    stage: String,
    end: f64,
    reached: Cell<bool>,
}

impl Observer for Window {
    fn on_progress(&self, event: &Progress) -> Result<()> {
        if event.kind == Kind::Accepted && event.stage == self.stage && event.progress >= self.end {
            self.reached.set(true);
        }
        Ok(())
    }

    fn stop_requested(&self) -> Result<bool> {
        Ok(self.reached.get())
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 6, "MESH GROUND CASE STAGE END NEW_RUN_JSON");
    let mesh: Mesh = serde_json::from_slice(&fs::read(&args[0])?)?;
    assert_eq!(mesh.tetrahedra.len(), 12);
    let ground: Ground = serde_json::from_slice(&fs::read(&args[1])?)?;
    let case: StageCase = serde_json::from_slice(&fs::read(&args[2])?)?;
    let window = Window {
        stage: args[3].clone(),
        end: args[4].parse().expect("numeric end fraction"),
        reached: Cell::new(false),
    };
    let run = solve_observed(&mesh, ground, case, &window)?;
    serde_json::to_writer(BufWriter::new(fs::File::create_new(&args[5])?), &run)?;
    assert!(
        matches!(
            run.failure_reason.as_deref(),
            Some(
                "stage attempt budget exceeded"
                    | "geostatic attempt budget exhausted"
                    | "checkpoint_requested"
            )
        ),
        "{:?}",
        run.failure_reason
    );
    Ok(())
}
