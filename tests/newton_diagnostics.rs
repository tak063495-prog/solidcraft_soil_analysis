use soil_fem_3d::{
    Error, Result,
    stage_control::{Observer, Progress},
    stages,
};
use std::cell::RefCell;
struct Diagnostic {
    values: RefCell<Vec<serde_json::Value>>,
    fail: bool,
}
impl Observer for Diagnostic {
    fn on_progress(&self, _: &Progress) -> Result<()> {
        Ok(())
    }
    fn diagnostic_enabled(&self, _: usize, _: f64) -> bool {
        true
    }
    fn on_diagnostic(&self, v: &serde_json::Value) -> Result<()> {
        if self.fail {
            return Err(Error::Input("diagnostic write failed".into()));
        }
        self.values.borrow_mut().push(v.clone());
        Ok(())
    }
}
fn normalize(v: &mut serde_json::Value) {
    match v {
        serde_json::Value::Object(m) => {
            m.remove("elapsed_seconds");
            for v in m.values_mut() {
                normalize(v);
            }
        }
        serde_json::Value::Array(a) => {
            for v in a {
                normalize(v);
            }
        }
        _ => (),
    }
}
#[test]
fn diagnostics_preserve_full_run_and_sink_failure_is_terminal() -> Result<()> {
    let mesh = serde_json::from_str(include_str!("../examples/joint-free-mesh.json"))?;
    let ground: soil_fem_3d::ground::Ground =
        serde_json::from_str(include_str!("../examples/ground-joint-mcc.json"))?;
    let mut case: stages::StageCase =
        serde_json::from_str(include_str!("../examples/stages-free-contact.json"))?;
    case.stages.truncate(1);
    case.verification_reject_once_stages.clear();
    let base = stages::solve(&mesh, ground.clone(), case.clone())?;
    let observer = Diagnostic {
        values: RefCell::new(vec![]),
        fail: false,
    };
    let observed = stages::solve_observed(&mesh, ground.clone(), case.clone(), &observer)?;
    let mut a = serde_json::to_value(&base)?;
    let mut b = serde_json::to_value(&observed)?;
    normalize(&mut a);
    normalize(&mut b);
    assert_eq!(a, b);
    let values = observer.values.borrow();
    assert!(values.iter().any(|v| v["kind"] == "iterate"));
    assert!(values.iter().any(|v| v["kind"] == "tangent"));
    assert!(values.iter().any(|v| v["kind"] == "line_search"));
    let broken = Diagnostic {
        values: RefCell::new(vec![]),
        fail: true,
    };
    assert!(matches!(
        stages::solve_observed(&mesh, ground, case, &broken),
        Err(Error::Control(_))
    ));
    Ok(())
}
