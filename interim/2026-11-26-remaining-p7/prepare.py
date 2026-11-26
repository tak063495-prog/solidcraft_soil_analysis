"""Install the opt-in experiment and passive diagnostics in a disposable copy of the sources."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
probe = Path(sys.argv[1])


def edit(path, pairs):
    s = path.read_text()
    for a,b in pairs:
        assert s.count(a)==1, a
        s = s.replace(a,b)
    path.write_text(s)


joint = probe/"src/stages_joint_path_error.rs"
edit(joint, [
    ("use super::*;", """use super::*;
thread_local! { static SPENT: std::cell::Cell<(usize, usize)> = const { std::cell::Cell::new((usize::MAX, 0)) }; }
pub(super) fn enabled() -> bool { std::env::var_os("P7_SWITCH_HALF").is_some() }
pub(super) fn charge(stage:usize, kind:&str) -> Result<()> {
    SPENT.with(|cell| {
        let (last,count)=cell.get();
        let used=if last==stage {count} else {0};
        if used>=1787 { return Err(Error::Control("P7 stage budget exhausted".into())); }
        cell.set((stage,used+1));
        if kind!="outer" { eprintln!("{}",serde_json::json!({"charge":kind,"stage":stage,"spent":used+1})); }
        Ok(())
    })
}"""),
    ("        let endpoint = |parent: &Committed, load: Target|", """        let p7 = enabled();
        let load = target.clone();
        let mid_load = midpoint.clone();
        let endpoint = |parent: &Committed, load: Target, extra: bool|"""),
    ("            let initial = self.trial", """            if extra { charge(load.stage, "quarter")?; }
            let initial = self.trial"""),
    ("            let mut corrected = self.trial", """            if extra { charge(load.stage, "quarter_rebalance")?; }
            let mut corrected = self.trial"""),
    ("endpoint(old, target.clone())?", "endpoint(old, target.clone(),false)?"),
    ("endpoint(old, midpoint)?", "endpoint(old, midpoint,false)?"),
    ("endpoint(&half.state, target)?", "endpoint(&half.state, target,false)?"),
    ("        let trial = if control.resolve_slip_path", (HERE/"switch_half.rs.inc").read_text()+
     "\n        let trial = if control.resolve_slip_path"),
])
stages = probe/"src/stages.rs"
edit(stages, [
    ("            count += 1;\n            let next = if matches!(stage,", """            if joint_path_error::enabled() && joint_path_error::charge(stage_index+1, "outer").is_err() {
                failure = Some("P7 stage budget exhausted".into());
                break;
            }
            count += 1;
            let next = if matches!(stage,"""),
    ("let attempt_budget = settings.attempt_budget_policy.stage_limit(model.mesh)?;",
     "let mut attempt_budget = settings.attempt_budget_policy.stage_limit(model.mesh)?; if joint_path_error::enabled() {attempt_budget = attempt_budget.min(1787);}"),
])
driver = probe/"examples/c01_methods.rs"
edit(driver, [("        let accepted: Vec<_>", """        if std::env::var_os("P7_SWITCH_HALF").is_some() {
            for a in &attempts {
                let path = a.joint_path_error.as_ref().map(|e| serde_json::json!({
                    "max_slip":e.maximum_slip_error_m,"slip_budget":e.slip_budget_m,
                    "max_work":e.maximum_dissipation_error_kpa_m,"work_budget":e.dissipation_budget_kpa_m,
                    "rebalance":e.endpoint_rebalance_from,"event":e.friction_event_target,
                    "proposed_increment":e.proposed_increment,"mid":e.midpoint_progress
                }));
                eprintln!("{}", serde_json::json!({"stage":name,"from":a.progress_from,"to":a.progress_to,
                    "accepted":a.accepted,"reason":a.reason,"parent_unchanged":a.committed_state_unchanged_on_rejection,
                    "force":a.relative_force_residual,"audit":path}));
            }
        }
        let accepted: Vec<_>""")])
