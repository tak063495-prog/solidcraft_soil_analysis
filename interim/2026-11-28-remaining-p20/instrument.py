"""Add a read-only, pre-trial local slip observer to the disposable p19 probe."""
from pathlib import Path
import sys

root = Path(sys.argv[1])
path = root/"src/stages.rs"
source = path.read_text()
anchor = '            let target: Vec<_> = if matches!(stage, Stage::ResetDisplay { .. }) {'
assert source.count(anchor) == 1
source = source.replace(anchor, '''            if p17 && (0.5625..0.625).contains(&progress)
                && std::env::var_os("P20_SLIP_OBSERVE").is_some()
                && let Some((previous, jumps)) = &p18_prior {
                use nalgebra::Vector3;
                let system = case.joint_system.as_ref().unwrap();
                let control = system.global_step_control.as_ref().unwrap().path_error.as_ref().unwrap();
                let mat = &system.joints[0].material;
                let ratio = (next-progress)/(progress-previous);
                let mut error = Vec::new();
                for (p,j) in committed.joints.as_ref().unwrap().points.iter().zip(jumps) {
                    let delta = (p.elastic_jump()-j)*ratio;
                    let start = Vector3::from(p.jump_m);
                    let mid = (start+delta*0.5).into();
                    let end = (start+delta).into();
                    let local = |old: &crate::joint_law::PointState, target| {
                        let n = mat.partition(old,target,system.maximum_increment_m)?;
                        mat.integrate_partition(old,target,n.max(control.minimum_contact_substeps.unwrap_or(1)))
                    };
                    let coarse = local(p,end)?;
                    let half = local(p,mid)?;
                    let fine = local(&half,end)?;
                    error.push(std::array::from_fn::<_,3,_>(|i|
                        coarse.plastic_slip_m[i]-fine.plastic_slip_m[i]));
                }
                eprintln!("{}",serde_json::json!({"p20_predict":count,"parent":progress,
                    "previous":previous,"to":next,"error_m":error,
                    "budget_m":stage.step(&case.settings).powi(2)*control.slip_scale_m}));
            }
''' + anchor)
path.write_text(source)
path = root/"src/stages_joint_path_error.rs"
source = path.read_text()
anchor = "        let trial = if control.resolve_slip_path"
assert source.count(anchor) == 1
path.write_text(source.replace(anchor, '''        if coarse.state.stage_index == 5 && (0.5625..0.625).contains(&start)
            && std::env::var_os("P20_SLIP_OBSERVE").is_some() {
            eprintln!("{}",serde_json::json!({"p20_actual":start,"to":fine.state.stage_progress,
                "error_m":audit.slip_errors_m,"budget_m":audit.slip_budget_m}));
        }
''' + anchor))
