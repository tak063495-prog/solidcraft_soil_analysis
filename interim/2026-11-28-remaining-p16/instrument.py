"""Install the accepted-parent reload predictor in a disposable p15 build."""
from pathlib import Path
import sys

p = Path(sys.argv[1]) / "src/stages.rs"
s = p.read_text()
key = '        let outer_gap ='
assert s.count(key) == 1
s = s.replace(key, '''        let p16 = p8_reload && std::env::var_os("P16_SLIP").is_some();
        let mut last_slip: Option<(f64,f64)> = None;
''' + key)
key = '            let nominal = if p14 &&'
assert s.count(key) == 1
s = s.replace(key, '''            let nominal = if p16 && (0.625..0.6875).contains(&progress) {
                let mut to = nominal;
                if let Some((last_width,ratio)) = last_slip {
                    // Step doubling estimates a quadratic local slip difference.
                    // Try the existing half-width cut before spending a trial.
                    let predicted = ratio*((nominal-progress)/last_width).powi(2);
                    if predicted > 1. {
                        to = progress + (nominal-progress)/2.;
                    }
                    eprintln!("{}",serde_json::json!({"p16_parent":progress,"attempt":count,
                        "last_width":last_width,"last_slip_ratio":ratio,"predicted_ratio":predicted,
                        "nominal":nominal,"to":to}));
                }
                to
            } else {nominal};
''' + key)
key = '                    progress = next;'
assert s.count(key) == 1
s = s.replace(key, '''                    if p16 {
                        last_slip = joint_path_error.as_ref().map(|a|
                            (next-from,a.maximum_slip_error_m/a.slip_budget_m));
                    }
''' + key)
p.write_text(s)

p = Path(sys.argv[1]) / "examples/c01_methods.rs"
s = p.read_text()
key = '        if self\n            .grid'
assert s.count(key) == 1
s = s.replace(key, '''        if stage == 3 && std::env::var_os("P16_FIELDS").is_some() {
            let pts = &state.joints.as_ref().unwrap().points;
            let vals: Vec<_> = pts.iter().map(|p|
                (p.traction_kpa,p.jump_m,p.plastic_slip_m,p.accumulated_plastic_slip_m)).collect();
            println!("{}",serde_json::json!({"t":progress,"u":state.displacement_m,
                "load":state.applied_nodal_load_kn,"points":vals}));
        }
''' + key)
p.write_text(s)
