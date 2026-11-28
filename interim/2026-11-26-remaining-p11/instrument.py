"""Apply opt-in p11 controls and accepted-candidate observation to the scratch probe."""
from pathlib import Path
import sys

def edit(file, needle, replacement):
    p = Path(sys.argv[1])/"src"/file
    s = p.read_text()
    assert s.count(needle)==1
    p.write_text(s.replace(needle,replacement))

needle = '            let next = if p8_predict || p9_window || local {'
edit("stages.rs",needle,'''
            static CROSS: std::sync::OnceLock<(Vec<[f64;2]>,Vec<([f64;2],f64)>)> = std::sync::OnceLock::new();
            let nominal = if p8_reload && let Some(path)=std::env::var_os("P11_CROSS") {
                let (bands, actions)=CROSS.get_or_init(|| {
                    let plan: serde_json::Value=serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
                    let bands=serde_json::from_value(plan["bands"].clone()).unwrap();
                    let actions=plan["controls"].as_array().unwrap().iter().map(|r|
                        (serde_json::from_value(r["span"].clone()).unwrap(),r["target"].as_f64().unwrap())).collect();
                    (bands,actions)
                });
                if bands.iter().any(|[a,b]| *a<=progress && progress<*b) { nominal } else {
                    actions.iter().filter(|([a,b],_)| progress<*a && *b<=nominal)
                        .map(|(_,x)| *x).fold(nominal,f64::min)
                }
            } else { nominal };
'''+needle)
needle='        if std::env::var_os("P10_KNOTS").is_some() && coarse.state.stage_index==3 {'
edit("stages_joint_path_error.rs",needle,'''
        if std::env::var_os("P11_OBSERVE").is_some() && coarse.state.stage_index==3 {
            let points=[a,b,d].map(|j| j.points.iter().map(|p|
                serde_json::json!([p.status,p.plastic_slip_m])).collect::<Vec<_>>());
            eprintln!("{}",serde_json::json!({"p11_span":[start,midpoint_progress,load.progress],"all":points}));
        }
'''+needle)
