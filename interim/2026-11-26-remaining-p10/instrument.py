"""Instrument the disposable P9 probe. Tables are read once per process."""
from pathlib import Path
import sys

src = Path(sys.argv[1]) / "src"
p = src / "stages.rs"
s = p.read_text()
key = "            let next = if p8_predict || p9_window {"
assert s.count(key)==1
s = s.replace(key, '''            static WINDOWS: std::sync::OnceLock<Vec<[f64;2]>> = std::sync::OnceLock::new();
            let local = p8_reload && std::env::var_os("P10_WINDOWS").is_some_and(|path| {
                let bands = WINDOWS.get_or_init(|| {
                    let data: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
                    serde_json::from_value(data["windows"].clone()).unwrap()
                });
                bands.iter().any(|[lo,hi]| *lo<=progress && progress<*hi)
            });
            static KNOTS: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
            let nominal = if stage.name()=="recontact" && let Some(path)=std::env::var_os("P10_KNOTS") {
                let knots = KNOTS.get_or_init(|| {
                    let data: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
                    serde_json::from_value(data["knots"].clone()).unwrap()
                });
                knots.iter().copied().find(|x| *x>progress).map_or(nominal,|x| x.min(nominal))
            } else {nominal};
            let next = if p8_predict || p9_window || local {''')
p.write_text(s)
p = src / "stages_joint_path_error.rs"
s = p.read_text()
key = "        let trial = if control.resolve_slip_path"
assert s.count(key)==1
p.write_text(s.replace(key,'''        if std::env::var_os("P10_LABELS").is_some() && coarse.state.stage_index==5 {
            let labels=[a,b,d].map(|j| j.points.iter().map(|p| &p.status).collect::<Vec<_>>());
            eprintln!("{}",serde_json::json!({"span":[start,midpoint_progress,load.progress],"labels":labels}));
        }
        if std::env::var_os("P10_KNOTS").is_some() && coarse.state.stage_index==3 {
            let points=[a,b,d].map(|j| [3,5,13].map(|i| {
                let p=&j.points[i];
                serde_json::json!([p.status,p.plastic_slip_m])
            }));
            eprintln!("{}",serde_json::json!({"contact_span":[start,midpoint_progress,load.progress],"points":points}));
        }
''' + key))
