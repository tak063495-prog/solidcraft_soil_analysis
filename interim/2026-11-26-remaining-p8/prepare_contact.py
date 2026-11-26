"""Install P8 diagnostics in the disposable P7 probe after running P7 prepare.py."""
import json
from pathlib import Path
import sys

from freeze import HERE

probe = Path(sys.argv[1])
table = json.loads((HERE/"c01/predictor.json").read_text())
caps = json.dumps(table["caps"])
casefile = probe/"src/stages.rs"
s = casefile.read_text()


def replace(a,b):
    global s
    assert s.count(a)==1,a
    s=s.replace(a,b)


replace("pub(crate) fn resume_observed(", "pub fn resume_observed(")
a = "                    if stop_after.is_some_and(|limit| accepted_this_run >= limit)"
# Only the construction-stage callback with a complete accepted reverse endpoint.
i = s.rindex(a)
s=s[:i]+s[i:].replace(a, """                    if (std::env::var_os("P8_PREFIX").is_some() && stage.name()=="reverse" && next==1.)
                        || stop_after.is_some_and(|limit| accepted_this_run >= limit)""",1)
a = "        while progress < 1. - 1e-12 {\n            if count >= attempt_budget"
replace(a,"""        let p8_reload = stage.name()=="reload";
        let p8_predict = p8_reload && std::env::var_os("P8_PREDICT").is_some();
        let mut p8_parents = Vec::new();
        if p8_reload && let Some(dir) = std::env::var_os("P8_REPLAY") {
            p8_parents = std::fs::read_dir(dir)?.map(|e| e.map(|x| x.path())).collect::<std::io::Result<Vec<_>>>()?;
            p8_parents.sort();
            assert_eq!(p8_parents.len(),163);
        }
        while progress < 1. - 1e-12 {
            if !p8_parents.is_empty() {
                if count == p8_parents.len() {
                    failure = Some("P8 paired observations finished".into());
                    break;
                }
                let bytes = std::fs::read(&p8_parents[count])?;
                let row: (usize, f64, f64, Committed) = serde_json::from_slice(&bytes)?;
                assert_eq!(bytes,serde_json::to_vec(&row)?);
                progress = row.1;
                increment = row.2 - row.1;
                committed = row.3;
            }
            if count >= attempt_budget""")
a="            let target: Vec<_> = if matches!(stage, Stage::ResetDisplay"
replace(a,"""            let nominal = if p8_parents.is_empty() { next } else {
                let bytes = std::fs::read(&p8_parents[count-1])?;
                let row: (usize,f64,f64,Committed) = serde_json::from_slice(&bytes)?;
                row.2
            };
            let next = if p8_predict {
                const CAPS: [f64;256] = """+caps+""";
                let cap=CAPS[((progress*256.).floor() as usize).min(255)];
                if nominal-progress > cap {
                    crate::joints::GlobalStepControl::next(progress,cap,stage.step(&case.settings)).min(nominal)
                } else { nominal }
            } else { nominal };
            let target: Vec<_> = if matches!(stage, Stage::ResetDisplay""")
a = "            let from = progress;\n            let linear_solves"
replace(a,"""            if p8_reload && let Some(dir) = std::env::var_os("P8_SAVE_PARENTS")
                && (trial.is_err() || force_reject || switch_cutback || path_error_cutback) {
                let path=std::path::Path::new(&dir).join(format!("{count:04}.json"));
                serde_json::to_writer(std::fs::File::create_new(path)?, &(count, progress, next, &committed))?;
            }
            if !p8_parents.is_empty() {
                let ok=trial.is_ok() && !force_reject && !switch_cutback && !path_error_cutback;
                let unchanged=before==serde_json::to_vec(&committed)?;
                assert!(unchanged);
                eprintln!("{}",serde_json::json!({"file":p8_parents[count-1].file_name().unwrap().to_str(),
                    "from":progress,"nominal":nominal,"to":next,"accepted":ok,"unchanged":unchanged,
                    "switch_reject":switch_cutback,"audit":joint_path_error,"error":trial.err().map(|e|e.to_string()),
                    "force":residual}));
                continue;
            }
            let from = progress;
            let linear_solves""")
casefile.write_text(s)
driver = probe/"examples/c01_methods.rs"
s = driver.read_text()
a = "    let run = stages::solve_observed(&mesh, ground, case, &study)?;"
replace(a,"""    let run = if let Some(path)=std::env::var_os("P8_RESUME") {
        stages::resume_observed(&mesh,ground,case,serde_json::from_reader(std::fs::File::open(path)?)?,&study)?
    } else { stages::solve_observed(&mesh, ground, case, &study)? };
    if let Some(path)=std::env::var_os("P8_PREFIX") {
        assert!(run.failure_reason.as_deref()==Some("checkpoint_requested"));
        serde_json::to_writer(std::fs::File::create_new(path)?, &run)?;
    }""")
driver.write_text(s)
