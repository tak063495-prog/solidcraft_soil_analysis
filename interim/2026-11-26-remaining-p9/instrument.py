"""Install P9-only switches in a disposable P8 probe."""
import json
from pathlib import Path
import sys

here = Path(__file__).resolve().parent
src = Path(sys.argv[1]) / "src"


def edit(file, old, new):
    path = src / file
    s = path.read_text()
    assert s.count(old) == 1, old
    path.write_text(s.replace(old,new))


enabled = [j in json.loads((here/"windows.json").read_text())["cells"] for j in range(256)]
edit("stages.rs","let next = if p8_predict {", """
            const P9_WINDOWS: [bool;256] = """ + json.dumps(enabled) + """;
            let p9_window = p8_reload && std::env::var_os("P9_WINDOWS").is_some()
                && P9_WINDOWS[((progress*256.).floor() as usize).min(255)];
            let next = if p8_predict || p9_window {""")
edit("stages.rs",'(std::env::var_os("P8_PREFIX").is_some()', '''(std::env::var_os("P9_HALVES").is_some() && stage.name()=="recontact" && next==1.)
                        || (std::env::var_os("P8_PREFIX").is_some()''')
edit("stages_joint_path_error.rs","        let trial = if control.resolve_slip_path", """
        if std::env::var_os("P9_HALVES").is_some() && coarse.state.stage_index <= 3 {
            let points: Vec<_> = [3,5,13].iter().map(|&i| {
                [a,b,d].map(|j| {
                    let p = &j.points[i];
                    serde_json::json!([p.status, p.plastic_slip_m])
                })
            }).collect();
            eprintln!("{}",serde_json::json!({"stage_index":coarse.state.stage_index,"accepted_candidate":[start,midpoint_progress,load.progress],
                                            "points":points}));
        }
        let trial = if control.resolve_slip_path""")
