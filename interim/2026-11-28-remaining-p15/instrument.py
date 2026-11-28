"""Install the two opt-in pretrial rules into the disposable p14 probe."""
from pathlib import Path
import sys

p = Path(sys.argv[1]) / "src/stages.rs"
s = p.read_text()
key = "if progress >= left && ![29,31].contains(&i) { continue }"
assert s.count(key) == 1
s = s.replace(key, """let added = std::env::var_os("P15_RELOAD").is_some() && [28,0].contains(&i);
                        if progress >= left && ![29,31].contains(&i) && !added { continue }""")
key = "        while progress < 1. - 1e-12 {\n            if !p8_parents.is_empty()"
assert s.count(key) == 1
s = s.replace(key, """        let outer_gap = stage.name()=="recontact" && std::env::var_os("P15_OUTER").is_some();
        let mut closing_parent: Option<(f64,[f64;32])> = None;
        let mut closing_prior = None;
""" + key)
key = "            let next = if p8_predict || p9_window || local {"
assert s.count(key) == 1
s = s.replace(key, """            let nominal = if outer_gap {
                // Grid cell with excess p14 rejections and quarter evaluations.
                let step = stage.step(&case.settings);
                let (left,right) = (952./1024.,953./1024.);
                let mut width_end = nominal;
                if progress >= left && progress < right {
                    let pts = &committed.joints.as_ref().unwrap().points;
                    let gaps: [f64;32] = std::array::from_fn(|i| {
                        let p = &pts[i];
                        p.normal.iter().zip(p.jump_m).map(|(a,b)| a*b).sum::<f64>()
                    });
                    if closing_parent.is_none_or(|(t,_)| t != progress) {
                        closing_prior = closing_parent;
                        closing_parent = Some((progress,gaps));
                    }
                    if let Some((t,old)) = closing_prior {
                        for j in 0..gaps.len() {
                            if old[j] > gaps[j] && gaps[j] > 0. {
                                let distance = gaps[j]*(progress-t)/(old[j]-gaps[j]);
                                width_end = width_end.min(progress + distance.max(step/512.));
                            }
                        }
                    }
                    eprintln!("{}",serde_json::json!({"outer_parent":progress,"attempt":count,
                        "prior":closing_prior,"gaps":gaps,"nominal":nominal,"to":width_end,"grid":step}));
                }
                width_end
            } else { nominal };
""" + key)
p.write_text(s)
